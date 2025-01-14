use crate::{
    atom::Atom,
    chirality::ChiralClass,
    molecular_formula::MolecularFormula,
    molecule::{
        bond::{BondOrder, BondTarget},
        node::Node,
        utils::{is_hueckel_satisfied,can_increase_bond},
        aromaticity::AromaticityType,
    },
};
use bitvec::vec::BitVec;
use bit_set::BitSet;
use chemistry_consts::ElementProperties;
use petgraph::{
    algo::subgraph_isomorphisms_iter,
    graph::{NodeIndex, UnGraph},
};
use tinyvec::{array_vec, ArrayVec};

use nohash_hasher::{IntMap, IntSet};
use std::collections::{HashMap, HashSet};

pub enum Chirality {
    None,
    CounterClockwise,
    Clockwise,
}


#[derive(Debug, Clone, Copy, Hash, PartialEq)]
pub struct AtomLabelImplicitHydrogens {
    atomic_number: u8,
    charge: i8,
    aromatic: bool,
    isotope: Option<u16>,
    chirality: ChiralClass,
    implicit_hydrogens: u8,
}

impl AtomLabelImplicitHydrogens {
    pub fn new(atomic_number: u8, charge: i8, aromatic: bool, isotope: Option<u16>, chirality: ChiralClass, implicit_hydrogens: u8) -> Self {
        Self { atomic_number, charge, aromatic, isotope, chirality, implicit_hydrogens }
    }
}


impl std::cmp::Eq for AtomLabelImplicitHydrogens {}

impl std::cmp::PartialOrd for AtomLabelImplicitHydrogens {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::cmp::Ord for AtomLabelImplicitHydrogens {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.atomic_number, self.charge, self.aromatic, self.isotope, self.chirality, self.implicit_hydrogens).cmp(&(other.atomic_number, other.charge, other.aromatic, other.isotope, other.chirality, other.implicit_hydrogens))
    }
}

#[derive(Debug, Clone, Copy, Hash, PartialEq)]
pub struct AtomLabel {
    atomic_number: u8, // This is also used for the bonds orders which are stored as u8, so we can use the same type for both
    charge: i8,
    aromatic: bool,
    isotope: Option<u16>,
    chirality: ChiralClass,
}


impl AtomLabel {
    pub fn new(atomic_number: u8, charge: i8, aromatic: bool, isotope: Option<u16>, chirality: ChiralClass) -> Self {
        Self { atomic_number, charge, aromatic, isotope, chirality }
    }
}

impl std::cmp::Eq for AtomLabel {}


impl std::cmp::PartialOrd for AtomLabel {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::cmp::Ord for AtomLabel {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.atomic_number, self.charge, self.aromatic,self.chirality, self.isotope).cmp(&(other.atomic_number, other.charge, self.aromatic, other.chirality, other.isotope))
    }
}


pub trait Molecule {
    fn atom_bonds(&self) -> &Vec<ArrayVec<[BondTarget; 10]>>;
    fn atom_bonds_mut(&mut self) -> &mut Vec<ArrayVec<[BondTarget; 10]>>;
    fn atomic_numbers(&self) -> &[u8];
    fn atomic_numbers_mut(&mut self) -> &mut Vec<u8>;
    fn chiral_classes(&self) -> Option<&[ChiralClass]>;
    fn chiral_classes_mut(&mut self) -> Option<&mut Vec<ChiralClass>>;
    fn charges(&self) -> &[i8];
    fn charges_mut(&mut self) -> &mut Vec<i8>;
    fn isotopes(&self) -> Option<&Vec<u16>>;
    fn atom_classes(&self) -> &Option<Vec<u8>>;
    fn atom_classes_mut(&mut self) -> &mut Option<Vec<u8>>;
    fn from_atoms(atoms: Vec<Atom>) -> Self;
    fn radical_states(&self) -> &[bool];
    fn radical_states_mut(&mut self) -> &mut Vec<bool>;
    fn pop_atom(&mut self) {
        self.atomic_numbers_mut().pop();
        self.charges_mut().pop();
        self.radical_states_mut().pop();
        self.atom_bonds_mut().pop();
        if let Some(classes) = self.atom_classes_mut() {
            classes.pop();
        }
        if let Some(chirals) = self.chiral_classes_mut() {
            chirals.pop();
        }
    }
    fn get_atomic_symbol(&self, atom_index: usize) -> Option<&str> {
        self.atomic_numbers().get(atom_index)?.atomic_symbol()
    }

    fn len(&self) -> usize {
        self.atomic_numbers().len()
    }
    fn is_empty(&self) -> bool {
        self.atomic_numbers().is_empty()
    }

    fn is_radical(&self) -> bool {
        self.radical_states().iter().any(|&is_radical| is_radical)
    }

    /// This function returns the valency_delta of an atom
    ///
    ///
    /// # Examples
    /// use molecules::prelude::*;
    /// let molecule = Molecule3D::from_smiles("C1CCCCC1CC2CCCCC2");
    /// assert_eq!(molecule.valency_delta(0), Some(0));
    fn valency_delta(&self, atom_index: usize) -> Option<i8> {
        let expected_valency = self.expected_valency(atom_index)?;
        let actual_valency = self.actual_valency(atom_index);
        Some(expected_valency - actual_valency)
    }

    fn get_atomic_number(&self, atom_index: usize) -> u8 {
        self.atomic_numbers()[atom_index]
    }

    fn set_atom_class(&mut self, atom_index: usize, class: u8) {
        if let Some(classes) = self.atom_classes_mut() {
            classes[atom_index] = class;
        } else {
            let size = self.atomic_numbers().len();
            let mut classes = vec![0; size];
            classes[atom_index] = class;
            *self.atom_classes_mut() = Some(classes);
        }
    }

    fn cmp_atom_charges(&self, atom_index1: usize, atom_index2: usize) -> std::cmp::Ordering {
        self.get_atom_charge(atom_index1)
            .cmp(&self.get_atom_charge(atom_index2))
    }
    /// Compares the electronegativities of two atoms based on the Pauling Scale, for detailed Values see the chemistry_consts crate
    ///
    /// # Arguments
    /// * 'atom1_index' - The index of the first atom.
    /// * 'atom2_index' - The index of the second atom.
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::molecule::{Molecule3D,Molecule};
    /// let molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    ///
    /// assert_eq!(molecule.cmp_electronegativities(0, 1), std::cmp::Ordering::Equal);
    /// ```
    fn cmp_electronegativities(
        &self,
        atom1_index: usize,
        atom2_index: usize,
    ) -> std::cmp::Ordering {
        self.electronegativity(atom1_index)
            .cmp(&self.electronegativity(atom2_index))
    }

    /// Returns the Pauling electronegativity of the atom at the specified index scaled by 100
    ///
    /// # Arguments
    /// * 'atom_index' - The index of the atom for which to retrieve the electronegativity.
    ///
    /// # Returns
    /// An `Option` containing the electronegativity value as a `u16` if the atom index is valid, `None` otherwise.
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::molecule::{Molecule3D,Molecule};
    /// let molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    /// assert_eq!(molecule.electronegativity(0), Some(254)); // Carbon's electronegativity
    /// assert_eq!(molecule.electronegativity(2), Some(220)); // Hydrogen's electronegativity
    /// ```
    fn electronegativity(&self, atom_index: usize) -> Option<u16> {
        self.atomic_numbers().get(atom_index)?.electronegativity()
    }

    /// This function returns the oxidation state of an atom
    ///
    /// # Arguments
    /// * 'atom_index' - The index of the atom
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::molecule::{Molecule3D,Molecule};
    /// let molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    /// assert_eq!(molecule.get_oxidation_state(0), -3);
    ///
    /// ```
    fn get_oxidation_state(&self, atom_index: usize) -> i8 {
        let mut state = self.get_atom_charge(atom_index);
        for bond in self.get_atom_bonds(atom_index).unwrap_or_default() {
            let target = bond.target();
            match self.cmp_electronegativities(atom_index, target) {
                std::cmp::Ordering::Less => state += 1,
                std::cmp::Ordering::Greater => state -= 1,
                _ => continue,
            }
        }
        state
    }
    /// This function returns the number of bonded elements of a specific type in case of invalid atom index it returns 0 to keep the api simple
    fn number_of_bonded_element(&self, atom_index: usize, element: u8) -> usize {
        let Some(bonds) = self.get_atom_bonds(atom_index) else {
            return 0;
        };
        bonds
            .iter()
            .filter(|bond| self.atomic_numbers()[bond.target()] == element)
            .count()
    }

    fn expected_valency(&self, atom_index: usize) -> Option<i8> {
        self.atomic_numbers()[atom_index].valencies()?.next()
    }

    /// This function returns the actual valency of an atom. In case of fractional valency it returns the rounded down value (e.g. fused aromatic systems are treated as aromatic) 
    ///
    /// # Arguments
    /// * 'atom_index' - The index of the atom
    ///
    /// # Returns
    /// The actual valency of the atom as an `i8`
    ///
    ///
    fn actual_valency(&self, atom_index: usize) -> i8 {
        (self.atom_bonds()[atom_index]
            .iter()
            .map(|bond| match bond.bond_order() {
                BondOrder::Single => 2,
                BondOrder::Double => 4,
                BondOrder::Triple => 6,
                BondOrder::Quadruple => 8,
                BondOrder::Aromatic => 3,
                BondOrder::Coordinate => 2,
            })
            .sum::<i8>()
            )
            / 2
    }
    fn degree(&self, atom_index: usize) -> Option<i8> {
        let expected_valency = self.expected_valency(atom_index)?;
        let actual_valency = self.actual_valency(atom_index);
        let charge = self.get_atom_charge(atom_index);
        let radical = self.is_atom_radical(atom_index) as i8;
        // May need to be changed for elements with unknown valencies
        Some(actual_valency - expected_valency + charge.abs() + radical)
    }

    fn degrees(&self) -> Vec<i8> {
        self.atomic_numbers()
            .iter()
            .enumerate()
            .map(|(index, _)| self.degree(index).unwrap_or(0))
            .collect()
    }

    fn monoisotopic_mass(&self) -> f64 {
        self.atomic_numbers()
            .iter()
            .map(|atom| atom.monoisotopic_mass().unwrap_or(0.0))
            .sum::<f64>()
    }

    fn number_of_atoms(&self) -> usize {
        self.atomic_numbers().len()
    }

    fn formal_charge(&self) -> i32 {
        self.charges().iter().sum::<i8>() as i32
    }

    fn get_atom_charge(&self, atom_index: usize) -> i8 {
        self.charges()[atom_index]
    }

    fn get_atom_charge_mut(&mut self, atom_index: usize) -> &mut i8 {
        &mut self.charges_mut()[atom_index]
    }

    fn get_isotope(&self, atom_index: usize) -> Option<u16> {
        self.isotopes()
            .and_then(|isotopes| isotopes.get(atom_index).copied())
    }

    fn get_chiral_class(&self, atom_index: usize) -> ChiralClass {
        self.chiral_classes()
            .and_then(|chirals| chirals.get(atom_index).copied())
            .unwrap_or(ChiralClass::None)
    }

    fn is_atom_radical(&self, atom_index: usize) -> bool;
    fn set_atom_radical(&mut self, atom_index: usize, is_radical: bool);
    fn set_atom_charge(&mut self, atom_index: usize, charge: i8) {
        if atom_index < self.charges().len() {
            self.charges_mut()[atom_index] = charge;
        } else {
            println!("Atom index out of bounds, could not set charge");
        }
    }
    fn get_atom_class(&self, atom_index: usize) -> u8 {
        let reference = self.atom_classes().as_ref();
        reference.map_or(0, |classes| classes[atom_index])
    }

    fn get_edges(&self) -> Vec<(usize, usize)> {
        self.atom_bonds()
            .iter()
            .enumerate()
            .flat_map(|(atom_index, bonds)| {
                bonds.iter().filter_map(move |bond| {
                    if atom_index < bond.target() {
                        Some((atom_index, bond.target()))
                    } else {
                        None
                    }
                })
            })
            .collect()
    }
    fn get_edges_with_type(&self) -> Vec<(usize, usize, BondOrder)> {
        self.atom_bonds()
            .iter()
            .enumerate()
            .flat_map(|(atom_index, bonds)| {
                bonds.iter().map(move |bond| {
                    if atom_index < bond.target() {
                        Some((atom_index, bond.target(), bond.bond_order()))
                    } else {
                        None
                    }
                })
            })
            .flatten()
            .collect()
    }

    fn to_ungraph(&self) -> UnGraph<u8, BondOrder> {
        let mut graph = UnGraph::<u8, BondOrder>::default();
        let mut node_indices = IntMap::<usize, NodeIndex<u32>>::default();

        for (index, &atomic_number) in self.atomic_numbers().iter().enumerate() {
            let node_index = graph.add_node(atomic_number);
            node_indices.insert(index, node_index);
        }

        for (index, bonds) in self.atom_bonds().iter().enumerate() {
            for bond in bonds {
                let source = node_indices[&index];
                let target = node_indices[&bond.target()];
                if source < target {
                    graph.add_edge(source, target, bond.bond_order());
                }
            }
        }
        graph
    }
    /// Returns a standardized string label for an atom containing all relevant atomic information
/// in a defined order: atomic_symbol, aromaticity, charge, hydrogens, etc.
///
/// # Arguments
/// * `mol` - The molecule containing the atom
/// * `atom_index` - The index of the atom to get the label for
///
/// # Returns
/// A string containing all atom information in a standardized format
///
/// # Example
/// ```
/// use molecules::prelude::*;
/// let mol = Molecule2D::from_smiles("[NH3+]").unwrap();
/// let mol = mol.first().unwrap();
/// assert_eq!(mol.get_atom_label(0), "N+1"); // Nitrogen, charge +1
/// ```
    fn get_atom_label(&self, atom_index: usize) -> AtomLabel {
        let atomic_number = self.atomic_numbers()[atom_index];
        // 2. Charge (c0 for neutral, c1 for +1, c-1 for -1, etc.)
        let charge = self.get_atom_charge(atom_index);
        let isotope = self.get_isotope(atom_index); 
        let is_aromatic = self.is_atom_aromatic(atom_index);
        let chirality = self.get_chiral_class(atom_index);

        AtomLabel::new(atomic_number,charge,is_aromatic,isotope,chirality)

    }

    fn get_atom_label_with_hydrogens(&self, atom_index: usize) -> AtomLabelImplicitHydrogens {
        let atomic_number = self.atomic_numbers()[atom_index];
        let charge = self.get_atom_charge(atom_index);
        let isotope = self.get_isotope(atom_index); 
        let is_aromatic = self.is_atom_aromatic(atom_index);
        let chirality = self.get_chiral_class(atom_index);
        let implicit_hydrogens = self.number_of_bonded_element(atom_index,1);
        AtomLabelImplicitHydrogens::new(atomic_number,charge,is_aromatic,isotope,chirality, implicit_hydrogens as u8)
    }


    fn to_ungraph_from_slice(&self, slice: &[usize]) -> UnGraph<AtomLabel, BondOrder> {
        let mut graph = UnGraph::<AtomLabel, BondOrder>::default();
        let mut node_indices = IntMap::<usize, NodeIndex<u32>>::default();

        for &index in slice {
            let node_index = graph.add_node(self.get_atom_label(index));
            node_indices.insert(index, node_index);
        }

        for &index in slice {
            let bonds = &self.atom_bonds()[index];
            for bond in bonds {
                let source = node_indices[&index];
                let target = node_indices[&bond.target()];
                if source < target && slice.contains(&bond.target()) {
                    graph.add_edge(source, target, bond.bond_order());
                }
            }
        }
        graph
    }

    fn match_submolecule(&self, other: &Self) -> Option<Vec<IntMap<usize, usize>>> {
        let mut self_components = self.get_components();
        let mut other_components = other.get_components();
        if self
            .atomic_numbers()
            .iter()
            .zip(other.atomic_numbers())
            .all(|(a, b)| a == b)
            && self
                .atom_bonds()
                .iter()
                .zip(other.atom_bonds())
                .all(|(a, b)| a == b)
            && self.len() == other.len()
        {
            return Some(vec![(0..self.number_of_atoms())
                .map(|index| (index, index))
                .collect()]);
        }

        self_components.retain(|component| component.len() > 1);
        other_components.retain(|component| component.len() > 1);

        for self_component in &mut self_components {
            let graph1 = self.to_ungraph_from_slice(self_component);
            let g_ref = &graph1;
            for other_component in &mut other_components {
                let graph2 = other.to_ungraph_from_slice(other_component);

                let h_ref = &graph2;

                if let Some(mappings) = subgraph_isomorphisms_iter(
                    &h_ref,
                    &g_ref,
                    &mut |node1, node2| node1 == node2,
                    &mut |edge1, edge2| edge1 == edge2,
                ) {
                    let mapping = mappings
                        .map(|mapping| {
                            mapping
                                .into_iter()
                                .zip(other_component.iter())
                                .map(|(index, &other_component_index)| {
                                    (self_component[index], other_component_index)
                                })
                                .collect::<IntMap<usize, usize>>()
                        })
                        .collect::<Vec<IntMap<usize, usize>>>();
                    if mapping.is_empty() {
                        continue;
                    }
                    return Some(mapping);
                };
            }
        }
        None
    }
    fn get_components(&self) -> Vec<Vec<usize>> {
        let mut connected_components = vec![];
        let mut visited_atoms = vec![false; self.atomic_numbers().len()];
        let mut stack = Vec::with_capacity(self.atomic_numbers().len());
        while let Some(index) = visited_atoms.iter().position(|&a| !a) {
            let mut component = self.traverse_component(&mut stack, index, &mut visited_atoms);
            component.sort();
            connected_components.push(component);
            stack.clear();
        }
        connected_components
    }

    fn get_atom_bonds(&self, atom_index: usize) -> Option<&[BondTarget]> {
        self.atom_bonds()
            .get(atom_index)
            .map(|bonds| bonds.as_slice())
    }
    fn traverse_component(
        &self,
        stack: &mut Vec<usize>,
        start: usize,
        visited_atoms: &mut [bool],
    ) -> Vec<usize> {
        let mut current_component = Vec::new();
        stack.push(start);
        visited_atoms[start] = true;

        while let Some(index) = stack.pop() {
            current_component.push(index);
            let Some(bonds) = self.get_atom_bonds(index) else {
                println!(
                    "No bonds found for atom {}, this means an out of bounds access",
                    index
                );
                continue;
            };
            for &bond in bonds {
                let target = bond.target();
                if !visited_atoms[target] {
                    stack.push(target);
                    visited_atoms[target] = true;
                }
            }
        }

        current_component
    }

    fn build_smiles_tree(
        &self,
        start_index: usize,
        bond_type: BondOrder,
        parent_index: Option<usize>,
        visited: &mut Vec<bool>,
        ring_closures: &mut Vec<(usize, usize, BondOrder)>,
        traversal_order: &mut Vec<usize>,
    ) -> Node {
        let mut root = Node::new(start_index, bond_type);
        visited[start_index] = true;
        traversal_order.push(start_index);  // Record traversal order
        
        let mut bonds = self.atom_bonds()[start_index];
        bonds.sort_by_key(|bond| bond.target());
        
        for neighbor in bonds {
            if Some(neighbor.target()) == parent_index {
                continue;
            }

            if self.atomic_numbers()[neighbor.target()] == 1 {
                continue;
            }

            if visited[neighbor.target()] {
                if neighbor.target() < start_index {
                    continue;
                }
                let bond_type = neighbor.bond_order();
                ring_closures.push((start_index, neighbor.target(), bond_type));
            } else {
                let child = self.build_smiles_tree(
                    neighbor.target(),
                    neighbor.bond_order(),
                    Some(start_index),
                    visited,
                    ring_closures,
                    traversal_order,
                );
                root.add_child(child);
            }
        }
        root
    }

    fn number_of_pi_electrons(&self, atom_index: usize) -> usize {
        let mut number_of_pi_electrons = 0;
        // Default to 0 if no bonds are found
        let Some(bonds) = self.get_atom_bonds(atom_index) else {
            return 0;
        };
        for bond in bonds {
            if bond.bond_order() == BondOrder::Double || bond.bond_order() == BondOrder::Aromatic {
                number_of_pi_electrons += 1;
            } else if bond.bond_order() == BondOrder::Triple {
                number_of_pi_electrons += 2;
            }
        }
        match self.get_atom_charge(atom_index) {
            2.. => number_of_pi_electrons += 1,
            ..=-2 => number_of_pi_electrons -= 1,
            _ => {}
        }
        match self.atomic_numbers()[atom_index] {
            7 => number_of_pi_electrons += 2,
            8 => number_of_pi_electrons += 2,
            _ => {}
        }
        number_of_pi_electrons
    }

    /// Checks if the molecule has any stereochemistry
    ///
    /// # Returns
    /// A boolean indicating whether the molecule has stereochemistry   
    ///
    /// # Example
    /// ```
    /// use molecules::prelude::*;
    /// let molecule = Molecule2D::from_smiles("C[C@H](N)O").unwrap();
    /// let molecule = molecule.first().unwrap();
    /// assert_eq!(molecule.has_stereochemistry(), true);
    /// ```
    fn has_stereochemistry(&self) -> bool {
        self.chiral_classes().map_or(false, |classes| classes.iter().any(|class| class != &ChiralClass::None))
    }

    fn molecular_formula(&self) -> MolecularFormula {
        MolecularFormula::from_molecule(self)
    }

    /// Converts the molecule to an undirected graph where both atoms and bonds are represented as nodes
    ///
    /// # Returns
    /// An undirected graph where:
    /// - Atom nodes are labeled with their atomic symbol (e.g., "C", "N", "O")
    /// - Bond nodes are labeled with the bond order (e.g., "-", "=", "#")
    /// - Edges connect atoms to their bonds (no direct atom-to-atom connections)
    ///
    /// # Example
    /// ```
    /// use molecules::prelude::*;
    /// let molecules = Molecule2D::from_smiles("C=C").unwrap();
    /// let molecule = molecules.first().unwrap();
    /// let graph = molecule.to_ungraph_with_edge_nodes();
    /// assert_eq!(graph.node_count(), 11); // 2 carbon atoms + 1 double bond
    /// assert_eq!(graph.edge_count(), 10); // 2 edges connecting atoms to bond
    /// ```
    fn to_ungraph_with_edge_nodes(&self) -> UnGraph<AtomLabel, ()> {
        let mut graph = UnGraph::<AtomLabel, ()>::default();
        let mut atom_indices = IntMap::<usize, NodeIndex<u32>>::default();
        let mut bond_indices = HashMap::<(usize, usize), NodeIndex<u32>>::default();

        // Add atom nodes
        for (index, _) in self.atomic_numbers().iter().enumerate() {
            let label = self.get_atom_label(index);
            let node_index = graph.add_node(label);
            atom_indices.insert(index, node_index);
        }

        // Add bond nodes and connect them to atoms
        for (index, bonds) in self.atom_bonds().iter().enumerate() {
            for bond in bonds {
                let target = bond.target();
                if index < target {
                    let bond_order = bond.bond_order();
                    let source_atom = atom_indices[&index];
                    let target_atom = atom_indices[&target];
                    if bond_order == BondOrder::Single || bond_order == BondOrder::Aromatic {
                        graph.add_edge(source_atom, target_atom, ());
                        continue;
                    }
                    // Only process each bond once
                    // Create bond node
                    let bond_label = AtomLabel::new(bond.bond_order().to_u8(), 0,false, None, ChiralClass::None);
                    let bond_node = graph.add_node(bond_label);
                    bond_indices.insert((index, target), bond_node);

                    // Connect atoms to bond
                    graph.add_edge(source_atom, bond_node, ());
                    graph.add_edge(target_atom, bond_node, ());
                }
            }
        }

        graph
    }

    fn to_ungraph_with_edge_nodes_and_implicit_hydrogens(&self) -> UnGraph<AtomLabelImplicitHydrogens, ()> {
        let mut graph = UnGraph::<AtomLabelImplicitHydrogens, ()>::default();
        let mut atom_indices = IntMap::<usize, NodeIndex<u32>>::default();
        let mut bond_indices = HashMap::<(usize, usize), NodeIndex<u32>>::default();

        // Add atom nodes
        for (index, _) in self.atomic_numbers().iter().enumerate() {
            if self.atomic_numbers()[index] == 1 {
                continue;
            }
            let label = self.get_atom_label_with_hydrogens(index);
            let node_index = graph.add_node(label);
            atom_indices.insert(index, node_index);
        }

        // Add bond nodes and connect them to atoms
        for (index, bonds) in self.atom_bonds().iter().enumerate() {
            for bond in bonds {
                let target = bond.target();
                if self.atomic_numbers()[index] == 1 || self.atomic_numbers()[target] == 1 {
                    continue;
                }
                if index < target {
                    let bond_order = bond.bond_order();
                    let source_atom = atom_indices[&index];
                    let target_atom = atom_indices[&target];
                    if bond_order == BondOrder::Single || bond_order == BondOrder::Aromatic {
                        graph.add_edge(source_atom, target_atom, ());
                        continue;
                    }
                    // Only process each bond once
                    // Create bond node
                    let bond_label = AtomLabelImplicitHydrogens::new(bond.bond_order().to_u8(), 0, false, None, ChiralClass::None, 0);
                    let bond_node = graph.add_node(bond_label);
                    bond_indices.insert((index, target), bond_node);

                    // Connect atoms to bond
                    graph.add_edge(source_atom, bond_node, ());
                    graph.add_edge(target_atom, bond_node, ());
                }
            }
        }

        graph
    }


    fn to_smiles(&self) -> String {
        if self.atomic_numbers().is_empty() {
            return String::new();
        }
        let mut smiles = String::new();
        let components = self.get_components();

        for component in components {
            let rings = &mut Vec::new();
            let mut visited = vec![false; self.atomic_numbers().len()];
            // Track traversal order of atoms
            let mut traversal_order = Vec::with_capacity(self.atomic_numbers().len());

            let root_index = component.iter().filter(|&&x| self.atomic_numbers()[x] != 1).next().unwrap_or(&component[0]);

            let root = self.build_smiles_tree(
                component[0],
                BondOrder::Single,
                None,
                &mut visited,
                rings,
                &mut traversal_order,
            );

            // Sort rings based on when their first atom was encountered in traversal
            let mut ring_closures: IntMap<usize, Vec<(usize, BondOrder)>> = IntMap::default();
            let mut ring_number = 1;
            
            // Process rings in order of first appearance in traversal
            for &atom_idx in &traversal_order {
                for (start, end, bond_order) in rings.iter() {
                    if *start == atom_idx {
                        ring_closures.entry(*start).or_default().push((ring_number, *bond_order));
                        ring_closures.entry(*end).or_default().push((ring_number, *bond_order));
                        ring_number += 1;
                    }
                }
            }
            //println!("Ring closures: {:?}", ring_closures);

            if !smiles.is_empty() {
                smiles.push('.');
            }

            smiles.push_str(&self.construct_smiles(&root, &ring_closures))
        }
        relabel_numbers(&smiles)
    }

    fn to_smiles_with_implicit_hydrogens(&self) -> String {
        if self.atomic_numbers().is_empty() {
            return String::new();
        }
        let mut smiles = String::new();
        let components = self.get_components();

        for component in components {
            let rings = &mut Vec::new();
            let mut visited = vec![false; self.atomic_numbers().len()];
            // Track traversal order of atoms
            let mut traversal_order = Vec::with_capacity(self.atomic_numbers().len());

            let root_index = component.iter().filter(|&&x| self.atomic_numbers()[x] != 1).next().unwrap_or(&component[0]);


            let root = self.build_smiles_tree(
                *root_index,
                BondOrder::Single,
                None,
                &mut visited,
                rings,
                &mut traversal_order,
            );

            // Sort rings based on when their first atom was encountered in traversal
            let mut ring_closures: IntMap<usize, Vec<(usize, BondOrder)>> = IntMap::default();
            let mut ring_number = 1;
            
            // Process rings in order of first appearance in traversal
            for &atom_idx in &traversal_order {
                for (start, end, bond_order) in rings.iter() {
                    if *start == atom_idx {
                        ring_closures.entry(*start).or_default().push((ring_number, *bond_order));
                        ring_closures.entry(*end).or_default().push((ring_number, *bond_order));
                        ring_number += 1;
                    }
                }
            }
            //println!("Ring closures: {:?}", ring_closures);

            if !smiles.is_empty() {
                smiles.push('.');
            }

            smiles.push_str(&self.construct_smiles(&root, &ring_closures))
        }
        relabel_numbers(&smiles)
    }

    /// Constructs the SMILES string recursively from the molecule's tree structure.
    ///
    /// # Arguments
    ///
    /// * `node` - The current node in the SMILES tree.
    /// * `ring_closures` - A map of ring closures for the molecule.
    ///
    /// # Returns
    ///
    /// A `String` representing the SMILES notation of the molecule.
    ///
    /// # Example
    ///
    /// ```rust
    /// use molecules::molecule::{Molecule2D, Molecule};
    ///
    /// let molecule = Molecule2D::from_smiles("c1ccccc1").unwrap();
    /// let molecule = molecule.first().unwrap();
    /// let smiles = molecule.to_smiles();
    /// assert_eq!(smiles, "c1ccccc1");
    /// ```
    fn construct_smiles(
        &self,
        node: &Node,
        ring_closures: &IntMap<usize, Vec<(usize, BondOrder)>>,
    ) -> String {
        let mut smiles = String::new();

        // Determine if the current bond is aromatic
        let is_current_bond_aromatic = matches!(node.bond_type(), BondOrder::Aromatic);

        // Add bond type symbol only if the bond is not aromatic
        if !is_current_bond_aromatic {
            match node.bond_type() {
                BondOrder::Single => {}
                BondOrder::Double => smiles.push('='),
                BondOrder::Triple => smiles.push('#'),
                BondOrder::Quadruple => smiles.push('$'),
                BondOrder::Aromatic => smiles.push(':'),
                BondOrder::Coordinate => {}
            }
        }

        let charge = self.charges()[node.index()];
        let number_of_hydrogens = self.number_of_bonded_element(node.index(), 1);
        let atomic_number = self.atomic_numbers()[node.index()];
        if atomic_number == 1 {

        // Add children
        if node.children().len() == 1 {
            smiles.push_str(&self.construct_smiles(&node.children()[0], ring_closures));
        } else {
            let length = node.children().len();
            for (index,child) in node.children().iter().enumerate() {
                if length - 1 == index {
                    // The terminal child does not need braces.
                    smiles.push_str(&self.construct_smiles(child, ring_closures));
                    continue;
                }
                smiles.push('(');
                smiles.push_str(&self.construct_smiles(child, ring_closures));
                smiles.push(')');
            }
        } return smiles;
        }

        // Get atomic symbol
        let Some(atomic_symbol) = atomic_number.atomic_symbol() else {
            return "".to_string();
        };
        let mut atomic_symbol = atomic_symbol.to_string();

        // Check for aromaticity: if any bond of the atom is aromatic, make the symbol lowercase
        let is_aromatic = self
            .get_atom_bonds(node.index())
            .unwrap()
            .iter()
            .any(|bond| bond.bond_order() == BondOrder::Aromatic);
        if is_aromatic {
            atomic_symbol = atomic_symbol.to_lowercase();
        }

        let is_charged = charge != 0;
        let is_organic_subset = matches!(atomic_number, 6 | 7 | 8 | 9 | 15| 16 | 17 | 35 | 53);

        let has_standard_hydrogens = match atomic_number {
            6 => number_of_hydrogens <= 4,       // Carbon
            7 => number_of_hydrogens <= 2,       // Nitrogen
            8 => number_of_hydrogens <= 2,       // Oxygen
            9 => number_of_hydrogens <= 1,       // Fluorine
            15 => number_of_hydrogens <= 2,      // Phosphorus
            16 => number_of_hydrogens <= 2,      // Sulfur
            17 | 35 | 53 => number_of_hydrogens <= 1, // Chlorine/Bromine/Iodine
            _ => false,
        };

        // Determine if we need square brackets
        let needs_brackets = !(!is_charged && is_organic_subset && has_standard_hydrogens);
        // Build atom string
        if needs_brackets {
            let hydrogen_str = match number_of_hydrogens {
                0 => String::new(),
                1 => "H".to_string(),
                _ => format!("H{}", number_of_hydrogens),
            };

            let charge_string = match charge {
                2.. => format!("+{}", charge.abs()),
                ..=-2 => format!("-{}", charge.abs()),
                1 => "+".to_string(),
                -1 => "-".to_string(),
                0 => String::new(),
            };

            smiles.push_str(&format!(
                "[{}{}{}",
                atomic_symbol, hydrogen_str, charge_string
            ));
            smiles.push(']');
        } else {
            smiles.push_str(&atomic_symbol);
        }

        // Add ring closures
        if let Some(closures) = ring_closures.get(&node.index()) {
            for closure in closures {
                if closure.1 == BondOrder::Single || closure.1 == BondOrder::Aromatic {
                    smiles.push_str(&format!("{}", closure.0));
                } else {
                    smiles.push_str(&format!("{}{}", closure.1, closure.0));
                }
            }
        }

        // Add children
        if node.children().len() == 1 {
            smiles.push_str(&self.construct_smiles(&node.children()[0], ring_closures));
        } else {
            let length = node.children().len();
            for (index,child) in node.children().iter().enumerate() {
                if length - 1 == index {
                    // The terminal child does not need braces.
                    smiles.push_str(&self.construct_smiles(child, ring_closures));
                    continue;
                }
                smiles.push('(');
                smiles.push_str(&self.construct_smiles(child, ring_closures));
                smiles.push(')');
            }
        }

        smiles
    }

    fn is_atom_aromatic(&self, atom_index: usize) -> bool {
        self.get_atom_bonds(atom_index).unwrap().iter().any(|bond| bond.bond_order() == BondOrder::Aromatic)
    }

    /// This function returns the rings in the molecule
    ///
    /// # Examples
    /// use molecules::prelude::*;
    /// let molecule = Molecule3D::from_smiles("C1CCCCC1CC2CCCCC2");
    ///
    /// println!("{molecule.find_rings()}");
    /// assert_eq!(molecule.find_rings(), vec![vec![0, 1, 2, 3, 4, 5],vec![6, 7, 8, 9, 10, 11]]);
    fn find_rings(&self) -> Vec<Vec<usize>> {
        let mut rings = Vec::new();
        let mut visited = vec![false; self.atomic_numbers().len()];
        let mut path = Vec::new();
        let mut path_set = HashSet::new();

        for start in 0..self.atomic_numbers().len() {
            if !visited[start] {
                self.dfs_find_rings(
                    start,
                    start,
                    &mut visited,
                    &mut path,
                    &mut path_set,
                    &mut rings,
                );
            }
        }

        rings
    }

    fn dfs_find_rings(
        &self,
        current: usize,
        parent: usize,
        visited: &mut Vec<bool>,
        path: &mut Vec<usize>,
        path_set: &mut HashSet<usize>,
        rings: &mut Vec<Vec<usize>>,
    ) {
        visited[current] = true;
        path.push(current);
        path_set.insert(current);

        for &BondTarget { target, .. } in &self.atom_bonds()[current] {
            if !visited[target] {
                self.dfs_find_rings(target, current, visited, path, path_set, rings);
            } else if target != parent && path_set.contains(&target) {
                // Found a cycle
                let cycle_start = path.iter().position(|&x| x == target).unwrap();
                let cycle = path[cycle_start..].to_vec();
                rings.push(cycle);
            }
        }

        path.pop();
        path_set.remove(&current);
    }

    fn add_hydrogens(&mut self) {
        let degrees = self.degrees();
        degrees.iter().enumerate().for_each(|(index, &degree)| {
            let number_of_hydrogens = -degree;
            let atomic_number = self.atomic_numbers()[index];
            (0..number_of_hydrogens).for_each(|_| {
                self.charges_mut().push(0);
                self.radical_states_mut().push(false)
            });

            if atomic_number == 1 {
                return;
            }

            //println!("Adding {} hydrogens to atom {}", number_of_hydrogens, index);
            for _ in 0..number_of_hydrogens {
                let hydrogen_index = self.atomic_numbers().len();
                self.atomic_numbers_mut().push(1);
                self.atom_bonds_mut()
                    .push(array_vec!([BondTarget;10] => BondTarget::single(index)));
                self.atom_bonds_mut()[index].push(BondTarget::single(hydrogen_index));
            }
        });
    }

    /// Returns true if the molecule contains explicitly defined hydrogen atoms
    fn has_explicit_hydrogens(&self) -> bool {
        self.atomic_numbers().iter().any(|&num| num == 1)
    }

    fn get_aromatic_atoms(&self) -> Vec<usize> {
        self.atom_bonds()
            .iter()
            .enumerate()
            .filter_map(|(index, bonds)| {
                if bonds
                    .iter()
                    .any(|bond| bond.bond_order == BondOrder::Aromatic)
                {
                    Some(index)
                } else {
                    None
                }
            })
            .collect()
    }

    fn change_bond_order(&mut self, atom1: usize, atom2: usize, new_order: BondOrder) {
        let bond1 = self.atom_bonds_mut()[atom1]
            .iter()
            .position(|bond| bond.target == atom2)
            .unwrap();
        let bond2 = self.atom_bonds_mut()[atom2]
            .iter()
            .position(|bond| bond.target == atom1)
            .unwrap();
        self.atom_bonds_mut()[atom1][bond1].bond_order = new_order;
        self.atom_bonds_mut()[atom2][bond2].bond_order = new_order;
    }

    /// Kekulizes the molecule by adding double bonds to aromatic atoms
    ///
    /// # Errors
    ///
    /// Returns an error if the molecule cannot be kekulized   
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::prelude::*;
    /// let molecule = Molecule2D::from_smiles("c1ccccc1").unwrap();
    /// let mut molecule = molecule.first().unwrap().clone();
    /// molecule.kekulize().unwrap();
    /// assert_eq!(molecule.to_smiles(), "C1=CC=CC=C1");
    /// let molecule = Molecule2D::from_smiles("c1ccccn1").unwrap();
    /// let mut molecule = molecule.first().unwrap().clone();
    /// molecule.kekulize().unwrap();
    /// assert_eq!(molecule.to_smiles(), "C1=CC=CC=N1");
    /// ```
    fn kekulize(&mut self) -> Result<(), String> {
        let mut aromatic_atoms: BitVec<usize> = BitVec::repeat(false, self.atom_bonds().len());
        let mut aromatic_rings: Vec<Vec<usize>> = Vec::new();
        let mut visited = BitVec::repeat(false, self.atomic_numbers().len());

        // Find aromatic atoms and rings
        for (index, bonds) in self.atom_bonds().iter().enumerate() {
            for bond in bonds {
                if bond.bond_order() == BondOrder::Aromatic {
                    aromatic_atoms.set(index, true);
                    break;
                }
            }
        }

        // Early return if no aromatic atoms
        if !aromatic_atoms.any() {
            return Ok(());
        }

        // Find aromatic rings using DFS
        let mut path = Vec::new();
        let mut path_set = HashSet::new();

        for start in 0..self.atomic_numbers().len() {
            if aromatic_atoms[start] && !visited[start] {
                self.dfs_find_aromatic_rings(
                    start,
                    start,
                    &mut visited,
                    &mut path,
                    &mut path_set,
                    &mut aromatic_rings,
                    &aromatic_atoms,
                );
            }
        }


        // Validate rings follow Hückel's rule
        //aromatic_rings.retain(|ring| self.is_valid_aromatic_cycle(ring));

        let mut is_subset = Vec::new();
        for outer_ring in &aromatic_rings {
            for (index, inner_ring) in aromatic_rings.iter().enumerate() {
                if inner_ring.len() < outer_ring.len() && inner_ring.iter().all(|atom_index| outer_ring.contains(atom_index)) {
                    is_subset.push(index)
                }
            }
        }

        for index in is_subset {
            aromatic_rings.remove(index);
        }

        // Convert aromatic bonds to alternating single/double bonds
        for ring in &aromatic_rings {
            let ring_len = ring.len();

            // First convert all bonds to single bonds
            for atom_index in ring {
                let atom_bonds = self.get_atom_bonds(*atom_index).unwrap().to_vec();  
                for bond in atom_bonds {
                    if bond.bond_order() == BondOrder::Aromatic {
                        self.change_bond_order(*atom_index, bond.target(), BondOrder::Single);
                    }
                }
            }

            let degrees = self.degrees();

            // Then add double bonds in alternating pattern
            let mut last_bond_was_double = false;
            for i in 0..ring_len {
                let atom1 = ring[i];
                let atom2 = ring[(i + 1) % ring_len];
                if can_increase_bond(atom1, atom2, &degrees) && !last_bond_was_double {
                    self.change_bond_order(atom1, atom2, BondOrder::Double);
                    last_bond_was_double = true;
                } else {
                    last_bond_was_double = false;
                }
            }
        }

        Ok(())
    }

    fn can_kekulize_bond(
        &self,
        atom1_idx: usize,
        atom2_idx: usize,
        aromatic_atoms: &BitVec<usize>,
    ) -> bool {
        aromatic_atoms[atom1_idx] && aromatic_atoms[atom2_idx]
    }

    // Helper methods
    fn get_bond_order(&self, atom1: usize, atom2: usize) -> BondOrder {
        self.atom_bonds()[atom1]
            .iter()
            .find(|bond| bond.target() == atom2)
            .map(|bond| bond.bond_order())
            .unwrap_or(BondOrder::Single)
    }

    fn dfs_find_aromatic_rings(
        &self,
        current: usize,
        parent: usize,
        visited: &mut BitVec<usize>,
        path: &mut Vec<usize>,
        path_set: &mut HashSet<usize>,
        rings: &mut Vec<Vec<usize>>,
        aromatic_atoms: &BitVec<usize>,
    ) {
        visited.set(current, true);
        path.push(current);
        path_set.insert(current);

        for &bond in &self.atom_bonds()[current] {
            let target = bond.target();
            if !aromatic_atoms[target] {
                continue;
            }

            if !visited[target] {
                self.dfs_find_aromatic_rings(
                    target,
                    current,
                    visited,
                    path,
                    path_set,
                    rings,
                    aromatic_atoms,
                );
            } else if target != parent && path_set.contains(&target) {
                // Found aromatic cycle
                let cycle_start = path.iter().position(|&x| x == target).unwrap();
                let cycle = path[cycle_start..].to_vec();
                if cycle.len() >= 5 {
                    // Common aromatic ring sizes
                    rings.push(cycle);
                }
            }
        }

        path.pop();
        path_set.remove(&current);
    }


    fn dfs_find_potential_aromatic_cycles(
        &self,
        current: usize,
        root: usize,
        visited: &mut Vec<bool>,
        current_path: &mut Vec<usize>,
        cycles: &mut Vec<Vec<usize>>,
        atoms_in_cycles: &mut IntSet<usize>,
        last_bond: BondOrder,
    ) {
        visited[current] = true;
        current_path.push(current);

        for &bond in &self.atom_bonds()[current] {
            let target = bond.target();

            let is_last_bond_double = last_bond == BondOrder::Double;
            let is_current_bond_double = bond.bond_order() == BondOrder::Double;

            let number_of_bonds = self.atom_bonds()[current].len();

            // Skip if this would create an invalid alternating pattern
            if (is_current_bond_double && is_last_bond_double) || number_of_bonds > 3 {
                continue;
            }

            // Found cycle back to root
            if target == root && current_path.len() >= 3 {
                // Only add cycle if it follows alternating pattern
                if bond.bond_order() != last_bond {
                    let cycle = current_path.clone();
                    // Add all atoms in cycle to tracked set
                    for &atom in &cycle {
                        atoms_in_cycles.insert(atom);
                    }
                    cycles.push(cycle);
                }
                continue;
            }

            // Continue DFS if unvisited
            if !visited[target] {
                self.dfs_find_potential_aromatic_cycles(
                    target,
                    root,
                    visited,
                    current_path,
                    cycles,
                    atoms_in_cycles,
                    bond.bond_order(),
                );
            }
        }

        current_path.pop();
        visited[current] = false;
    }

    fn count_pi_electrons_in_cycle(&self, cycle: &[usize]) -> usize {
        cycle
            .iter()
            .map(|&atom_idx| {
                let atomic_number = self.atomic_numbers()[atom_idx];
                let charge = self.get_atom_charge(atom_idx);
                match atomic_number {
                    7 => {
                        // Nitrogen
                        if charge == 0 {
                            2 // Neutral N contributes 2 electrons
                        } else if charge == 1 {
                            1 // N+ contributes 1 electron
                        } else {
                            3 // N- contributes 3 electrons
                        }
                    }
                    8 => 2,  // Oxygen contributes 2 electrons
                    16 => 2, // Sulfur contributes 2 electrons
                    15 => {
                        // Phosphorus
                        if charge == 0 {
                            2
                        } else {
                            0
                        }
                    }
                    6 => 1, // Carbon contributes 1 electron
                    _ => 0, // Other atoms contribute 0 electrons
                }
            })
            .sum()
    }

    fn is_valid_aromatic_cycle(&self, cycle: &[usize]) -> bool {
        // Check cycle size (most common aromatic rings are 5 or 6 membered)
        if cycle.len() < 5 {
            return false;
        }

        // Count pi electrons
        let pi_electrons = self.count_pi_electrons_in_cycle(cycle);

        // Check Hückel's rule (4n + 2)
        is_hueckel_satisfied(pi_electrons)
    }

    fn can_be_aromatic(&self, atom_idx: usize) -> bool {
        let atomic_number = self.atomic_numbers()[atom_idx];
        match atomic_number {
            6 => {
                let valence = self.actual_valency(atom_idx);
                valence <= 3 // C should have 3 or fewer bonds
            }
            7 => {
                // Nitrogen
                let valence = self.actual_valency(atom_idx);
                valence <= 3 // N should have 3 or fewer bonds
            }
            8 => {
                // Oxygen
                let valence = self.actual_valency(atom_idx);
                valence <= 2 // O should have 2 or fewer bonds
            }
            16 => {
                // Sulfur
                let valence = self.actual_valency(atom_idx);
                valence <= 2 // S should have 2 or fewer bonds in aromatic systems
            }
            15 => {
                // Phosphorus
                let valence = self.actual_valency(atom_idx);
                valence <= 3 // P should have 3 or fewer bonds
            }
            _ => false,
        }
    }

    fn aromatize(&mut self) { 
        let mut cycles: Vec<BitSet> = Vec::new();

        let mut atoms_in_cycles = IntSet::default();

        // First pass: Find all potential cycles
        self.find_all_conjugated_cycles(&mut cycles, &mut atoms_in_cycles);

        // Group and analyze systems
        let mut fused_systems = self.group_fused_cycles(&cycles);

        // Order based on size in descending order
        fused_systems.sort_by(|a, b| b.len().cmp(&a.len()));


        // Process each system
        for system in fused_systems {
            let aromaticity_type = self.determine_aromaticity_type(&system);
            match aromaticity_type {
                AromaticityType::Aromatic => self.aromatize_system(&system),
                AromaticityType::MoebiusAromatic => (), //self.aromatize_moebius_system(&system),
                AromaticityType::AntiAromatic => continue,
                AromaticityType::Metallocenic => (), //self.aromatize_metallocenic_system(&system),
                AromaticityType::NonAromatic => continue,
            }
        }
    }

    fn determine_aromaticity_type(&self, system: &BitSet) -> AromaticityType {
        let pi_electrons = self.count_system_pi_electrons(system);
        let has_metal_center = self.check_for_metal_center(system);
        let is_moebius = self.check_moebius_topology(system);

        match (pi_electrons, has_metal_center, is_moebius) {
            (e, false, false) if is_hueckel_satisfied(e) => AromaticityType::Aromatic,
            (e, false, true) if (e % 4 == 2) => AromaticityType::MoebiusAromatic,
            (e, false, false) if e % 4 == 0 => AromaticityType::AntiAromatic,
            (_, true, _) => AromaticityType::Metallocenic,
            _ => AromaticityType::NonAromatic,
        }
    }

    fn count_system_pi_electrons(&self, system: &BitSet) -> usize {
        let mut total = 0;
        let mut counted_atoms = IntSet::default();

        for atom_idx in system.iter() {
            if counted_atoms.insert(atom_idx) {
                total += self.count_atom_pi_electrons_extended_within_system(atom_idx, system);
            }
        }

        total
    }

    fn count_atom_pi_electrons_extended_within_system(&self, atom_idx: usize, system: &BitSet) -> usize {
        let mut count = 0;
        let atomic_number = self.atomic_numbers()[atom_idx];
        let charge = self.get_atom_charge(atom_idx);
        let is_radical = self.is_atom_radical(atom_idx);

        if self.is_atom_aromatic(atom_idx) {
            return 1;
        }

        // Add bond contributions
        for bond in self.get_atom_bonds(atom_idx).unwrap_or_default() {
            if !system.contains(bond.target()) {
                if bond.bond_order() == BondOrder::Double {
                    return 0;
                }
                continue;
            }
            count += match bond.bond_order() {
                BondOrder::Double => 2,
                BondOrder::Aromatic => if atomic_number == 6 {1} else {0},
                _ => 0,
            };
        }

        count /= 2;

        // Base pi electrons from atomic configuration, considering charge
        count += match atomic_number {
            6 if count == 0 => match charge {  // Carbon
                0 => 0,
                1 => 0,  // C+ has no free electrons
                -1 => 2, // C- has two extra electrons
                _ if charge > 1 => 0,

                _ => 2,  // C2- or lower has two extra electrons
            },
            7 if count == 0 => match charge {  // Nitrogen
                1 => 1,  // N+ has one free electron
                0 => {
                    if self.get_atom_bonds(atom_idx).unwrap_or_default().len() > 2 {
                        2
                    } else {
                        1
                    }
                },  // Neutral N has two free electrons
                -1 => 3, // N- has three free electrons
                _ if charge > 1 => 0,
                _ => 4,  // N2- or lower has four free electrons
            },
            8 | 16 => match charge {  // Oxygen, Sulfur
                1 => 1,  // O+/S+ has one free electron
                0 => 2,  // Neutral O/S has two free electrons
                -1 => 3, // O-/S- has three free electrons
                _ if charge > 1 => 0,
                _ => 4,  // O2-/S2- or lower has four free electrons
            },
            15 => match charge {  // Phosphorus
                1 => 0,  // P+ has no free electrons
                0 => 2,  // Neutral P has two free electrons
                -1 => 3, // P- has three free electrons
                _ if charge > 1 => 0,
                _ => 4,  // P2- or lower has four free electrons
            },
            _ => 0,
        };

        // Add radical contribution
        if is_radical {
            count += 1;
        }

        count
    }

    fn check_moebius_topology(&self, system: &BitSet) -> bool {
        false // Placeholder
    }

    fn find_all_conjugated_cycles(&self, cycles: &mut Vec<BitSet>, atoms_in_cycles: &mut IntSet<usize>) {
        let mut visited = vec![false; self.atomic_numbers().len()];
        let mut path = Vec::new();
        let mut path_set = HashSet::new();

        for start in 0..self.atomic_numbers().len() {
            if !atoms_in_cycles.contains(&start) && self.can_be_conjugated(start) {
                self.dfs_find_conjugated_cycles(
                    start,
                    start,
                    &mut visited,
                    &mut path,
                    &mut path_set,
                    cycles,
                    atoms_in_cycles,
                );
            }
        }
    }

    fn dfs_find_conjugated_cycles(
        &self,
        current: usize,
        parent: usize,
        visited: &mut Vec<bool>,
        path: &mut Vec<usize>,
        path_set: &mut HashSet<usize>,
        cycles: &mut Vec<BitSet>,
        atoms_in_cycles: &mut IntSet<usize>,
    ) {
        visited[current] = true;
        path.push(current);
        path_set.insert(current);

        for &bond in &self.atom_bonds()[current] {
            let target = bond.target();
            
            // Skip parent and non-conjugated atoms
            if target == parent || !self.can_be_conjugated(target) {
                continue;
            }

            // Found a cycle
            if path_set.contains(&target) {
                if target == path[0] && path.len() >= 5 {
                    // Valid cycle found (minimum size 5)
                    let cycle = path.clone();
                    if self.is_potentially_aromatic_cycle(&cycle) {
                        let cycle_bitset = cycle_to_bitset(&cycle);
                        if !cycles.contains(&cycle_bitset) {
                            cycles.push(cycle_bitset);
                        }

                        for atom in path.iter() {
                            atoms_in_cycles.insert(*atom);
                        }
                    }
                }
            } else if !visited[target] {
                self.dfs_find_conjugated_cycles(
                    target,
                    current,
                    visited,
                    path,
                    path_set,
                    cycles,
                    atoms_in_cycles,
                );
            }
        }

        path.pop();
        path_set.remove(&current);
        visited[current] = false;
    }

    fn is_valid_aromatic_system(&self, system: &BitSet) -> bool {
        // Count total pi electrons for the system
        let pi_electrons = system
            .iter()
            .map(|atom_idx| self.count_atom_pi_electrons_within_system(atom_idx, &system))
            .sum();

        // Check if the total system follows Hückel's rule
        is_hueckel_satisfied(pi_electrons)
    }

    fn are_all_atoms_conjugated(&self, system: &BitSet) -> bool {
        for atom_idx in system.iter() {
            if !self.can_be_conjugated_within_system(atom_idx, system) {
                return false;
            }
        }
        true
    }

    fn aromatize_system(&mut self, system: &BitSet) {
        // Get all bonds in the system
        let mut system_bonds = Vec::new();
        for atom_idx in system.iter() {
            for &bond in &self.atom_bonds()[atom_idx] {
                if system.contains(bond.target()) {

                    system_bonds.push((atom_idx, bond.target()));
                }
            }
        }

        // Convert all bonds in the system to aromatic
        for (atom1, atom2) in system_bonds {
            if atom1 < atom2 {
                self.change_bond_order(atom1, atom2, BondOrder::Aromatic);
            }
        }

    }

    fn count_atom_pi_electrons(&self, atom_idx: usize) -> usize {
        // Count pi electrons from multiple bonds
        let pi_electrons = self.atom_bonds()[atom_idx]
            .iter()
            .map(|bond| match bond.bond_order() {
                BondOrder::Double => 1, // One pi electron per double bond
                BondOrder::Triple => 2, // Two pi electrons per triple bond
                _ => 0,
            })
            .sum();

        // Add electrons from lone pairs for certain atoms
        let atomic_number = self.atomic_numbers()[atom_idx];
        match atomic_number {
            7 => {
                // Nitrogen
                if self.get_atom_charge(atom_idx) == 0 && pi_electrons == 0 {
                    2 // Neutral N with no pi bonds contributes lone pair
                } else {
                    pi_electrons
                }
            }
            8 | 16 => {
                // Oxygen or Sulfur
                if self.get_atom_charge(atom_idx) == 0 && pi_electrons == 0 {
                    2 // Neutral O/S with no pi bonds contributes lone pair
                } else {
                    pi_electrons
                }
            }
            _ => pi_electrons,
        }
    }


    fn count_atom_pi_electrons_within_system(&self, atom_idx: usize, system: &BitSet) -> usize {
        // Count pi electrons from multiple bonds
        let pi_electrons = self.atom_bonds()[atom_idx]
            .iter()
            .map(|bond| 
                if system.contains(bond.target()) {
                    match bond.bond_order() {
                        BondOrder::Double => 1, // One pi electron per double bond
                        BondOrder::Triple => 2, // Two pi electrons per triple bond
                        _ => 0,
                    }
                }else {
                    0
                }
            )
            .sum();

        // Add electrons from lone pairs for certain atoms
        let atomic_number = self.atomic_numbers()[atom_idx];
        match atomic_number {
            7 => {
                // Nitrogen
                if self.get_atom_charge(atom_idx) == 0 && pi_electrons == 0 {
                    2 // Neutral N with no pi bonds contributes lone pair
                } else {
                    pi_electrons
                }
            }
            8 | 16 => {
                // Oxygen or Sulfur
                if self.get_atom_charge(atom_idx) == 0 && pi_electrons == 0 {
                    2 // Neutral O/S with no pi bonds contributes lone pair
                } else {
                    pi_electrons
                }
            }
            _ => pi_electrons,
        }
    }

    fn check_for_metal_center(&self, system: &BitSet) -> bool {
        // Check if any atom is a transition metal
        system.iter().any(|atom_idx| {
            let atomic_number = self.atomic_numbers()[atom_idx];
            // Transition metals: 21-30, 39-48, 57-80, 89-112
            matches!(atomic_number,
                21..=30 | 39..=48 | 57..=80 | 89..=112
            )
        })
    }

    fn aromatize_metallocenic_system(&mut self, system: &BitSet) {
        // // First, identify the metal center
        // let metal_center = system.iter()
        //     .flat_map(|cycle| cycle.iter())
        //     .find(|&&atom_idx| {
        //         let atomic_number = self.atomic_numbers()[atom_idx];
        //         matches!(atomic_number, 21..=30 | 39..=48 | 57..=80 | 89..=112)
        //     });

        // if let Some(&metal_idx) = metal_center {
        //     // Convert all bonds in cyclopentadienyl rings to aromatic
        //     for cycle in system {
        //         if cycle.len() == 5 { // Cp rings are typically 5-membered
        //             for i in 0..cycle.len() {
        //                 let atom1 = cycle[i];
        //                 let atom2 = cycle[(i + 1) % cycle.len()];
        //                 self.change_bond_order(atom1, atom2, BondOrder::Aromatic);
                        
        //                 // Add coordinate bond from metal to ring atoms
        //                 if atom1 != metal_idx {
        //                     self.change_bond_order(metal_idx, atom1, BondOrder::Coordinate);
        //                 }
        //             }
        //         }
        //     }
        // }
    }

    fn can_be_conjugated(&self, atom_idx: usize) -> bool {
        let atomic_number = self.atomic_numbers()[atom_idx];
        let number_of_bonds = self.atom_bonds()[atom_idx].len();
        let number_of_double_or_triple_bonds = self.atom_bonds()[atom_idx].iter().filter(|bond| bond.bond_order() == BondOrder::Double || bond.bond_order() == BondOrder::Triple).count();

        if number_of_double_or_triple_bonds > 1 {
            return false;
        }
        
        match atomic_number {
            6 => {

                let number_of_aromatic_bonds = self.atom_bonds()[atom_idx].iter().filter(|bond| bond.bond_order() == BondOrder::Aromatic).count();
                let charge = self.get_atom_charge(atom_idx);
                number_of_bonds <= 3 && (number_of_double_or_triple_bonds == 1 || number_of_aromatic_bonds > 1 || charge != 0) }, // Carbon
            7 => number_of_bonds <= 4, // Nitrogen
            8 => number_of_bonds <= 2, // Oxygen
            15 => number_of_bonds <= 3, // Phosphorus
            16 => number_of_bonds <= 2, // Sulfur
            21..=30 | 39..=48 | 57..=80 | 89..=112 => true, // Transition metals
            _ => false,
        }
    }

    fn can_be_conjugated_within_system(&self, atom_idx: usize, system: &BitSet) -> bool {
        if self.get_atom_bonds(atom_idx).unwrap_or_default().iter().any(|bond| bond.bond_order() == BondOrder::Aromatic) {
            return true;
        }

        let atomic_number = self.atomic_numbers()[atom_idx];
        let number_of_bonds = self.atom_bonds()[atom_idx].len();
        let number_of_double_or_triple_bonds = self.atom_bonds()[atom_idx].iter().filter(|bond| bond.bond_order() == BondOrder::Double || bond.bond_order() == BondOrder::Triple).count();
        if number_of_double_or_triple_bonds > 1 {
            return false;
        }
        
        match atomic_number {
            6 => {
                let number_of_double_or_triple_bonds_within_system = self.atom_bonds()[atom_idx].iter().filter(|bond| system.contains(bond.target()) && (bond.bond_order() == BondOrder::Double || bond.bond_order() == BondOrder::Triple) ).count();
                let has_exocyclic_double_bond = self.atom_bonds()[atom_idx].iter().any(|bond| !system.contains(bond.target()) && bond.bond_order() == BondOrder::Double);
                //let is_bonded_to_exocyclic_more_electronegative_atom_via_double_bond = self.atom_bonds()[atom_idx].iter().filter(|bond| !system.contains(bond.target())).any(|bond| (self.cmp_electronegativities(atom_idx, bond.target()) == std::cmp::Ordering::Less) && bond.bond_order() == BondOrder::Double);
                number_of_bonds <= 3 && (
                 number_of_double_or_triple_bonds == 1 &&
                 number_of_double_or_triple_bonds_within_system == 1
                 ||
                 has_exocyclic_double_bond
                 // is_bonded_to_exocyclic_more_electronegative_atom_via_double_bond
                )
                 
                 }, // Carbon
            7 => number_of_bonds <= 4, // Nitrogen
            8 => number_of_bonds <= 2, // Oxygen
            15 => number_of_bonds <= 3, // Phosphorus
            16 => number_of_bonds <= 2, // Sulfur
            21..=30 | 39..=48 | 57..=80 | 89..=112 => true, // Transition metals
            _ => false,
        }
    }

    fn number_of_bond_order(&self, atom_idx: usize, bond_order: BondOrder) -> usize {
        self.atom_bonds()[atom_idx].iter().filter(|bond| bond.bond_order() == bond_order).count()
    }

    fn is_potentially_aromatic_cycle(&self, cycle: &[usize]) -> bool {
        cycle.len() > 4
    }


    fn count_cycle_pi_electrons(&self, cycle: &[usize]) -> usize {
        let mut total = 0;
        
        // Count electrons from atoms
        for &atom_idx in cycle {
            total += self.count_atom_pi_electrons_extended_within_system(atom_idx, &cycle_to_bitset(cycle));
        }

        total
    }


    /// Checks if two cycles share any atoms
    /// 
    /// # Arguments
    /// * `cycle1_bits` - BitSet representation of first cycle
    /// * `cycle2_bits` - BitSet representation of second cycle
    /// 
    /// # Returns
    /// true if the cycles share any atoms
    fn cycles_share_atoms(cycle1_bits: &BitSet, cycle2_bits: &BitSet) -> bool {
        cycle1_bits.intersection(cycle2_bits).next().is_some()
    }

    /// Groups cycles into fused systems using BitSet operations
    fn group_fused_cycles(&self, cycles: &[BitSet]) -> Vec<BitSet> {
        let mut fused_systems = Vec::new();
        let mut used = vec![false; cycles.len()];

        for (i, _) in cycles.iter().enumerate() {
            if used[i] {
                continue;
            }

            let mut current_system = cycles[i].clone();
            used[i] = true;

            let mut changed = true;
            while changed {
                changed = false;
                for (j, _) in cycles.iter().enumerate() {
                    if used[j] {
                        continue;
                    }

                    // Check if cycle j shares atoms with current system
                    if Self::cycles_share_atoms(&current_system, &cycles[j]) {
                        // Union of bits
                        current_system.union_with(&cycles[j]);
                        used[j] = true;
                        changed = true;
                    }
                }
            }

            if self.is_rdkit_aromatic_system(&current_system) && !fused_systems.contains(&current_system) && !cycles.contains(&current_system) {
                fused_systems.push(current_system);
            }
        }

        fused_systems.extend_from_slice(&cycles);

        fused_systems
    }

    /// Gets the union of atoms in a set of cycles
    /// 
    /// # Arguments
    /// * `cycles` - Vector of cycles to get the union of
    /// 
    /// # Returns
    /// A BitSet representing all atoms present in any of the cycles
    fn get_system_atoms(&self, cycles: &[Vec<usize>]) -> BitSet {
        let mut system_atoms = BitSet::new();
        for cycle in cycles {
            for &atom_idx in cycle {
                system_atoms.insert(atom_idx);
            }
        }
        system_atoms
    }

    /// Performs a simple aromatization following RDKit's model
    ///
    /// # Examples
    /// ```
    /// use molecules::prelude::*;
    /// let mut mol = Molecule2D::from_smiles("c1ccccc1").unwrap();
    /// let mut mol = &mut mol[0];
    /// mol.simple_aromatization();
    /// assert!(mol.get_aromatic_atoms().len() == 6);
    /// ```
    fn simple_aromatization(&mut self) {
        let mut cycles = Vec::new();
        let mut atoms_in_cycles = IntSet::default();
        
        // Find all potential cycles
        self.find_all_conjugated_cycles(&mut cycles, &mut atoms_in_cycles);
        self.find_and_add_subcycles(&mut cycles, &mut atoms_in_cycles);
        
        // Group cycles into fused systems
        let mut fused_systems = self.group_fused_cycles(&cycles);
        fused_systems.sort_by(|a, b| b.len().cmp(&a.len()));
        let mut aromatic_systems: Vec<BitSet> = Vec::new();


        for system in fused_systems {
            if aromatic_systems.iter().any(|s| s.is_superset(&system)) {
                continue;
            }
            self.aromatize_system(&system);
            aromatic_systems.push(system);
        }

        self.add_aromatic_bonds();
    }

    fn is_rdkit_aromatic_system(&self, system: &BitSet) -> bool {
        for atom_idx in system.iter() {
            if !self.can_be_conjugated_within_system(atom_idx, system) {
                return false;
            }
        }


        let pi_electrons = self.count_system_pi_electrons(system);
        is_hueckel_satisfied(pi_electrons)
    }

    fn is_electronegative(&self, atom_idx: usize) -> bool {
        let atomic_number = self.atomic_numbers()[atom_idx];
        matches!(atomic_number, 7 | 8 | 9 | 17 | 35 | 53)  // N, O, F, Cl, Br, I
    }

    fn is_amine_nitrogen(&self, atom_idx: usize) -> bool {
        let double_bonds = self.has_double_bonds(atom_idx);
        let hydrogen_count = self.number_of_bonded_element(atom_idx, 1);
        // Amine nitrogen has at least one hydrogen
        !double_bonds || hydrogen_count > 0
    }

    fn has_double_bonds(&self, atom_idx: usize) -> bool {
        self.atom_bonds()[atom_idx].iter().any(|bond| bond.bond_order() == BondOrder::Double)
    }

    /// Finds all subcycles within a given cycle and adds them to the cycles list
    /// 
    /// # Arguments
    /// * `cycles` - Vector of BitSets representing cycles to check for subcycles
    /// * `atoms_in_cycles` - Set of atoms that are part of any cycle
    /// 
    /// # Returns
    /// A new vector containing all original cycles plus any valid subcycles found
    fn find_and_add_subcycles(&self, cycles: &mut Vec<BitSet>, atoms_in_cycles: &mut IntSet<usize>) {
        let mut subcycles = Vec::new();

        // For each cycle, try to find subcycles by removing one atom at a time
        for cycle in cycles.clone() {
            let cycle_vec: Vec<_> = cycle.iter().collect();
            if cycle_vec.len() <= 5 {
                continue; // Skip cycles that are too small to contain subcycles
            }

            // Try removing each atom to find potential subcycles
            for skip_atom in cycle_vec.iter() {
                let mut visited = vec![false; self.atomic_numbers().len()];
                let mut path = Vec::new();
                let mut path_set = HashSet::new();

                // Start DFS from each remaining atom in the cycle
                for &start in cycle_vec.iter().filter(|&&atom| atom != *skip_atom) {
                    if visited[start] {
                        continue;
                    }

                    self.dfs_find_subcycles(
                        start,
                        start,
                        *skip_atom,
                        &cycle,
                        &mut visited,
                        &mut path,
                        &mut path_set,
                        &mut subcycles,
                        atoms_in_cycles,
                    );
                }
            }
        }

        // Add valid subcycles to the result
        for subcycle in subcycles {
            if !cycles.contains(&subcycle) && self.is_potentially_aromatic_cycle(&subcycle.iter().collect::<Vec<_>>()) {
                cycles.push(subcycle);
            }
        }
    }

    /// Helper function for DFS subcycle finding
    fn dfs_find_subcycles(
        &self,
        current: usize,
        parent: usize,
        skip_atom: usize,
        original_cycle: &BitSet,
        visited: &mut Vec<bool>,
        path: &mut Vec<usize>,
        path_set: &mut HashSet<usize>,
        subcycles: &mut Vec<BitSet>,
        atoms_in_cycles: &mut IntSet<usize>,
    ) {
        visited[current] = true;
        path.push(current);
        path_set.insert(current);

        for &bond in &self.atom_bonds()[current] {
            let target = bond.target();
            
            // Skip the atom we're excluding and atoms not in the original cycle
            if target == skip_atom || !original_cycle.contains(target) {
                continue;
            }

            // Found a cycle
            if path_set.contains(&target) {
                if target == path[0] && path.len() >= 5 {
                    let subcycle = cycle_to_bitset(path);
                    if !subcycles.contains(&subcycle) {
                        subcycles.push(subcycle);
                        for atom in path.clone() {
                            atoms_in_cycles.insert(atom);
                        }
                    }
                }
            } else if !visited[target] {
                self.dfs_find_subcycles(
                    target,
                    current,
                    skip_atom,
                    original_cycle,
                    visited,
                    path,
                    path_set,
                    subcycles,
                    atoms_in_cycles,
                );
            }
        }

        path.pop();
        path_set.remove(&current);
        visited[current] = false;
    }
    
    /// Searches for aromatic cycles within the molecule.
    ///
    /// An aromatic cycle is defined as a ring where all participating bonds are aromatic.
    ///
    /// # Returns
    ///
    /// A vector of aromatic cycles, each represented as a vector of atom indices.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::prelude::*;
    /// 
    /// // Benzene molecule with a single aromatic cycle
    /// let molecules = Molecule2D::from_smiles("c1ccccc1").unwrap();
    /// let molecule = molecules.first().unwrap();
    /// let cycles = molecule.find_aromatic_cycles();
    /// assert_eq!(cycles.len(), 1);
    /// assert_eq!(cycles[0], vec![0, 1, 2, 3, 4, 5]);
    /// 
    /// // Naphthalene with two aromatic cycles
    /// let molecules = Molecule2D::from_smiles("c1ccc2ccccc2c1").unwrap();
    /// let molecule = molecules.first().unwrap();
    /// let cycles = molecule.find_aromatic_cycles();
    /// assert_eq!(cycles.len(), 2);
    /// assert!(cycles.contains(&vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]));
    /// 
    /// // Cyclohexane with no aromatic cycles
    /// let molecules = Molecule2D::from_smiles("C1CCCCCC1").unwrap();
    /// let molecule = molecules.first().unwrap();
    /// let cycles = molecule.find_aromatic_cycles();
    /// assert_eq!(cycles.len(), 0);
    /// ```
    fn find_aromatic_cycles(&self) -> Vec<Vec<usize>> {
        self.find_rings().into_iter().filter(|ring| {
            // Check if all consecutive bonds in the ring are aromatic
            let all_aromatic = ring.windows(2).all(|w| {
                let atom1 = w[0];
                let atom2 = w[1];
                self.atom_bonds()[atom1].iter().any(|bond| 
                    bond.target() == atom2 && bond.bond_order() == BondOrder::Aromatic
                )
            }) && {
                // Ensure the bond closing the ring is aromatic
                let last_atom = ring.last().unwrap();
                let first_atom = ring[0];
                self.atom_bonds()[*last_atom].iter().any(|bond| 
                    bond.target() == first_atom && bond.bond_order() == BondOrder::Aromatic
                )
            };
            all_aromatic
        }).collect()
    }

    fn add_aromatic_bonds(&mut self) {
        for atom_idx in 0..self.len() {
            let is_atom_aromatic = self.is_atom_aromatic(atom_idx);
            if !is_atom_aromatic {
                continue;
            }
            let bonds = self.get_atom_bonds(atom_idx).unwrap_or_default().to_vec();
            for neighbor_idx in bonds {
                let is_neighbor_aromatic = self.is_atom_aromatic(neighbor_idx.target());
                if !is_neighbor_aromatic || neighbor_idx.bond_order() == BondOrder::Aromatic {
                    continue;
                }
                self.change_bond_order(atom_idx, neighbor_idx.target(), BondOrder::Aromatic);
            }
        }
    }

    /// Fixes any valence inconsistencies in the molecule.
    ///
    /// This function ensures that each atom adheres to its expected valency by adjusting
    /// charges and bond orders as necessary.
    ///
    /// # Errors
    ///
    /// Returns an error if valency cannot be satisfied for any atom.
    fn fix_valence_inconsistencies(&mut self) -> Result<(), String> {
        // Vector to hold all bond changes to be made
        let mut bond_changes: Vec<(usize, usize, BondOrder)> = Vec::new();
        // Vector to hold charge adjustments (atom_index, charge_change)
        let mut charge_changes: Vec<(usize, i8)> = Vec::new();

        // Iterate over all atoms to identify valence issues
        for atom_index in 0..self.len() {
            let expected_valency = match self.expected_valency(atom_index) {
                Some(val) => val,
                None => continue, // Skip atoms with undefined valency
            };
            let actual_valency = self.actual_valency(atom_index);
            let charge = self.get_atom_charge(atom_index);

            if actual_valency > expected_valency {
                // Valency exceeded: need to reduce valency by converting single bonds to double bonds
                let excess = actual_valency - expected_valency;
                let mut bonds_to_modify = 0;

                for bond in self.atom_bonds()[atom_index].iter() {
                    if bond.bond_order == BondOrder::Single && bonds_to_modify < excess {
                        bond_changes.push((bond.target(), atom_index, BondOrder::Double));
                        bonds_to_modify += 1;
                    }
                    if bonds_to_modify >= excess {
                        break;
                    }
                }

                if bonds_to_modify < excess {
                    return Err(format!(
                        "Cannot fix valency for atom {}: expected {}, got {}",
                        atom_index, expected_valency, actual_valency
                    ));
                }
            } else if actual_valency < expected_valency {
                // Valency deficit: need to increase valency by converting single bonds to double bonds
                let deficit = expected_valency - actual_valency;
                let mut bonds_to_modify = 0;

                for bond in self.atom_bonds()[atom_index].iter() {
                    if bond.bond_order == BondOrder::Single && bonds_to_modify < deficit {
                        bond_changes.push((bond.target(), atom_index, BondOrder::Double));
                        bonds_to_modify += 1;
                    }
                    if bonds_to_modify >= deficit {
                        break;
                    }
                }

                if bonds_to_modify < deficit {
                    // If unable to fix valency by bond modifications, adjust charge
                    charge_changes.push((atom_index, -(deficit as i8)));
                }
            }
        }

        // Apply all bond changes
        for (atom1, atom2, new_order) in bond_changes {
            self.change_bond_order(atom1, atom2, new_order);
        }

        // Apply all charge adjustments
        for (atom_index, charge_change) in charge_changes {
            self.set_atom_charge(atom_index, self.get_atom_charge(atom_index) + charge_change);
        }

        Ok(())
    }

    /// Sanitizes the molecule by removing invalid aromatic rings and correcting inconsistencies.
    ///
    /// This function performs the following steps:
    /// 1. Identifies all aromatic rings in the molecule.
    /// 2. Validates each aromatic ring against Hückel's rule and other aromaticity criteria.
    /// 3. Removes aromaticity from rings that do not satisfy aromaticity requirements.
    /// 4. Fixes any valence inconsistencies resulting from the removal.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::prelude::*;
    ///
    /// // Example with a valid aromatic ring
    /// let mut benzene = Molecule2D::from_smiles("c1ccccc1").unwrap().first().unwrap().clone();
    /// benzene.sanitize_mol().unwrap();
    /// assert_eq!(benzene.find_aromatic_cycles().len(), 1);
    ///
    /// // Example with an invalid aromatic ring
    /// let mut invalid_ring = Molecule2D::from_smiles("c1ccccc1C").unwrap().first().unwrap().clone();
    /// // Manually corrupt the aromaticity
    /// invalid_ring.change_bond_order(0, 1, BondOrder::Single);
    /// invalid_ring.sanitize_mol().unwrap();
    /// assert_eq!(invalid_ring.find_aromatic_cycles().len(), 0);
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an error if the molecule cannot be sanitized due to internal inconsistencies.
    /// 
    
    fn sanitize_mol(&mut self) -> Result<(), String> {
        // Step 1: Find all aromatic cycles
        let aromatic_cycles = self.find_aromatic_cycles();

        // Step 2: Validate each aromatic ring
        let mut invalid_cycles = Vec::new();
        for cycle in aromatic_cycles {
            if !self.is_valid_aromatic_cycle(&cycle) {
                invalid_cycles.push(cycle);
            }
        }

        // Step 3: Remove aromaticity from invalid cycles
        for cycle in invalid_cycles {
            for window in cycle.windows(2) {
                let atom1 = window[0];
                let atom2 = window[1];
                self.change_bond_order(atom1, atom2, BondOrder::Single);
            }
            // Close the ring by setting the bond between the last and first atom to single
            let last_atom = cycle.last().unwrap();
            let first_atom = cycle[0];
            self.change_bond_order(*last_atom, first_atom, BondOrder::Single);
        }

        // Step 4: Fix valence inconsistencies
        //self.fix_valence_inconsistencies()?;

        Ok(())
    }


}


fn relabel_numbers(input: &str) -> String {
    let mut result = String::with_capacity(input.len());
    let mut number_map = HashMap::new();
    let mut next_number = 1;
    let mut chars = input.chars().peekable();

    while let Some(c) = chars.next() {
        if c == '%' {
            // Handle multi-digit numbers
            let mut number = String::new();
            while let Some(&next_char) = chars.peek() {
                if next_char.is_ascii_digit() {
                    number.push(chars.next().unwrap());
                } else {
                    break;
                }
            }
            
            // Map the number if we haven't seen it before
            let mapped = number_map
                .entry(number.clone())
                .or_insert_with(|| {
                    let current = next_number;
                    next_number += 1;
                    current.to_string()
                });
            
            result.push('%');
            result.push_str(mapped);
        } else if c.is_ascii_digit() {
            // Handle single digit numbers
            let number = c.to_string();
            let mapped = number_map
                .entry(number)
                .or_insert_with(|| {
                    let current = next_number;
                    next_number += 1;
                    current.to_string()
                });
            
            result.push_str(mapped);
        } else {
            result.push(c);
        }
    }
    
    result
}

// Converts a cycle (represented as a Vec<usize>) to a BitSet for efficient comparison
/// 
/// # Arguments
/// * `cycle` - Vector of atom indices representing a cycle
/// 
/// # Returns
/// A BitSet where each set bit represents the presence of an atom in the cycle
/// 
/// # Example
/// ```
/// use molecules::molecule::base::cycle_to_bitset;
/// 
/// let cycle = vec![0, 2, 4];
/// let bitset = cycle_to_bitset(&cycle);
/// assert_eq!(bitset.len(), 3);
/// assert!(bitset.contains(0));
/// assert!(bitset.contains(2));
/// assert!(bitset.contains(4));
/// ```
pub fn cycle_to_bitset(cycle: &[usize]) -> BitSet {
    let mut bitset = BitSet::new();
    for &atom_idx in cycle {
        bitset.insert(atom_idx);
    }
    bitset
}

/// Converts a list of cycles to a vector of BitSets for efficient comparison
/// 
/// # Returns
/// A vector of BitSets, each representing a cycle in the molecule
pub fn cycles_to_bitsets(cycles: &[Vec<usize>]) -> Vec<BitSet> {
    cycles.iter()
        .map(|cycle| cycle_to_bitset(cycle))
        .collect()
}
