use crate::atom::Atom;
use crate::consts::{BOND_SEARCH_THRESHOLD, BOND_TOLERANCE};
use crate::vector::Vector;
use chemistry_consts::ElementProperties;
use core::fmt::{Display, Formatter};
use itertools::Itertools;
use kiddo::{KdTree, SquaredEuclidean};
use nohash_hasher::{IntMap, IntSet};
use petgraph::algo::subgraph_isomorphisms_iter;
use petgraph::prelude::*;
use rayon::prelude::*;

use std::{
    collections::{HashSet, VecDeque},
    default::Default,
    fs::File,
    io::{self, BufRead},
    num::ParseFloatError,
    path::Path,
};
#[derive(Debug, Default, Clone, Copy, Eq, PartialEq)]
pub enum ChiralClass {
    S,
    R,
    AL(u8),
    SP(u8),
    TB(u8),
    OH(u8),
    #[default]
    None,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub struct BondTarget {
    target: usize,
    bond_type: BondType,
}

impl BondTarget {
    pub fn new(target: usize, bond_type: BondType) -> BondTarget {
        BondTarget { target, bond_type }
    }
    fn single(target: usize) -> BondTarget {
        BondTarget {
            target,
            bond_type: BondType::Single,
        }
    }
    fn double(target: usize) -> BondTarget {
        BondTarget {
            target,
            bond_type: BondType::Double,
        }
    }
    fn triple(target: usize) -> BondTarget {
        BondTarget {
            target,
            bond_type: BondType::Triple,
        }
    }
    fn aromatic(target: usize) -> BondTarget {
        BondTarget {
            target,
            bond_type: BondType::Aromatic,
        }
    }
    pub fn bond_type(&self) -> BondType {
        self.bond_type
    }
    pub fn target(&self) -> usize {
        self.target
    }
    pub fn order(&self) -> u8 {
        match self.bond_type {
            BondType::Single => 1,
            BondType::Double => 2,
            BondType::Triple => 3,
            BondType::Quadruple => 4,
            // TODO: Implement aromatic bond order
            BondType::Aromatic => 1,
        }
    }
}

#[derive(Debug, Default, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub enum BondType {
    #[default]
    Single,
    Double,
    Triple,
    Aromatic,
    Quadruple,
}

impl Display for BondType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            BondType::Single => write!(f, "-"),
            BondType::Double => write!(f, "="),
            BondType::Triple => write!(f, "#"),
            BondType::Aromatic => write!(f, ":"),
            BondType::Quadruple => write!(f, "$"),
        }
    }
}

/// This function   
#[derive(Debug, Default, Clone)]
pub struct Molecule2D {
    pub atomic_numbers: Vec<u8>,
    pub atom_classes: Option<Vec<u8>>,
    pub charge: i32,
    pub charges: Vec<i8>,
    pub chiral_classes: Option<Vec<ChiralClass>>,
    pub isotopes: Option<Vec<u16>>,
    radical_states: Vec<bool>,
    atom_bonds: Vec<Vec<BondTarget>>,
}

#[derive(Debug, Default, Clone)]
pub struct Molecule3D {
    pub atomic_numbers: Vec<u8>,
    pub atom_classes: Option<Vec<u8>>,
    pub charge: i32,
    pub charges: Vec<i8>,
    pub chiral_classes: Option<Vec<ChiralClass>>,
    pub isotopes: Option<Vec<u16>>,
    pub positions: Vec<Vector>,
    radical_states: Vec<bool>,
    atom_bonds: Vec<Vec<BondTarget>>,
}

pub trait Molecule {
    fn atom_bonds(&self) -> &Vec<Vec<BondTarget>>;
    fn atom_bonds_mut(&mut self) -> &mut Vec<Vec<BondTarget>>;
    fn atomic_numbers(&self) -> &[u8];
    fn atomic_numbers_mut(&mut self) -> &mut Vec<u8>;
    fn chiral_classes(&self) -> Option<&[ChiralClass]>;
    fn chiral_classes_mut(&mut self) -> Option<&mut Vec<ChiralClass>>;
    fn charges(&self) -> &[i8];
    fn charges_mut(&mut self) -> &mut Vec<i8>;
    fn isotopes(&self) -> Option<&Vec<u16>>;
    fn atom_classes(&self) -> &Option<Vec<u8>>;
    fn atom_classes_mut(&mut self) -> &mut Option<Vec<u8>>;
    fn add_hydrogens(&mut self);
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

    fn valency_delta(&self, atom_index: usize) -> Option<i8> {
        let expected_valency = self.expected_valency(atom_index)?;
        let actual_valency = self.actual_valency(atom_index);
        Some(actual_valency - expected_valency)
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
    fn actual_valency(&self, atom_index: usize) -> i8 {
        self.atom_bonds()[atom_index]
            .iter()
            .map(|bond| match bond.bond_type() {
                BondType::Single => 2,
                BondType::Double => 4,
                BondType::Triple => 6,
                BondType::Quadruple => 8,
                BondType::Aromatic => 3,
            })
            .sum::<i8>()
            / 2
            + if self.is_atom_radical(atom_index) {
                1
            } else {
                0
            }
            + self.get_atom_charge(atom_index).abs()
    }
    fn degree(&self, atom_index: usize) -> Option<i8> {
        let expected_valency = self.expected_valency(atom_index)?;
        let actual_valency = self.actual_valency(atom_index);
        // May need to be changed for elements with unknown valencies
        Some(actual_valency - expected_valency)
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
    fn get_charge(&self) -> i32;
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
    fn get_edges_with_type(&self) -> Vec<(usize, usize, BondType)> {
        self.atom_bonds()
            .iter()
            .enumerate()
            .flat_map(|(atom_index, bonds)| {
                bonds.iter().map(move |bond| {
                    if atom_index < bond.target() {
                        Some((atom_index, bond.target(), bond.bond_type()))
                    } else {
                        None
                    }
                })
            })
            .flatten()
            .collect()
    }

    fn to_ungraph(&self) -> UnGraph<u8, ()> {
        let mut graph = UnGraph::<u8, ()>::default();
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
                    graph.add_edge(source, target, ());
                }
            }
        }
        graph
    }

    fn to_ungraph_from_slice(&self, slice: &[usize]) -> UnGraph<u8, ()> {
        let mut graph = UnGraph::<u8, ()>::default();
        let mut node_indices = IntMap::<usize, NodeIndex<u32>>::default();

        for &index in slice {
            let atomic_number = self.atomic_numbers()[index];
            let node_index = graph.add_node(atomic_number);
            node_indices.insert(index, node_index);
        }

        for &index in slice {
            let bonds = &self.atom_bonds()[index];
            for bond in bonds {
                let source = node_indices[&index];
                let target = node_indices[&bond.target()];
                if source < target && slice.contains(&bond.target()) {
                    graph.add_edge(source, target, ());
                }
            }
        }
        graph
    }
    fn match_submolecule(&self, other: &Self) -> Option<Vec<IntMap<usize, usize>>> {
        let mut self_components = self.get_components();
        let mut other_components = other.get_components();

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
        self.atom_bonds().get(atom_index).map(|bonds| bonds.as_slice())
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
                println!("No bonds found for atom {}, this means an out of bounds access", index);
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
        bond_type: BondType,
        parent_index: Option<usize>,
        visited: &mut Vec<bool>,
        ring_closures: &mut IntMap<usize, Vec<(usize, BondType)>>,
        ring_counter: &mut usize,
    ) -> Node {
        let mut root = Node::new(start_index, bond_type);
        visited[start_index] = true;
        let bonds = &self.atom_bonds()[start_index];
        for &neighbor in bonds {
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
                let bond_type = neighbor.bond_type();
                ring_closures
                    .entry(start_index)
                    .or_default()
                    .push((*ring_counter, bond_type));
                ring_closures
                    .entry(neighbor.target)
                    .or_default()
                    .push((*ring_counter, bond_type));
                *ring_counter += 1;
            } else {
                let child = self.build_smiles_tree(
                    neighbor.target,
                    neighbor.bond_type(),
                    Some(start_index),
                    visited,
                    ring_closures,
                    ring_counter,
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
            if bond.bond_type() == BondType::Double || bond.bond_type() == BondType::Aromatic {
                number_of_pi_electrons += 2;
            } else if bond.bond_type() == BondType::Triple {
                number_of_pi_electrons += 4;
            }
        }
        number_of_pi_electrons
    }

    fn molecular_formula(&self) -> MolecularFormula {
        MolecularFormula::from_molecule(self)
    }

    fn to_smiles(&self) -> String {
        let mut smiles = String::new();
        let components = self.get_components();

        for component in components {
            let ring_closures = &mut IntMap::default();
            let mut ring_counter = 1;
            let mut visited = vec![false; self.atomic_numbers().len()];
            let root = self.build_smiles_tree(
                component[0],
                BondType::Single,
                None,
                &mut visited,
                ring_closures,
                &mut ring_counter,
            );
            if !smiles.is_empty() {
                smiles.push('.');
            }
            smiles.push_str(&self.construct_smiles(&root, ring_closures))
        }
        smiles
    }

    fn construct_smiles(
        &self,
        node: &Node,
        ring_closures: &IntMap<usize, Vec<(usize, BondType)>>,
    ) -> String {
        let mut smiles = String::new();
        match node.bond_type {
            BondType::Single => {}
            BondType::Double => smiles.push('='),
            BondType::Triple => smiles.push('#'),
            BondType::Quadruple => smiles.push('$'),
            BondType::Aromatic => smiles.push('~'),
        }
        let charge = self.charges()[node.index];
        let number_of_hydrogens = self.number_of_bonded_element(node.index, 1);

        let hydrogen_str = match number_of_hydrogens {
            0 => String::new(),
            1 => "H".to_string(),
            _ => format!("H{}", number_of_hydrogens),
        };

        // TODO: Introduce error handling here
        let Some(atomic_symbol) = self.atomic_numbers()[node.index].atomic_symbol() else {
            return "".to_string();
        };

        let charge_string = match charge {
            2.. => format!("[{}{}+{}]", atomic_symbol, hydrogen_str, charge.abs()),
            ..=-2 => format!("[{}{}-{}]", atomic_symbol, hydrogen_str, charge.abs()),
            1 => format!("[{}{}+1]", atomic_symbol, hydrogen_str),
            -1 => format!("[{}{}-1]", atomic_symbol, hydrogen_str),
            0 => {
                if number_of_hydrogens > 0 {
                    format!("[{}{}]", atomic_symbol, hydrogen_str)
                } else {
                    atomic_symbol.to_string()
                }
            }
        };
        smiles.push_str(&charge_string);

        if let Some(closures) = ring_closures.get(&node.index) {
            for closure in closures {
                if closure.1 == BondType::Single {
                    smiles.push_str(&format!("{}", closure.0));
                } else {
                    smiles.push_str(&format!("{}{}", closure.1, closure.0));
                }
            }
        }

        if node.children.len() == 1 {
            smiles.push_str(&self.construct_smiles(&node.children[0], ring_closures));
        } else {
            for child in &node.children {
                smiles.push('(');
                smiles.push_str(&self.construct_smiles(child, ring_closures));
                smiles.push(')');
            }
        }
        smiles
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
                self.dfs_find_rings(start, start, &mut visited, &mut path, &mut path_set, &mut rings);
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
}

impl Molecule for Molecule3D {
    fn atom_bonds(&self) -> &Vec<Vec<BondTarget>> {
        &self.atom_bonds
    }
    fn atom_bonds_mut(&mut self) -> &mut Vec<Vec<BondTarget>> {
        &mut self.atom_bonds
    }
    fn atom_classes(&self) -> &Option<Vec<u8>> {
        &self.atom_classes
    }
    fn atom_classes_mut(&mut self) -> &mut Option<Vec<u8>> {
        &mut self.atom_classes
    }
    fn charges(&self) -> &[i8] {
        &self.charges
    }
    fn charges_mut(&mut self) -> &mut Vec<i8> {
        &mut self.charges
    }
    fn atomic_numbers(&self) -> &[u8] {
        &self.atomic_numbers
    }
    fn atomic_numbers_mut(&mut self) -> &mut Vec<u8> {
        &mut self.atomic_numbers
    }
    
    fn chiral_classes(&self) -> Option<&[ChiralClass]> {
        self.chiral_classes.as_deref()
    }

    fn chiral_classes_mut(&mut self) -> Option<&mut Vec<ChiralClass>> {
        self.chiral_classes.as_mut()
    }

    fn get_charge(&self) -> i32 {
        self.charge
    }

    fn isotopes(&self) -> Option<&Vec<u16>> {
        self.isotopes.as_ref()
    }

    fn radical_states(&self) -> &[bool] {
        &self.radical_states
    }

    fn radical_states_mut(&mut self) -> &mut Vec<bool> {
        &mut self.radical_states
    }

    fn is_atom_radical(&self, atom_index: usize) -> bool {
        self.radical_states[atom_index]
    }
    fn set_atom_radical(&mut self, atom_index: usize, is_radical: bool) {
        if atom_index < self.radical_states.len() {
            self.radical_states[atom_index] = is_radical;
        } else {
            println!("Atom index out of bounds, could not set radical state");
        }
    }

    fn add_hydrogens(&mut self) {
        let degrees = self.degrees();
        degrees.iter().enumerate().for_each(|(index, &degree)| {
            let number_of_hydrogens = -degree;
            for _ in 0..number_of_hydrogens {
                self.atomic_numbers.push(1);
                self.atom_bonds.push(Vec::from([BondTarget::single(index)]));
                self.atom_bonds[index].push(BondTarget::single(self.atomic_numbers.len() - 1));
            }
        });
    }
    fn from_atoms(atoms: Vec<Atom>) -> Self {
        let atomic_numbers = atoms.iter().map(|atom| atom.atomic_number).collect();
        let charges = atoms.iter().map(|atom| atom.charge).collect();
        let radical_states = atoms.iter().map(|atom| atom.is_radical).collect();
        let mut isotopes = None;
        if atoms.iter().any(|atom| atom.isotope.is_some()) {
            isotopes = Some(
                atoms
                    .iter()
                    .map(|atom| {
                        if atom.isotope.is_some() {
                            atom.isotope().unwrap()
                        } else {
                            atom.atomic_number()
                                .isotopes()
                                .unwrap()
                                .next()
                                .unwrap()
                                .mass
                                .round() as u16
                        }
                    })
                    .collect(),
            );
        }
        let atom_bonds = atoms.iter().map(|atom| atom.bonds().clone()).collect();
        let positions: Vec<Vector> = atoms
            .iter()
            .map(|atom| atom.position_vector.unwrap_or_default())
            .collect();

        let mut chirals =None;
        if atoms
            .iter()
            .any(|atom| atom.chiral_class != ChiralClass::None)
        {
            chirals = Some(atoms.into_iter().map(|atom| atom.chiral_class).collect());
        }
        let mut molecule = Molecule3D {
            atomic_numbers,
            charges,
            radical_states,
            atom_bonds,
            positions,
            isotopes,
            chiral_classes: chirals,
            ..Default::default()
        };
        molecule.identify_bonds(BOND_SEARCH_THRESHOLD);
        molecule
    }
}

impl Molecule for Molecule2D {
    fn atom_bonds(&self) -> &Vec<Vec<BondTarget>> {
        &self.atom_bonds
    }
    fn atom_bonds_mut(&mut self) -> &mut Vec<Vec<BondTarget>> {
        &mut self.atom_bonds
    }
    fn atom_classes(&self) -> &Option<Vec<u8>> {
        &self.atom_classes
    }
    fn atom_classes_mut(&mut self) -> &mut Option<Vec<u8>> {
        &mut self.atom_classes
    }
    fn charges(&self) -> &[i8] {
        &self.charges
    }
    fn charges_mut(&mut self) -> &mut Vec<i8> {
        &mut self.charges
    }
    fn atomic_numbers(&self) -> &[u8] {
        &self.atomic_numbers
    }
    fn atomic_numbers_mut(&mut self) -> &mut Vec<u8> {
        &mut self.atomic_numbers
    }
    fn chiral_classes(&self) -> Option<&[ChiralClass]> {
        self.chiral_classes.as_deref()
    }
    fn chiral_classes_mut(&mut self) -> Option<&mut Vec<ChiralClass>> {
        self.chiral_classes.as_mut()
    }
    fn get_charge(&self) -> i32 {
        self.charge
    }
    fn isotopes(&self) -> Option<&Vec<u16>> {
        self.isotopes.as_ref()
    }
    fn radical_states(&self) -> &[bool] {
        &self.radical_states
    }
    fn radical_states_mut(&mut self) -> &mut Vec<bool> {
        &mut self.radical_states
    }
    fn is_atom_radical(&self, atom_index: usize) -> bool {
        self.radical_states[atom_index]
    }
    fn set_atom_radical(&mut self, atom_index: usize, is_radical: bool) {
        self.radical_states[atom_index] = is_radical;
    }
    fn add_hydrogens(&mut self) {
        let degrees = self.degrees();
        degrees.iter().enumerate().for_each(|(index, &degree)| {
            let number_of_hydrogens = -degree;
            for _ in 0..number_of_hydrogens {
                let hydrogen_index = self.atomic_numbers.len();
                self.atomic_numbers_mut().push(1);
                self.atom_bonds_mut()
                    .push(Vec::from([BondTarget::single(index)]));
                self.atom_bonds_mut()[index].push(BondTarget::single(hydrogen_index));
            }
        });
    }

    fn from_atoms(atoms: Vec<Atom>) -> Self {
        let atomic_numbers = atoms.iter().map(|atom| atom.atomic_number).collect();
        let charges = atoms.iter().map(|atom| atom.charge).collect();
        let radical_states = atoms.iter().map(|atom| atom.is_radical).collect();

        let atom_bonds = atoms.iter().map(|atom| atom.bonds.clone()).collect();
        let mut isotopes = None;
        let mut atom_classes = None;
        let mut chiral_classes = None;

        if atoms.iter().any(|atom| atom.isotope.is_some()) {
            isotopes = Some(
                atoms
                    .iter()
                    .map(|atom| {
                        if atom.isotope.is_some() {
                            atom.isotope.unwrap()
                        } else {
                            atom.atomic_number
                                .isotopes()
                                .unwrap()
                                .next()
                                .unwrap()
                                .mass
                                .round() as u16
                        }
                    })
                    .collect(),
            );
        }
        if atoms.iter().any(|atom| atom.atom_class.is_some()) {
            atom_classes = Some(
                atoms
                    .iter()
                    .map(|atom| atom.atom_class.unwrap_or(0))
                    .collect(),
            );
        }
        if atoms
            .iter()
            .any(|atom| atom.chiral_class != ChiralClass::None)
        {
            chiral_classes = Some(atoms.into_iter().map(|atom| atom.chiral_class).collect());
        }

        Molecule2D {
            atomic_numbers,
            charges,
            radical_states,
            atom_bonds,
            isotopes,
            chiral_classes,
            atom_classes,
            ..Default::default()
        }
    }
}

pub struct MolecularSystem {
    pub molecules: Vec<Molecule3D>,
}

impl MolecularSystem {
    pub fn new(molecules: Vec<Molecule3D>) -> Self {
        MolecularSystem { molecules }
    }
    pub fn from_xyz<P: AsRef<Path>>(path: P) -> MolecularSystem {
        let molecule = Molecule3D::from_xyz(&path);
        MolecularSystem::new(vec![molecule])
    }
    pub fn from_pdb<P: AsRef<Path>>(path: P) -> MolecularSystem {
        let file = File::open(&path).expect("Could not open file");
        let reader = io::BufReader::new(&file);

        let mut atom_groups: Vec<Vec<String>> = vec![Vec::new()];

        for line in reader.lines().map_while(Result::ok) {
            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                atom_groups.last_mut().unwrap().push(line);
            } else if line.starts_with("TER") {
                atom_groups.push(Vec::new());
            } else if line.starts_with("ENDMDL") {
                break;
            }
        }
        let path = path.as_ref().to_str().unwrap();
        if path.ends_with("pdb") {
            MolecularSystem {
                molecules: atom_groups
                    .into_iter()
                    .map(|group| {
                        let atoms: Vec<Atom> = group
                            .par_iter()
                            .map(|line| extract_atom_pdb(line).expect("Could not parse Atom"))
                            .collect();
                        Molecule3D::from_atoms(atoms)
                    })
                    .collect(),
            }
        } else if path.ends_with("cif") {
            MolecularSystem {
                molecules: atom_groups
                    .into_iter()
                    .map(|group| {
                        let atoms: Vec<Atom> = group
                            .par_iter()
                            .map(|line| extract_atom_cif(line).expect("Could not parse Atom"))
                            .collect();
                        Molecule3D::from_atoms(atoms)
                    })
                    .collect(),
            }
        } else {
            MolecularSystem {
                molecules: vec![Molecule3D::default()],
            }
        }
    }

    pub fn center(&self) -> Vector {
        let mut sum = Vector::new(0.0, 0.0, 0.0);
        let mut count = 0.0;
        for molecule in &self.molecules {
            for &position in &molecule.positions {
                sum += position;
            }
            count += molecule.atomic_numbers.len() as f64
        }
        sum / count
    }

    pub fn find_bonds(&mut self, threshold: f64) {
        self.molecules
            .iter_mut()
            .for_each(|molecule| molecule.identify_bonds(threshold))
    }

    pub fn number_of_atoms(&self) -> usize {
        self.molecules
            .iter()
            .map(|molecule| molecule.atomic_numbers.len())
            .sum()
    }

    pub fn find_angles(&self) -> Vec<Vec<BondAngle>> {
        self.molecules
            .iter()
            .map(|molecule| molecule.find_angles())
            .collect()
    }

    pub fn add_hydrogens(&mut self) {
        self.molecules
            .iter_mut()
            .for_each(|molecule| molecule.add_hydrogens())
    }

    // pub fn find_dihedrals(&mut self) {
    //     self.molecules
    //         .iter_mut()
    //         .for_each(|molecule| molecule.find_dihedrals())
    // }
    // pub fn par_find_dihedrals(&mut self) {
    //     self.molecules
    //         .par_iter_mut()
    //         .for_each(|molecule| molecule.par_find_dihedrals())
    // }
}

impl Molecule3D {
    pub fn atom_bonds_mut(&mut self, atom_index: usize) -> &mut Vec<BondTarget> {
        &mut self.atom_bonds[atom_index]
    }

    pub fn class_of_atom(&self, atom_index: usize) -> u8 {
        if let Some(classes) = self.atom_classes() {
            *classes.get(atom_index).unwrap_or(&0)
        } else {
            0
        }
    }

    pub fn get_charge(&self, atom_index: usize) -> i8 {
        self.charges.get(atom_index).copied().unwrap_or(0)
    }

    pub fn get_isotope(&self, atom_index: usize) -> Option<u16> {
        let Some(isotopes) = &self.isotopes else {
            return None;
        };
        isotopes.get(atom_index).copied()
    }

    pub fn get_chiral_class(&self, atom_index: usize) -> ChiralClass {
        self.chiral_classes
            .as_ref()
            .and_then(|chirals| chirals.get(atom_index).copied())
            .unwrap_or(ChiralClass::None)
    }

    pub fn atom_charge_mut(&mut self, atom_index: usize) -> &mut i8 {
        &mut self.charges[atom_index]
    }

    pub fn is_atom_radical(&self, atom_index: usize) -> bool {
        self.radical_states[atom_index]
    }

    pub fn add_atom(&mut self, atom: Atom) {
        self.atomic_numbers.push(atom.atomic_number);
        self.positions.push(atom.position_vector.unwrap_or_default());
        self.charges.push(atom.charge);
        self.radical_states.push(atom.is_radical);
        self.atom_bonds.push(atom.bonds);
        self.isotopes
            .get_or_insert_with(Vec::new)
            .push(atom.isotope.unwrap_or_default());
        self.chiral_classes
            .get_or_insert_with(Vec::new)
            .push(atom.chiral_class);
    }

    // This function returns a reference to outgoing bonds of an atom
    pub fn atom_bonds(&self) -> &Vec<Vec<BondTarget>> {
        &self.atom_bonds
    }
    /// Creates a new molecule from an xyz file
    ///
    /// # Panics
    /// Function panics if the file cannot be opened or the file is empty
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::{Molecule3D,Molecule};
    ///
    /// let molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    /// assert_eq!(molecule.atomic_numbers.len(), 8);
    /// ```
    pub fn from_xyz<P: AsRef<Path>>(filepath: P) -> Molecule3D {
        let file = File::open(filepath).expect("Could not open file");
        let reader = io::BufReader::new(file);
        let atoms = reader
            .lines() // .skip(2) is possible here but could cause trouble for not much gain in malformatted files
            .filter_map(|line| {
                if let Ok(line) = line {
                    Atom::from_xyz_line(&line).ok()
                } else {
                    None
                }
            })
            .collect::<Vec<Atom>>();
        Molecule3D::from_atoms(atoms)
    }

    pub fn build_tree(&self) -> KdTree<f64, 3> {
        let mut tree: KdTree<f64, 3> = KdTree::with_capacity(self.atomic_numbers.len());
        for (index, position_vector) in self.positions.iter().enumerate() {
            tree.add(&position_vector.as_array(), index as u64);
        }
        tree
    }

    // PROTOTYPE
    pub fn extract_transformation_rule(&self, other: &Self) -> String {
        let _transformation_rule = String::new();
        let _atom_changes: Vec<(usize, usize)> = Vec::new();
        let charge_changes = self.charge_changes(other);
        let changed_nodes: Vec<usize> = charge_changes.iter().map(|changes| changes.0).collect();
        let _bond_changes = self.identify_bond_changes(other);

        let mut possible_paths = Vec::new();

        changed_nodes.iter().combinations(2).for_each(|window| {
            let atom1 = window[0];
            let atom2 = window[1];
            if let Some(path) = self.shortest_path(*atom1, *atom2) {
                possible_paths.push(path)
            }
        });

        possible_paths
            .iter()
            .filter(|path| changed_nodes.iter().all(|node| path.contains(node)))
            .for_each(|path| println!("{path:?}"));

        "".to_string()
    }

    /// Compares the electronegativities of two atoms
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
    /// assert_eq!(molecule.cmp_electronegativities(0, 1), std::cmp::Ordering::Equal);
    /// ```
    pub fn cmp_electronegativities(
        &self,
        atom1_index: usize,
        atom2_index: usize,
    ) -> std::cmp::Ordering {
        self.electronegativity(atom1_index)
            .cmp(&self.electronegativity(atom2_index))
    }

    pub fn electronegativity(&self, atom_index: usize) -> Option<u16> {
        self.atomic_numbers.get(atom_index)?.electronegativity()
    }

    pub fn identify_bond_changes(&self, other: &Self) -> Vec<BondChange> {
        let mut bond_changes: Vec<BondChange> = Vec::new();
        self.atom_bonds
            .iter()
            .zip(other.atom_bonds.iter())
            .enumerate()
            .for_each(|(atom_index, (self_bonds, other_bonds))| {
                // Check for broken bonds
                // If the other molecule does not contain the bond, it must have been broken
                for &self_bond in self_bonds {
                    if !other_bonds
                        .iter()
                        .map(|other_bond| other_bond.target())
                        .contains(&self_bond.target())
                        && self_bond.target() > atom_index
                    {
                        bond_changes.push(BondChange::broken(atom_index, self_bond.target()));
                    }
                }

                // Check for formed bonds
                // If the self molecule does not contain the bond, it must have been formed
                for &other_bond in other_bonds {
                    if !self_bonds
                        .iter()
                        .map(|self_bond| self_bond.target())
                        .contains(&other_bond.target())
                        && other_bond.target() > atom_index
                    {
                        bond_changes.push(BondChange::formed(atom_index, other_bond.target()));
                    }
                }
            });
        bond_changes
    }

    pub fn charge_changes(&self, other: &Molecule3D) -> Vec<(usize, i8, i8)> {
        self.charges
            .iter()
            .zip(&other.charges)
            .enumerate()
            .filter_map(|(index, (&charge1, &charge2))| {
                if charge1 != charge2 {
                    Some((index, charge1, charge2))
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn update_atom_charge(&mut self, atom_index: usize, charge: i8) {
        if atom_index < self.charges().len() {
            self.charge += charge as i32;
            self.charges_mut()[atom_index] += charge;
        } else {
            println!("Atom index out of bounds, could not update charge");
        }
    }

    pub fn identify_difference(&self, other: &Molecule3D) {
        let _charge_changes = self.charge_changes(other);
        let _bond_changes = self.identify_bond_changes(other);
        let _bond_type_changes = self.identify_bond_type_changes(other);
    }

    fn identify_bond_type_changes(&self, other: &Molecule3D) -> Vec<BondTypeChange> {
        self.atom_bonds
            .iter()
            .zip(&other.atom_bonds)
            .enumerate()
            .flat_map(|(index, (self_bonds, other_bonds))| {
                self_bonds
                    .iter()
                    .filter_map(|self_bond| {
                        other_bonds
                            .iter()
                            .find(|other_bond| other_bond.target() == self_bond.target())
                            .and_then(|other_bond| {
                                if self_bond.bond_type() != other_bond.bond_type() {
                                    Some(BondTypeChange {
                                        atom_index: index,
                                        target: self_bond.target(),
                                        from: self_bond.bond_type(),
                                        to: other_bond.bond_type(),
                                    })
                                } else {
                                    None
                                }
                            })
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    }
    pub fn identify_bonds(&mut self, threshold: f64) {
        let threshold_squared = threshold.powi(2);
        let kdtree = self.build_tree();
        let bonds = self
            .positions
            .par_iter()
            .enumerate()
            .map(|(index, position)| {
                let mut bonds = kdtree
                    .within::<SquaredEuclidean>(
                        &[position.x, position.y, position.z],
                        threshold_squared,
                    )
                    .iter()
                    .filter_map(|neighbor| {
                        if neighbor.item == index as u64 {
                            return None;
                        }
                        let distance = neighbor.distance;
                        let is_bonded =
                            self.is_bonded(index, neighbor.item as usize, distance, BOND_TOLERANCE);
                        if is_bonded {
                            Some(BondTarget::single(neighbor.item as usize))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<BondTarget>>();
                bonds.sort_by_key(|bond| bond.target);
                bonds
            })
            .collect::<Vec<Vec<BondTarget>>>();
        self.atom_bonds = bonds;
    }
    /// This function checks if two atoms are bonded based on their atomic numbers and distance
    fn is_bonded(
        &self,
        index1: usize,
        index2: usize,
        squared_distance: f64,
        tolerance: f64,
    ) -> bool {
        let Some(atomic_number1) = self.atomic_numbers.get(index1) else {
            return false;
        };
        let Some(atomic_number2) = self.atomic_numbers.get(index2) else {
            return false;
        };
        if let (Some(covalent_radius1), Some(covalent_radius2)) = (
            atomic_number1.covalent_radius(),
            atomic_number2.covalent_radius(),
        ) {
            squared_distance < ((covalent_radius1 + covalent_radius2) * tolerance).powi(2)
        } else {
            false
        }
    }

    pub fn charge(&self) -> i32 {
        if self.atomic_numbers.len() > 1000 {
            self.charges.par_iter().map(|&charge| charge as i32).sum()
        } else {
            self.charges.iter().map(|&charge| charge as i32).sum()
        }
    }

    /// Returns a vector to the center of the Molecule
    ///
    /// # Examples
    /// ```
    /// use molecules::prelude::*;
    /// let atom1 = Atom::new(6).with_position((0.0, 0.0, 0.0));
    /// let atom2 = Atom::new(6).with_position((1.0, 0.0, 0.0));
    /// let atom3 = Atom::new(6).with_position((2.0, 0.0, 0.0));
    /// let molecule = Molecule3D::from_atoms(vec![atom1, atom2, atom3]);
    /// assert_eq!(molecule.center(), Vector::new(1.0, 0.0, 0.0));
    /// ```
    pub fn center(&self) -> Vector {
        let mut sum = Vector::new(0.0, 0.0, 0.0);

        for &position_vector in &self.positions {
            sum += position_vector;
        }

        let count = self.atomic_numbers.len() as f64;
        sum / count
    }

    pub fn number_of_charged_atoms(&self) -> usize {
        self.charges.iter().filter(|&&charge| charge != 0).count()
    }
    pub fn find_angles(&self) -> Vec<BondAngle> {
        let bond_angles: Vec<BondAngle> = self
            .atom_bonds
            .par_iter()
            .enumerate()
            .flat_map(|(atom_index, bonds)| {
                let mut local_bond_angles = Vec::with_capacity(12);
                bonds.iter().combinations(2).for_each(|slice| {
                    let neighbor1_index = slice[0].target();
                    let neighbor2_index = slice[1].target();
                    let position_vector1 = &self.positions[neighbor1_index];
                    let position_vector2 = &self.positions[atom_index];
                    let position_vector3 = &self.positions[neighbor2_index];

                    let angle =
                        position_vector1.angle_between_points(position_vector2, position_vector3);
                    local_bond_angles.push(BondAngle::new(
                        angle,
                        (neighbor1_index, atom_index, neighbor2_index),
                    ));
                });
                local_bond_angles
            })
            .collect();
        bond_angles
    }
    /// This function identifies all dihedral angles in a molecule
    /// It is currently broken due to some internal changes
    ///
    /// # Panics
    /// This function does not panic
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::prelude::*;
    ///
    /// let mut molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    /// let dihedrals = molecule.dihedrals();
    /// assert_eq!(dihedrals.len(), 9);
    ///
    ///```
    // PROTOTYPE
    pub fn dihedrals(&mut self) -> Vec<((usize, usize, usize, usize), f64)> {
        let bond_angles = self.find_angles();
        let dihedrals: Vec<((usize,usize,usize,usize),f64)>= bond_angles
            .iter()
            .flat_map(|angle| {
                let (atom1, atom2, atom3) = angle.atoms;

                let Some(bonds) = self.atom_bonds.get(atom1) else {
                    return Vec::new();
                };

                bonds
                    .iter()
                    .filter_map(|&bond| {
                        if bond.target() != atom3 && bond.target() != atom2 && atom3 < bond.target()
                        {
                            let dihedral = (bond.target(), atom1, atom2, atom3);
                            self.dihedral_angle(&dihedral)
                                .map(|dihedral_angle| (dihedral, dihedral_angle))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        dihedrals
    }
    /// This function calculates the dihedral angle for all atoms
    ///
    /// # Panics
    /// This function does not panic.
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::molecule::{Molecule3D,Molecule};
    /// let mut molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    /// let dihedral_angle = molecule.dihedral_angle(&(0, 1, 2, 3));
    /// assert_eq!(dihedral_angle, Some(0.5807210503904102));
    /// ```
    pub fn dihedral_angle(&self, dihedral: &(usize, usize, usize, usize)) -> Option<f64> {
        if [dihedral.0, dihedral.1, dihedral.2, dihedral.3]
            .iter()
            .any(|&index| index >= self.positions.len())
        {
            return None;
        }
        let (a, b, c, d) = (
            self.positions[dihedral.0],
            self.positions[dihedral.1],
            self.positions[dihedral.2],
            self.positions[dihedral.3],
        );

        let v3 = d - c;
        let v2 = c - b;
        let v1 = b - a;
        let normal1 = v1.cross(&v2);
        let normal2 = v2.cross(&v3);
        let angle = normal1.angle_between(&normal2);
        let sign = normal1.cross(&normal2).dot(&v2);
        if sign < 0.0 {
            Some(-angle)
        } else {
            Some(angle)
        }
    }
    /// This function calculates the degrees of all atoms
    ///
    /// # Example
    /// ```
    /// use molecules::molecule::{Molecule3D,Molecule};
    /// let mut molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    /// let degrees = molecule.degrees();
    /// assert_eq!(degrees, vec![0, 0, 0, 0, 0, 0, 0, 0]);
    /// ```
    pub fn degrees(&self) -> Vec<i8> {
        (0..self.atomic_numbers.len())
            .map(|index| self.degree(index).unwrap_or(0))
            .collect::<Vec<i8>>()
    }

    pub fn degree(&self, atom_index: usize) -> Option<i8> {
        let atomic_number = self.atomic_numbers[atom_index];
        let expected_valency = atomic_number.valencies()?.next()?;
        let actual_valency = self.actual_valency(atom_index);
        // May need to be changed for elements with unknown valencies
        Some(actual_valency - expected_valency)
    }
    /// Return the actual valency of the atom based on the number of bonds
    ///
    /// # Examples
    /// ```
    /// use molecules::prelude::*;
    /// let atom = Atom::new(6);
    /// assert_eq!(atom.actual_valency(), 0);
    /// ```
    pub fn actual_valency(&self, atom_index: usize) -> i8 {
        self.atom_bonds[atom_index]
            .iter()
            .map(|bond| match bond.bond_type() {
                BondType::Single => 2,
                BondType::Double => 4,
                BondType::Triple => 6,
                BondType::Quadruple => 8,
                BondType::Aromatic => 3,
            })
            .sum::<i8>()
            / 2
            + if self.radical_states[atom_index] { 1 } else { 0 }
            + self.get_charge(atom_index).abs()
    }
    pub fn shift_charge(&mut self, index: usize, neighbor_index: usize, charge: i8) {
        self.update_atom_charge(index, -charge);
        self.update_atom_charge(neighbor_index, charge);
    }

    pub fn update_charge(&mut self, charge: i8) {
        self.charge += charge as i32;
    }

    pub fn identify_bond_types(&mut self) {
        //, previous_state: Option<IntMap<usize,Vec<Bond>>>) {
        let mut degrees = self.degrees();

        //TODO: Implement a way to keep track of previous state

        if degrees.iter().all(|&a| a == 0) {
            return;
        }

        // I am not sure if this is the right place for this.
        // filter for both carbon and a positive degree
        let oversaturated_carbons: Vec<(usize, i8)> = degrees
            .iter()
            .enumerate()
            .flat_map(|(index, degree)| {
                if self.atomic_numbers[index] == 6 && *degree > 0 {
                    Some((index, *degree))
                } else {
                    None
                }
            })
            .collect();

        oversaturated_carbons.iter().for_each(|(index, degree)| {
            let bonds = &self.atom_bonds[*index];
            let mut max_degree = 0;
            let mut max_degree_index = 0;

            for bond in bonds {
                let neighbor_index = bond.target();
                let neighbor_number_of_carbons = self.number_of_bonded_element(neighbor_index, 6);
                if neighbor_number_of_carbons > max_degree {
                    max_degree = neighbor_number_of_carbons;
                    max_degree_index = neighbor_index;
                }
            }

            if max_degree != 0 {
                self.shift_charge(*index, max_degree_index, *degree);
            }
        });

        if !backtrack_bonding(self, &mut degrees) {
            relaxed_backtrack_bonding(self, &mut degrees);
        }
    }

    pub fn split_components(&self) -> Vec<Molecule3D> {
        let components = self.get_components();
        if components.len() == 1 {
            return vec![self.clone()];
        }
        let mut new_atom_indices: IntMap<usize, usize> = IntMap::default();
        for component in &components {
            for (new_index, &old_index) in component.iter().enumerate() {
                new_atom_indices.insert(old_index, new_index);
            }
        }

        let mut molecules = Vec::with_capacity(components.len());
        for component in components {
            let atomic_numbers: Vec<u8> = component
                .iter()
                .map(|&index| self.atomic_numbers[index])
                .collect();
            let mut new_bonds = Vec::with_capacity(component.len());
            let new_positions = component
                .iter()
                .map(|&index| self.positions[index])
                .collect();
            let new_charges = component.iter().map(|&index| self.charges[index]).collect();
            let new_is_radical = component
                .iter()
                .map(|&index| self.radical_states[index])
                .collect();

            // TODO: Properly implement thi
            for &index in component.iter() {
                let mut bonds = self.atom_bonds[index].clone();
                bonds
                    .iter_mut()
                    .for_each(|bond| bond.target = new_atom_indices[&bond.target]);
                new_bonds.push(bonds);
            }

            // TODO: Not sure if this is the right way to handle this
            if component.len() == 2 && atomic_numbers.iter().all(|&atom| atom == 7) {
                continue;
            }

            let mut molecule = Molecule3D {
                atomic_numbers,
                charges: new_charges,
                radical_states: new_is_radical,
                atom_bonds: new_bonds,
                positions: new_positions,
                ..Default::default()
            };
            molecule.identify_bonds(BOND_SEARCH_THRESHOLD);
            molecules.push(molecule);
        }
        molecules
    }
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct MolecularFormula {
    elements: IntMap<u8, usize>,
}

impl MolecularFormula {
    pub fn new(elements: IntMap<u8, usize>) -> Self {
        Self { elements }
    }
    pub fn from_molecule(molecule: &(impl Molecule + ?Sized)) -> Self {
        let mut formula = MolecularFormula::default();
        for &atomic_number in molecule.atomic_numbers() {
            *formula.elements.entry(atomic_number).or_insert(0) += 1;
        }
        formula
    }
    pub fn elements(&self) -> &IntMap<u8, usize> {
        &self.elements
    }
    /// Merges two molecular formulas.
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// let mut water_methane = water.clone();
    /// water_methane.merge(&methane);
    /// assert_eq!(water_methane, "CH4H2O".parse::<MolecularFormula>().unwrap());
    /// assert_eq!(water, "H2O".parse::<MolecularFormula>().unwrap());
    /// assert_eq!(methane, "CH4".parse::<MolecularFormula>().unwrap());
    /// assert_eq!(water_methane, "CH4H2O".parse::<MolecularFormula>().unwrap());
    ///
    /// ```

    pub fn merge(&mut self, other: &MolecularFormula) {
        for (atom, count) in &other.elements {
            *self.elements.entry(*atom).or_insert(0) += count;
        }
    }

    /// Combines two molecular formulas.
    /// # Examples
    /// ```
    /// use molecules::molecule::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// let water_methane = water.combine(&methane);
    /// assert_eq!(water_methane, "CH4H2O".parse::<MolecularFormula>().unwrap());
    /// ```
    pub fn combine(&self, other: &MolecularFormula) -> Self {
        let mut combined = MolecularFormula::default();
        for (atom, count) in &self.elements {
            *combined.elements.entry(*atom).or_insert(0) += count;
        }
        for (atom, count) in &other.elements {
            *combined.elements.entry(*atom).or_insert(0) += count;
        }
        combined
    }
    /// Calculates the monoisotopic mass of the molecular formula.
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// println!("{:?}", water);
    /// assert_eq!(water.monoisotopic_mass(), 18.01056468403);
    /// assert_eq!(methane.monoisotopic_mass(), 16.03130012892);
    ///
    /// ```
    pub fn monoisotopic_mass(&self) -> f64 {
        self.elements.iter().fold(0.0, |acc, (atom, count)| {
            // TODO check if this is safe
            acc + atom.monoisotopic_mass().unwrap() * *count as f64
        })
    }

    /// Calculates the molecular mass of the molecular formula.
    /// # Examples
    /// ```
    /// use molecules::molecule::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// assert_eq!(water.molecular_mass().round(), 18.0);
    /// assert_eq!(methane.molecular_mass().round(), 16.0);
    pub fn molecular_mass(&self) -> f64 {
        self.elements.iter().fold(0.0, |acc, (atom, count)| {
            // TODO check if this is safe
            acc + atom.standard_atomic_weight().unwrap() * *count as f64
        })
    }
}

impl Display for MolecularFormula {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut elements = self
            .elements
            .iter()
            .map(|(&atom, &count)| (atom.atomic_symbol().unwrap().to_string(), count))
            .collect::<Vec<_>>();
        elements.sort_by(|(atom1, _), (atom2, _)| atom1.cmp(atom2));
        for (atom, count) in elements {
            write!(f, "{}{}", atom.clone(), count)?;
        }
        Ok(())
    }
}

use std::str::FromStr;

#[derive(Debug)]
pub struct ParseFormulaError;

impl std::fmt::Display for ParseFormulaError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "could not parse the provided string as a molecular formula"
        )
    }
}

impl std::error::Error for ParseFormulaError {}

impl FromStr for MolecularFormula {
    type Err = ParseFormulaError;

    /// Parses a molecular formula from a string.
    ///
    /// # Panics
    ///
    /// This function does not panic. But it will ignore everything that is neither a number nor a
    /// letter
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// assert_eq!(water.monoisotopic_mass(), 18.01056468403);
    /// assert_eq!(methane.monoisotopic_mass(), 16.03130012892);
    ///
    /// ```

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(MolecularFormula::default());
        }
        let mut number_buffer = String::new();
        let mut element_buffer = String::new();
        let mut elements: IntMap<u8, usize> = IntMap::default();
        for character in s.chars() {
            if character.is_alphabetic() {
                if !element_buffer.is_empty() && character.is_uppercase() {
                    let Some(atomic_number) = element_buffer.as_str().atomic_number() else {
                        return Err(ParseFormulaError);
                    };
                    let count = if !number_buffer.is_empty() {
                        number_buffer.parse::<usize>().unwrap_or(1)
                    } else {
                        1
                    };
                    *elements.entry(atomic_number).or_insert(0) += count;
                    number_buffer.clear();
                    element_buffer.clear();
                }
                element_buffer.push(character);
            } else if character.is_numeric() {
                number_buffer.push(character);
            }
        }
        if !element_buffer.is_empty() {
            let Some(atomic_number) = element_buffer.as_str().atomic_number() else {
                return Err(ParseFormulaError);
            };
            let count = if !number_buffer.is_empty() {
                number_buffer.parse::<usize>().unwrap_or(1)
            } else {
                1
            };
            *elements.entry(atomic_number).or_insert(0) += count;
        }
        Ok(MolecularFormula { elements })
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondAngle {
    angle: f64,
    atoms: (usize, usize, usize),
}

impl BondAngle {
    fn new(angle: f64, atoms: (usize, usize, usize)) -> Self {
        Self { angle, atoms }
    }
}
#[derive(Debug, Clone, PartialEq, Default)]
pub enum BondState {
    Formed,
    #[default]
    Broken,
}

#[derive(Debug, Default)]
pub struct BondChange {
    atom1_index: usize,
    atom2_index: usize,
    bond_state: BondState,
}

impl Display for BondChange {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self.bond_state {
            BondState::Formed => write!(
                f,
                "Bond formed: {} and {}",
                self.atom1_index, self.atom2_index
            ),
            BondState::Broken => write!(
                f,
                "Bond broken: {} and {}",
                self.atom1_index, self.atom2_index
            ),
        }
    }
}

impl BondChange {
    pub fn broken(atom1_index: usize, atom2_index: usize) -> Self {
        BondChange {
            atom1_index,
            atom2_index,
            bond_state: BondState::Broken,
        }
    }
    pub fn formed(atom1_index: usize, atom2_index: usize) -> Self {
        BondChange {
            atom1_index,
            atom2_index,
            bond_state: BondState::Formed,
        }
    }
    pub fn is_broken(&self) -> bool {
        matches!(self.bond_state, BondState::Broken)
    }
    pub fn is_formed(&self) -> bool {
        matches!(self.bond_state, BondState::Formed)
    }
    pub fn atom1_index(&self) -> usize {
        self.atom1_index
    }
    pub fn atom2_index(&self) -> usize {
        self.atom2_index
    }
    pub fn bond_state(&self) -> &BondState {
        &self.bond_state
    }
}

struct SmilesParameters {
    n_bonds: u8,
    n_non_hydrogen: u8,
    atomic_number: u8,
    sign_of_charge: u8,
    n_hydrogen: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    index: usize,
    bond_type: BondType,
    children: Vec<Node>,
}

impl Node {
    fn new(index: usize, bond_type: BondType) -> Self {
        Node {
            index,
            bond_type,
            children: Vec::new(),
        }
    }

    fn add_child(&mut self, child: Node) {
        self.children.push(child);
    }
}

fn backtrack_bonding(molecule: &mut Molecule3D, degrees: &mut [i8]) -> bool {
    // We can include the expected charge as an argument here that we have derived from the
    // simulation summary file

    if degrees.iter().all(|&degree| degree >= 0) {
        return true;
    }

    // Find the first unsatisfied atom (negative degree)
    let mut unsatisfied_atom_index = degrees.iter().position(|&degree| degree < 0).unwrap();

    let number_of_double_bonds = molecule.atom_bonds[unsatisfied_atom_index]
        .iter()
        .filter(|&bond| bond.bond_type == BondType::Double)
        .count();

    // If there are more than two unsatisfied atoms and an atom already has a double bond then we need to find the next unsatisfied atom
    if number_of_double_bonds > 0 && degrees.iter().filter(|&degree| *degree < 0).count() > 2 {
        for (index, degree) in degrees.iter().enumerate() {
            if *degree < 0 && index > unsatisfied_atom_index {
                unsatisfied_atom_index = index;
                break;
            }
        }
    }

    for &neighbor in &molecule.atom_bonds[unsatisfied_atom_index].clone() {
        let number_of_neighbor_atom_double_bonds = molecule.atom_bonds[neighbor.target]
            .iter()
            .filter(|bond| bond.bond_type == BondType::Double)
            .count();

        if can_increase_bond(unsatisfied_atom_index, neighbor.target, degrees)
            && number_of_neighbor_atom_double_bonds == 0
        {
            increase_bonds(molecule, unsatisfied_atom_index, neighbor.target);
            update_degrees(degrees, unsatisfied_atom_index, neighbor.target, false);
            if backtrack_bonding(molecule, degrees) {
                return true;
            }
            decrease_bonds(molecule, unsatisfied_atom_index, neighbor.target);
            update_degrees(degrees, unsatisfied_atom_index, neighbor.target, true);
        }
    }
    false
}

pub fn relaxed_backtrack_bonding(molecule: &mut Molecule3D, degrees: &mut [i8]) -> bool {
    // If all degrees are zero or only one degree is negative (positive charge) then we are done
    if degrees.iter().all(|&degree| degree >= 0) {
        return true;
    }

    // Find the first unsatisfied atom (negative degree)
    let unsatisfied_atom_index = degrees.iter().position(|&degree| degree < 0).unwrap();

    for &neighbor in &molecule.atom_bonds[unsatisfied_atom_index].clone() {
        let number_of_neighbor_atom_double_bonds = molecule.atom_bonds[neighbor.target]
            .iter()
            .filter(|bond| bond.bond_type == BondType::Double)
            .count();

        if can_increase_bond(unsatisfied_atom_index, neighbor.target, degrees)
            && number_of_neighbor_atom_double_bonds == 0
        {
            increase_bonds(molecule, unsatisfied_atom_index, neighbor.target);
            update_degrees(degrees, unsatisfied_atom_index, neighbor.target, false);
            if relaxed_backtrack_bonding(molecule, degrees) {
                return true;
            }
            decrease_bonds(molecule, unsatisfied_atom_index, neighbor.target);
            update_degrees(degrees, unsatisfied_atom_index, neighbor.target, true);
        }
    }

    false
}

pub fn can_increase_bond(atom_index: usize, neighbor_index: usize, degrees: &[i8]) -> bool {
    degrees[atom_index] < 0 && degrees[neighbor_index] < 0
}

pub fn increase_bonds(molecule: &mut Molecule3D, atom_index: usize, neighbor_index: usize) {
    for bond in &mut molecule.atom_bonds[atom_index] {
        if bond.target == neighbor_index {
            increase_bond(bond);
        }
    }
    for bond in &mut molecule.atom_bonds[neighbor_index] {
        if bond.target == atom_index {
            increase_bond(bond);
        }
    }
}

pub fn increase_bond(bond: &mut BondTarget) {
    match bond.bond_type {
        BondType::Single => bond.bond_type = BondType::Double,
        BondType::Double => bond.bond_type = BondType::Triple,
        // TODO handle this in the case of metals
        BondType::Triple => bond.bond_type = BondType::Quadruple,
        BondType::Quadruple => panic!("Cannot increase bond beyond quadruple bond"),
        BondType::Aromatic => panic!("Cannot increase bond beyond aromatic bond"),
    }
}

pub fn decrease_bonds(molecule: &mut Molecule3D, atom_index: usize, neighbor_index: usize) {
    for bond in &mut molecule.atom_bonds[atom_index] {
        if bond.target == neighbor_index {
            decrease_bond(bond);
        }
    }
    for bond in &mut molecule.atom_bonds[neighbor_index] {
        if bond.target == atom_index {
            decrease_bond(bond);
        }
    }
}

pub fn decrease_bond(bond: &mut BondTarget) {
    match bond.bond_type {
        BondType::Single => panic!("Cannot decrease bond beyond single bond"),
        BondType::Double => bond.bond_type = BondType::Single,
        BondType::Triple => bond.bond_type = BondType::Double,
        BondType::Quadruple => bond.bond_type = BondType::Triple,
        BondType::Aromatic => bond.bond_type = BondType::Single,
    }
}

pub fn update_degrees(
    degrees: &mut [i8],
    atom_index: usize,
    neighbor_index: usize,
    decrease: bool,
) {
    let factor = if decrease { -1 } else { 1 };
    degrees[atom_index] += factor;
    degrees[neighbor_index] += factor;
}

/// This function reads a pdb file line and extracts the atom information
///
/// # Panics
/// This function returns a ParseFloatError if the line has errors in the float region, it does not check for other errors.
///
/// ```
/// use molecules::molecule::extract_atom_pdb;
///
/// let line =  "ATOM   2073  CB  ALA B 128      11.390 -11.172  71.797  1.00 16.79           C";
/// let atom = extract_atom_pdb(line).unwrap();
/// let position = atom.position_vector.unwrap();
/// assert_eq!(position.x, 11.390, "Incorrect x coordinate");
/// assert_eq!(position.y, -11.172, "Incorrect y coordinate");
/// assert_eq!(position.z, 71.797, "Incorrect z coordinate");
/// assert_eq!(atom.atomic_symbol(), Some("C"), "Incorrect atom name");
/// assert_eq!(atom.atomic_number(), 6, "Incorrect atom type");
/// ```
pub fn extract_atom_pdb(line: &str) -> Result<Atom, ParseFloatError> {
    // This is a very ugly function, but it is also very fast
    // TODO; Make this function more readable and robust
    let position = Vector {
        x: line[31..=37].trim().parse::<f64>()?,
        y: line[38..=45].trim().parse::<f64>()?,
        z: line[46..=53].trim().parse::<f64>()?,
    };
    let name = line[76..=77].trim().to_string();
    let atomic_number = name
        .as_str()
        .atomic_number()
        .unwrap_or_else(|| panic!("Could not find atom with symbol {}", name));
    Ok(Atom {
        position_vector: Some(position),
        atomic_number,
        ..Default::default()
    })
}
/// This function reads a cif file line and extracts the atom information
///
/// # Panics
///
/// This function returns a ParseFloatError if the line has errors in the float region, it does not check for other errors.
///
/// # Example
///
/// ```
/// use molecules::prelude::*;
/// use molecules::molecule::extract_atom_cif;
///
/// let line = "ATOM   1    N N   . GLN A 1 1   ? 201.310 198.892 131.429 1.00 70.25  ? 1   GLN A N   1";
/// let atom = extract_atom_cif(line).unwrap();
/// let position = atom.position_vector.unwrap();
/// assert_eq!(position.x, 201.310, "Incorrect x coordinate");
/// assert_eq!(position.y, 198.892, "Incorrect y coordinate");
/// assert_eq!(position.z, 131.429, "Incorrect z coordinate");
/// assert_eq!(atom.atomic_symbol(), Some("N"), "Incorrect atom name");
/// assert_eq!(atom.atomic_number(), 7, "Incorrect atom type");
/// ```
pub fn extract_atom_cif(line: &str) -> Result<Atom, ParseFloatError> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    let atomic_number = fields[2]
        .atomic_number()
        .unwrap_or_else(|| panic!("Could not find atom with symbol {}", fields[2]));

    let position = Vector {
        x: fields[10].parse::<f64>()?,
        y: fields[11].parse::<f64>()?,
        z: fields[12].parse::<f64>()?,
    };
    Ok(Atom {
        position_vector: Some(position),
        atomic_number,
        ..Default::default()
    })
}

impl Molecule3D {
    /// Finds the shortest path of a set of atoms
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::{Molecule3D,Molecule};
    ///
    /// let mut molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    /// let atom_set = vec![0,1,2].into_iter().collect();
    ///
    /// let mut shortest_path = molecule.find_shortest_path_of_set(&atom_set).unwrap();
    /// // HashSets are unordered, so we need to sort the path
    /// if !(shortest_path[0] < shortest_path[shortest_path.len()-1]) {
    ///     shortest_path.reverse();
    /// }
    /// assert_eq!(shortest_path, vec![1,0,2]);
    /// ```
    ///
    pub fn find_shortest_path_of_set(&self, atom_set: &IntSet<usize>) -> Option<Vec<usize>> {
        let shortest_paths: Vec<Vec<usize>> = atom_set
            .iter()
            .combinations(2)
            .par_bridge() // This does not preserve order!
            .flat_map(|combination| self.shortest_path(*combination[0], *combination[1]))
            .collect();

        if shortest_paths.is_empty() {
            return None;
        }

        let shortest_path = shortest_paths
            .iter()
            .filter(|path| {
                path.len() >= atom_set.len() && atom_set.iter().all(|&atom| path.contains(&atom))
            })
            .sorted_by_key(|path| path.len())
            .collect::<Vec<_>>();

        shortest_path
            .first()
            .map(|shortest_path| shortest_path.to_vec())
    }

    /// Finds the shortest path between two atoms
    ///
    /// # Arguments
    ///
    /// * `start` - The index of the starting atom
    /// * `end` - The index of the ending atom
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::{Molecule3D,Molecule};
    ///
    /// let mut molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    /// let shortest_path = molecule.shortest_path(0,1).unwrap();
    /// assert_eq!(shortest_path, vec![0,1]);
    /// ```
    ///
    pub fn shortest_path(&self, start: usize, end: usize) -> Option<Vec<usize>> {
        let mut queue = VecDeque::new();
        let mut visited = vec![false; self.atomic_numbers.len()];
        let mut parent = vec![None; self.atomic_numbers.len()];

        visited[start] = true;
        queue.push_back(start);

        while let Some(current_node) = queue.pop_front() {
            if current_node == end {
                return Some(reconstruct_path(current_node, &parent));
            }
            for &neighbor in &self.atom_bonds[current_node] {
                let target = neighbor.target();
                if !visited[target] {
                    parent[target] = Some(current_node);
                    queue.push_back(target);
                    visited[target] = true;
                }
            }
        }
        None
    }

    fn find_minimum_cycle(&self, atom_set: &HashSet<usize>) -> Option<Vec<usize>> {
        let mut minimum_cycle = None;
        let mut minimum_cycle_length = usize::MAX;

        for &atom1 in atom_set {
            for &atom2 in atom_set {
                if atom1 != atom2 {
                    if let Some(cycle) = self.bfs_for_cycle(atom1, atom2) {
                        let cycle_length = cycle.len();
                        if cycle_length < minimum_cycle_length {
                            minimum_cycle = Some(cycle);
                            minimum_cycle_length = cycle_length;
                        }
                    }
                }
            }
        }

        minimum_cycle
    }

    fn bfs_for_cycle(&self, start: usize, end: usize) -> Option<Vec<usize>> {
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();
        let mut parent = vec![None; self.atomic_numbers.len()];

        queue.push_back(start);

        while let Some(current_node) = queue.pop_front() {
            if current_node == end {
                return Some(reconstruct_path(current_node, &parent));
            }

            for &neighbor in &self.atom_bonds[current_node] {
                if !visited.contains(&neighbor.target()) {
                    parent[neighbor.target()] = Some(current_node);
                    queue.push_back(neighbor.target());
                    visited.insert(neighbor.target());
                }
            }
        }
        None
    }

    fn canonicalize(&self) -> IntMap<usize, usize> {
        // Find the node with the highest class
        let ec_classes = self.morgans_algorithm(None);
        let start_node = ec_classes
            .iter()
            .enumerate()
            .max_by_key(|&(_, class)| class)
            .map(|(idx, _)| idx)
            .unwrap();

        // Initialize BFS structures
        let mut queue = VecDeque::new();
        queue.push_back(start_node);
        let mut assigned = IntMap::default();
        let mut assignment_order = 1;

        // Perform BFS, ordering neighbors based on their classes
        while let Some(current_node) = queue.pop_front() {
            if let Some(bonds) = self.atom_bonds.get(current_node) {
                // Extract necessary information about neighbors
                let mut neighbors_info: Vec<(usize, u8, BondType)> = bonds
                    .iter()
                    .map(|bond| {
                        let target = bond.target();
                        let atomic_number = self.atomic_numbers[target];
                        let bond_type = bond.bond_type();
                        (target, atomic_number, bond_type)
                    })
                    .collect();

                // Sort neighbors based on class, then atomic number, then bond type
                neighbors_info.sort_by(
                    |&(node_a, atomic_num_a, bond_type_a), &(node_b, atomic_num_b, bond_type_b)| {
                        match ec_classes[node_a].cmp(&ec_classes[node_b]) {
                            std::cmp::Ordering::Equal => match atomic_num_a.cmp(&atomic_num_b) {
                                std::cmp::Ordering::Equal => bond_type_a.cmp(&bond_type_b),
                                other => other,
                            },
                            other => other,
                        }
                    },
                );

                for &(neighbor, _, _) in &neighbors_info {
                    if let std::collections::hash_map::Entry::Vacant(entry) =
                        assigned.entry(neighbor)
                    {
                        entry.insert(assignment_order);
                        assignment_order += 1;
                        queue.push_back(neighbor);
                    }
                }
            }
        }

        assigned
    }

    pub fn charges(&self) -> &[i8] {
        self.charges.as_slice()
    }

    pub fn names(&self) -> Vec<&str> {
        self.atomic_numbers
            .iter()
            .map(|atom| atom.atomic_symbol().unwrap())
            .collect()
    }

    pub fn morgans_algorithm(&self, max_depth: Option<usize>) -> Vec<usize> {
        let mut vertex_degrees = self
            .atom_bonds
            .iter()
            .map(|bonds| bonds.len())
            .collect::<Vec<_>>();
        let mut buffer = vec![0; self.atomic_numbers.len()];
        let mut last_count = 0;
        let mut depth = 0;
        let max_depth = max_depth.unwrap_or(100);

        loop {
            depth += 1;
            if depth > max_depth {
                println!("Morgan's algorithm did not converge after 100 iterations");
                break;
            }

            for (index, bonds) in self.atom_bonds.iter().enumerate() {
                let sum = bonds.iter().map(|bond| vertex_degrees[bond.target()]).sum();
                buffer[index] = sum;
            }
            for (index, degree) in buffer.iter().enumerate() {
                vertex_degrees[index] = *degree;
            }

            let count = vertex_degrees.iter().unique().count();
            if count == last_count {
                break;
            }
            last_count = count;
        }
        vertex_degrees
    }
}

fn reconstruct_path(mut current_node: usize, parents: &[Option<usize>]) -> Vec<usize> {
    let mut path = Vec::new();

    while let Some(parent_atom) = parents[current_node] {
        path.push(current_node);
        current_node = parent_atom;
    }
    path.push(current_node);

    path.reverse();
    path
}

// Helper function to determine the order of indices based on their degrees
fn order_by_degree(degrees: &[usize]) -> Vec<usize> {
    let mut degree_indices = (0..degrees.len()).collect::<Vec<_>>();
    degree_indices.sort_by(|&a, &b| degrees[a].cmp(&degrees[b]).then_with(|| a.cmp(&b)));
    degree_indices
}

struct BondTypeChange {
    atom_index: usize,
    target: usize,
    from: BondType,
    to: BondType,
}
// Unit tests (Still needs to be improved)
#[cfg(test)]
mod tests {
    #[test]
    fn test_extract_atom() {
        let line =
            "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N  ";
        let atom = super::extract_atom_pdb(line).unwrap();
        let position = atom.position_vector.unwrap();
        assert_eq!(atom.atomic_number, 7);
        assert_eq!(position.x, 10.0);
        assert_eq!(position.y, 10.0);
        assert_eq!(position.z, 10.0);
    }
    #[test]
    fn test_vector_angle() {
        let v1 = super::Vector {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        };
        let v2 = super::Vector {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        };
        let angle = v1.angle_between(&v2);
        assert_eq!(angle, std::f64::consts::FRAC_PI_2);
    }
    #[test]
    fn test_cross_product() {
        let v1 = super::Vector {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        };
        let v2 = super::Vector {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        };
        let v3 = v1.cross(&v2);
        assert_eq!(v3.x, 0.0);
        assert_eq!(v3.y, 0.0);
        assert_eq!(v3.z, 1.0);
    }
    #[test]
    fn test_bond_angle() {
        use super::*;
        let atom1 = Atom::new(6).with_position((1.5, 0.0, 0.0));
        let atom2 = Atom::new(6).with_position((0.0, 0.0, 0.0));
        let atom3 = Atom::new(6).with_position((0.0, 0.0, 1.5));
        let molecule = Molecule3D::from_atoms(vec![atom1, atom2, atom3]);
        let angles = molecule.find_angles();
        println!("{:?}", angles);
        assert_eq!(angles.len(), 1);
        assert_eq!(angles[0].angle, std::f64::consts::FRAC_PI_2);
    }
}

use std::fmt;
#[derive(Debug)]
pub enum MoleculeError {
    BondAlreadyExists,
    BondNotFound,
}

impl fmt::Display for MoleculeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MoleculeError::BondAlreadyExists => write!(f, "Bond already exists"),
            MoleculeError::BondNotFound => write!(f, "Bond not found"),
        }
    }
}

impl std::error::Error for MoleculeError {}

fn is_hueckel_satisfied(number_of_pi_electrons: usize) -> bool {
    number_of_pi_electrons % 4 == 2
}
