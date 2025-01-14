use crate::{
    atom::Atom,
    chirality::{ChiralClass, ChiralClassifier},
    consts::BOND_TOLERANCE,
    molecule::base::Molecule,
    molecule::bond::{BondAngle, BondOrder, BondTarget},
    molecule::utils::{can_increase_bond, decrease_bond, reconstruct_path, update_degrees},
    vector::Vector,
};
use chemistry_consts::ElementProperties;
use itertools::Itertools;
use kiddo::{KdTree, SquaredEuclidean};
use nohash_hasher::{IntMap, IntSet};
use rayon::prelude::*;
use tinyvec::ArrayVec;

use std::{
    collections::{HashSet, VecDeque},
    default::Default,
    fs::File,
    io::{self, BufRead},
    path::Path,
};

#[derive(Debug, Default, Clone)]
pub struct Molecule3D {
    pub atomic_numbers: Vec<u8>,
    pub atom_classes: Option<Vec<u8>>,
    pub charges: Vec<i8>,
    pub chiral_classes: Option<Vec<ChiralClass>>,
    pub isotopes: Option<Vec<u16>>,
    pub positions: Vec<Vector>,
    pub radical_states: Vec<bool>,
    pub atom_bonds: Vec<ArrayVec<[BondTarget; 10]>>,
}

impl Molecule3D {
    pub fn atom_bonds_mut(&mut self, atom_index: usize) -> &mut ArrayVec<[BondTarget; 10]> {
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
        self.positions
            .push(atom.position_vector.unwrap_or_default());
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
    pub fn atom_bonds(&self) -> &Vec<ArrayVec<[BondTarget; 10]>> {
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

    pub fn update_atom_charge(&mut self, atom_index: usize, charge: i8) {
        if atom_index < self.charges().len() {
            self.charges_mut()[atom_index] += charge;
        } else {
            println!("Atom index out of bounds, could not update charge");
        }
    }

    pub fn identify_bonds(&mut self, tolerance: f64) {
        // The threshold is dynamically determined by the largest covalent radius in the molecule
        let threshold_squared = (self
            .atomic_numbers
            .iter()
            .fold(f64::NEG_INFINITY, |prev, &atomic_number| {
                prev.max(atomic_number.covalent_radius().unwrap_or_default())
            })
            * 3.0)
            .powi(2);
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
                            self.is_bonded(index, neighbor.item as usize, distance, tolerance);
                        if is_bonded {
                            Some(BondTarget::single(neighbor.item as usize))
                        } else {
                            None
                        }
                    })
                    .collect::<ArrayVec<[BondTarget; 10]>>();
                bonds.sort_by_key(|bond| bond.target());
                bonds
            })
            .collect::<Vec<ArrayVec<[BondTarget; 10]>>>();
        self.atom_bonds = bonds;
    }

    pub fn identify_bonds_alternate_covalent_radii(
        &mut self,
        covalent_radii: &[f64],
        tolerance: f64,
    ) {
        // The threshold is dynamically determined by the largest covalent radius in the molecule
        let threshold_squared =
            (self
                .atomic_numbers
                .iter()
                .fold(f64::NEG_INFINITY, |prev, &atomic_number| {
                    let covalent_radius = covalent_radii.get(atomic_number as usize);
                    prev.max(*covalent_radius.unwrap_or(&0.0))
                })
                * 3.0)
                .powi(2);
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
                        let atomic_number1 = self.get_atomic_number(index);
                        let atomic_number2 = self.get_atomic_number(neighbor.item as usize);
                        let covalent_radius1 = covalent_radii.get(atomic_number1 as usize);
                        let covalent_radius2 = covalent_radii.get(atomic_number2 as usize);
                        if covalent_radius1.is_none() || covalent_radius2.is_none() {
                            return None;
                        }

                        let is_bonded = distance
                            < ((covalent_radius1.unwrap() + covalent_radius2.unwrap()) * tolerance)
                                .powi(2);
                        if is_bonded {
                            Some(BondTarget::single(neighbor.item as usize))
                        } else {
                            None
                        }
                    })
                    .collect::<ArrayVec<[BondTarget; 10]>>();
                bonds.sort_by_key(|bond| bond.target());
                bonds
            })
            .collect::<Vec<ArrayVec<[BondTarget; 10]>>>();
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

        let count = self.len() as f64;
        sum / count
    }

    pub fn get_atom_position(&self, atom_index: usize) -> Option<Vector> {
        self.positions.get(atom_index).copied()
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

                    let Some(angle) =
                        position_vector1.angle_between_points(position_vector2, position_vector3)
                    else {
                        // This should never happen
                        return;
                    };

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
    pub fn dihedrals(&mut self) -> Vec<((usize, usize, usize, usize), f64)> {
        let bond_angles = self.find_angles();
        let dihedrals: Vec<((usize, usize, usize, usize), f64)> = bond_angles
            .iter()
            .flat_map(|angle| {
                let (atom1, atom2, atom3) = angle.atoms();

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

        let v1 = b - a;
        let v2 = c - b;
        let v3 = d - c;
        let normal1 = v1.cross(&v2);
        let normal2 = v2.cross(&v3);
        let angle = normal1.angle_between(&normal2)?;

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
        Some(
            actual_valency - expected_valency
                + self.get_charge(atom_index).abs()
                + self.is_atom_radical(atom_index) as i8,
        )
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
            .map(|bond| match bond.bond_order() {
                BondOrder::Single => 2,
                BondOrder::Double => 4,
                BondOrder::Triple => 6,
                BondOrder::Quadruple => 8,
                BondOrder::Aromatic => 3,
                BondOrder::Coordinate => 2,
            })
            .sum::<i8>()
            / 2
    }
    pub fn shift_charge(&mut self, index: usize, neighbor_index: usize, charge: i8) {
        self.update_atom_charge(index, -charge);
        self.update_atom_charge(neighbor_index, charge);
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
                let mut bonds = self.atom_bonds[index];
                bonds
                    .iter_mut()
                    .for_each(|bond| bond.target = new_atom_indices[&bond.target]);
                new_bonds.push(bonds);
            }

            // TODO: Not sure if this is the right way to handle this
            if component.len() == 2 && atomic_numbers.iter().all(|&atom| atom == 7) {
                continue;
            }

            let molecule = Molecule3D {
                atomic_numbers,
                charges: new_charges,
                radical_states: new_is_radical,
                atom_bonds: new_bonds,
                positions: new_positions,
                ..Default::default()
            };
            molecules.push(molecule);
        }
        molecules
    }

    pub fn rotate_around_center(&mut self, x_angle: f64, y_angle: f64, z_angle: f64) {
        let center = self.center();
        for position in self.positions.iter_mut() {
            *position = position.rotate_around(center, x_angle, y_angle, z_angle);
        }
    }

    fn to_smiles_with_preprocessing(&mut self) -> String {
        if self.chiral_classes().is_none() {
            self.identify_chiral_classes();
        }
        self.to_smiles()
    }
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
                let mut neighbors_info: Vec<(usize, u8, BondOrder)> = bonds
                    .iter()
                    .map(|bond| {
                        let target = bond.target();
                        let atomic_number = self.atomic_numbers[target];
                        let bond_type = bond.bond_order();
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
    // Molecule3D specific implementation
    fn from_atoms_raw(atoms: Vec<Atom>) -> Self {
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

        let mut chirals = None;
        if atoms
            .iter()
            .any(|atom| atom.chiral_class != ChiralClass::None)
        {
            chirals = Some(atoms.into_iter().map(|atom| atom.chiral_class).collect());
        }
        Molecule3D {
            atomic_numbers,
            charges,
            radical_states,
            atom_bonds,
            positions,
            isotopes,
            chiral_classes: chirals,
            ..Default::default()
        }
    }

    pub fn from_atoms_alternate_covalent_radii(atoms: Vec<Atom>, covalent_radii: &[f64]) -> Self {
        let mut molecule = Molecule3D::from_atoms_raw(atoms);
        molecule.identify_bonds_alternate_covalent_radii(covalent_radii, BOND_TOLERANCE);
        molecule
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
/// This function tries to increase bonds in a greedy way to satisfy all bonds to saturate all atoms valence shells
fn backtrack_bonding(molecule: &mut Molecule3D, degrees: &mut [i8]) -> bool {
    if degrees.iter().all(|&degree| degree >= 0) {
        return true;
    }

    // Find the first unsatisfied atom (negative degree)
    let mut unsatisfied_atom_index = degrees.iter().position(|&degree| degree < 0).unwrap();

    let number_of_double_bonds = molecule.atom_bonds[unsatisfied_atom_index]
        .iter()
        .filter(|&bond| bond.bond_order == BondOrder::Double)
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
            .filter(|bond| bond.bond_order == BondOrder::Double)
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
    match bond.bond_order {
        BondOrder::Single => bond.bond_order = BondOrder::Double,
        BondOrder::Double => bond.bond_order = BondOrder::Triple,
        // TODO handle this in the case of metals
        BondOrder::Triple => bond.bond_order = BondOrder::Quadruple,
        BondOrder::Quadruple => panic!("Cannot increase bond beyond quadruple bond"),
        BondOrder::Aromatic => panic!("Cannot increase bond beyond aromatic bond"),
        // TODO handle this in the case of metals
        BondOrder::Coordinate => panic!("Cannot increase bond beyond coordinate bond"),
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
            .filter(|bond| bond.bond_order == BondOrder::Double)
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

impl Molecule for Molecule3D {
    fn atom_bonds(&self) -> &Vec<ArrayVec<[BondTarget; 10]>> {
        &self.atom_bonds
    }
    fn atom_bonds_mut(&mut self) -> &mut Vec<ArrayVec<[BondTarget; 10]>> {
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

    // Molecule3D specific implementation
    fn from_atoms(atoms: Vec<Atom>) -> Self {
        let mut molecule = Molecule3D::from_atoms_raw(atoms);
        molecule.identify_bonds(BOND_TOLERANCE);
        molecule
    }
}
