use crate::consts::{
    ATOMIC_NUMBERS, ATOMIC_RADII, ATOMIC_SYMBOLS, ELECTRONEGATIVITIES, MONOISOTOPIC_MASSES, PRIMES,
    STANDARD_ATOMIC_WEIGHTS, VALENCIES, BOND_SEARCH_THRESHOLD, BOND_TOLERANCE,
};
use crate::vector::Vector;
use core::fmt::{Display, Formatter};
use itertools::Itertools;
use kiddo::{KdTree, SquaredEuclidean};
use nohash_hasher::IntMap;
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

#[derive(Debug, Clone)]
pub struct Atom {
    pub atomic_number: u8,
    pub charge: i8,
    pub is_radical: bool,
    pub position_vector: Vector,
    bonds: Vec<Bond>,
}

impl Default for Atom {
    fn default() -> Atom {
        Atom {
            atomic_number: 1,
            charge: 0,
            is_radical: false,
            position_vector: Vector::default(),
            bonds: Vec::new(),
        }
    }
}

impl Atom {
    pub fn new(atomic_number: u8, position: (f64, f64, f64)) -> Atom {
        let position = Vector::new(position.0, position.1, position.2);
        Atom {
            atomic_number,
            position_vector: position,
            ..Default::default()
        }
    }

    pub fn degree(&self) -> i8 {
        self.actual_valency() - self.expected_valency() + self.charge().abs()
        //+ if self.is_radical() { 1 } else { 0 }
    }

    pub fn from_xyz_line(line: &str) -> Result<Atom, Box<dyn std::error::Error>> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 4 {
            let atomic_number = ATOMIC_NUMBERS[parts[0].to_uppercase().as_str()];
            let x = parts[1].parse::<f64>()?;
            let y = parts[2].parse::<f64>()?;
            let z = parts[3].parse::<f64>()?;
            Ok(Atom::new(atomic_number, (x, y, z)))
        } else {
            Err("Could not parse Line".into())
        }
    }
    /// Return the actual valency of the atom based on the number of bonds
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::Atom;
    /// let atom = Atom::new(6, (0.0, 0.0, 0.0));
    /// assert_eq!(atom.actual_valency(), 0);
    /// ```
    pub fn actual_valency(&self) -> i8 {
        self.bonds
            .iter()
            .map(|bond| match bond.bond_type() {
                BondType::Single => 1,
                BondType::Double => 2,
                BondType::Triple => 3,
                BondType::Aromatic => 1,
            })
            .sum::<i8>()
            + if self.is_radical { 1 } else { 0 }
    }
    /// Returns the expected valency of the atom based on its atomic number
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::Atom;
    /// let atom = Atom::new(6, (0.0, 0.0, 0.0));
    /// assert_eq!(atom.expected_valency(), 4);
    /// ```
    pub fn expected_valency(&self) -> i8 {
        match VALENCIES.get(self.atomic_number as usize) {
            Some(valency) => *valency,
            None => {
                println!("Could not find valency for {}", self.atomic_number);
                0
            }
        }
    }

    pub fn cmp_electronegativity(&self, other: &Self) -> std::cmp::Ordering {
        self.electronegativity().cmp(&other.electronegativity())
    }

    pub fn bonds(&self) -> &Vec<Bond> {
        &self.bonds
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn set_charge(&mut self, charge: i8) {
        self.charge = charge;
    }

    pub fn update_charge(&mut self, charge: i8) {
        self.charge += charge;
    }

    pub fn is_radical(&self) -> bool {
        self.is_radical
    }

    pub fn set_radical(&mut self, is_radical: bool) {
        self.is_radical = is_radical;
    }

    pub fn name(&self) -> &str {
        match ATOMIC_SYMBOLS.get(self.atomic_number as usize) {
            Some(atomic_symbol) => atomic_symbol,
            None => {
                println!("Could not find atomic symbol for {}", self.atomic_number);
                ""
            }
        }
    }

    pub fn electronegativity(&self) -> Option<u32> {
        ELECTRONEGATIVITIES
            .get(self.atomic_number as usize)
            .copied()
    }

    pub fn is_more_electronegative(&self, other: &Self) -> bool {
        self.electronegativity() > other.electronegativity()
    }

    fn monoisotopic_mass(&self) -> f64 {
        match MONOISOTOPIC_MASSES.get(self.atomic_number as usize) {
            Some(monoisotopic_mass) => *monoisotopic_mass,
            None => {
                println!(
                    "Could not find monoisotopic mass for {}",
                    self.atomic_number
                );
                0.0
            }
        }
    }

    /// Determines if two atoms are bonded based on their distance to each other, their expected covalent radii and a constant tolerance factor
    ///
    /// # Parameters
    /// * 'self' - A reference to the current atom.
    /// * 'other' - A reference to the other atom.
    /// * 'threshold' - The threshold for the distance between the two atoms.
    /// * 'tolerance' - The tolerance factor for the distance between the two atoms.
    ///
    /// # Returns
    /// A boolean value indicating whether the two atoms are bonded.
    ///
    /// # Example
    /// ```
    /// use molecules::molecule::Atom;
    /// use molecules::vector::Vector;
    /// use molecules::consts::BOND_TOLERANCE;
    /// let mut atom1 = Atom::default();
    /// let mut atom2 = Atom::default();
    /// atom1.position_vector = Vector::new(1.0,0.0,0.0);
    /// atom2.position_vector = Vector::new(0.0,0.0,0.0);
    /// assert!(atom1.is_bonded(&atom2, 0.5, BOND_TOLERANCE));
    /// ```
    pub fn is_bonded(&self, other: &Self, distance: f64, tolerance: f64) -> bool {
        distance < ((self.get_atomic_radius() + other.get_atomic_radius()) * tolerance).powi(2)
    }

    /// Calculates the distance between two atoms.
    ///
    /// # Arguments
    /// * 'self' - A reference to the current atom.
    /// * 'other' - A reference to the other atom.
    ///
    /// # Returns
    /// The distance between the two atoms as a 64-bit floating-point number.
    ///
    /// # Example
    /// ```
    /// use molecules::molecule::Atom;
    /// use molecules::vector::Vector;
    /// let mut atom1 = Atom::default();
    /// let mut atom2 = Atom::default();
    /// atom1.position_vector = Vector::new(1.0,0.0,0.0);
    /// atom2.position_vector = Vector::new(0.0,0.0,0.0);
    /// assert_eq!(atom1.distance(&atom2),1.0);
    /// ```
    pub fn distance(&self, other: &Atom) -> f64 {
        self.position_vector.distance(&other.position_vector)
    }

    /// Calculates the squared distance between two atoms saving computation time.
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::molecule::Atom;
    /// use molecules::vector::Vector;
    /// let mut atom1 = Atom::default();
    /// let mut atom2 = Atom::default();
    /// atom1.position_vector = Vector::new(1.0,0.0,0.0);
    /// atom2.position_vector = Vector::new(0.0,0.0,0.0);
    /// assert_eq!(atom1.distance_squared(&atom2),1.0);
    /// ```
    ///
    pub fn distance_squared(&self, other: &Atom) -> f64 {
        self.position_vector
            .distance_squared(&other.position_vector)
    }

    /// Calculates the potential energy of a diatomic system using the harmonic oscillator model.
    ///
    /// # Parameters
    /// * 'other' - A reference to the other atom.
    /// * 'force_constant' - The force constant (k) of the diatomic system.
    /// * 'equilibrium_distance' - The equilibrium distance between the two atoms.
    ///
    /// # Returns
    /// The potential energy of the diatomic system as a 64-bit floating-point number.
    pub fn potential_energy(
        &self,
        other: &Atom,
        force_constant: f64,
        equilibrium_distance: f64,
    ) -> f64 {
        0.5 * force_constant * (self.distance(other) - equilibrium_distance).powi(2)
    }

    pub fn standard_atomic_weight(&self) -> f64 {
        let atomic_weight = match STANDARD_ATOMIC_WEIGHTS.get(self.atomic_number as usize) {
            Some(atomic_weight) => *atomic_weight,
            None => {
                println!("Could not find atomic weight for {}", self.atomic_number);
                0.0
            }
        };
        atomic_weight
    }

    pub fn get_atomic_radius(&self) -> f64 {
        let atomic_radius = match ATOMIC_RADII.get(self.atomic_number as usize) {
            Some(atomic_radius) => *atomic_radius,
            None => {
                println!("Could not find atomic radius for {}", self.atomic_number);
                0.0
            }
        };
        atomic_radius
    }

    pub fn atomic_number(&self) -> u8 {
        self.atomic_number
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub struct Bond {
    target: usize,
    bond_type: BondType,
}

impl Bond {
    fn new(target: usize, bond_type: BondType) -> Bond {
        Bond { target, bond_type }
    }
    fn single(target: usize) -> Bond {
        Bond {
            target,
            bond_type: BondType::Single,
        }
    }
    fn double(target: usize) -> Bond {
        Bond {
            target,
            bond_type: BondType::Double,
        }
    }
    fn triple(target: usize) -> Bond {
        Bond {
            target,
            bond_type: BondType::Triple,
        }
    }
    fn aromatic(target: usize) -> Bond {
        Bond {
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
            // TODO: Implement aromatic bond order
            BondType::Aromatic => 1,
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub enum BondType {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl Display for BondType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            BondType::Single => write!(f, "-"),
            BondType::Double => write!(f, "="),
            BondType::Triple => write!(f, "#"),
            BondType::Aromatic => write!(f, ":"),
        }
    }
}

/// This function   
#[derive(Debug, Default, Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub charge: i32,
}

pub struct MolecularSystem {
    pub molecules: Vec<Molecule>,
}

impl MolecularSystem {
    pub fn new(molecules: Vec<Molecule>) -> Self {
        MolecularSystem { molecules }
    }
    pub fn from_xyz<P: AsRef<Path>>(path: P) -> MolecularSystem {
        let molecule = Molecule::from_xyz(&path);
        MolecularSystem::new(vec![molecule])
    }
    pub fn from_pdb<P: AsRef<Path>>(path: P) -> MolecularSystem {
        let file = File::open(&path).expect("Could not open file");
        let reader = io::BufReader::new(&file);

        let mut atom_groups: Vec<Vec<String>> = vec![Vec::new()];

        for line in reader.lines().flatten() {
            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                atom_groups.last_mut().unwrap().push(line);
            } else if line.starts_with("TER") {
                atom_groups.push(Vec::new());
            } else if line.starts_with("ENDMDL") {
                break;
            } else {
                continue;
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
                        Molecule::from_atoms(atoms)
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
                        Molecule::from_atoms(atoms)
                    })
                    .collect(),
            }
        } else {
            MolecularSystem {
                molecules: vec![Molecule::default()],
            }
        }
    }

    pub fn center(&self) -> Vector {
        let mut sum = Vector::new(0.0, 0.0, 0.0);
        let mut count = 0.0;
        for molecule in &self.molecules {
            for atom in &molecule.atoms {
                sum += atom.position_vector
            }
            count += molecule.atoms.len() as f64
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
            .map(|molecule| molecule.atoms.len())
            .sum()
    }

    pub fn find_angles(&self) -> Vec<Vec<BondAngle>> {
        self.molecules
            .iter()
            .map(|molecule| molecule.find_angles())
            .collect()
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

fn push_atom(atomic_number: u8, position: (f64, f64, f64), atoms: &mut Vec<Atom>) {
    atoms.push(Atom::new(atomic_number, position));
}

impl Molecule {

    pub fn from_atoms(atoms: Vec<Atom>) -> Self {
        Molecule {
            atoms,
            ..Default::default()
        }
    }
    
    pub fn from_smiles(smiles: &str) -> Self {
        let mut atoms: Vec<Atom> = Vec::new();
        let mut ring_number: Option<usize> = None;
        let mut current_atom: Option<Atom> = None;
        let mut current_atom_index: usize = 0;
        let mut is_negative = false;
        let mut element_buffer = String::new();
        
        for (index, c) in smiles.chars().enumerate() {
            match c {
               'A'..='Z' => {
                    if let Some(atom) = current_atom {
                        atoms.push(atom);
                        current_atom = None;
                    }
                    element_buffer.push(c);
                    current_atom_index = index;
                }, 

                'a'..='z' => {
                    element_buffer.push(c);
                },

                '-' => {
                    is_negative = true;
                    if let Some(atom) = current_atom {
                        atoms.push(atom);
                        current_atom = None;
                    }}
                    ,
                '1'..='9' => {
                    if is_negative {
                        atoms[current_atom_index].charge = -(c.to_digit(10).unwrap() as i8);
                    } else {
                        atoms[current_atom_index].charge = c.to_digit(10).unwrap() as i8;
                    }
                },
                _ => (),
            }
        }

        
        
        Molecule::from_atoms(atoms)


    }
    
    pub fn atoms(&self) -> &Vec<Atom> {
        &self.atoms
    }

    pub fn atoms_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.atoms
    }

    pub fn atom_bonds_mut(&mut self, atom_index: usize) -> &mut Vec<Bond> {
        &mut self.atoms[atom_index].bonds
    }

    pub fn nth_atom(&self, index: usize) -> &Atom {
        &self.atoms[index]
    }

    pub fn get_charge(&self, atom_index: usize) -> i8 {
        self.atoms[atom_index].charge
    }

    pub fn atom_charge_mut(&mut self, atom_index: usize) -> &mut i8 {
        &mut self.atoms[atom_index].charge
    }

    pub fn is_atom_radical(&self, atom_index: usize) -> bool {
        self.atoms[atom_index].is_radical
    }

    pub fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom);
    }

    pub fn cmp_atom_charges(&self, atom_index1: usize, atom_index2: usize) -> std::cmp::Ordering {
        self.get_charge(atom_index1)
            .cmp(&self.get_charge(atom_index2))
    }

    pub fn number_of_bonded_element(&self, atom_index: usize, element: u8) -> usize {
        self.atoms[atom_index]
            .bonds()
            .iter()
            .filter(|bond| self.atoms[bond.target()].atomic_number == element)
            .count()
    }

    /// Creates a new molecule from an xyz file
    ///
    /// # Panics
    /// Function panics if the file cannot be opened or the file is empty
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::Molecule;
    ///
    /// let molecule = Molecule::from_xyz("tests/ethane.xyz");
    /// assert_eq!(molecule.atoms.len(), 8);
    /// ```
    pub fn from_xyz<P: AsRef<Path>>(filepath: P) -> Molecule {
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
        let mut molecule = Molecule::from_atoms(atoms);
        molecule.identify_bonds(BOND_SEARCH_THRESHOLD);
        molecule
    }

    pub fn build_tree(&self) -> KdTree<f64, 3> {
        let mut tree: KdTree<f64, 3> = KdTree::with_capacity(self.atoms.len());
        for (index, atom) in self.atoms.iter().enumerate() {
            tree.add(&atom.position_vector.as_array(), index as u64);
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

    pub fn identify_bond_changes(&self, other: &Self) -> Vec<BondChange> {
        let mut bond_changes: Vec<BondChange> = Vec::new();
        self.atoms
            .iter()
            .map(|atom| atom.bonds())
            .zip(other.atoms.iter().map(|atom| atom.bonds()))
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

    pub fn charge_changes(&self, other: &Molecule) -> Vec<(usize, i8, i8)> {
        self.atoms
            .iter()
            .zip(&other.atoms)
            .enumerate()
            .filter(|(_, (atom1, atom2))| atom1.charge() != atom2.charge())
            .map(|(index, (atom1, atom2))| (index, atom1.charge(), atom2.charge()))
            .collect()
    }

    pub fn update_atom_charge(&mut self, atom_index: usize, charge: i8) {
        self.atoms[atom_index].update_charge(charge);
    }

    pub fn identify_difference(&self, other: &Molecule) {
        let _charge_changes = self.charge_changes(other);
        let _bond_changes = self.identify_bond_changes(other);
        let _bond_type_changes: Vec<(usize, Vec<&Bond>, Vec<&Bond>)> = self
            .atoms
            .iter()
            .zip(&other.atoms)
            .enumerate()
            .filter(|(_, (atom1, atom2))| {
                let mut atom1_bonds = atom1.bonds.clone();
                let mut atom2_bonds = atom2.bonds.clone();
                atom1_bonds.sort_by_key(|a| a.target());
                atom2_bonds.sort_by_key(|a| a.target());
                atom1_bonds != atom2_bonds
            })
            .map(|(index, (atom1, atom2))| {
                let mut pre_transform = Vec::new();
                let mut post_transform = Vec::new();
                for bond1 in atom1.bonds() {
                    for bond2 in atom2.bonds() {
                        if bond1.target() == bond2.target()
                            && bond1.bond_type() != bond2.bond_type()
                        {
                            pre_transform.push(bond1);
                            post_transform.push(bond2);
                        }
                    }
                }
                (index, pre_transform, post_transform)
            })
            .collect();
    }

    pub fn identify_bonds(&mut self, threshold: f64) {
        let threshold_squared = threshold.powi(2);
        let kdtree = self.build_tree();

        let bonds = self
            .atoms
            .par_iter()
            .enumerate()
            .map(|(index, atom)| {
                let mut bonds = kdtree
                    .within::<SquaredEuclidean>(
                        &[
                            atom.position_vector.x,
                            atom.position_vector.y,
                            atom.position_vector.z,
                        ],
                        threshold_squared,
                    )
                    .iter()
                    .filter_map(|neighbor| {
                        if neighbor.item == index as u64 {
                            return None;
                        }

                        let distance = neighbor.distance;

                        let is_bonded = atom.is_bonded(
                            &self.atoms[neighbor.item as usize],
                            distance,
                            BOND_TOLERANCE,
                        );

                        if is_bonded {
                            Some(Bond {
                                target: neighbor.item as usize,
                                bond_type: BondType::Single,
                            })
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<Bond>>();
                bonds.sort_by_key(|bond| bond.target);
                bonds
            })
            .collect::<Vec<Vec<Bond>>>();
        self.atoms
            .par_iter_mut()
            .zip(bonds)
            .for_each(|(atom, bonds)| {
                atom.bonds = bonds;
            });
    }

    /// This function returns the bonds for a given atom
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::{Molecule,Bond,Atom};
    /// let mut molecule = Molecule::default();
    /// molecule.atoms.push(Atom::new(6, (0.0, 0.0, 0.0)));
    /// molecule.atoms.push(Atom::new(6, (1.0, 0.0, 0.0)));
    /// molecule.atoms.push(Atom::new(6, (2.0, 0.0, 0.0)));
    /// molecule.identify_bonds(1.5);
    /// assert_eq!(molecule.get_bonds(1).len(), 2);
    /// ```
    pub fn get_bonds(&self, atom_index: usize) -> &Vec<Bond> {
        &self.atoms[atom_index].bonds
    }

    pub fn get_edges(&self) -> Vec<(u32, u32)> {
        self.atoms()
            .iter()
            .enumerate()
            .flat_map(|(atom_index, atom)| {
                atom.bonds().iter().map(move |bond| {
                    if atom_index < bond.target() {
                        Some((atom_index as u32, bond.target() as u32))
                    } else {
                        None
                    }
                })
            })
            .flatten()
            .collect()
    }

    pub fn charge(&self) -> i32 {
        if self.atoms.len() > 1000 {
            self.atoms.par_iter().map(|atom| atom.charge() as i32).sum()
        } else {
            self.atoms.iter().map(|atom| atom.charge() as i32).sum()
        }
    }

    /// Returns a vector to the center of the Molecule
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::{Molecule,Atom};
    /// use molecules::vector::Vector;
    /// let mut molecule = Molecule::default();
    /// molecule.atoms.push(Atom::new(6, (0.0, 0.0, 0.0)));
    /// molecule.atoms.push(Atom::new(6, (1.0, 0.0, 0.0)));
    /// molecule.atoms.push(Atom::new(6, (2.0, 0.0, 0.0)));
    /// assert_eq!(molecule.center(), Vector::new(1.0, 0.0, 0.0));
    /// ```
    pub fn center(&self) -> Vector {
        let mut sum = Vector::new(0.0, 0.0, 0.0);

        for atom in &self.atoms {
            sum += atom.position_vector
        }

        let count = self.atoms.len() as f64;
        sum / count
    }

    pub fn number_of_charged_atoms(&self) -> usize {
        self.atoms.iter().filter(|atom| atom.charge() != 0).count()
    }

    pub fn find_angles(&self) -> Vec<BondAngle> {
        let bond_angles: Vec<BondAngle> = self
            .atoms
            .par_iter()
            .enumerate()
            .flat_map(|(atom_index, atom)| {
                let mut local_bond_angles = Vec::with_capacity(12);
                atom.bonds.iter().combinations(2).for_each(|slice| {
                    let neighbor1_index = slice[0].target;
                    let neighbor2_index = slice[1].target;
                    let angle = self.atoms[neighbor1_index]
                        .position_vector
                        .angle_between_points(
                            &self.atoms[atom_index].position_vector,
                            &self.atoms[neighbor2_index].position_vector,
                        );
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
    /// use molecules::molecule::Molecule;
    ///
    /// let mut molecule = Molecule::from_xyz("tests/ethane.xyz");
    /// let dihedrals = molecule.dihedrals();
    /// assert_eq!(dihedrals.len(), 9);
    ///
    ///```
    // PROTOTYPE
    pub fn dihedrals(&mut self) -> Vec<((usize, usize, usize, usize), f64)> {
        let bond_angles = self.find_angles();
        let dihedrals = bond_angles
            .par_iter()
            .flat_map(|angle| {
                let (atom1, atom2, atom3) = angle.atoms;

                self.get_bonds(atom1)
                    .iter()
                    .filter_map(|&bond| {
                        if bond.target() != atom3 && bond.target() != atom2 && atom3 < bond.target()
                        {
                            let dihedral = (bond.target(), atom1, atom2, atom3);
                            Some((dihedral, self.dihedral_angle(&(dihedral))))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        dihedrals
    }

    // pub fn find_dihedrals(&mut self) {
    //     if self.bond_angles.is_empty() {
    //         self.find_angles()
    //     }
    //     let dihedrals = self
    //         .bond_angles
    //         .iter()
    //         .flat_map(|angle| {
    //             let ((atom1, atom2, atom3), _) = *angle;
    //             self.get_bonds(atom3)
    //                 .iter()
    //                 .filter_map(|&bond| {
    //                     // This logic needs to be updated
    //                     if bond.target != atom2 && atom1 < bond.target {
    //                         let dihedral = (atom1, atom2, atom3, bond.target);
    //                         Some((dihedral, self.dihedral_angle(&(dihedral))))
    //                     } else {
    //                         None
    //                     }
    //                 })
    //                 .collect::<Vec<_>>()
    //         })
    //         .collect();
    //     self.dihedral_angles = dihedrals;
    // }

    // pub fn par_find_dihedrals(&mut self) {
    //     if self.bond_angles.is_empty() {
    //         self.find_angles()
    //     }
    //     if !self.bonds_calculated {
    //         self.identify_bonds(BOND_SEARCH_THRESHOLD)
    //     }

    //     let dihedrals = self
    //         .bond_angles
    //         .par_iter()
    //         .flat_map(|angle| {
    //             self.bonds
    //                 .iter()
    //                 .filter_map(|&(bond0, bond1)| {
    //                     let ((atom1, atom2, atom3), _) = *angle;
    //                     if bond0 == atom2 || bond1 == atom2 {
    //                         return None;
    //                     }
    //                     let dihedral = match (atom1, atom2) {
    //                         _ if atom1 == bond0 => (bond1, atom1, atom2, atom3),
    //                         _ if atom1 == bond1 => (bond0, atom1, atom2, atom3),
    //                         _ if atom3 == bond0 => (atom1, atom2, atom3, bond1),
    //                         _ if atom3 == bond1 => (atom1, atom2, atom3, bond0),
    //                         _ => return None,
    //                     };
    //                     if dihedral.0 < dihedral.3 {
    //                         Some((dihedral, self.dihedral_angle(&dihedral)))
    //                     } else {
    //                         None
    //                     }
    //                 })
    //                 .collect::<Vec<((usize, usize, usize, usize), f64)>>()
    //         })
    //         .collect::<Vec<((usize, usize, usize, usize), f64)>>();
    //     self.dihedral_angles = dihedrals;
    // }
    /// This function calculates the dihedral angle for all atoms
    ///
    /// # Panics
    /// This function does not panic.
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::molecule::Molecule;
    /// let mut molecule = Molecule::from_xyz("tests/ethane.xyz");
    /// let dihedral_angle = molecule.dihedral_angle(&(0, 1, 2, 3));
    /// assert_eq!(dihedral_angle, 0.5807210503904102);
    /// ```
    pub fn dihedral_angle(&self, dihedral: &(usize, usize, usize, usize)) -> f64 {
        let (a, b, c, d) = (
            self.atoms[dihedral.0].position_vector,
            self.atoms[dihedral.1].position_vector,
            self.atoms[dihedral.2].position_vector,
            self.atoms[dihedral.3].position_vector,
        );
        let v3 = d - c;
        let v2 = c - b;
        let v1 = b - a;

        let normal1 = v1.cross(&v2);
        let normal2 = v2.cross(&v3);

        let angle = normal1.angle_between(&normal2);
        let sign = normal1.cross(&normal2).dot(&v2);
        if sign < 0.0 {
            -angle
        } else {
            angle
        }
    }
    /// This function calculates the degrees of all atoms
    ///
    /// # Example
    /// ```
    /// use molecules::molecule::Molecule;
    /// let mut molecule = Molecule::from_xyz("tests/ethane.xyz");
    /// molecule.identify_bonds(2.0);
    /// let degrees = molecule.degrees();
    /// assert_eq!(degrees, vec![0, 0, 0, 0, 0, 0, 0, 0]);
    /// ```
    pub fn degrees(&self) -> Vec<i8> {
        self.atoms
            .iter()
            .map(|atom| atom.degree())
            .collect::<Vec<i8>>()
    }

    pub fn shift_charge(&mut self, index: usize, neighbor_index: usize, charge: i8) {
        self.atoms[index].update_charge(-charge);
        self.atoms[neighbor_index].update_charge(charge);
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
                if self.atoms()[index].atomic_number() == 6 && *degree > 0 {
                    Some((index, *degree))
                } else {
                    None
                }
            })
            .collect();

        oversaturated_carbons.iter().for_each(|(index, degree)| {
            let bonds = self.atoms[*index].bonds();
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

    pub fn split_components(&self) -> Vec<Molecule> {
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
            let mut atoms = Vec::with_capacity(component.len());
            for index in component {
                let mut atom = self.atoms[index].clone();
                atom.bonds
                    .iter_mut()
                    .for_each(|bond| bond.target = new_atom_indices[&bond.target]);
                atoms.push(self.atoms[index].clone());
            }

            if atoms.len() == 2 && atoms.iter().all(|atom| atom.atomic_number == 7) {
                continue;
            }
            let mut molecule = Molecule::from_atoms(atoms);
            molecule.identify_bonds(BOND_SEARCH_THRESHOLD);
            molecules.push(molecule);
        }
        molecules
    }

    pub fn get_components(&self) -> Vec<Vec<usize>> {
        let mut connected_components = vec![];
        let mut visited_atoms = vec![false; self.atoms.len()];
        let mut stack = Vec::with_capacity(self.atoms.len());
        while let Some(index) = visited_atoms.iter().position(|&a| !a) {
            let mut component = self.traverse_component(&mut stack, index, &mut visited_atoms);
            component.sort();
            connected_components.push(component);
            stack.clear();
        }
        connected_components
    }

    pub fn traverse_component(
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
            for &bond in self.get_bonds(index) {
                let target = bond.target();
                if !visited_atoms[target] {
                    stack.push(target);
                    visited_atoms[target] = true;
                }
            }
        }
        current_component
    }

    pub fn to_ungraph(&self) -> UnGraph<u8, ()> {
        let mut graph = UnGraph::<u8, ()>::default();
        let mut node_indices = IntMap::<usize, NodeIndex<u32>>::default();

        for (index, atom) in self.atoms.iter().enumerate() {
            let node_index = graph.add_node(atom.atomic_number());
            node_indices.insert(index, node_index);
        }

        for (index, atom) in self.atoms.iter().enumerate() {
            for bond in atom.bonds() {
                let source = node_indices[&index];
                let target = node_indices[&bond.target()];
                if source < target {
                    graph.add_edge(source, target, ());
                }
            }
        }
        graph
    }

    pub fn to_ungraph_from_slice(&self, slice: &[usize]) -> UnGraph<u8, ()> {
        let mut graph = UnGraph::<u8, ()>::default();
        let mut node_indices = IntMap::<usize, NodeIndex<u32>>::default();

        for &index in slice {
            let atom = &self.atoms[index];
            let node_index = graph.add_node(atom.atomic_number());
            node_indices.insert(index, node_index);
        }

        for &index in slice {
            let atom = &self.atoms[index];
            for bond in atom.bonds() {
                let source = node_indices[&index];
                let target = node_indices[&bond.target()];
                if source < target && slice.contains(&bond.target()) {
                    graph.add_edge(source, target, ());
                }
            }
        }
        graph
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
        let atom = &self.atoms[start_index];
        for &neighbor in atom.bonds() {
            if Some(neighbor.target()) == parent_index {
                continue;
            }

            if self.atoms[neighbor.target()].atomic_number == 1 {
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

    // PROTOTYPE
    pub fn to_cannonical_smiles(&self) -> String {
        let _primes = &PRIMES;

        let _smiles_parameters = self.atoms.iter().map(|atom| {
            let n_bonds = atom.bonds.len() as u8;
            let n_non_hydrogen = atom
                .bonds
                .iter()
                .filter(|bond| self.atoms[bond.target].atomic_number != 1)
                .count() as u8;
            let n_hydrogen = n_bonds - n_non_hydrogen;
            let atomic_number = atom.atomic_number;
            let sign_of_charge = atom.charge().signum() as u8;
            SmilesParameters {
                n_bonds,
                n_non_hydrogen,
                atomic_number,
                sign_of_charge,
                n_hydrogen,
            }
        });

        "".to_string()
    }
    pub fn monoisotopic_mass(&self) -> f64 {
        self.atoms
            .iter()
            .map(|atom| atom.monoisotopic_mass())
            .sum::<f64>()
    }

    pub fn molecular_formula(&self) -> MolecularFormula {
        MolecularFormula::from_molecule(self)
    }

    pub fn to_smiles(&self) -> String {
        let mut smiles = String::new();
        let components = self.get_components();

        for component in components {
            let ring_closures = &mut IntMap::default();
            let mut ring_counter = 1;
            let mut visited = vec![false; self.atoms.len()];
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

    pub fn format_atom(&self, atom_index: usize) -> String {
        let atom = match self.atoms.get(atom_index) {
            Some(atom) => atom,
            None => return "".to_string(),
        };
        let hydrogen_str = if atom.atomic_number() == 1 {
            "H".to_string()
        } else {
            "".to_string()
        };

        let number_of_hydrogens = self.number_of_bonded_element(atom_index, 1);

        let charge = atom.charge();

        match charge {
            2.. => format!("[{}{}+{}]", atom.name(), hydrogen_str, atom.charge.abs()),
            ..=-2 => format!("[{}{}-{}]", atom.name(), hydrogen_str, atom.charge.abs()),
            1 => format!("[{}{}+1]", atom.name(), hydrogen_str),
            -1 => format!("[{}{}-1]", atom.name(), hydrogen_str),
            0 => {
                if number_of_hydrogens > 0 {
                    format!("[{}{}]", atom.name(), hydrogen_str)
                } else {
                    format!("[{}]", atom.name())
                }
            }
        }
    }
    // This function works but it does not provide the desired output
    // pub fn to_canon_smiles(&self) -> String {
    //     use nauty_pet::prelude::*;
    //     let molecules = self.split_components();

    //     for molecule in &molecules {
    //         println!("Morgan\n{:?}\n", molecule.canonicalize());
    //         let graph = molecule.to_ungraph();
    //         println!("Graph\n {:?}\n", graph);

    //         let canon_labeling = graph.clone().into_canon();
    //         println!("Canon labeling\n{:?}\n", canon_labeling);

    //         println!("Is the same: {:?}", canon_labeling.is_identical(&graph));
    //     }

    //     println!("====================\n");
    //     String::new()
    // }

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
            BondType::Aromatic => smiles.push('~'),
        }
        let charge = self.atoms[node.index].charge;
        let number_of_hydrogens = self.number_of_bonded_element(node.index, 1);

        let hydrogen_str = if number_of_hydrogens > 0 {
            format!("H{}", number_of_hydrogens)
        } else {
            String::new()
        };

        let charge_string = match charge {
            2.. => format!(
                "[{}{}+{}]",
                self.atoms[node.index].name(),
                hydrogen_str,
                charge.abs()
            ),
            ..=-2 => format!(
                "[{}{}-{}]",
                self.atoms[node.index].name(),
                hydrogen_str,
                charge.abs()
            ),
            1 => format!("[{}{}+1]", self.atoms[node.index].name(), hydrogen_str),
            -1 => format!("[{}{}-1]", self.atoms[node.index].name(), hydrogen_str),
            0 => {
                if number_of_hydrogens > 0 {
                    format!("[{}{}]", self.atoms[node.index].name(), hydrogen_str)
                } else {
                    format!("[{}]", self.atoms[node.index].name())
                }
            }
        };
        smiles.push_str(&charge_string);

        if let Some(closures) = ring_closures.get(&node.index) {
            for closure in closures {
                smiles.push_str(&format!("{}{}", closure.1, closure.0));
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
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct MolecularFormula {
    elements: IntMap<u8, usize>,
}

impl MolecularFormula {
    pub fn new(elements: IntMap<u8, usize>) -> Self {
        Self { elements }
    }
    pub fn from_molecule(molecule: &Molecule) -> Self {
        let mut formula = MolecularFormula::default();
        for atom in &molecule.atoms {
            *formula.elements.entry(atom.atomic_number()).or_insert(0) += 1;
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
            acc + MONOISOTOPIC_MASSES[*atom as usize] * *count as f64
        })
    }

    /// Calculates the molecular mass of the molecular formula.
    /// # Examples
    /// ```
    /// use molecules::molecule::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// assert_eq!(water.molecular_mass(), 18.015);
    /// assert_eq!(methane.molecular_mass(), 16.043);
    pub fn molecular_mass(&self) -> f64 {
        self.elements.iter().fold(0.0, |acc, (atom, count)| {
            acc + STANDARD_ATOMIC_WEIGHTS[*atom as usize] * *count as f64
        })
    }
}

impl Display for MolecularFormula {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut elements = self
            .elements
            .iter()
            .map(|(&atom, &count)| (ATOMIC_SYMBOLS[atom as usize], count))
            .collect::<Vec<_>>();
        elements.sort_by_key(|&(atom, _)| atom);
        for (atom, count) in elements {
            write!(f, "{}{}", atom, count)?;
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
        s.chars().for_each(|c| {
            if c.is_alphabetic() {
                if !element_buffer.is_empty() && c.is_uppercase() {
                    let atomic_number = ATOMIC_NUMBERS[&element_buffer];
                    let count = if !number_buffer.is_empty() {
                        number_buffer.parse::<usize>().unwrap_or(1)
                    } else {
                        1
                    };
                    *elements.entry(atomic_number).or_insert(0) += count;
                    number_buffer.clear();
                    element_buffer.clear();
                }
                element_buffer.push(c);
            } else if c.is_numeric() {
                number_buffer.push(c);
            }
        });
        if !element_buffer.is_empty() {
            let atomic_number = ATOMIC_NUMBERS[&element_buffer];
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
struct Node {
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

fn backtrack_bonding(molecule: &mut Molecule, degrees: &mut [i8]) -> bool {
    // We can include the expected charge as an argument here that we have derived from the
    // simulation summary file

    if degrees.iter().all(|&degree| degree >= 0) {
        return true;
    }

    // Find the first unsatisfied atom (negative degree)
    let mut unsatisfied_atom_index = degrees.iter().position(|&degree| degree < 0).unwrap();

    let number_of_double_bonds = molecule.atoms[unsatisfied_atom_index]
        .bonds
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

    for &neighbor in &molecule.atoms[unsatisfied_atom_index].bonds.clone() {
        let number_of_neighbor_atom_double_bonds = molecule.atoms[neighbor.target]
            .bonds
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

pub fn relaxed_backtrack_bonding(molecule: &mut Molecule, degrees: &mut [i8]) -> bool {
    // If all degrees are zero or only one degree is negative (positive charge) then we are done
    // TODO check what is going on here
    // We can include the expected charge as an argument here that we have derived from the
    // simulation summary file

    if degrees.iter().all(|&degree| degree >= 0) {
        println!("All degrees are positive in relaxed backtrack bonding");
        return true;
    }

    // Find the first unsatisfied atom (negative degree)
    let unsatisfied_atom_index = degrees.iter().position(|&degree| degree < 0).unwrap();

    for &neighbor in &molecule.atoms[unsatisfied_atom_index].bonds.clone() {
        let number_of_neighbor_atom_double_bonds = molecule.atoms[neighbor.target]
            .bonds
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
    //  let mut stack = Vec::with_capacity(degrees.len());
    //  let mut backtrack_stack = Vec::with_capacity(degrees.len());

    //  for degree_index in 0..degrees.len() {
    //      if degrees[degree_index] < 0 {
    //          stack.push(degree_index);
    //      }
    //      while let Some(index) = stack.pop() {
    //          if degrees.iter().filter(|&degree| *degree < 0).sum::<i8>() == 0   {
    //              return true;
    //          }

    //          for neighbor in molecule.atoms[index].bonds.clone() {
    //              if can_increase_bond(index, neighbor.target(), degrees)
    //              {
    //                  // #[cfg(debug_assertions)]
    //                  increase_bonds(molecule, index, neighbor.target());
    //                  update_degrees(degrees, index, neighbor.target(), false);
    //                  backtrack_stack.push((index, neighbor.target()));
    //                  stack.push(degrees.iter().position(|&degree| degree < 0).unwrap());
    //                  break;
    //              }
    //          }
    //          if let Some((index1, index2)) = backtrack_stack.pop() {
    //              decrease_bonds(molecule, index1, index2);
    //              update_degrees(degrees, index1, index2, true);
    //          }
    //      }
    //  }
    //  false
}

pub fn can_increase_bond(atom_index: usize, neighbor_index: usize, degrees: &[i8]) -> bool {
    degrees[atom_index] < 0 && degrees[neighbor_index] < 0
}

pub fn increase_bonds(molecule: &mut Molecule, atom_index: usize, neighbor_index: usize) {
    for bond in &mut molecule.atoms[atom_index].bonds {
        if bond.target == neighbor_index {
            increase_bond(bond);
        }
    }
    for bond in &mut molecule.atoms[neighbor_index].bonds {
        if bond.target == atom_index {
            increase_bond(bond);
        }
    }
}

pub fn increase_bond(bond: &mut Bond) {
    match bond.bond_type {
        BondType::Single => bond.bond_type = BondType::Double,
        BondType::Double => bond.bond_type = BondType::Triple,
        BondType::Triple => panic!("Cannot increase bond beyond triple bond"),
        BondType::Aromatic => panic!("Cannot increase bond beyond aromatic bond"),
    }
}

pub fn decrease_bonds(molecule: &mut Molecule, atom_index: usize, neighbor_index: usize) {
    for bond in &mut molecule.atoms[atom_index].bonds {
        if bond.target == neighbor_index {
            decrease_bond(bond);
        }
    }
    for bond in &mut molecule.atoms[neighbor_index].bonds {
        if bond.target == atom_index {
            decrease_bond(bond);
        }
    }
}

pub fn decrease_bond(bond: &mut Bond) {
    match bond.bond_type {
        BondType::Single => panic!("Cannot decrease bond beyond single bond"),
        BondType::Double => bond.bond_type = BondType::Single,
        BondType::Triple => bond.bond_type = BondType::Double,
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
/// assert_eq!(atom.position_vector.x, 11.390, "Incorrect x coordinate");
/// assert_eq!(atom.position_vector.y, -11.172, "Incorrect y coordinate");
/// assert_eq!(atom.position_vector.z, 71.797, "Incorrect z coordinate");
/// assert_eq!(atom.name(), "C", "Incorrect atom name");
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
    let atomic_number = *ATOMIC_NUMBERS.get(&name).unwrap_or(&0);
    Ok(Atom {
        position_vector: position,
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
/// use molecules::molecule::{Atom,extract_atom_cif};
///
/// let line = "ATOM   1    N N   . GLN A 1 1   ? 201.310 198.892 131.429 1.00 70.25  ? 1   GLN A N   1";
/// let atom = extract_atom_cif(line).unwrap();
/// assert_eq!(atom.position_vector.x, 201.310, "Incorrect x coordinate");
/// assert_eq!(atom.position_vector.y, 198.892, "Incorrect y coordinate");
/// assert_eq!(atom.position_vector.z, 131.429, "Incorrect z coordinate");
/// assert_eq!(atom.name(), "N", "Incorrect atom name");
/// assert_eq!(atom.atomic_number(), 7, "Incorrect atom type");
/// ```
pub fn extract_atom_cif(line: &str) -> Result<Atom, ParseFloatError> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    let atomic_number = ATOMIC_NUMBERS[fields[2]];

    let position = Vector {
        x: fields[10].parse::<f64>()?,
        y: fields[11].parse::<f64>()?,
        z: fields[12].parse::<f64>()?,
    };
    Ok(Atom {
        position_vector: position,
        atomic_number,
        ..Default::default()
    })
}

// Unit tests (Still needs to be improved)
#[cfg(test)]
mod tests {
    #[test]
    fn test_extract_atom() {
        let line =
            "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N  ";
        let atom = super::extract_atom_pdb(line).unwrap();
        assert_eq!(atom.atomic_number, 7);
        assert_eq!(atom.position_vector.x, 10.0);
        assert_eq!(atom.position_vector.y, 10.0);
        assert_eq!(atom.position_vector.z, 10.0);
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
        let atom1 = Atom {
            position_vector: Vector {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },

            ..Default::default()
        };
        let atom2 = Atom {
            position_vector: Vector {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
            ..Default::default()
        };
        let atom3 = Atom {
            position_vector: Vector {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            atomic_number: 6,
            ..Default::default()
        };
        let mut molecule = Molecule::from_atoms(vec![atom1, atom2, atom3]);
        molecule.identify_bonds(2.0);
        let angles = molecule.find_angles();
        assert_eq!(angles.len(), 1);
        assert_eq!(angles[0].angle, std::f64::consts::FRAC_PI_2);
    }
}

struct GMLReactionRule {
    before: Molecule,
    context: Molecule,
    after: Molecule,
}

impl Molecule {
    /// Finds the shortest path of a set of atoms
    ///
    /// # Examples
    /// ```
    /// use molecules::molecule::Molecule;
    ///
    /// let mut molecule = Molecule::from_xyz("tests/ethane.xyz");
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
    pub fn find_shortest_path_of_set(&self, atom_set: &HashSet<usize>) -> Option<Vec<usize>> {
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
    /// use molecules::molecule::Molecule;
    ///
    /// let mut molecule = Molecule::from_xyz("tests/ethane.xyz");
    /// let shortest_path = molecule.shortest_path(0,1).unwrap();
    /// assert_eq!(shortest_path, vec![0,1]);
    /// ```
    ///
    pub fn shortest_path(&self, start: usize, end: usize) -> Option<Vec<usize>> {
        let mut queue = VecDeque::new();
        let mut visited = vec![false; self.atoms.len()];
        let mut parent = vec![None; self.atoms.len()];

        visited[start] = true;
        queue.push_back(start);

        while let Some(current_node) = queue.pop_front() {
            if current_node == end {
                return Some(reconstruct_path(current_node, &parent));
            }
            for &neighbor in self.atoms[current_node].bonds() {
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
        let mut parent = vec![None; self.atoms.len()];

        queue.push_back(start);

        while let Some(current_node) = queue.pop_front() {
            if current_node == end {
                return Some(reconstruct_path(current_node, &parent));
            }

            for &neighbor in self.atoms[current_node].bonds() {
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
            if let Some(atom) = self.atoms().get(current_node) {
                let bonds = atom.bonds();
                // Extract necessary information about neighbors
                let mut neighbors_info: Vec<(usize, u8, BondType)> = bonds
                    .iter()
                    .map(|bond| {
                        let target = bond.target();
                        let atomic_number = self
                            .atoms()
                            .get(target)
                            .map_or(0, |atom| atom.atomic_number());
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

    pub fn charges(&self) -> Vec<i8> {
        self.atoms().iter().map(|atom| atom.charge()).collect()
    }

    pub fn names(&self) -> Vec<&str> {
        self.atoms().iter().map(|atom| atom.name()).collect()
    }

    
    pub fn morgans_algorithm(&self, max_depth: Option<usize>) -> Vec<usize> {
        let mut vertex_degrees = self
            .atoms()
            .iter()
            .map(|atom| atom.bonds().len())
            .collect::<Vec<_>>();
        let mut buffer = vec![0; self.atoms.len()];
        let mut last_count = 0;
        let mut depth = 0;
        let max_depth = max_depth.unwrap_or(100);

        loop {
            depth += 1;
            if depth > max_depth {
                println!("Morgan's algorithm did not converge after 100 iterations");
                break;
            }
            println!("{:?}", vertex_degrees);
            for (index, atom) in self.atoms().iter().enumerate() {
                let sum = atom
                    .bonds()
                    .iter()
                    .map(|bond| vertex_degrees[bond.target()])
                    .sum();
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
