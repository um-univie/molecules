use std::path::Path;

use crate::{
    atom::Atom,
    molecule::bond::BondAngle,
    molecule::{base::Molecule, molecule3d::Molecule3D},
    vector::Vector,
};
use chemistry_consts::ElementProperties;
use std::fs::File;
use std::io::{self, BufRead};
use std::num::ParseFloatError;

use rayon::prelude::*;

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

    pub fn identify_bonds(&mut self, tolerance: f64) {
        self.molecules
            .iter_mut()
            .for_each(|molecule| molecule.identify_bonds(tolerance))
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

/// This function reads a pdb file line and extracts the atom information
///
/// # Panics
/// This function returns a ParseFloatError if the line has errors in the float region, it does not check for other errors.
///
/// ```
/// use molecules::molecule::molecular_system::extract_atom_pdb;
/// 
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
/// use molecules::molecule::molecular_system::extract_atom_cif;
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
