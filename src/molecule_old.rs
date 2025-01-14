struct SmilesParameters {
    n_bonds: u8,
    n_non_hydrogen: u8,
    atomic_number: u8,
    sign_of_charge: u8,
    n_hydrogen: u8,
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
    from: BondOrder,
    to: BondOrder,
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

