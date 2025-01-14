use core::num;
use std::os;

use crate::molecule::base::{AtomLabel, AtomLabelImplicitHydrogens, Molecule};
use crate::prelude::*;
use canonizer::canonization::{first_non_trivial, graph_canon};

pub fn mol_canonization<T: Molecule>(mol: &T) -> T where T: Clone {
    let reorder_labeling: Vec<usize>;

    if mol.isotopes().is_none() {
        let mut sorted_mol = sort_mol_by_atom_labels_implicit_hydrogens(mol);
        let graph = sorted_mol.to_ungraph_with_edge_nodes_and_implicit_hydrogens();

        //println!("Graph: {:#?}", graph);
        let (hash_labeling, automorphisms, trace_impact) = graph_canon(&graph, first_non_trivial, true);

        let number_of_hydrogens = sorted_mol.atomic_numbers().iter().filter(|&&x| x == 1).count();

        let mut labeling = vec![0;sorted_mol.atomic_numbers().len()];

        for (key, value) in hash_labeling {
            if key < (sorted_mol.atomic_numbers().len() - number_of_hydrogens) {
                labeling[value+number_of_hydrogens] = key + number_of_hydrogens
            }
        }



        for (index, _atomic_number) in sorted_mol.atomic_numbers().iter().enumerate().filter(|(_index, atomic_number)| **atomic_number == 1) {
            labeling[index] = index;
        }

        reorder_labeling = labeling;

        reorder_molecule(&mut sorted_mol, &reorder_labeling);
        sorted_mol
    } else {
        let mut sorted_mol = sort_mol_by_atom_labels(mol);
        let graph = sorted_mol.to_ungraph_with_edge_nodes();

        let (hash_labeling, automorphisms, trace_impact) = graph_canon(&graph, first_non_trivial, true);
    
        let mut labeling = vec![0;mol.atomic_numbers().len()];

        for (key, value) in hash_labeling {
            if key < mol.atomic_numbers().len() {
                labeling[value] = key
            }
        }


        reorder_labeling = labeling;

        reorder_molecule(&mut sorted_mol, &reorder_labeling);

        sorted_mol
    }

}

/// Converts a canonical labeling from a graph with edge nodes back to the original molecule labeling
/// by removing the edge node labels.
///
/// # Arguments
/// * `labeling` - The canonical labeling from the graph with edge nodes
/// * `num_atoms` - The number of atoms in the original molecule
///
/// # Returns
/// A vector containing the canonical labeling for just the atom nodes
///
/// # Example
/// ```
/// use molecules::graph_algorithms::mol_canonization::contract_edge_nodes;
/// let edge_node_labeling = vec![0, 4, 1, 5, 2, 6, 3];  // labeling with edge nodes
/// let num_atoms = 4;
/// let atom_labeling = contract_edge_nodes(&edge_node_labeling, num_atoms);
/// assert_eq!(atom_labeling, vec![0, 1, 2, 3]);  // just atom nodes
/// ```
pub fn contract_edge_nodes(labeling: &[usize], num_atoms: usize) -> Vec<usize> {
    labeling
        .iter()
        .filter(|&&x| x < num_atoms)
        .cloned()
        .collect()
}

pub fn contract_edge_nodes_in_place(labeling: &mut Vec<usize>, num_atoms: usize) {
    let mut write_idx = 0;
    for read_idx in 0..labeling.len() {
        if labeling[read_idx] < num_atoms {
            labeling.swap(write_idx, read_idx);
            write_idx += 1;
        }
    }
    labeling.truncate(write_idx);
}

/// Generates a canonical SMILES string for the molecule by first canonizing the atoms
/// and then using the canonical ordering to generate the SMILES string.
///
/// # Arguments
/// * `mol` - The molecule to generate a canonical SMILES string for
///
/// # Returns
/// A canonical SMILES string representation of the molecule
///
/// # Example
/// ```
/// use molecules::prelude::*;
/// use molecules::graph_algorithms::mol_canonization::to_canonical_smiles;
/// let mols = Molecule2D::from_smiles("CC=CC.C(C=CC)").unwrap();
/// let canonical_smiles = to_canonical_smiles(&mols[0]);
/// let canonical_smiles2 = to_canonical_smiles(&mols[1]);
///
/// assert_eq!(canonical_smiles, "C(=CC)C");  
/// assert_eq!(canonical_smiles2, "C(=CC)C");  
/// ```
pub fn to_canonical_smiles<T: Molecule + Clone>(mol: &T) -> String {
    let canonical_mol = mol_canonization(mol);

    //let inverse_relabeling = relabeling.iter().map(|&x| canonical_labeling.iter().position(|&y| y == x).unwrap()).collect::<Vec<usize>>();

    // Create a new molecule with reordered atoms according to canonical labeling
    //println!("Reordered atomic numbers: {:?}", reordered_mol.atom_bonds().iter().map(|bonds| bonds.len()).collect::<Vec<usize>>());
    

    // Generate SMILES from reordered molecule
    canonical_mol.to_smiles_with_implicit_hydrogens()
}

fn sort_mol_by_atom_labels<T: Molecule + Clone>(mol: &T) -> T {
    // Get canonical labeling
    let mut atom_indices_and_labels: Vec<(usize, AtomLabel, usize)> = (0..mol.atomic_numbers().len()).map(|index| { 
            let label = mol.get_atom_label(index);    
            let number_of_bonds = mol.get_atom_bonds(index).iter().count();
            (index,label,number_of_bonds)
            }
    ).collect();
    atom_indices_and_labels.sort_by_key(|(_index, label, number_of_bonds)| (*label, *number_of_bonds) );

    let relabeling = atom_indices_and_labels
        .into_iter()
        .map(|(index,_label,_number_of_bonds)| index)
        .collect::<Vec<usize>>();

    let mut reordered_mol = mol.clone();
    reorder_molecule(&mut reordered_mol, &relabeling);
    reordered_mol
}
fn sort_mol_by_atom_labels_implicit_hydrogens<T: Molecule + Clone>(mol: &T) -> T {
    // Get canonical labeling
    let mut atom_indices_and_labels: Vec<(usize, AtomLabelImplicitHydrogens)> = (0..mol.atomic_numbers().len()).map(|index| { 
            let label = mol.get_atom_label_with_hydrogens(index);    
            (index,label)
            }
    ).collect();

    atom_indices_and_labels.sort_by_key(|(_index, label)| *label );

    let relabeling = atom_indices_and_labels
        .into_iter()
        .map(|(index,_label)| index)
        .collect::<Vec<usize>>();

    let mut reordered_mol = mol.clone();
    reorder_molecule(&mut reordered_mol, &relabeling);
    reordered_mol
}



fn canonicalize_mol<T: Molecule + Clone>(mol: &T) -> T {

   
    let canonical_mol = mol_canonization(mol);
    
    // Create a new molecule with reordered atoms according to canonical labeling
    canonical_mol
}

/// Reorders the atoms in a molecule according to a given labeling
fn reorder_molecule<T: Molecule>(mol: &mut T, labeling: &[usize]) {
    // Create inverse mapping
    let mut inverse_mapping = vec![0; labeling.len()];
    for (new_idx, &old_idx) in labeling.iter().enumerate() {
        inverse_mapping[old_idx] = new_idx;
    }

    // Reorder atomic numbers
    let mut new_atomic_numbers = vec![0; mol.atomic_numbers().len()];
    for (new_idx, &old_idx) in labeling.iter().enumerate() {
        new_atomic_numbers[new_idx] = mol.atomic_numbers()[old_idx];
    }
    *mol.atomic_numbers_mut() = new_atomic_numbers;

    // Reorder charges
    let mut new_charges = vec![0; mol.atomic_numbers().len()];
    for (new_idx, &old_idx) in labeling.iter().enumerate() {
        new_charges[new_idx] = mol.charges().get(old_idx).copied().unwrap_or(0);
    }
    *mol.charges_mut() = new_charges;

    // Reorder bonds and update bond targets
    let mut new_bonds = vec![mol.atom_bonds()[0]; mol.atom_bonds().len()];
    for (new_idx, &old_idx) in labeling.iter().enumerate() {
        let atom_bonds = mol.atom_bonds();
        let mut bonds = atom_bonds[old_idx];
        // Update bond targets according to new ordering
        for bond in &mut bonds {
            bond.target = inverse_mapping[bond.target()];
        }
        new_bonds[new_idx] = bonds;
    }
    *mol.atom_bonds_mut() = new_bonds;
    mol.atom_bonds_mut()
        .iter_mut()
        .for_each(|bonds| bonds.sort_by_key(|bond| bond.target()));
}

/// Compares two strings and returns a formatted comparison using ANSI color codes.
/// Matching characters are shown in green, differences in red.
///
/// # Arguments
/// * `s1` - First string to compare
/// * `s2` - Second string to compare
///
/// # Returns
/// String containing a formatted comparison with colored differences
///
/// # Example
/// ```
/// use molecules::graph_algorithms::mol_canonization::compare_strings_formatted;
/// let comparison = compare_strings_formatted("hello", "hallo");
/// // Will print both strings with the 'e' and 'a' in red, other chars in green
/// ```
pub fn compare_strings_formatted(s1: &str, s2: &str) -> String {
    const GREEN: &str = "\x1b[32m";
    const RED: &str = "\x1b[31m";
    const RESET: &str = "\x1b[0m";
    
    let mut result = String::new();
    let s1_chars: Vec<char> = s1.chars().collect();
    let s2_chars: Vec<char> = s2.chars().collect();
    let max_len = s1.len().max(s2.len());
    
    // Format first string with colors
    for i in 0..max_len {
        match (s1_chars.get(i), s2_chars.get(i)) {
            (Some(&c1), Some(&c2)) if c1 == c2 => {
                result.push_str(GREEN);
                result.push(c1);
                result.push_str(RESET);
            }
            (Some(&c1), _) => {
                result.push_str(RED);
                result.push(c1);
                result.push_str(RESET);
            }
            (None, _) => break,
        }
    }
    result.push('\n');
    
    // Format second string with colors
    for i in 0..max_len {
        match (s1_chars.get(i), s2_chars.get(i)) {
            (Some(&c1), Some(&c2)) if c1 == c2 => {
                result.push_str(GREEN);
                result.push(c2);
                result.push_str(RESET);
            }
            (_, Some(&c2)) => {
                result.push_str(RED);
                result.push(c2);
                result.push_str(RESET);
            }
            (_, None) => break,
        }
    }
    
    result
}

/// Canonicalizes and sorts a SMILES string that may contain multiple molecules separated by dots.
///
/// # Arguments
/// * `smiles` - A SMILES string potentially containing multiple molecules separated by dots
///
/// # Returns
/// A canonical SMILES string with the individual molecules sorted and recombined
///
/// # Example
/// ```
/// use molecules::graph_algorithms::mol_canonization::canonicalize_and_sort_smiles;
/// let smiles = "CC=CC.C(C=CC)";
/// let canonical_smiles = canonicalize_and_sort_smiles(smiles);
/// assert_eq!(canonical_smiles, "C(=CC)C.C(=CC)C");
/// ```
pub fn canonicalize_and_sort_smiles(smiles: &str) -> String {
    let mut mols = Molecule2D::from_smiles_with_sanitization(smiles).unwrap();
    

    // Aromatize molecules
    mols.iter_mut().for_each(|mol| mol.simple_aromatization());

    // Canonicalize each molecule
    let mut canonical_smiles: Vec<String> = mols.iter()
        .map(to_canonical_smiles)
        .collect();

    // Sort the canonical SMILES strings
    canonical_smiles.sort();

    // Recombine the sorted canonical SMILES strings
    canonical_smiles.join(".")
}



 #[cfg(test)]
 mod tests {
     use crate::prelude::*;
     use std::fs::OpenOptions;
     use std::io::Write;
     use crate::graph_algorithms::mol_canonization::{to_canonical_smiles, mol_canonization, compare_strings_formatted, canonicalize_and_sort_smiles, canonicalize_mol};
     use rayon::prelude::*;
     use std::sync::Mutex;
     use std::fs::File;
     use std::io::{BufRead, BufReader};
     use std::time::Instant;

     #[test]
     fn test_mol_canonization() {
         let mol = Molecule3D::from_sdf("tests/cis_2_butene.sdf").unwrap()[0].clone();
         let labeling = mol_canonization(&mol);
         println!("{:?}", labeling);
     }

     #[test]
     fn test_canonical_smiles() {
         let mol = Molecule3D::from_sdf("tests/cis_2_butene.sdf").unwrap()[0].clone();
         let canonical_smiles = to_canonical_smiles(&mol);

         // Test that different input orderings give same canonical SMILES
         let mol2 = Molecule2D::from_smiles(&mol.to_smiles()).unwrap();
         let mol2 = &mol2[0];

         let canonical_smiles2 = to_canonical_smiles(mol2);
         assert_eq!(canonical_smiles, canonical_smiles2);

         let mol2 = Molecule2D::from_smiles("c12c3c(ccc2)c(ccc4)c5c4cccc5c3ccc1").unwrap();
         let mol2 = &mol2[0];

         let canonical_smiles2 = to_canonical_smiles(mol2);


         assert_eq!(canonical_smiles2, "c1c2c3c(c4c5c2cccc5ccc4)cccc3cc1");
     }


#[test]
fn benchmark_canonization_from_file() {
    let smiles_file_path = "tests/smiles_tree_permutations_filtered_c1(=N)ccccc1=C,N.txt";

    // Open the SMILES file
    let file = File::open(smiles_file_path)
        .unwrap_or_else(|_| panic!("Failed to open SMILES file: {}", smiles_file_path));
    let reader = BufReader::new(file);
    let lines = reader.lines().map_while(Result::ok).collect::<Vec<_>>();
    //let mut wrong_smiles_file = Mutex::new(File::create("tests/wrong_smiles.txt").unwrap());

    // Record the start time
    let start_time = Instant::now();

    // Process each line in parallel and accumulate counts
    let (total_smiles, correct_smiles) = lines
        .par_iter()
        .map(|line| {
            let mut iterator = line.split_whitespace();

            // Extract the comparison SMILES string
            let comparison_smiles = match iterator.next() {
                Some(smiles) => smiles,
                None => return (0, 0),
            };

            // Canonicalize the comparison SMILES
            let canonical_comparison_smiles = canonicalize_and_sort_smiles(comparison_smiles);

            // Assuming `canonicalize_mol` and related functions are defined elsewhere
            // Initialize local counters
            let mut local_total = 0;
            let mut local_correct = 0;

            // Iterate through each SMILES in the line
            for smiles in iterator {
                local_total += 1;

                // Canonicalize the current SMILES
                let canonical_smiles = canonicalize_and_sort_smiles(smiles);

                // Update the correct SMILES counter if canonicalizations match
                if canonical_comparison_smiles == canonical_smiles {
                    local_correct += 1;
                } else {
                    //wrong_smiles_file.lock().unwrap().write_all(format!("{} vs {}\n", canonical_comparison_smiles, canonical_smiles).as_bytes()).unwrap();
                }
            }

            (local_total, local_correct)
        })
        .fold(
            || (0, 0),
            |acc, counts| (acc.0 + counts.0, acc.1 + counts.1)
        )
        .reduce(
            || (0, 0),
            |a, b| (a.0 + b.0, a.1 + b.1)
        );

    // Record the end time
    let duration = start_time.elapsed();

    // Calculate success rate
    let success_rate = if total_smiles > 0 {
        (correct_smiles as f64 / total_smiles as f64) * 100.0
    } else {
        0.0
    };

    // Print benchmark results
    println!("=== Benchmark Results ===");
    println!("Total SMILES processed: {}", total_smiles);
    println!("Correct SMILES: {}", correct_smiles);
    println!("Success rate: {:.2}%", success_rate);
    println!("Total time for canonization: {:?}", duration);
}

     #[test]
     fn test_canonization_from_file() {
         use std::io::BufRead;
         use std::sync::atomic::{AtomicUsize, Ordering};
         
         // Initialize the log file with append mode
         let log_file = OpenOptions::new()
             .create(true)
             .append(true)
             .open("tests/canonization_log.txt")
             .expect("Failed to open log file");
         
         let log_file = Mutex::new(log_file);
         
         let file = File::open("tests/smiles_tree_permutations_filtered_c1(=N)ccccc1=C,N.txt")
             .expect("Failed to open SMILES file");
         let reader = BufReader::new(file);
         
         // Atomic counters for thread-safe counting
         let total_smiles = AtomicUsize::new(0);
         let correct_smiles = AtomicUsize::new(0);


         reader.lines().enumerate().par_bridge().for_each(|(line_counter, line)| {
             let line = line.expect("Failed to read line");
             let mut iterator = line.split_whitespace(); 
             let comparison_smiles = match iterator.next() {
                 Some(smiles) => smiles,
                 None => return,
             };

             let mut mols = Molecule2D::from_smiles_with_sanitization(comparison_smiles).unwrap();
             let mol = &mut mols[0];
             mol.simple_aromatization();

             // Create a debug info collector
             let mut debug_info: Vec<String> = Vec::new();
             let written_smiles = mol.to_smiles();

             let comparison_canonical_smiles = canonicalize_and_sort_smiles(comparison_smiles);
             let comparison_canonical_mol = canonicalize_mol(mol);

             let comparison_total_charge: i8 = comparison_canonical_mol.charges().iter().sum();
             let comparison_number_of_atoms = mol.atomic_numbers().len();
             let comparison_number_of_hydrogens = mol.atomic_numbers().iter().filter(|&&num| num == 1).count();
             let comparison_number_of_single_bonds = mol.atom_bonds().iter().filter(|bonds| bonds.iter().all(|bond| bond.bond_order() == BondOrder::Single)).count();
             let comparison_number_of_double_bonds = mol.atom_bonds().iter().filter(|bonds| bonds.iter().any(|bond| bond.bond_order() == BondOrder::Double)).count();
             let comparison_number_of_aromatic_bonds = mol.atom_bonds().iter().filter(|bonds| bonds.iter().any(|bond| bond.bond_order() == BondOrder::Aromatic)).count();
             let comparison_number_of_rings = mol.find_rings().len();
             let comparison_bond_pattern = comparison_canonical_mol.atom_bonds()
                 .iter()
                 .map(|bonds| bonds.iter().map(|bond| bond.bond_order()).collect::<Vec<BondOrder>>())
                 .collect::<Vec<Vec<BondOrder>>>();

             for smiles in iterator {
                 total_smiles.fetch_add(1, Ordering::Relaxed);
                 let canonical_smiles = canonicalize_and_sort_smiles(smiles);
                 let mut mols = Molecule2D::from_smiles_with_sanitization(smiles).unwrap();
                 let mol = &mut mols[0];
                 mol.simple_aromatization();
                 let canonical_mol = canonicalize_mol(mol);

                 let number_of_atoms = canonical_mol.atomic_numbers().len();
                 let number_of_hydrogens = canonical_mol.atomic_numbers().iter().filter(|&&num| num == 1).count();
                 let number_of_single_bonds = canonical_mol.atom_bonds().iter().filter(|bonds| bonds.iter().all(|bond| bond.bond_order() == BondOrder::Single)).count();
                 let number_of_double_bonds = canonical_mol.atom_bonds().iter().filter(|bonds| bonds.iter().any(|bond| bond.bond_order() == BondOrder::Double)).count();
                 let number_of_aromatic_bonds = canonical_mol.atom_bonds().iter().filter(|bonds| bonds.iter().any(|bond| bond.bond_order() == BondOrder::Aromatic)).count();
                 let bond_pattern = canonical_mol.atom_bonds()
                     .iter()
                     .map(|bonds| bonds.iter().map(|bond| bond.bond_order()).collect::<Vec<BondOrder>>())
                     .collect::<Vec<Vec<BondOrder>>>();

                 if number_of_hydrogens != comparison_number_of_hydrogens {
                     debug_info.push(format!(
                         "Hydrogens mismatch: got {} expected {}",
                         number_of_hydrogens, comparison_number_of_hydrogens
                     ));
                 }
                 if number_of_aromatic_bonds != comparison_number_of_aromatic_bonds {
                     debug_info.push(format!(
                         "Aromatic bonds mismatch: got {} expected {}",
                         number_of_aromatic_bonds, comparison_number_of_aromatic_bonds
                     ));
                 }
                 if number_of_single_bonds != comparison_number_of_single_bonds {
                     debug_info.push(format!(
                         "Single bonds mismatch: got {} expected {}",
                         number_of_single_bonds, comparison_number_of_single_bonds
                     ));
                 }
                 if number_of_double_bonds != comparison_number_of_double_bonds {
                     debug_info.push(format!(
                         "Double bonds mismatch: got {} expected {}",
                         number_of_double_bonds, comparison_number_of_double_bonds
                     ));
                 }
                 if number_of_atoms != comparison_number_of_atoms {
                     debug_info.push(format!(
                         "Atoms mismatch: got {} expected {}",
                         number_of_atoms, comparison_number_of_atoms
                     ));
                 }
                 if bond_pattern != comparison_bond_pattern {
                     let mismatches = bond_pattern.iter().zip(comparison_bond_pattern.iter())
                         .enumerate()
                         .filter(|(_, (got, expected))| got != expected)
                         .map(|(idx, (got, expected))| {
                             format!(
                                 "Atom {} ({}): {:?}, {:?}",
                                 idx,
                                 canonical_mol.get_atomic_symbol(idx).unwrap_or("unknown"),
                                 got,
                                 expected
                             )
                         })
                         .collect::<Vec<_>>();
                     
                     debug_info.push(format!(
                         "Bond pattern mismatches:\n{}",
                         mismatches.join("\n")
                     ));
                 }
                 // Check stereochemistry if relevant
                 if canonical_mol.has_stereochemistry() {
                     for (chiral_class, comparison_chiral_class) in canonical_mol.chiral_classes().unwrap()
                         .iter()
                         .zip(comparison_canonical_mol.chiral_classes().unwrap())
                     {
                         if chiral_class != comparison_chiral_class {
                             debug_info.push(format!(
                                 "Stereochemistry differs for SMILES {}",
                                 canonical_mol.to_smiles()
                             ));
                         }
                     }
                 }
                 if comparison_canonical_smiles == canonical_smiles {
                     correct_smiles.fetch_add(1, Ordering::Relaxed);
                 } else {
                     debug_info.push("Checking SMILES reader and writer".to_string());
                     debug_info.push(
                         compare_strings_formatted(comparison_smiles, &written_smiles)
                     .to_string());
                     debug_info.push(format!("Line number: {}", line_counter));
                     debug_info.push(format!("Original SMILES: {}", comparison_smiles));
                     debug_info.push(format!(
                         "Canonical comparison: {}",
                         comparison_canonical_smiles
                     ));
                     
                     let mut mols = Molecule2D::from_smiles_with_sanitization(smiles).unwrap();
                     let mol = &mut mols[0];
                     
                     debug_info.push(format!(
                         "Original vs raw:\n{}",
                         compare_strings_formatted(smiles, &mol.to_smiles())
                     ));
                     debug_info.push(format!(
                         "Expected vs actual canonical:\n{}",
                         compare_strings_formatted(&comparison_canonical_smiles, &canonical_smiles)
                     ));
                     let debug_message = format!(
                         "\n=== Debug Information for SMILES at line {} ===\n{}\n===============================================\n",
                         line_counter,
                         debug_info.join("\n")
                     );
                     // Write debug information to the log file
                     let mut file = log_file.lock().unwrap();
                     file.write_all(debug_message.as_bytes()).unwrap();
                     
                     // Exit early on mismatch
                     return;
                 }
                 if (
                     number_of_atoms,
                     number_of_hydrogens,
                     number_of_single_bonds,
                     number_of_double_bonds,
                     number_of_aromatic_bonds
                 ) != (
                     comparison_number_of_atoms,
                     comparison_number_of_hydrogens,
                     comparison_number_of_single_bonds,
                     comparison_number_of_double_bonds,
                     comparison_number_of_aromatic_bonds
                 ) {
                     return;
                 }
             }
         });

         // Write final statistics to the log file
         let final_statistics = format!(
             "\n=== Final Statistics ===\nTotal SMILES processed: {}\nCorrect SMILES: {}\nSuccess rate: {:.2}%",
             total_smiles.load(Ordering::Relaxed),
             correct_smiles.load(Ordering::Relaxed),
             if total_smiles.load(Ordering::Relaxed) > 0 {
                 (correct_smiles.load(Ordering::Relaxed) as f64 / total_smiles.load(Ordering::Relaxed) as f64) * 100.0
             } else {
                 0.0
             }
         );

         let mut file = log_file.lock().unwrap();
         file.write_all(final_statistics.as_bytes()).unwrap();
         
         assert_eq!(
             correct_smiles.load(Ordering::Relaxed),
             total_smiles.load(Ordering::Relaxed),
             "Not all SMILES were correct!"
         );
     }

 }

