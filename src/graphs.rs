use crate::molecule::{BondChange, Molecule};
use nohash_hasher::IntMap;
use std::collections::{HashMap, HashSet};
use std::fmt::Write;

impl Molecule {
    /// This function converts a `Molecule` to a gml file
    ///
    /// # Panics
    /// This function does not panic
    ///
    /// # Example
    /// ```
    /// use molecules::molecule::Molecule;
    /// let mut molecule = Molecule::from_xyz("tests/ethane.xyz");
    /// let gml = molecule.to_gml();
    /// println!("{}", &gml);
    /// let test = std::fs::read_to_string("tests/ethane.gml").unwrap();
    /// assert_eq!(gml, test);
    /// // Dummy to fail
    /// ```
    ///
    pub fn to_gml(&self) -> String {
        let mut bonds = String::new();
        let mut atoms = String::new();
        let components = self.get_components();
        let mut result_string = String::new();

        for component in components {
            for index in component.iter() {
                let atom = &self.atoms[*index];
                let mut charge = match atom.charge {
                    0 => String::new(),
                    1 => "+".to_string(),
                    -1 => "-".to_string(),
                    2..=i8::MAX => format!("{}+", atom.charge),
                    i8::MIN..=-2 => format!("{}-", atom.charge.abs()),
                };

                if atom.is_radical {
                    charge.push('.');
                }

                let _ = writeln!(
                    atoms,
                    "node [id {index} label \"{}{}\"]",
                    atom.name(),
                    charge
                );

                atom.bonds().iter().for_each(|bond| {
                    if *index < bond.target() {
                        let label = bond.bond_type().to_string();
                        let _ = writeln!(
                            bonds,
                            "edge [source {index} target {} label \"{}\"]",
                            bond.target(),
                            label
                        );
                    }
                });
            }
            result_string.push_str(&format!("graph [\n{}{}]\n", atoms, bonds));
            atoms.clear();
            bonds.clear();
        }
        result_string
    }

    pub fn to_gml_with_slice(&self, filter: &[usize], outer_filter: &[usize]) -> String {
        let mut bonds = String::new();
        let mut atoms = String::new();
        let mut result_string = String::new();

        for index in filter.iter() {
            let atom = &self.atoms[*index];
            let mut charge = match atom.charge {
                1..=i8::MAX => format!("{}+", atom.charge),
                i8::MIN..=-1 => format!("{}-", atom.charge.abs()),
                0 => String::new(),
            };
            if atom.is_radical {
                charge.push('.');
            }
            atoms.push_str(&format!(
                "node [id {index} label \"{}{}\"]\n",
                atom.name(),
                charge
            ));
            for &bond in atom.bonds().iter() {
                if *index < bond.target()
                    && (filter.contains(&bond.target()) || outer_filter.contains(&bond.target()))
                {
                    let label = bond.bond_type().to_string();
                    bonds.push_str(&format!(
                        "edge [source {index} target {} label \"{}\"]\n",
                        bond.target(),
                        label
                    ));
                }
            }
        }
        result_string.push_str("graph [");
        if !atoms.is_empty() {
            result_string.push('\n');
            result_string.push_str(&atoms);
            result_string.push('\n');
        }
        if !bonds.is_empty() {
            if atoms.is_empty() {
                result_string.push('\n');
            }
            result_string.push_str(&bonds);
            result_string.push('\n');
        }
        result_string.push(']');
        result_string
    }

    pub fn extract_reaction_rule(&self, other: &Molecule, radius: usize) {
        let charge_changes = self.charge_changes(other);

        let changed_nodes: HashSet<usize> =
            charge_changes.iter().map(|changes| changes.0).collect();

        let bond_changes = self.identify_bond_changes(other);

        let broken_bonds = bond_changes
            .iter()
            .filter(|bond_change| bond_change.is_broken())
            .map(|bond_change| (bond_change.atom1_index(), bond_change.atom2_index()))
            .collect::<Vec<(usize, usize)>>();

        println!("New Rule:\n ==========\n");

        if broken_bonds.len() > 1 {
            println!("\n\nBroken bonds: {:?}", broken_bonds);
            println!("Charge changes: {:?}", charge_changes);

            // Here we still need to implement the radius method.
            let minimum_path = self.find_shortest_path_of_set(&changed_nodes);
            println!("Minimum path: {:?}\n\n", minimum_path);
        } else if !broken_bonds.is_empty() {
            println!("\n\nBroken bonds: {:?}", broken_bonds);
            println!("Charge changes: {:?}", charge_changes);
            let minimum_path = match self.find_shortest_path_of_set(&changed_nodes) {
                Some(path) => path,
                None => {
                    println!("No path found");
                    return;
                }
            };
            // TODO: make this work for radius > 1 as well and change to Option
            let additional_nodes = if radius > 0 {
                minimum_path
                    .iter()
                    .flat_map(|index| {
                        self.atoms[*index]
                            .bonds()
                            .iter()
                            .filter(|bond| !minimum_path.contains(&bond.target()))
                            .map(|bond| bond.target())
                            .collect::<Vec<usize>>()
                    })
                    .collect::<Vec<usize>>()
            } else {
                Vec::new()
            };

            let rule = assemble_reaction_rule(
                self,
                other,
                &charge_changes,
                &bond_changes,
                &minimum_path,
                &additional_nodes,
            );
            println!("\n\n{}", rule);
            println!("Charge changes: {:?}", charge_changes);
            println!("Minimum path {:?}\n\n", broken_bonds);
        }
    }
}

fn assemble_reaction_rule(
    prior_molecule: &Molecule,
    post_molecule: &Molecule,
    charge_changes: &[(usize, i8, i8)],
    bond_changes: &[BondChange],
    minimum_path: &[usize],
    additional_nodes: &[usize],
) -> String {
    let charge_change_indices = charge_changes
        .iter()
        .map(|charge_change| charge_change.0)
        .collect::<Vec<usize>>();

    let broken_bonds = bond_changes
        .iter()
        .filter(|bond_change| bond_change.is_broken())
        .map(|bond_change| (bond_change.atom1_index(), bond_change.atom2_index()))
        .collect::<Vec<(usize, usize)>>();

    let formed_bonds = bond_changes
        .iter()
        .filter(|bond_change| bond_change.is_formed())
        .map(|bond_change| (bond_change.atom1_index(), bond_change.atom2_index()))
        .collect::<Vec<(usize, usize)>>();

    let all_bonds = minimum_path
        .iter()
        .flat_map(|index| {
            prior_molecule.atoms[*index]
                .bonds()
                .iter()
                .map(|bond| (*index, bond.target()))
                .collect::<Vec<(usize, usize)>>()
        })
        .collect::<Vec<(usize, usize)>>();

    let left_nodes = charge_changes
        .iter()
        .map(|charge_change| {
            format!(
                "node [id {} label \"{}{}\"]",
                charge_change.0,
                prior_molecule.atoms[charge_change.0].name(),
                {
                    let charge = prior_molecule.atoms[charge_change.0].charge();
                    match charge {
                        0 => String::new(),
                        1.. => format!("{}+", charge),
                        ..=-1 => format!("{}-", charge.abs()),
                    }
                },
            )
        })
        .collect::<Vec<String>>()
        .join("\n");

    let left_edges = broken_bonds
        .iter()
        .map(|(index1, index2)| {
            format!(
                "edge [source {} target {} label \"{}\"]",
                index1,
                index2,
                    prior_molecule.atoms[*index1]
                        .bonds()
                        .iter()
                        .find(|bond| bond.target() == *index2)
                        .unwrap()
                        .bond_type()
            )
        })
        .collect::<Vec<String>>()
        .join("\n");

    let context_nodes = minimum_path
        .iter()
        .filter(|index| !charge_change_indices.contains(index))
        .fold(String::new(), |mut acc, index| {
            let _ = writeln!(
                acc,
                "node [id {} label \"{}\"]",
                index,
                prior_molecule.atoms[*index].name()
            );
            acc
        });

    let context_edges = all_bonds
        .iter()
        .filter(|(index1, index2)| {
            !charge_change_indices.contains(index1) && !charge_change_indices.contains(index2)
        })
        .map(|(index1, index2)| {
            format!(
                "edge [source {} target {} label \"{}\"]",
                index1,
                index2,
                prior_molecule.atoms[*index1]
                        .bonds()
                        .iter()
                        .find(|bond| bond.target() == *index2)
                        .unwrap()
                        .bond_type()
                )
                .to_string()
        })
        .collect::<Vec<String>>()
        .join("\n");

    let right_nodes = charge_changes
        .iter()
        .map(|charge_change| {
            format!(
                "node [id {} label \"{}{}\"]",
                charge_change.0,
                post_molecule.atoms[charge_change.0].name(),
                {
                    let charge = post_molecule.atoms[charge_change.0].charge();
                    match charge {
                        0 => String::new(),
                        1.. => format!("{}+", charge),
                        ..=-1 => format!("{}-", charge.abs()),
                    }
                },
            )
        })
        .collect::<Vec<String>>()
        .join("\n");
    let right_edges = formed_bonds
        .iter()
        .map(|(index1, index2)| {
            format!(
                "edge [source {} target {} label \"{}\"]",
                index1,
                index2,
                post_molecule.atoms[*index1]
                        .bonds()
                        .iter()
                        .find(|bond| bond.target() == *index2)
                        .unwrap()
                        .bond_type()
            )
        })
        .collect::<Vec<String>>()
        .join("\n");

    let left_hand = format!("{}\n{}", left_nodes, left_edges);
    let context = format!("{}\n{}", context_nodes, context_edges);
    let right_hand = format!("{}\n{}", right_nodes, right_edges);
    #[cfg(debug_assertions)]
    {
        println!("Left nodes: {}", left_nodes);
        println!("Left edges: {}", left_edges);
        println!("Context nodes: {}", context_nodes);
        println!("Context edges: {}", context_edges);
        println!("Right nodes: {}", right_nodes);
        println!("Right edges: {}", right_edges);
    }

    let rule = format!(
        "rule[\nleft[\n{}\n]\ncontext[\n{}\n]\nright[\n{}\n]\n]",
        left_hand, context, right_hand
    );
    reindex_numbers_in_string(&rule, 1)
}

fn reindex_numbers_in_string(input: &str, start_index: usize) -> String {
    let mut index = start_index;
    let mut number_mapping = HashMap::new();

    // First pass: identify numbers and assign new indices
    for word in input.split_whitespace() {
        if word.chars().all(char::is_numeric) {
            number_mapping.entry(word.to_string()).or_insert_with(|| {
                let current_index = index.to_string();
                index += 1;
                current_index
            });
        }
    }

    // Second pass: replace numbers with their new indices
    input
        .split_whitespace()
        .map(|word| {
            number_mapping
                .get(word)
                .unwrap_or(&word.to_string())
                .clone()
        })
        .collect::<Vec<String>>()
        .join(" ")
}

#[allow(dead_code)]
fn reindex_string(input: &str) -> String {
    let mut new_index = 1;
    let mut index_map = IntMap::default();
    let mut output = String::new();
    let keywords = ["id", "source", "target"];

    input.split_whitespace().enumerate().for_each(|(i, word)| {
        if keywords.contains(&word) {
            if let Some(next_word) = input.split_whitespace().nth(i + 1) {
                if let Ok(id) = next_word.parse::<usize>() {
                    index_map.entry(id).or_insert_with(|| {
                        let current = new_index;
                        new_index += 1;
                        current
                    });
                }
            }
        }
    });

    let mut is_number_next = false;
    for word in input.split_whitespace() {
        if is_number_next {
            if let Ok(id) = word.parse::<usize>() {
                if let Some(&new_id) = index_map.get(&id) {
                    output.push_str(&new_id.to_string());
                } else {
                    output.push_str(word);
                }
            } else {
                output.push_str(word);
            }
            is_number_next = false;
        } else {
            if keywords.contains(&word) {
                is_number_next = true;
            }
            output.push_str(word);
        }
        output.push(' ');
    }

    output.trim_end().to_string()
}
