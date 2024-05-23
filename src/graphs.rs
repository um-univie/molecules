use crate::molecule::{BondChange, Molecule};
use nohash_hasher::{IntMap,IntSet};
use std::collections::{HashMap};
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
                    atom.atomic_symbol().unwrap_or("Unknown"),
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
                atom.atomic_symbol().unwrap_or("Unknown"),
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

}
