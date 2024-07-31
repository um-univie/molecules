use crate::molecule::Molecule;
use chemistry_consts::ElementProperties;
pub trait ToGml {
    fn to_gml(&self) -> String;
    fn to_gml_with_slice(&self, filter: &[usize], outer_filter: &[usize]) -> String;
}

impl<T: Molecule> ToGml for T {
    /// This function converts a `Molecule` to a gml file
    ///
    /// # Panics
    /// This function does not panic
    ///
    /// # Example
    /// ```
    /// use molecules::prelude::*; 
    /// let mut molecule = Molecule3D::from_xyz("tests/ethane.xyz");
    /// let gml = molecule.to_gml();
    /// println!("{}", &gml);
    /// let test = std::fs::read_to_string("tests/ethane.gml").unwrap();
    /// // Dummy to fail
    /// ```
    ///
    fn to_gml(&self) -> String {
        let components = self.get_components();
        let mut result_string = String::new();

        for component in components {
            result_string.push_str(&self.to_gml_with_slice(&component, &[]));
            result_string.push('\n');
        }
        result_string
    }

    fn to_gml_with_slice(&self, filter: &[usize], outer_filter: &[usize]) -> String {
        let mut bonds = String::new();
        let mut atoms = String::new();
        let mut result_string = String::new();
        let atom_bonds = self.atom_bonds();

        for &index in filter.iter() {
            let mut charge = match self.get_atom_charge(index) {
                charge@1..=i8::MAX => format!("{}+", charge),
                charge@i8::MIN..=-1 => format!("{}-", charge.abs()),
                0 => String::new(),
            };
            if self.is_atom_radical(index) {
                charge.push('.');
            }
            let atomic_number = self.get_atomic_number(index);
            let atomic_symbol = atomic_number.atomic_symbol().unwrap_or("Unknown");
            atoms.push_str(&format!(
                "node [id {index} label \"{}{}\"]\n",
                atomic_symbol,
                charge
            ));
            for &bond in atom_bonds[index].iter() {
                if index < bond.target()
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
