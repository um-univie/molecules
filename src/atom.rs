use crate::molecule::{BondTarget, BondType, ChiralClass};
use crate::vector::Vector;
use chemistry_consts::ElementProperties;

#[derive(Default, Debug, Clone, PartialEq)]
pub struct Atom {
    pub atomic_number: u8,
    pub charge: i8,
    pub is_radical: bool,
    pub isotope: Option<u16>,
    pub chiral_class: ChiralClass,
    pub atom_class: Option<u8>,
    pub position_vector: Option<Vector>,
    pub bonds: Vec<BondTarget>,
}

impl Atom {
    pub fn new(atomic_number: u8) -> Atom {
        Atom {
            atomic_number,
            ..Default::default()
        }
    }

    pub fn isotope(&self) -> Option<u16> {
        self.isotope
    }

    pub fn with_chiral_class(mut self, chiral_class: ChiralClass) -> Self {
        self.chiral_class = chiral_class;
        self
    }

    pub fn with_position(mut self, position: (f64, f64, f64)) -> Self {
        self.position_vector = Some(Vector::new(position.0, position.1, position.2));
        self
    }

    pub fn with_charge(mut self, charge: i8) -> Self {
        self.charge = charge;
        self
    }

    pub fn with_isotope(mut self, isotope: u16) -> Self {
        self.isotope = Some(isotope);
        self
    }

    pub fn with_atom_class(mut self, atom_class: Option<u8>) -> Self {
        self.atom_class = atom_class;
        self
    }

    pub fn from_xyz_line(line: &str) -> Result<Atom, Box<dyn std::error::Error>> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 4 {
            let atomic_number = parts[0]
                .to_uppercase()
                .as_str()
                .atomic_number()
                .ok_or("Could not find atomic number")?;
            let x = parts[1].parse::<f64>()?;
            let y = parts[2].parse::<f64>()?;
            let z = parts[3].parse::<f64>()?;
            Ok(Atom::new(atomic_number).with_position((x, y, z)))
        } else {
            Err("Could not parse Line".into())
        }
    }

    /// Returns the expected valency of the atom based on its atomic number
    ///
    /// # Examples
    /// ```
    /// use molecules::prelude::*;
    /// let atom = Atom::new(6);
    /// assert_eq!(atom.expected_valency(), Some(4));
    /// ```
    pub fn expected_valency(&self) -> Option<i8> {
        // TODO: Implement for all elements
        self.atomic_number.valencies()?.next()
    }

    pub fn cmp_electronegativity(&self, other: &Self) -> std::cmp::Ordering {
        self.electronegativity().cmp(&other.electronegativity())
    }

    pub fn bonds(&self) -> &Vec<BondTarget> {
        &self.bonds
    }

    pub fn add_bond(&mut self, bond: BondTarget) {
        self.bonds.push(bond);
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn is_aromatic(&self) -> bool {
        self.bonds
            .iter()
            .any(|bond| bond.bond_type() == BondType::Aromatic)
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

    pub fn atomic_symbol(&self) -> Option<&str> {
        self.atomic_number.atomic_symbol()
    }

    pub fn electronegativity(&self) -> Option<u16> {
        self.atomic_number.electronegativity()
    }

    pub fn is_more_electronegative(&self, other: &Self) -> bool {
        self.electronegativity() > other.electronegativity()
    }

    pub fn monoisotopic_mass(&self) -> Option<f64> {
        self.atomic_number.monoisotopic_mass()
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
    /// use molecules::prelude::*;
    /// let atom1 = Atom::new(6).with_position((0.0,0.0,0.0));
    /// let atom2 = Atom::new(6).with_position((1.0,0.0,0.0));
    /// assert_eq!(atom1.distance(&atom2),1.0);
    /// ```
    pub fn distance(&self, other: &Atom) -> f64 {
        let Some(self_position_vector) = self.position_vector else {
            return 0.0;
        };
        let Some(other_position_vector) = other.position_vector else {
            return 0.0;
        };
        self_position_vector.distance(&other_position_vector)
    }

    /// Calculates the squared distance between two atoms saving computation time.
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::prelude::*;
    /// let atom1 = Atom::new(6).with_position((0.0,0.0,0.0));
    /// let atom2 = Atom::new(6).with_position((1.0,0.0,0.0));
    /// assert_eq!(atom1.distance_squared(&atom2),1.0);
    /// ```
    ///
    pub fn distance_squared(&self, other: &Atom) -> f64 {
        let Some(self_position_vector) = self.position_vector else {
            return 0.0;
        };
        let Some(other_position_vector) = other.position_vector else {
            return 0.0;
        };
        self_position_vector.distance_squared(&other_position_vector)
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
    ///
    /// # Example
    /// ```
    /// use molecules::prelude::*;
    /// let atom1 = Atom::new(6).with_position((0.0,0.0,0.0));
    /// let atom2 = Atom::new(6).with_position((1.0,0.0,0.0));
    /// assert_eq!(atom1.potential_energy(&atom2, 1.0, 1.0),0.0);
    /// ```
    pub fn potential_energy(
        &self,
        other: &Atom,
        force_constant: f64,
        equilibrium_distance: f64,
    ) -> f64 {
        0.5 * force_constant * (self.distance(other) - equilibrium_distance).powi(2)
    }

    pub fn standard_atomic_weight(&self) -> Option<f64> {
        self.atomic_number.standard_atomic_weight()
    }

    pub fn atomic_number(&self) -> u8 {
        self.atomic_number
    }
    pub fn degree(&self) -> Option<i8> {
        let atomic_number = self.atomic_number;
        let expected_valency = atomic_number.valencies()?.next()?;
        let actual_valency = self.actual_valency();
        // May need to be changed for elements with unknown valencies
        Some(actual_valency - expected_valency)
    }
    pub fn actual_valency(&self) -> i8 {
        self.bonds
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
            + if self.is_radical { 1 } else { 0 }
            + self.charge.abs()
    }
}
