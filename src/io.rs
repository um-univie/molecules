use crate::atom::Atom;
pub use crate::{
    chirality::ChiralClass,
    molecule::base::Molecule,
    molecule::bond::{BondOrder, BondTarget},
    molecule::molecule2d::Molecule2D,
};
pub use chemistry_consts::ElementProperties;
pub use nohash_hasher::IntMap;

#[derive(Debug, Clone)]
pub enum ParseError {
    ElementNotFound(String),
    BondNotFound,
    RingIndexError,
    InvalidBranch,
    InvalidAromatic(String),
    SMILESComplexError,
    Charge,
    ChiralClass(String),
    AtomClassMissing,
    EOL,
}

#[derive(Default)]
pub struct SMILESParser {
    atoms: Vec<Atom>,
    bonds: Vec<(usize, usize, BondOrder)>,
    current_atom_index: usize,
    current_atom_charge: Option<i8>,
    element_buffer: String,
    isotope: Option<u16>,
    chiral_class: ChiralClass,
    current_atom_class: Option<u8>,
    last_bond_type: BondOrder,
    ring_number: Option<usize>,
    ring_bonds: IntMap<usize, (Option<usize>, Option<usize>, Option<BondOrder>)>,
    rings: Vec<(usize,usize,BondOrder)>,
    hydrogens: IntMap<usize, u8>,
    branch_stack: Vec<usize>,
    branch_exits: usize,
    is_multiple_branch: bool,
    is_double_digit: bool,
}

impl SMILESParser {
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
    /// Parses a SMILES string and returns a Molecule
    /// # Arguments
    /// * `smiles` - A string slice that holds the SMILES string
    ///
    /// # Example
    /// ```
    /// use molecules::prelude::*;
    /// let molecules = Molecule2D::from_smiles("C(C(C))COCCl").unwrap();
    /// assert_eq!(molecules[0].atomic_numbers().len(), 18);
    /// assert_eq!(molecules[0].get_edges().len(), 17);
    /// ```
    pub fn parse_smiles(smiles: &str) -> Result<Vec<Molecule2D>, ParseError> {
        let mut molecules = SMILESParser::parse_smiles_raw(smiles)?;
        for molecule in molecules.iter_mut() {
            molecule.add_hydrogens();
            molecule.add_aromatic_bonds();
        }
        Ok(molecules)
    }
    pub fn parse_smiles_raw(smiles: &str) -> Result<Vec<Molecule2D>, ParseError> {
        let mut parser = SMILESParser::default();
        let mut molecules = Vec::new();

        // SMILES should only be ASCII so we could bypass the UTF-8 checks using bytes()
        let bytes = smiles.as_bytes();
        let mut pointer = 0;
        while pointer < bytes.len() {
            // This is kind of redundant maybe remove
            let Some(&byte) = bytes.get(pointer) else {
                break;
            };
            pointer += 1;
            match byte {
                b'A'..=b'Z' => {
                    // Handle the previous element
                    parser.handle_atom(Some(byte))?
                }
                b'a'..=b'z' => {
                    // According to the SMILES specification, the lowercase letters are used to denote aromatic atoms if they are not complex
                    match byte {
                        b'b' | b'c' | b'n' | b'o' | b's' | b'p' => {
                            parser.handle_atom(Some(byte))?
                        }
                        b'r' => {
                            if parser.element_buffer == "B" {
                                parser.element_buffer.push('R')
                            } else {
                                return Err(ParseError::InvalidAromatic(
                                    parser.element_buffer.clone(),
                                ));
                            }
                        }

                        b'l' => {
                            if parser.element_buffer == "C" {
                                parser.element_buffer.push('L')
                            } else {
                                return Err(ParseError::InvalidAromatic(
                                    parser.element_buffer.clone(),
                                ));
                            }
                        }
                        anything_else => {
                            parser.element_buffer.push(anything_else as char);
                            return Err(ParseError::InvalidAromatic(parser.element_buffer.clone()));
                        }
                    }
                }

                b'.' => {
                    parser.handle_atom(None)?;
                    parser.add_all_bonds();
                    let molecule = Molecule2D::from_atoms(parser.atoms);
                    molecules.push(molecule);
                    parser = SMILESParser::default();
                }

                b'1'..=b'9' => {
                    parser.handle_number(byte)?;
                }
                b'-' => {
                    parser.last_bond_type = BondOrder::Single;
                }
                b'=' => {
                    parser.last_bond_type = BondOrder::Double;
                }
                b'#' => {
                    parser.last_bond_type = BondOrder::Triple;
                }
                b'$' => {
                    parser.last_bond_type = BondOrder::Quadruple;
                }
                b':' => {
                    parser.last_bond_type = BondOrder::Aromatic;
                }
                b'(' => {
                    if parser.branch_exits > 0 {
                        parser.is_multiple_branch = true;
                    } else {
                        parser.branch_stack.push(parser.current_atom_index);
                    }
                }
                b')' => parser.branch_exits += 1,
                b'[' => {
                    parser.handle_complex_atom(bytes, &mut pointer)?;
                }

                b'%' => {
                    parser.is_double_digit = true;
                }
                _ => (),
            }
        }

        if !parser.element_buffer.is_empty() {
            parser.handle_atom(None)?
        }

        parser.add_all_bonds();
        let molecule = Molecule2D::from_atoms(parser.atoms);
        molecules.push(molecule);
        Ok(molecules)
    }

    /// Adds all bonds to the current molecule, ensuring aromatic bonds are correctly assigned.
    ///
    /// This method processes ring bonds and explicit hydrogen bonds, assigning the appropriate
    /// bond orders based on the aromaticity of the connected atoms.
    ///
    /// Aromatic bonds between two aromatic atoms are marked as `BondOrder::Aromatic`.
    /// All other bonds default to their parsed bond orders.
    fn add_all_bonds(&mut self) {
        for entry in self.ring_bonds.iter() {
            match (entry.1.0, entry.1.1) {
                (Some(start_idx), Some(end_idx)) => self.rings.push((start_idx, end_idx, entry.1.2.unwrap_or(BondOrder::Single))),
                (None, None) => (),
                _ => ()
            }
        }


        for (start_idx, end_idx, order) in self.rings.iter() {
                if *start_idx >= self.atoms.len() || *end_idx >= self.atoms.len() {
                    continue; // Skip invalid bonds
                }

                // This is false if they are not in the same system
                let bond_order = if self.atoms[*start_idx].aromatic && self.atoms[*end_idx].aromatic
                {
                    BondOrder::Aromatic
                } else {
                    *order
                };

                self.bonds.push((*start_idx, *end_idx, bond_order));
        }

        for (atom_index, hydrogen_count) in self.hydrogens.iter() {
            for _ in 0..*hydrogen_count {
                let hydrogen_index = self.atoms.len();
                self.atoms.push(Atom::new(1)); // Hydrogen
                self.bonds
                    .push((*atom_index, hydrogen_index, BondOrder::Single));
                //println!(
                //    "Bond added between Atom {} and Hydrogen {} with bond order {:?} in add_all_bonds",
                //    atom_index, hydrogen_index, BondOrder::Single
                //);
            }
        }

        // Assign bonds to the molecule
        for bond in self.bonds.iter_mut() {
            // Ensure bond indices are within bounds
            if bond.0 >= self.atoms.len() || bond.1 >= self.atoms.len() {
                println!(
                    "Error: Bond indices out of bounds: bond.0={}, bond.1={}, atoms.len={}",
                    bond.0,
                    bond.1,
                    self.atoms.len()
                );
                continue; // Skip invalid bonds
            }

            let atom1 = &mut self.atoms[bond.0];
            atom1.add_bond(BondTarget::new(bond.1, bond.2));
            let atom2 = &mut self.atoms[bond.1];
            atom2.add_bond(BondTarget::new(bond.0, bond.2));
        }

        for index in 0..self.atoms.len() {
            let is_atom_aromatic = self.atoms[index].aromatic;
            for bond_index in 0..self.atoms[index].bonds.len() {
                let bond = &self.atoms[index].bonds[bond_index];
                let is_bond_target_aromatic = self.atoms[bond.target].aromatic;
                if is_atom_aromatic && is_bond_target_aromatic {
                    self.atoms[index].bonds[bond_index].bond_order = BondOrder::Aromatic;
                }
            }
        }
    }

    fn handle_branch(&mut self) -> Result<(), ParseError> {
        // Pop the stack until we find the last branch atom
        let mut branch_atom = None;
        while let Some(branch) = self.branch_stack.pop() {
            branch_atom = Some(branch);
            if self.branch_exits > 0 {
                self.branch_exits -= 1;
            } else {
                return Err(ParseError::InvalidBranch);
            }
            if self.branch_exits == 0 {
                break;
            }
        }

        if self.branch_exits > 0 {
            return Err(ParseError::InvalidBranch);
        }

        if let Some(branch_atom) = branch_atom {
            self.bonds
                .push((branch_atom, self.current_atom_index, self.last_bond_type));
            if self.is_multiple_branch {
                // If we have multiple branches, we need to push the current atom back on the stack
                self.branch_stack.push(branch_atom);
                self.is_multiple_branch = false;
            }
        } else {
            return Err(ParseError::BondNotFound);
        }
        Ok(())
    }

    /// Handles the creation of a bond between the current atom and the previous atom.
    ///
    /// This method determines the bond order based on the last bond type parsed
    /// and the aromaticity of the connected atoms. If both atoms are aromatic,
    /// the bond is set to `BondOrder::Aromatic`, overriding any previously assigned bond order.
    ///
    fn handle_bond(&mut self) -> Result<(), ParseError> {
        if self.current_atom_index == 0 {
            return Ok(()); // No bond to handle if it's the first atom
        }

        if self.branch_exits > 0 && !self.branch_stack.is_empty() {
            self.handle_branch()?;
        } else {
            let previous_atom = self.current_atom_index - 1;
            let current_atom = self.current_atom_index;

            let bond_order = if let (Some(prev_atom), Some(curr_atom)) =
                (self.atoms.get(previous_atom), self.atoms.get(current_atom))
            {
                // Check if both connected atoms are aromatic
                if prev_atom.aromatic && curr_atom.aromatic {
                    BondOrder::Aromatic
                } else {
                    // Otherwise, use the last bond type parsed
                    self.last_bond_type
                }
            } else {
                self.last_bond_type
            };

            // Add the bond and print debug information
            self.bonds.push((previous_atom, current_atom, bond_order));
        }
        self.last_bond_type = BondOrder::Single; // Reset to default bond type
        Ok(())
    }

    /// Handles the parsing of an atom in the SMILES string.
    ///
    /// This method processes the current element buffer, determines the atomic number,
    /// and updates the atom's properties such as aromaticity, isotope, chiral class, and charge.
    ///
    /// If the `byte` parameter is provided, it indicates the next character in the SMILES string
    /// and may influence the bond type.
    ///
    /// # Arguments
    ///
    /// * `byte` - An optional byte representing the next character in the SMILES string.
    ///
    fn handle_atom(&mut self, byte: Option<u8>) -> Result<(), ParseError> {
        if !self.element_buffer.is_empty() {
            let atomic_number = self
                .element_buffer
                .to_uppercase()
                .as_str()
                .atomic_number()
                .ok_or(ParseError::ElementNotFound(self.element_buffer.clone()))?;

            let mut atom = Atom::new(atomic_number).with_atom_class(self.current_atom_class);

            self.current_atom_index += 1;

            // Determine aromaticity based on lowercase symbol or aromatic flag
            if self.element_buffer.chars().next().unwrap().is_lowercase() {
                atom.aromatic = true;
            }

            // Handle bond if additional byte is present
            if byte.is_some() {
                self.handle_bond()?
            }

            if let Some(isotope) = self.isotope {
                if is_valid_isotope(atomic_number, isotope) {
                    atom = atom.with_isotope(isotope);
                    self.isotope = None;
                } else {
                    println!(
                        "Isotope {} is not valid for atomic number {}, skipping it",
                        isotope, atomic_number
                    );
                }
            }

            if self.chiral_class != ChiralClass::None {
                atom = atom.with_chiral_class(self.chiral_class);
                self.chiral_class = ChiralClass::None;
            }

            if let Some(charge) = self.current_atom_charge {
                atom = atom.with_charge(charge);
                self.current_atom_charge = None;
            }

            self.atoms.push(atom);
            self.element_buffer.clear();
        }

        if let Some(character_byte) = byte {
            self.element_buffer.push(character_byte as char);
        }

        Ok(())
    }

    /// Handles numeric characters in SMILES strings, particularly for ring closures
    ///
    /// # Arguments
    /// * `byte` - ASCII byte representing a digit (1-9)
    ///
    /// # Returns
    /// * `Result<(), ParseError>` - Ok if number handled successfully, Err if invalid ring format
    fn handle_number(&mut self, byte: u8) -> Result<(), ParseError> {
        // Handle double-digit ring numbers (after %)
        if self.is_double_digit {
            let digit = byte_to_number(byte) as usize;
            self.ring_number = Some(match self.ring_number {
                Some(existing) => existing * 10 + digit,
                None => digit,
            });
            self.is_double_digit = false;
            return Ok(());
        }

        // Handle single-digit ring numbers
        if self.ring_number.is_none() {
            self.ring_number = Some(byte_to_number(byte) as usize);
        }

        // Process the complete ring number
        let ring = self.ring_number.take().unwrap();
        let entry = self.ring_bonds.entry(ring).or_insert((None, None, None));

        match (entry.0, entry.1) {
            (None, _) => {
                entry.0 = Some(self.current_atom_index);
                entry.2 = Some(self.last_bond_type); // Store bond order when first encountered
            },
            (Some(_), None) => {
                entry.1 = Some(self.current_atom_index);
                entry.2 = Some(self.last_bond_type);
            },
            _ => {
                // Ring closure is complete, create new entry
                self.rings.push((
                    entry.0.unwrap(),
                    entry.1.unwrap(),
                    entry.2.unwrap_or(BondOrder::Single)
                ));
                *entry = (Some(self.current_atom_index), None, Some(self.last_bond_type));
            }
        }

        self.last_bond_type = BondOrder::Single; // Reset bond type after handling

        // Validate that the ring closure isn't connecting an atom to itself
        if entry.0 == entry.1 && entry.1.is_some() {
            return Err(ParseError::RingIndexError);
        }

        Ok(())
    }

    pub fn handle_complex_atom(
        &mut self,
        bytes: &[u8],
        position: &mut usize,
    ) -> Result<(), ParseError> {
        self.handle_atom(None)?;

        if bytes[*position].is_ascii_digit() {
            let mut temp_number: u16 = 0;
            while bytes[*position].is_ascii_digit() {
                let number = byte_to_number(bytes[*position]);
                temp_number = temp_number * 10 + number as u16;
                *position += 1;
            }
            self.isotope = Some(temp_number);
        }

        let mut is_se_or_as = false;
        let test = [bytes[*position], bytes[*position + 1]];
        match &test {
            b"se" => {
                is_se_or_as = true;
                self.element_buffer.push_str("SE");
            }
            b"as" => {
                is_se_or_as = true;
                self.element_buffer.push_str("AS");
            }
            _ => (),
        }

        if bytes[*position].is_ascii() && !is_se_or_as {
            self.element_buffer.push(bytes[*position] as char);
            if bytes[*position + 1].is_ascii_lowercase() {
                *position += 1;
                self.element_buffer.push(bytes[*position] as char);
            }
        } 

        *position += 1;
        let start_position = *position;

        if bytes[*position] == b'@' {
            let mut temp_counter = 0;
            *position += 1;
            while bytes[*position].is_ascii_uppercase() {
                if bytes[*position] == b'H' && temp_counter == 0 {
                    break;
                }
                *position += 1;
                temp_counter += 1;
                if temp_counter > 2 {
                    return Err(ParseError::ChiralClass(
                        "Chiral class is too long".to_string(),
                    ));
                }
            }
            while bytes[*position].is_ascii_digit() {
                *position += 1;
            }
            self.chiral_class = parse_chiral_class(&bytes[start_position..*position])?;

        }

        if bytes[*position] == b'H' {
            *position += 1;
            // Theoretically we could have more than 9 hydrogen in extreme cases but it is not accepted in SMILES
            if bytes[*position].is_ascii_digit() {
                self.hydrogens
                    .insert(self.current_atom_index, byte_to_number(bytes[*position]));
                *position += 1;
            } else {
                self.hydrogens.insert(self.current_atom_index, 1);
            }
        }

        if bytes[*position] == b'+' || bytes[*position] == b'-' {
            let sign = bytes[*position];
            *position += 1;
            if bytes[*position].is_ascii_digit() {
                match sign {
                    b'+' => self.current_atom_charge = Some(byte_to_number(bytes[*position]) as i8),
                    b'-' => {
                        self.current_atom_charge = Some(-(byte_to_number(bytes[*position]) as i8))
                    }
                    _ => (), // This cant happen
                }
            } else {
                match sign {
                    b'+' => self.current_atom_charge = Some(1),
                    b'-' => self.current_atom_charge = Some(-1),
                    _ => (), // This cant happen
                }
            }
        }

        if bytes[*position] == b':' {
            *position += 1;
            if bytes[*position].is_ascii_digit() {
                self.current_atom_class = Some(byte_to_number(bytes[*position]));
            } else {
                return Err(ParseError::AtomClassMissing);
            }
        }
        self.handle_bond()?;
        Ok(())
    }
}

fn byte_to_number(byte: u8) -> u8 {
    byte - b'0'
}

fn parse_chiral_class(slice: &[u8]) -> Result<ChiralClass, ParseError> {
    match slice {
        s if s.starts_with(b"@@") => Ok(ChiralClass::TH(2)),
        s if s.starts_with(b"@AL") => {
            let number = parse_number_on_end_of_chiral_class(s);
            Ok(ChiralClass::AL(number))
        }
        s if s.starts_with(b"@SP") => {
            let number = parse_number_on_end_of_chiral_class(s);
            Ok(ChiralClass::SP(number))
        }
        s if s.starts_with(b"@TB") => {
            let number = parse_number_on_end_of_chiral_class(s);
            Ok(ChiralClass::TB(number))
        }
        s if s.starts_with(b"@OH") => {
            let number = parse_number_on_end_of_chiral_class(s);
            Ok(ChiralClass::OH(number))
        }
        s if s.starts_with(b"@") => Ok(ChiralClass::TH(1)),
        _ => Err(ParseError::ChiralClass(
            String::from_utf8_lossy(slice).to_string(),
        )),
    }
}

fn parse_number_on_end_of_chiral_class(chiral_class: &[u8]) -> u8 {
    let mut number = 0;
    for &byte in chiral_class {
        if byte.is_ascii_digit() {
            number = number * 10 + byte_to_number(byte);
        }
    }
    number
}

fn is_valid_isotope(atomic_number: u8, isotope: u16) -> bool {
    let Some(isotopes) = atomic_number.isotopes() else {
        // If the atomic number is not found, we assume that the isotope is not valid
        // This is a bit of a hack but it should work
        return false;
    };

    for iso in isotopes {
        // TODO check if this is correct for all cases
        if iso.mass.round() as u16 == isotope {
            return true;
        }
    }
    false
}

pub trait ToSMILES {
    fn to_smiles(&self) -> String;
}

pub trait FromSMILES {
    fn from_smiles(smiles: &str) -> Result<Vec<Self>, ParseError>
    where
        Self: Sized;
    fn from_smiles_with_sanitization(smiles: &str) -> Result<Vec<Self>, ParseError>
    where 
        Self: Sized;
}

impl FromSMILES for Molecule2D {
    fn from_smiles(smiles: &str) -> Result<Vec<Molecule2D>, ParseError> {
        let molecules = SMILESParser::parse_smiles(smiles)?;
        Ok(molecules)
    }
    fn from_smiles_with_sanitization(smiles: &str) -> Result<Vec<Molecule2D>, ParseError> {
        let mut molecules = SMILESParser::parse_smiles(smiles)?;
        for molecule in molecules.iter_mut() {
            molecule.add_aromatic_bonds();
        }
        Ok(molecules)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_smiles() {
        let molecules = SMILESParser::parse_smiles("C(C(C))COCCl").unwrap();
        assert_eq!(molecules[0].atomic_numbers[0], 6);
        assert_eq!(molecules[0].atomic_numbers.len(), 18);
    }

    #[test]
    fn test_parse_smiles_with_isotope() {
        let molecules = SMILESParser::parse_smiles("CCCC[13C]").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 17);
        assert_eq!(molecules[0].get_isotope(4).unwrap(), 13);
    }

    #[test]
    fn test_parse_smiles_with_chiral_class() {
        let molecules = SMILESParser::parse_smiles("C[C@](F)(Cl)Br").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 8);
        assert_eq!(molecules[0].get_edges().len(), 7);
        assert_eq!(molecules[0].get_chiral_class(1), ChiralClass::TH(1));
    }

    #[test]
    fn test_parse_smiles_with_complex_atom_and_isotope() {
        let molecules = SMILESParser::parse_smiles("C[C@](F)(Cl)[13C](Br)").unwrap();
        println!("{:?}", molecules);
        println!("{:?}", molecules[0].get_edges());
        assert_eq!(molecules[0].atomic_numbers.len(), 11);
        assert_eq!(molecules[0].get_isotope(4).unwrap(), 13);
        assert_eq!(molecules[0].get_chiral_class(1), ChiralClass::TH(1));
        assert_eq!(molecules[0].get_edges().len(), 10);
    }

    #[test]
    fn test_parse_smiles_with_complex_atom_and_hydrogen() {
        let molecules = SMILESParser::parse_smiles("C[C@H](Cl)C").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 11);
        assert_eq!(molecules[0].get_edges().len(), 10);
    }

    #[test]
    fn test_bond_parsing() {
        let molecules = SMILESParser::parse_smiles("C-C-C#C").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 10);
        assert_eq!(molecules[0].get_edges().len(), 9);
    }

    #[test]
    fn parse_test_file() {
        let smiles = std::fs::read_to_string("tests/smiles.txt").unwrap();
        for smile in smiles.lines() {
            println!("Parsing: {}", smile);
            match SMILESParser::parse_smiles(smile) {
                Ok(molecules) => {
                    println!("{:?}", molecules);
                }
                Err(e) => {
                    println!("Error: {:?}", e);
                    panic!();
                }
            }
        }
    }

    #[test]
    fn test_submolecule() {
        let molecules = SMILESParser::parse_smiles("C(C(C))COCCl.C(C(C))").unwrap();
        let submolecule = molecules[0].match_submolecule(&molecules[1]);
        assert_eq!(submolecule, None);
    }

    #[test]
    fn test_ring_closure_with_charge() {
        let molecules = SMILESParser::parse_smiles("c1ccc2c(c1)=[NH+]C(=O)C=2Nc1ccc(C(OCCC)=O)cc1").unwrap();

        println!("Molecule: {}", molecules[0].to_smiles());
        assert_eq!(molecules[0].atomic_numbers.len(), 33);
    }

    #[test]
    fn test_parse_smiles_with_charged_atoms() {
        // Test positive charges
        let molecules = SMILESParser::parse_smiles("[NH4+]").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 5);  // N + 4H
        assert_eq!(molecules[0].get_atom_charge(0), 1);

        // Test negative charges
        let molecules = SMILESParser::parse_smiles("[O-]").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 2); // O + 1H automatically added
        assert_eq!(molecules[0].get_atom_charge(0), -1); 

        // Test multiple charges
        let molecules = SMILESParser::parse_smiles("[NH3+][O-]").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 5);  // N + 3H + O
        assert_eq!(molecules[0].get_atom_charge(0), 1);
        assert_eq!(molecules[0].get_atom_charge(1), -1); // Hydrogens are added last
    }

    #[test]
    fn test_parse_smiles_with_explicit_hydrogens() {
        // Test NH2 group
        let molecules = SMILESParser::parse_smiles("[NH2]").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 4);  // N + 2H + 1H automatically added
        
        // Test CH3 group
        let molecules = SMILESParser::parse_smiles("[CH3]").unwrap();
        println!("{:?}", molecules[0].atomic_numbers);
        assert_eq!(molecules[0].atomic_numbers.len(), 5);  // C + 3H + 1H automatically added
        
        // Test complex molecule with explicit hydrogens
        let molecules = SMILESParser::parse_smiles("[CH3][NH2]").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 7);  // C + 3H + N + 2H
        assert_eq!(molecules[0].get_edges().len(), 6);     // All single bonds
    }

    #[test]
    fn test_parse_smiles_with_charged_complex_molecules() {
        // Test zwitterion (amino acid-like structure)
        let molecules = SMILESParser::parse_smiles("[NH3+]CC[O-]").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 11);  // N + 3H + 2C + O + 4H automatically added
        assert_eq!(molecules[0].get_atom_charge(0), 1);
        assert_eq!(molecules[0].get_atom_charge(6), -1);
        
        // Test multiple charged centers
        let molecules = SMILESParser::parse_smiles("[NH4+]CC[NH3+]").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 10); // 2N + 7H + 2C
        assert_eq!(molecules[0].get_atom_charge(0), 1);
        assert_eq!(molecules[0].get_atom_charge(6), 1);
    }
}

#[cfg(test)]
mod aromatic_bond_tests {
    use super::*;

    #[test]
    fn test_parse_aromatic_bonds() {
        let molecules = SMILESParser::parse_smiles("c1ccccc1").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 12);
        assert_eq!(molecules[0].get_edges().len(), 12);
        let mut aromatic_bonds = 0;
        for bond in molecules[0].get_edges_with_type().into_iter() {
            if bond.2 == BondOrder::Aromatic {
                aromatic_bonds += 1;
            }
        }
        assert_eq!(aromatic_bonds, 6);
    }

    #[test]
    fn test_parse_mixed_bonds_with_aromatic_atoms() {
        let molecules = SMILESParser::parse_smiles("C1=CC=CC=C1").unwrap();

        assert_eq!(molecules[0].atomic_numbers.len(), 12);
        assert_eq!(molecules[0].get_edges().len(), 12);
    }
}
