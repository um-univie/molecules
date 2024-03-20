use crate::consts::{ATOMIC_NUMBERS, ISOTOPES};
use crate::molecule::{Atom, Bond, BondType, ChiralClass, Molecule};
use nohash_hasher::IntMap;

#[derive(Debug)]
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
    bonds: Vec<(usize, usize, BondType)>,
    current_atom_index: usize,
    current_atom_charge: Option<i8>,
    element_buffer: String,
    isotope: Option<u16>,
    chiral_class: ChiralClass,
    atom_class: Option<u8>,
    atom_classes: IntMap<usize, u8>,
    last_bond_type: BondType,
    ring_number: Option<usize>,
    ring_bonds: IntMap<usize, (Option<usize>, Option<usize>)>,
    hydrogens: IntMap<usize, u8>,
    branch_stack: Vec<usize>,
    branch_exits: usize,
    is_multiple_branch: bool,
    is_double_digit: bool,
    is_aromatic: bool,
}


impl SMILESParser {
    /// Parses a SMILES string and returns a Molecule
    /// # Arguments
    /// * `smiles` - A string slice that holds the SMILES string
    ///
    /// # Example
    /// ```
    /// use molecules::io::SMILESParser;
    /// let molecules = SMILESParser::parse_smiles("C(C(C))COcCl").unwrap();
    /// assert_eq!(molecules[0].atoms.len(), 7);
    /// ```
    pub fn parse_smiles(smiles: &str) -> Result<Vec<Molecule>, ParseError> {
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
            println!("Byte: {}", byte as char);
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
                                return Err(ParseError::InvalidAromatic(parser.element_buffer.clone()));
                            }
                        }

                        b'l' => {
                            if parser.element_buffer == "C" {
                                parser.element_buffer.push('L')
                            } else {
                                return Err(ParseError::InvalidAromatic(parser.element_buffer.clone()));
                            }
                        }
                        anything_else => {parser.element_buffer.push(anything_else as char);return Err(ParseError::InvalidAromatic(parser.element_buffer.clone()))},
                    }
                }

                b'.' => {
                    molecules.push(Molecule::from_atoms(parser.atoms).with_classes(parser.atom_classes));
                    parser = SMILESParser::default();
                }

                b'1'..=b'9' => {
                    parser.handle_number(byte)?;
                }
                b'-' => {
                    parser.last_bond_type = BondType::Single;
                }
                b'=' => {
                    parser.last_bond_type = BondType::Double;
                }
                b'#' => {
                    parser.last_bond_type = BondType::Triple;
                }
                b'$' => {
                    parser.last_bond_type = BondType::Quadruple;
                }
                b':' => {
                    parser.last_bond_type = BondType::Aromatic;
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

        // Add ring bonds
        for (start, end) in parser.ring_bonds.values() {
            if start.is_some() && end.is_some() {
                parser
                    .bonds
                    .push((start.unwrap(), end.unwrap(), BondType::Aromatic));
            }
        }

        for (atom, number) in parser.hydrogens.iter() {
            let hydrogen_index = parser.current_atom_index;
            for index in 0..*number {
                parser.atoms.push(Atom::new(1));
                parser.bonds.push((*atom, hydrogen_index + index as usize, BondType::Single));
            }
        }

        // Add bonds to the molecules
        for bond in parser.bonds.iter() {
            let atom1 = &mut parser.atoms[bond.0];
            atom1.add_bond(Bond::new(bond.1, bond.2));
            let atom2 = &mut parser.atoms[bond.1];
            atom2.add_bond(Bond::new(bond.0, bond.2));
        }


        molecules.push(Molecule::from_atoms(parser.atoms).with_classes(parser.atom_classes));
        Ok(molecules)
    }

    fn handle_branch(&mut self) -> Result<(), ParseError> {
        let Some(mut branch_atom) = self.branch_stack.pop() else {
            return Err(ParseError::InvalidBranch);
        };
        self.branch_exits -= 1;
        // Pop the stack until we find the last branch atom
        while self.branch_exits > 0 {
            branch_atom = self.branch_stack.pop().unwrap();
            self.branch_exits -= 1;
        }
        
        self.bonds
            .push((branch_atom, self.current_atom_index, self.last_bond_type));
        if self.is_multiple_branch {
            // If we have multiple branches, we need to push the current atom back on the stack 
            self.branch_stack.push(branch_atom);
            self.is_multiple_branch = false;
        }
        Ok(())
    }
    fn handle_bond(&mut self) -> Result<(), ParseError> {
        if self.current_atom_index == 0 {
            return Ok(());
        }
        if self.branch_exits > 0 && !self.branch_stack.is_empty() {
            self.handle_branch()?
        } else {
            self.bonds.push((
                self.current_atom_index - 1,
                self.current_atom_index,
                self.last_bond_type,
            ));
        }
        self.last_bond_type = BondType::Single;
        Ok(())
    }

    fn handle_atom(&mut self, byte: Option<u8>) -> Result<(), ParseError> {
        if !self.element_buffer.is_empty() {
            let atomic_number = ATOMIC_NUMBERS
                .get(&self.element_buffer.to_uppercase())
                .ok_or(ParseError::ElementNotFound(self.element_buffer.to_owned()))?;

            let mut atom = Atom::new(*atomic_number);

            self.current_atom_index += 1;

            // TODO check if this is correct
            if byte.is_some() {
                self.handle_bond()?
            }

            if self.element_buffer.chars().next().unwrap().is_lowercase() || self.is_aromatic {
                atom = Atom::new(*atomic_number).aromatic();
            }

            if let Some(isotope) = self.isotope {
                if is_valid_isotope(*atomic_number, isotope) {
                    atom = atom.with_isotope(isotope);
                    self.isotope = None;
                } else {
                    println!(
                        "Isotope {} is not valid for atomic number {}, skipping it",
                        isotope, *atomic_number
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

            if let Some(atom_class) = self.atom_class {
                self.atom_classes.insert(self.current_atom_index, atom_class);
                self.atom_class = None;
            }

            self.atoms.push(atom);
            self.element_buffer.clear();
        }
        let Some(character_byte) = byte else {
            return Ok(());
        };
        self.element_buffer.push(character_byte as char);
        Ok(())
    }

    fn handle_number(&mut self, byte: u8) -> Result<(), ParseError> {
        if self.is_double_digit {
            if let Some(ring_number_value) = self.ring_number {
                // TODO make this generic for any size
                self.ring_number = Some(ring_number_value * 10 + byte_to_number(byte) as usize);
                self.is_double_digit = false;
            } else {
                self.ring_number = Some(byte_to_number(byte) as usize);
                return Ok(());
            }
        }

        let Some(ring) = self.ring_number else {return Ok(())};

        let (start, end) = self.ring_bonds.entry(ring).or_insert((None, None));
        // If start is None, then we are at the start of the bond
        if start.is_none() {
            *start = Some(self.current_atom_index);
        } else if end.is_none() {
            *end = Some(self.current_atom_index);
        } else {
            return Err(ParseError::RingIndexError);
        }
        self.ring_number = None;
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
                self.is_aromatic = true;
            }
            b"as" => {
                is_se_or_as = true;
                self.element_buffer.push_str("AS");
                self.is_aromatic = true;
            }
            _ => (),
        }

        if bytes[*position].is_ascii_uppercase() && !is_se_or_as {
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
                self.hydrogens.insert(self.current_atom_index, byte_to_number(bytes[*position]));
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
                    b'-' => self.current_atom_charge = Some(-(byte_to_number(bytes[*position]) as i8)),
                    _ => (), // This should never happen
                }
            } else {
                match sign {
                    b'+' => self.current_atom_charge = Some(1),
                    b'-' => self.current_atom_charge = Some(-1),
                    _ => (), // This should never happen
                }
            }
        }

        if bytes[*position] == b':' {
            *position += 1;
            if bytes[*position].is_ascii_digit() {
                self.atom_class = Some(byte_to_number(bytes[*position]));
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
    println!("Slice: {:?}", slice);
    match slice {
        s if s.starts_with(b"@@") => Ok(ChiralClass::R),
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
        s if s.starts_with(b"@") => Ok(ChiralClass::S),
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
    let isotopes = ISOTOPES.get(atomic_number as usize).unwrap();
    for iso in isotopes.iter().flatten() {
        // TODO check if this is correct
        if iso.mass.round() as u16 == isotope {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_smiles() {
        let molecules = SMILESParser::parse_smiles("C(C(C))COcCl").unwrap();
        assert_eq!(molecules[0].atoms.len(), 7);
    }

    #[test]
    fn test_parse_smiles_with_isotope() {
        let molecules = SMILESParser::parse_smiles("CCCC[13C]").unwrap();
        assert_eq!(molecules[0].atoms.len(), 5);
        assert_eq!(molecules[0].atoms[4].isotope.unwrap(), 13);
    }

    #[test]
    fn test_parse_smiles_with_chiral_class() {
        let molecules = SMILESParser::parse_smiles("C[C@](F)(Cl)Br").unwrap();
        assert_eq!(molecules[0].atoms.len(), 5); assert_eq!(molecules[0].atoms[1].chiral_class, ChiralClass::S);
    }

    #[test]
    fn test_parse_smiles_with_complex_atom_and_isotope() {
        let molecules = SMILESParser::parse_smiles("C[C@](F)(Cl)Br[13C]").unwrap();
        assert_eq!(molecules[0].atoms.len(), 6);
        assert_eq!(molecules[0].atoms[5].isotope.unwrap(), 13);
    }
    
    #[test]
    fn test_parse_smiles_with_complex_atom_and_hydrogen() {
        let molecules  = SMILESParser::parse_smiles("C[C@H](Cl)Br").unwrap();
        assert_eq!(molecules[0].atoms.len(), 5);
        assert_eq!(molecules[0].atoms[1].bonds().len(), 4);
    }

    #[test]
    fn test_bond_parsing() {
        let molecules = SMILESParser::parse_smiles("C-C=C#C").unwrap();
        assert_eq!(molecules[0].atoms.len(), 4);
        assert_eq!(molecules[0].get_edges().len(), 3);
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
}
