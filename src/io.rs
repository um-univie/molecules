use crate::consts::{ATOMIC_NUMBERS,ISOTOPES};
use crate::molecule::{Atom, Bond, BondType, Molecule, ChiralClass};

use nohash_hasher::IntMap;

#[derive(Debug)]
pub enum ParseError { ElementNotFound,
    BondNotFound,
    RingIndexError,
    SMILESComplexError,
    Charge,
    ChiralClass(String),
    EOL,
}

struct SMILESParser {
    atoms: Vec<Atom>,
    bonds: Vec<(usize, usize, BondType)>,
    current_atom_index: usize,
    element_buffer: String,
    last_bond_type: BondType,
    chiral_class: ChiralClass,
    ring_number: Option<usize>,
    ring_bonds: IntMap<usize, (Option<usize>, Option<usize>)>,
    branch_stack: Vec<usize>,
    branch_exits: usize,
    isotope: Option<u16>,
    is_double_digit: bool,
    is_explicit_hydrogen: bool,
}


impl SMILESParser {
    pub fn new() -> SMILESParser {
        SMILESParser {
            atoms: Vec::new(),
            bonds: Vec::new(),
            ring_number: None,
            ring_bonds: IntMap::default(),
            element_buffer: String::new(),
            chiral_class: ChiralClass::None,
            branch_stack: Vec::new(),
            isotope: None,
            last_bond_type: BondType::Single,
            branch_exits: 0,
            is_double_digit: false,
            is_explicit_hydrogen: false,
            current_atom_index: 0,
        }
    }

    fn handle_bond(&mut self) {
        if self.branch_exits > 0 && !self.branch_stack.is_empty() {
            let mut branch_atom = self.branch_stack.pop().unwrap();
            self.branch_exits -= 1;
            while self.branch_exits > 0 {
                branch_atom = self.branch_stack.pop().unwrap();
                self.branch_exits -= 1;
            }
            self.bonds
                .push((branch_atom, self.current_atom_index, self.last_bond_type));
            self.branch_exits = 0;
        } else {
            self.bonds.push((
                self.current_atom_index - 1,
                self.current_atom_index,
                self.last_bond_type,
            ));
        }
        self.last_bond_type = BondType::Single;
    }
    fn handle_atom(&mut self, byte: Option<u8>) -> Result<(), ParseError> {
        if !self.element_buffer.is_empty() {
            // The map is case sensitive, so we need to convert to uppercase
            let atomic_number = ATOMIC_NUMBERS
                .get(&self.element_buffer.to_uppercase())
                .ok_or(ParseError::ElementNotFound)?;
            // TODO: Handle the case where we have only one hydrogen
            if *atomic_number == 1 {
                self.is_explicit_hydrogen = true;
                return Ok(());
            }

            self.current_atom_index += 1;

            if self.current_atom_index > 0 && byte.is_some() {
                self.handle_bond();
            }
            
            let mut atom = Atom::new(*atomic_number);

            if self.element_buffer.chars().next().unwrap().is_lowercase() {
                atom = Atom::new(*atomic_number).aromatic();
            }

            if let Some(isotope) = self.isotope {
                if is_valid_isotope(*atomic_number, isotope) {
                    atom = atom.with_isotope(isotope);
                    self.isotope = None;
                } else {
                    println!("Isotope {} is not valid for atomic number {}, skipping it", isotope, *atomic_number);
                } 
            }

            if self.chiral_class != ChiralClass::None {
                atom = atom.with_chiral_class(self.chiral_class);
                self.chiral_class = ChiralClass::None;
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
                self.ring_number =
                    Some(ring_number_value * 10 + byte_to_number(byte) as usize);
                self.is_double_digit = false;
            } else {
                self.ring_number = Some(byte_to_number(byte) as usize);
                return Ok(());
            }
        }
        let Some(ring) = self.ring_number else {
            return Err(ParseError::RingIndexError);
        };
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

    pub fn handle_complex_atom(&mut self, bytes: &[u8], position: &mut usize) -> Result<(), ParseError> {
        let mut atomic_number = None;
        let mut chiral_class_buffer = String::new();

        if bytes[*position].is_ascii_digit() {
            let mut temp_number: u16 = 0;
            while bytes[*position].is_ascii_digit() {
                let number = byte_to_number(bytes[*position]);
                temp_number = temp_number * 10 + number as u16;
                *position += 1;
            }
            self.isotope = Some(temp_number);
        }
        
        if bytes[*position].is_ascii_uppercase() {
            let mut temp_string = String::new();
            while bytes[*position].is_ascii_lowercase() {
                temp_string.push(bytes[*position] as char);
                *position += 1;
            }
            let temp_atomic_number = ATOMIC_NUMBERS.get(&temp_string).ok_or(ParseError::ElementNotFound)?;
            atomic_number = Some(*temp_atomic_number);
        }

        let start_position = *position;
        if bytes[*position] == b'@' {
            let mut temp_counter = 0;
            *position += 1;
            while bytes[*position].is_ascii_uppercase() {
                chiral_class_buffer.push(bytes[*position] as char);
                *position += 1;
                temp_counter += 1;
                if temp_counter > 2 {
                    return Err(ParseError::ChiralClass("Chiral class is too long".to_string()));
                }
            }
            while bytes[*position].is_ascii_digit() {
                *position += 1;
            }
            self.chiral_class = parse_chiral_class(&bytes,start_position, position)?;
        }

        self.handle_atom(atomic_number)?;

        Ok(())
    }
}

/// Parses a SMILES string and returns a Molecule
///
/// # Arguments
/// * `smiles` - A string slice that holds the SMILES string
///
/// # Example
/// ```
/// use molecules::io::parse_smiles;
/// let molecule = parse_smiles("C(C(C))COcCl").unwrap();
/// assert!(false);
/// ```
pub fn parse_smiles(smiles: &str) -> Result<Molecule, ParseError> {
    let mut parser = SMILESParser::new();

    // SMILES should only be ASCII so we could bypass the UTF-8 checks using bytes()
    let bytes = smiles.as_bytes();
    let mut pointer = 0;
    while pointer < bytes.len() { 
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
                       b'b' | b'c' | b'n' | b'o' | b's' | b'p' => parser.handle_atom(Some(byte))?,
                       b'r' => {
                            if parser.element_buffer == "B" {
                                parser.element_buffer.push('R')
                            } else {
                                return Err(ParseError::SMILESComplexError);
                            }
                        }

                       b'l' => {
                            if parser.element_buffer == "C" {
                                parser.element_buffer.push('L')
                            } else {
                                return Err(ParseError::SMILESComplexError);
                            }
                        }
                        _ => return Err(ParseError::SMILESComplexError),
                    }
            }

           b'1'..=b'9' => {
                parser.handle_number(byte)?;
            }
           b'-' => {
                    parser.bonds.push((
                        parser.current_atom_index,
                        parser.current_atom_index + 1,
                        BondType::Single,
                    ));
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
                parser.branch_stack.push(parser.current_atom_index);
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

    // Add bonds to the molecules
    for bond in parser.bonds.iter() {
        let atom1 = &mut parser.atoms[bond.0];
        atom1.add_bond(Bond::new(bond.1, bond.2));
        let atom2 = &mut parser.atoms[bond.1];
        atom2.add_bond(Bond::new(bond.0, bond.2));
    }

    Ok(Molecule::from_atoms(parser.atoms))
}


fn byte_to_number(byte: u8) -> u8 {
    byte - b'0'
}

fn parse_chiral_class(bytes: &[u8],start_position: usize, position: &mut usize) -> Result<ChiralClass, ParseError> {
    // TODO make this safe
    let slice = &bytes[start_position..*position];

    match slice {
        s if s.starts_with(b"@") => {
            Ok(ChiralClass::R)
        }
        s if s.starts_with(b"AL") => {
            let number = parse_number_on_end_of_chiral_class(s);
            Ok(ChiralClass::AL(number))
        }
        s if s.starts_with(b"SP") => {
            let number = parse_number_on_end_of_chiral_class(s);
            Ok(ChiralClass::SP(number))
        },
        s if s.starts_with(b"TB") => {
            let number = parse_number_on_end_of_chiral_class(s);
            Ok(ChiralClass::TB(number))
        },
        s if s.starts_with(b"OH") => {
            let number = parse_number_on_end_of_chiral_class(s);
            Ok(ChiralClass::OH(number))
        },
        _ => Err(ParseError::ChiralClass(String::from_utf8_lossy(slice).to_string())),
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

fn is_valid_isotope(atomic_number: u8,isotope: u16) -> bool {
    let isotopes = ISOTOPES.get(atomic_number as usize).unwrap();
    for iso in isotopes.iter().flatten() {
        // TODO check if this is correct
        if iso.mass.round() as u16 == isotope {
            return true;
        }
    }
    false

}
