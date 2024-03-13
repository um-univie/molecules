use crate::consts::ATOMIC_NUMBERS;
use crate::molecule::{Atom, Bond, BondType, Molecule};
use nohash_hasher::IntMap;

#[derive(Debug)]
pub enum ParseError {
    ElementNotFound,
    BondNotFound,
    RingIndexError,
    SMILESComplexError,
    Charge,
}

struct SMILESParser {
    atoms: Vec<Atom>,
    bonds: Vec<(usize, usize, BondType)>,
    ring_number: Option<usize>,
    ring_bonds: IntMap<usize, (Option<usize>, Option<usize>)>,
    element_buffer: String,
    branch_stack: Vec<usize>,
    last_bond_type: BondType,
    branch_exits: usize,
    is_complex: bool,
    is_negative_charge: bool,
    is_double_digit: bool,
    is_explicit_hydrogen: bool,
    current_atom_index: usize,
}

impl SMILESParser {
    pub fn new() -> SMILESParser {
        SMILESParser {
            atoms: Vec::new(),
            bonds: Vec::new(),
            ring_number: None,
            ring_bonds: IntMap::default(),
            element_buffer: String::new(),
            branch_stack: Vec::new(),
            last_bond_type: BondType::Single,
            branch_exits: 0,
            is_complex: false,
            is_negative_charge: false,
            is_double_digit: false,
            is_explicit_hydrogen: false,
            current_atom_index: 0,
        }
    }
    fn handle_bond(&mut self) {
        if self.branch_exits > 0 && !self.branch_stack.is_empty() {
            println!("Branch exit");
            let mut branch_atom = self.branch_stack.pop().unwrap();
            self.branch_exits -= 1;
            while self.branch_exits > 0 {
                branch_atom = self.branch_stack.pop().unwrap();
                self.branch_exits -= 1;
            }
            println!("{} {}", branch_atom, self.current_atom_index);
            println!("Branch stack {:?}", self.branch_stack);
            self.bonds
                .push((branch_atom, self.current_atom_index, self.last_bond_type));
            self.branch_exits = 0;
        } else {
            println!("Normal bond");
            println!(
                "{} {}",
                self.current_atom_index - 1,
                self.current_atom_index
            );
            self.bonds.push((
                self.current_atom_index - 1,
                self.current_atom_index,
                self.last_bond_type,
            ));
        }
        self.last_bond_type = BondType::Single;
    }

    fn handle_atom(&mut self, character: Option<char>) -> Result<(), ParseError> {
        if !self.element_buffer.is_empty() {
            // The map is case sensitive, so we need to convert to uppercase
            let upper = self.element_buffer.to_uppercase().to_string();
            println!("Upper: {}", upper);
            let atomic_number = ATOMIC_NUMBERS
                .get(&upper)
                .ok_or(ParseError::ElementNotFound)?;
            if *atomic_number == 1 {
                self.is_explicit_hydrogen = true;
                return Ok(());
            }

            self.current_atom_index += 1;

            if self.current_atom_index > 0 && character.is_some() {
                self.handle_bond();
            }

            let mut atom = Atom::new(*atomic_number, None);
            if self.element_buffer.chars().next().unwrap().is_lowercase() {
                atom = Atom::new(*atomic_number, None).aromatic();
            }
            self.atoms.push(atom);
            self.element_buffer.clear();
        }
        let Some(character) = character else {
            return Ok(());
        };
        self.element_buffer.push(character);
        Ok(())
    }

    fn handle_number(&mut self, character: char) -> Result<(), ParseError> {
        if self.is_complex {
            if self.is_explicit_hydrogen {
                let number_of_hydrogens = character.to_digit(10).unwrap() as usize;
                for ith_hydrogen in 1..=number_of_hydrogens {
                    self.atoms.push(Atom::new(1, None));
                    self.bonds.push((
                        self.current_atom_index,
                        self.current_atom_index + ith_hydrogen,
                        self.last_bond_type,
                    ));
                }
            }
            if self.is_negative_charge {
                self.atoms[self.current_atom_index].charge =
                    -(character.to_digit(10).unwrap() as i8);
            } else {
                self.atoms[self.current_atom_index].charge = character.to_digit(10).unwrap() as i8;
            }
        } else {
            if self.is_double_digit {
                if let Some(ring_number_value) = self.ring_number {
                    // TODO make this generic for any size
                    self.ring_number =
                        Some(ring_number_value * 10 + character.to_digit(10).unwrap() as usize);
                    self.is_double_digit = false;
                } else {
                    self.ring_number = Some(character.to_digit(10).unwrap() as usize);
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
        }
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

    for character in smiles.chars() {
        match character {
            'A'..='Z' => {
                // Handle the previous element
                parser.handle_atom(Some(character))?
            }

            'a'..='z' => {
                if parser.is_complex {
                    // If we are in a complex, then we need to check if the last character was a
                    // capital letter, if so, then we need to add the element to the buffer
                    if parser.element_buffer.chars().last().unwrap().is_uppercase() {
                        parser.handle_atom(None)?;
                    }
                } else {
                    // According to the SMILES specification, the lowercase letters are used to denote aromatic atoms if they are not complex
                    match character {
                        'b' | 'c' | 'n' | 'o' | 's' | 'p' => parser.handle_atom(Some(character))?,
                        'r' => {
                            if parser.element_buffer == "B" {
                                parser.element_buffer.push('R')
                            } else {
                                return Err(ParseError::SMILESComplexError);
                            }
                        }

                        'l' => {
                            if parser.element_buffer == "C" {
                                parser.element_buffer.push('L')
                            } else {
                                return Err(ParseError::SMILESComplexError);
                            }
                        }
                        _ => return Err(ParseError::SMILESComplexError),
                    }
                }
            }

            '1'..='9' => {
                parser.handle_number(character)?;
            }
            '-' => {
                if parser.is_complex {
                    parser.is_negative_charge = true;
                } else {
                    parser.bonds.push((
                        parser.current_atom_index,
                        parser.current_atom_index + 1,
                        BondType::Single,
                    ));
                }
            }

            '=' => {
                parser.last_bond_type = BondType::Double;
            }
            '#' => {
                parser.last_bond_type = BondType::Triple;
            }
            ':' => {
                parser.last_bond_type = BondType::Aromatic;
            }
            '(' => {
                parser.branch_stack.push(parser.current_atom_index);
            }
            ')' => parser.branch_exits += 1,
            '[' => {
                parser.is_complex = true;
            }
            ']' => {
                parser.is_complex = false;
            }
            '%' => {
                parser.is_double_digit = true;
            }
            _ => (),
        }
    }

    println!("bonds: {:?}", parser.bonds);
    if !parser.element_buffer.is_empty() {
        println!("Element buffer {:?}", parser.element_buffer);
        parser.handle_atom(None)?
    }

    for (start, end) in parser.ring_bonds.values() {
        if start.is_some() && end.is_some() {
            parser
                .bonds
                .push((start.unwrap(), end.unwrap(), BondType::Aromatic));
        }
    }

    println!("{:?}", parser.ring_bonds);
    println!("{:?}", parser.atoms);
    println!("{:?}", parser.bonds);

    for bond in parser.bonds.iter() {
        let atom1 = &mut parser.atoms[bond.0];
        atom1.add_bond(Bond::new(bond.1, bond.2));
        let atom2 = &mut parser.atoms[bond.1];
        atom2.add_bond(Bond::new(bond.0, bond.2));
    }

    Ok(Molecule::from_atoms(parser.atoms))
}
