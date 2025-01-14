impl SMILESParser {
    /// Adds all bonds to the atoms, correctly handling aromatic bonds.
    fn add_all_bonds(&mut self) {
        for (start, end) in self.ring_bonds.values() {
            // Debugging information
            println!("{:?}, {:?}", start, end);
            if let (Some(start_idx), Some(end_idx)) = (start, end) {
                let bond_order = if self.atoms[start_idx].aromatic && self.atoms[end_idx].aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                };
                self.bonds.push((start_idx, end_idx, bond_order));
            }
        }

        for (atom, number) in self.hydrogens.iter() {
            for _ in 0..*number {
                self.atoms.push(Atom::new(1));
                self.bonds
                    .push((*atom, self.atoms.len() - 1, BondOrder::Single));
            }
        }

        // Add bonds to the molecule
        for bond in self.bonds.iter() {
            let atom1 = &mut self.atoms[bond.0];
            atom1.add_bond(BondTarget::new(bond.1, bond.2));
            let atom2 = &mut self.atoms[bond.1];
            atom2.add_bond(BondTarget::new(bond.0, bond.2));
        }
    }

    fn handle_atom(&mut self, byte: Option<u8>) -> Result<(), ParseError> {
        if !self.element_buffer.is_empty() {
            let atomic_number = &self
                .element_buffer
                .to_uppercase()
                .as_str()
                .atomic_number()
                .ok_or(ParseError::ElementNotFound(self.element_buffer.to_owned()))?;

            let mut atom = Atom::new(*atomic_number).with_atom_class(self.current_atom_class);

            self.current_atom_index += 1;

            if byte.is_some() {
                self.handle_bond()?
            }

            if self.element_buffer.chars().next().unwrap().is_lowercase() || self.is_aromatic {
                atom.aromatic = true;
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

            self.atoms.push(atom);
            self.element_buffer.clear();
        }

        if let Some(character_byte) = byte {
            self.element_buffer.push(character_byte as char);
        }

        // Reset the aromatic flag after processing the atom
        self.is_aromatic = false;

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
            if bytes[*position].is_ascii_digit() {
                self.hydrogens
                    .insert(self.current_atom_index, byte_to_number(bytes[*position]));
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
                self.current_atom_class = Some(byte_to_number(bytes[*position]));
            } else {
                return Err(ParseError::AtomClassMissing);
            }
        }
        self.handle_bond()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_aromatic_bonds() {
        let molecules = SMILESParser::parse_smiles("c1ccccc1").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 6);
        assert_eq!(molecules[0].get_edges().len(), 6);
        for bond in molecules[0].get_edges() {
            assert_eq!(bond.order, BondOrder::Aromatic);
        }
    }

    #[test]
    fn test_parse_mixed_bonds_with_aromatic_atoms() {
        let molecules = SMILESParser::parse_smiles("C1=CC=CC=C1").unwrap();
        assert_eq!(molecules[0].atomic_numbers.len(), 6);
        assert_eq!(molecules[0].get_edges().len(), 6);
        // Depending on implementation, could assert specific bond orders
        // Here, ensure at least one aromatic bond exists
        assert!(molecules[0].get_edges().iter().any(|bond| bond.order == BondOrder::Aromatic));
    }
} 