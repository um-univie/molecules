
            if bytes[position] == b'H' {
                let number_of_hydrogens = (byte as char).to_digit(10).unwrap() as usize;
                for ith_hydrogen in 1..=number_of_hydrogens {
                    self.atoms.push(Atom::new(1, None));
                    self.bonds.push((
                        self.current_atom_index,
                        self.current_atom_index + ith_hydrogen,
                        self.last_bond_type,
                    ));
            }

            if self.is_negative_charge {
                self.atoms[self.current_atom_index].charge = -(byte_to_number(byte) as i8);
            } else {
                self.atoms[self.current_atom_index].charge = byte_to_number(byte) as i8;
            }
                if self.is_explicit_hydrogen {
                    self.atoms.push(Atom::new(1, None));
                    self.bonds.push((
                        self.current_atom_index,
                        self.current_atom_index + 1,
                        self.last_bond_type,
                    ));
                }
            }
