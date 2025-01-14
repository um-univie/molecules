use crate::io::{ParseError, SMILESParser};
use crate::{
    atom::Atom,
    chirality::ChiralClass,
    molecule::bond::{BondOrder, BondTarget},
    molecule::Molecule,
};
use chemistry_consts::ElementProperties;
use tinyvec::ArrayVec;

/// This function   
#[derive(Debug, Default, Clone)]
pub struct Molecule2D {
    pub atomic_numbers: Vec<u8>,
    pub atom_classes: Option<Vec<u8>>,
    pub charges: Vec<i8>,
    pub chiral_classes: Option<Vec<ChiralClass>>,
    pub isotopes: Option<Vec<u16>>,
    radical_states: Vec<bool>,
    atom_bonds: Vec<ArrayVec<[BondTarget; 10]>>,
}

impl Molecule for Molecule2D {
    fn atom_bonds(&self) -> &Vec<ArrayVec<[BondTarget; 10]>> {
        &self.atom_bonds
    }
    fn atom_bonds_mut(&mut self) -> &mut Vec<ArrayVec<[BondTarget; 10]>> {
        &mut self.atom_bonds
    }
    fn atom_classes(&self) -> &Option<Vec<u8>> {
        &self.atom_classes
    }
    fn atom_classes_mut(&mut self) -> &mut Option<Vec<u8>> {
        &mut self.atom_classes
    }
    fn charges(&self) -> &[i8] {
        &self.charges
    }
    fn charges_mut(&mut self) -> &mut Vec<i8> {
        &mut self.charges
    }
    fn atomic_numbers(&self) -> &[u8] {
        &self.atomic_numbers
    }
    fn atomic_numbers_mut(&mut self) -> &mut Vec<u8> {
        &mut self.atomic_numbers
    }
    fn chiral_classes(&self) -> Option<&[ChiralClass]> {
        self.chiral_classes.as_deref()
    }
    fn chiral_classes_mut(&mut self) -> Option<&mut Vec<ChiralClass>> {
        self.chiral_classes.as_mut()
    }
    fn isotopes(&self) -> Option<&Vec<u16>> {
        self.isotopes.as_ref()
    }
    fn radical_states(&self) -> &[bool] {
        &self.radical_states
    }
    fn radical_states_mut(&mut self) -> &mut Vec<bool> {
        &mut self.radical_states
    }

    fn is_atom_radical(&self, atom_index: usize) -> bool {
        self.radical_states[atom_index]
    }

    fn set_atom_radical(&mut self, atom_index: usize, is_radical: bool) {
        self.radical_states[atom_index] = is_radical;
    }

    fn from_atoms(atoms: Vec<Atom>) -> Self {
        let atomic_numbers = atoms.iter().map(|atom| atom.atomic_number).collect();
        let charges = atoms.iter().map(|atom| atom.charge).collect();
        let radical_states = atoms.iter().map(|atom| atom.is_radical).collect();
        let atom_bonds = atoms.iter().map(|atom| atom.bonds).collect();
        let mut isotopes = None;
        let mut atom_classes = None;
        let mut chiral_classes = None;

        if atoms.iter().any(|atom| atom.isotope.is_some()) {
            isotopes = Some(
                atoms
                    .iter()
                    .map(|atom| {
                        if atom.isotope.is_some() {
                            atom.isotope.unwrap()
                        } else {
                            atom.atomic_number
                                .isotopes()
                                .unwrap()
                                .next()
                                .unwrap()
                                .mass
                                .round() as u16
                        }
                    })
                    .collect(),
            );
        }
        if atoms.iter().any(|atom| atom.atom_class.is_some()) {
            atom_classes = Some(
                atoms
                    .iter()
                    .map(|atom| atom.atom_class.unwrap_or(0))
                    .collect(),
            );
        }
        if atoms
            .iter()
            .any(|atom| atom.chiral_class != ChiralClass::None)
        {
            chiral_classes = Some(atoms.into_iter().map(|atom| atom.chiral_class).collect());
        }

        Molecule2D {
            atomic_numbers,
            charges,
            radical_states,
            atom_bonds,
            isotopes,
            chiral_classes,
            atom_classes,
        }
    }
}

impl Molecule2D {
    pub fn from_smiles(smiles: &str) -> Result<Vec<Self>, ParseError> {
        let parser = SMILESParser::parse_smiles(smiles)?;
        Ok(parser)
    }

    fn count_atom_pi_electrons(&self, atom_idx: usize) -> usize {
        // Count pi electrons from multiple bonds
        let pi_electrons = self.atom_bonds()[atom_idx]
            .iter()
            .map(|bond| match bond.bond_order() {
                BondOrder::Double => 1, // One pi electron per double bond
                BondOrder::Triple => 2, // Two pi electrons per triple bond
                _ => 0,
            })
            .sum();

        // Add electrons from lone pairs for certain atoms
        let atomic_number = self.atomic_numbers()[atom_idx];
        match atomic_number {
            7 => {
                // Nitrogen
                if self.get_atom_charge(atom_idx) == 0 && pi_electrons == 0 {
                    2 // Neutral N with no pi bonds contributes lone pair
                } else {
                    pi_electrons
                }
            }
            8 | 16 => {
                // Oxygen or Sulfur
                if self.get_atom_charge(atom_idx) == 0 && pi_electrons == 0 {
                    2 // Neutral O/S with no pi bonds contributes lone pair
                } else {
                    pi_electrons
                }
            }
            _ => pi_electrons,
        }
    }
}
