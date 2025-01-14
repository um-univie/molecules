use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondAngle {
    angle: f64,
    atoms: (usize, usize, usize),
}

impl BondAngle {
    pub fn new(angle: f64, atoms: (usize, usize, usize)) -> Self {
        Self { angle, atoms }
    }
    pub fn atoms(&self) -> (usize, usize, usize) {
        self.atoms
    }
    pub fn angle(&self) -> f64 {
        self.angle
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum BondState {
    Formed,
    #[default]
    Broken,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct BondChange {
    atom1_index: usize,
    atom2_index: usize,
    bond_state: BondState,
}

impl Display for BondChange {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self.bond_state {
            BondState::Formed => write!(
                f,
                "Bond formed: {} and {}",
                self.atom1_index, self.atom2_index
            ),
            BondState::Broken => write!(
                f,
                "Bond broken: {} and {}",
                self.atom1_index, self.atom2_index
            ),
        }
    }
}

impl BondChange {
    pub fn broken(atom1_index: usize, atom2_index: usize) -> Self {
        BondChange {
            atom1_index,
            atom2_index,
            bond_state: BondState::Broken,
        }
    }
    pub fn formed(atom1_index: usize, atom2_index: usize) -> Self {
        BondChange {
            atom1_index,
            atom2_index,
            bond_state: BondState::Formed,
        }
    }
    pub fn is_broken(&self) -> bool {
        matches!(self.bond_state, BondState::Broken)
    }
    pub fn is_formed(&self) -> bool {
        matches!(self.bond_state, BondState::Formed)
    }
    pub fn atom1_index(&self) -> usize {
        self.atom1_index
    }
    pub fn atom2_index(&self) -> usize {
        self.atom2_index
    }
    pub fn bond_state(&self) -> &BondState {
        &self.bond_state
    }
}

#[derive(Debug, Default, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub struct BondTarget {
    pub target: usize,
    pub bond_order: BondOrder,
}

impl BondTarget {
    pub fn new(target: usize, bond_type: BondOrder) -> BondTarget {
        BondTarget {
            target,
            bond_order: bond_type,
        }
    }
    pub fn single(target: usize) -> BondTarget {
        BondTarget {
            target,
            bond_order: BondOrder::Single,
        }
    }
    pub fn double(target: usize) -> BondTarget {
        BondTarget {
            target,
            bond_order: BondOrder::Double,
        }
    }
    pub fn triple(target: usize) -> BondTarget {
        BondTarget {
            target,
            bond_order: BondOrder::Triple,
        }
    }
    fn aromatic(target: usize) -> BondTarget {
        BondTarget {
            target,
            bond_order: BondOrder::Aromatic,
        }
    }
    pub fn bond_order(&self) -> BondOrder {
        self.bond_order
    }
    pub fn target(&self) -> usize {
        self.target
    }
    pub fn order(&self) -> u8 {
        match self.bond_order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            // TODO: Implement aromatic bond order
            BondOrder::Aromatic => 1,
            BondOrder::Coordinate => 1,
        }
    }
}

#[derive(Debug, Default, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub enum BondOrder {
    #[default]
    Single,
    Double,
    Triple,
    Aromatic,
    Quadruple,
    Coordinate,
}

impl BondOrder {
    pub fn from_sdf_str(slice: &str) -> Result<BondOrder, String> {
        match slice.trim().to_lowercase().as_str() {
            "1" => Ok(BondOrder::Single),
            "2" => Ok(BondOrder::Double),
            "3" => Ok(BondOrder::Triple),
            "4" => Ok(BondOrder::Quadruple),
            _ => Err(format!("Invalid bond order: {}", slice)),
        }
    }
    pub fn to_u8(&self) -> u8 {
        match self {
            BondOrder::Single => 119,
            BondOrder::Double => 120,
            BondOrder::Triple => 121,
            BondOrder::Aromatic => 122,
            BondOrder::Quadruple => 123,
            BondOrder::Coordinate => 124,
        }
    }   
}

impl Display for BondOrder {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            BondOrder::Single => write!(f, "-"),
            BondOrder::Double => write!(f, "="),
            BondOrder::Triple => write!(f, "#"),
            BondOrder::Aromatic => write!(f, ":"),
            BondOrder::Quadruple => write!(f, "$"),
            // To be implemented
            BondOrder::Coordinate => write!(f, ""),
        }
    }
}

#[derive(Debug)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub bond_order: BondOrder,
}
