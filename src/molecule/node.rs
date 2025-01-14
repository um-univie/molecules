use crate::molecule::bond::BondOrder;

#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    index: usize,
    bond_type: BondOrder,
    children: Vec<Node>,
}

impl Node {
    pub fn new(index: usize, bond_type: BondOrder) -> Self {
        Node {
            index,
            bond_type,
            children: Vec::new(),
        }
    }

    pub fn add_child(&mut self, child: Node) {
        self.children.push(child);
    }

    pub fn bond_type(&self) -> BondOrder {
        self.bond_type
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn children(&self) -> &Vec<Node> {
        &self.children
    }
}
