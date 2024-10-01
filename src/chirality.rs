use crate::{
    molecule::{
        Molecule, Molecule3D
    },
    bond::{BondOrder,BondTarget},
    vector::Vector
};

/// Represents the chirality of a molecule.
/// 
/// The definitions of the chiral classes are based on the following source:
/// 
/// http://opensmiles.org/opensmiles.html
#[derive(Debug, Default, Clone, Copy, Eq, PartialEq)]
pub enum ChiralClass {
    TH(u8), // @ 1..=2 1 = anti-clockwise, 2 = clockwise
    AL(u8), // @ 1..=2 1 = anti-clockwise, 2 = clockwise
    SP(u8), // @ 1..=3 1 = U, 2 = Z, 3 = 4
    TB(u8), // @ 1..=20 
    OH(u8), // @ 1..=30
    Clockwise,
    Counterclockwise,
    #[default]
    None,
}
/// Represents the shape of an octahedral configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
enum OctahedralShape {
    U,   // U-shaped configuration
    Z,   // Z-shaped configuration
    Four, // 4-shaped configuration
    None,
}

/// Represents the winding direction of ligands in an octahedral configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
enum WindingDirection {
    Clockwise,
    Counterclockwise,
    None
}

impl WindingDirection {
    pub fn reverse(&self) -> WindingDirection {
        match self {
            WindingDirection::Clockwise => WindingDirection::Counterclockwise,
            WindingDirection::Counterclockwise => WindingDirection::Clockwise,
            WindingDirection::None => WindingDirection::None,
        }
    }
}

/// Represents the direction of the viewing axis in an octahedral configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
enum AxisDirection {
    AtoF,
    AtoE,
    AtoD,
    AtoC,
    AtoB,
    None
}

/// Represents the orientation of the viewing axis in an octahedral configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
enum AxisOrientation {
    Normal,
    Reversed,
    None
}

impl Molecule3D {
    pub fn identify_chiral_classes(&mut self) {
        if self.chiral_classes.is_some() {
            return;
        }
        let mut chiral_classes = vec![ChiralClass::None; self.len()];

        for (index, &atomic_number) in self.atomic_numbers().iter().enumerate() {
            // If we already know the chirality, skip it
            if chiral_classes[index] != ChiralClass::None {
                continue;
            }
            let Some(neighbors) = self.get_atom_bonds(index) else {
                continue;
            };

            chiral_classes[index] = if self.is_potential_tetrahedral_center(atomic_number, neighbors) {
                let tetrahedral_chirality = self.determine_tetrahedral_chirality(index, neighbors);
                if tetrahedral_chirality != ChiralClass::None {
                    tetrahedral_chirality
                } else {
                    self.determine_square_planar_chirality(index, neighbors)
                }
            } else if self.is_tetrahedral_allene(atomic_number, neighbors) {
                self.determine_allene_chirality(neighbors).unwrap_or(ChiralClass::None)
            } 
            else if self.is_potential_cis_trans(neighbors) {
                self.determine_cis_trans_chirality(index, neighbors, &mut chiral_classes)
            }
            else if self.is_potential_trigonal_bipyramidal(neighbors) {
                self.determine_trigonal_bipyramidal_chirality(index, neighbors)
            } 
            else if self.is_potential_octahedral(neighbors) {
                self.determine_octahedral_chirality(index, neighbors)
            } else {
                ChiralClass::None
            };
        }

        self.chiral_classes = Some(chiral_classes);
    }

    /// Determines the cis/trans chirality of a double bond between two atoms and assigns
    /// the appropriate chiral classes using atom-based specification (`@` and `@@`).
    ///
    /// This method updates the `chiral_classes` array for both the central atom and its
    /// double bond neighbor based on their spatial configuration. It ensures that each
    /// chirality assignment is atom-specific, avoiding the complexities of bond-based
    /// directional labels.
    ///
    /// # Arguments
    ///
    /// * `center_atom_index` - The index of the central atom in the double bond.
    /// * `center_neighbors` - A slice of `BondTarget` representing the neighboring atoms
    ///   bonded to the central atom.
    /// * `chiral_classes` - A mutable slice of `ChiralClass` where the determined
    ///   chirality will be stored.
    ///
    /// # Returns
    ///
    /// Returns the `ChiralClass` of the central atom after determination.
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::molecule::{Molecule3D, Molecule};
    /// use molecules::bond::BondTarget;
    /// use molecules::chirality::ChiralClass;
    /// 
    /// let mut molecule = Molecule3D::from_sdf("tests/cis_2_butene.sdf").unwrap()[0].clone();
    /// molecule.identify_chiral_classes();
    /// let chiral_classes = molecule.chiral_classes().unwrap();
    /// println!("{:?}", chiral_classes);
    /// assert!(chiral_classes[0] == ChiralClass::Counterclockwise);
    /// assert!(chiral_classes[1] == ChiralClass::Clockwise);
    /// ```
    pub fn determine_cis_trans_chirality(
        &self,
        center_atom_index: usize,
        center_neighbors: &[BondTarget],
        chiral_classes: &mut [ChiralClass],
    ) -> ChiralClass {
        // Identify the neighbor atom involved in the double bond
        let double_bond_neighbor = center_neighbors
            .iter()
            .find(|&bond| bond.bond_order() == BondOrder::Double)
            .map(|bond| bond.target());

        let double_bond_neighbor_index = match double_bond_neighbor {
            Some(index) => index,
            None => return ChiralClass::None,
        };

        // Retrieve the bonds of the double bond neighbor
        let double_bond_neighbor_bonds = self.get_atom_bonds(double_bond_neighbor_index)
            .unwrap_or_default();

        // Validate that both the central atom and its double bond neighbor have exactly two single bonds
        if center_neighbors.len() != 3 || double_bond_neighbor_bonds.len() != 3 {
            return ChiralClass::None;
        }

        // Obtain the spatial positions of the central atom and its double bond neighbor
        let center_position = match self.get_atom_position(center_atom_index) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };
        let double_bond_neighbor_position = match self.get_atom_position(double_bond_neighbor_index) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };

        // Identify single bond neighbors excluding the double bond connection
        let center_single_bond_neighbors: Vec<_> = center_neighbors.iter()
            .filter(|&bond| bond.bond_order() == BondOrder::Single)
            .collect();
        let double_bond_single_neighbors: Vec<_> = double_bond_neighbor_bonds.iter()
            .filter(|&bond| bond.bond_order() == BondOrder::Single && bond.target() != center_atom_index)
            .collect();

        if center_single_bond_neighbors.len() != 2 || double_bond_single_neighbors.len() != 2 {
            return ChiralClass::None;
        }

        // Retrieve positions of all single bond neighbors
        let center_neighbor1_position = match self.get_atom_position(center_single_bond_neighbors[0].target()) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };
        let center_neighbor2_position = match self.get_atom_position(center_single_bond_neighbors[1].target()) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };
        let double_bond_neighbor1_position = match self.get_atom_position(double_bond_single_neighbors[0].target()) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };
        let double_bond_neighbor2_position = match self.get_atom_position(double_bond_single_neighbors[1].target()) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };

        // Compute vectors from central atoms to their single bond neighbors
        let center_vector1 = center_neighbor1_position - center_position;
        let center_vector2 = center_neighbor2_position - center_position;
        let double_bond_vector1 = double_bond_neighbor1_position - double_bond_neighbor_position;
        let double_bond_vector2 = double_bond_neighbor2_position - double_bond_neighbor_position;
        let double_bond_axis = double_bond_neighbor_position - center_position;

        // Calculate the perpendicular axis to the double bond plane
        let perpendicular_axis = double_bond_axis.cross(&center_vector1);

        // Determine chirality of the central atom based on the orientation of its neighbors
        let central_cross = center_vector1.cross(&center_vector2);
        let central_chirality = if central_cross.dot(&perpendicular_axis) > 0.0 {
            ChiralClass::Counterclockwise
        } else {
            ChiralClass::Clockwise
        };

        // Determine chirality of the double bond neighbor atom based on the orientation of its neighbors
        let neighbor_cross = double_bond_vector1.cross(&double_bond_vector2);
        let neighbor_chirality = if neighbor_cross.dot(&perpendicular_axis) > 0.0 {
            ChiralClass::Counterclockwise
        } else {
            ChiralClass::Clockwise
        };

        // Assign the determined chirality to both atoms
        chiral_classes[center_atom_index] = central_chirality;
        chiral_classes[double_bond_neighbor_index] = neighbor_chirality;

        // Return the chirality of the central atom
        central_chirality
    }
    
    /// Determines the chirality of a square planar structure.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the central atom.
    /// * `neighbors` - A slice of `BondTarget` representing the neighboring atoms.
    ///
    /// # Returns
    ///
    /// A `ChiralClass` indicating the chirality of the square planar center.
    ///
    /// This geometric analysis allows us to classify the square planar structure into
    /// one of three chiral arrangements: U-shaped, Z-shaped, or 4-shaped.
    fn determine_square_planar_chirality(&self, index: usize, neighbors: &[BondTarget]) -> ChiralClass {
        if neighbors.len() != 4 {
            return ChiralClass::None;
        }

        let self_position = match self.get_atom_position(index) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };

        let positions: Vec<Vector> = neighbors
            .iter()
            .flat_map(|&neighbor| self.get_atom_position(neighbor.target()))
            .collect();

        if positions.len() != 4 {
            return ChiralClass::None;
        }

        // Calculate vectors from central atom to neighbors
        let vectors: Vec<Vector> = positions.iter().map(|&pos| pos - self_position).collect();

        // Calculate cross products to determine chirality
        let cross1: Vector = vectors[0].cross(&vectors[1]);
        let cross2: Vector = vectors[1].cross(&vectors[2]);
        let cross3: Vector = vectors[2].cross(&vectors[3]);
        let cross4: Vector = vectors[3].cross(&vectors[0]);
        
        // This is used for Z and 4
        let magnitude1 = cross1.magnitude();

        // Check if all cross products are pointing in the same direction
        let dot1 = cross1.dot(&cross2);
        let dot2 = cross2.dot(&cross3);
        let dot3 = cross3.dot(&cross4);
        let dot4 = cross4.dot(&cross1);

        match (dot1.abs(), dot2.abs(), dot3.abs(), dot4.abs()) {
            (0.1..f64::MAX, 0.1..f64::MAX, 0.1..f64::MAX, 0.1..f64::MAX) => ChiralClass::SP(1), // U
            (-0.1..0.1, -0.1..0.1, -0.1..0.1, -0.1..0.1) => {
                if magnitude1 > 0.1 {
                    ChiralClass::SP(3) // Z
                } else {
                    ChiralClass::SP(2) // 4
                }
            }
            _ => ChiralClass::None,
        }
    }

    /// Determines the chirality of a trigonal bipyramidal structure.  /// /// # Arguments ///
    /// * `index` - The index of the central atom.
    /// * `neighbors` - A slice of `BondTarget` representing the neighboring atoms.
    ///
    /// # Returns
    ///
    /// A `ChiralClass` indicating the chirality of the trigonal bipyramidal center.
    fn determine_trigonal_bipyramidal_chirality(&self, index: usize, neighbors: &[BondTarget]) -> ChiralClass {
        if neighbors.len() != 5 {
            return ChiralClass::None;
        }

        let central_position = match self.get_atom_position(index) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };

        let neighbor_positions: Vec<(Vector, usize)> = neighbors
            .iter()
            .map(|&neighbor| (self.get_atom_position(neighbor.target()).unwrap(), neighbor.target()))
            .collect();

        if neighbor_positions.len() != 5 {
            return ChiralClass::None;
        }

        // Calculate vectors from central atom to neighbors
        let vectors_and_indices: Vec<(Vector, usize)> = neighbor_positions.iter()
            .map(|&(pos, index)| (pos - central_position, index))
            .collect();

        // Identify axial and equatorial ligands
        let (axial_indices, equatorial_indices) = self.identify_axial_equatorial(&vectors_and_indices);
        let is_ordered = axial_indices[0].1 < axial_indices[1].1;
        if axial_indices.len() != 2 || equatorial_indices.len() != 3 {
            return ChiralClass::None;
        }
        let axial_indices = [axial_indices[0].0, axial_indices[1].0];
        let equatorial_indices = [equatorial_indices[0], equatorial_indices[1], equatorial_indices[2]];
        // Is ordering clockwise?
        let mut is_clockwise = equatorial_indices[0].1 < equatorial_indices[1].1 && equatorial_indices[1].1 < equatorial_indices[2].1;
        if is_ordered {
            is_clockwise = !is_clockwise;
        }

        // Determine chirality based on the viewing axis and order
        match (&axial_indices, is_clockwise) {
            ([0, 4], false) => ChiralClass::TB(1),  // a -> e, (b, c, d) anti-clockwise
            ([0, 4], true) => ChiralClass::TB(2),  // a -> e, (b, c, d) clockwise
            ([0, 3], false) => ChiralClass::TB(3),  // a -> d, (b, c, e) anti-clockwise
            ([0, 3], true) => ChiralClass::TB(4),  // a -> d, (b, c, e) clockwise
            ([0, 2], false) => ChiralClass::TB(5),  // a -> c, (b, d, e) anti-clockwise
            ([0, 2], true) => ChiralClass::TB(6),  // a -> c, (b, d, e) clockwise
            ([0, 1], false) => ChiralClass::TB(7),  // a -> b, (c, d, e) anti-clockwise
            ([0, 1], true) => ChiralClass::TB(8),  // a -> b, (c, d, e) clockwise
            ([1, 4], false) => ChiralClass::TB(9),  // b -> e, (a, c, d) anti-clockwise
            ([1, 4], true) => ChiralClass::TB(11),  // b -> e, (a, c, d) clockwise
            ([1, 3], false) => ChiralClass::TB(10),  // b -> d, (a, c, e) anti-clockwise
            ([1, 3], true) => ChiralClass::TB(12),  // b -> d, (a, c, e) clockwise
            ([1, 2], false) => ChiralClass::TB(13),  // b -> c, (a, d, e) anti-clockwise
            ([1, 2], true) => ChiralClass::TB(14),  // b -> c, (a, d, e) clockwise
            ([2, 4], false) => ChiralClass::TB(15),  // c -> e, (a, b, d) anti-clockwise
            ([2, 4], true) => ChiralClass::TB(20),  // c -> e, (a, b, d) clockwise
            ([2, 3], false) => ChiralClass::TB(16),  // c -> d, (a, b, e) anti-clockwise
            ([2, 3], true) => ChiralClass::TB(19),  // c -> d, (a, b, e) clockwise
            ([3, 4], false) => ChiralClass::TB(17),  // d -> e, (a, b, c) anti-clockwise
            ([3, 4], true) => ChiralClass::TB(18),  // d -> e, (a, b, c) clockwise
            _ => ChiralClass::None,
        }
    }

    /// Identifies axial and equatorial ligands in a trigonal bipyramidal structure.
    ///
    /// # Arguments
    ///
    /// * `vectors` - A slice of `Vector` representing the vectors from the central atom to its neighbors.
    ///
    /// # Returns
    ///
    /// A tuple containing two vectors of indices:
    /// - The first vector contains indices of the two axial ligands.
    /// - The second vector contains indices of the three equatorial ligands.
    fn identify_axial_equatorial(&self, vectors_and_indices: &[(Vector, usize)]) -> (Vec<(usize, usize)>, Vec<(usize, usize)>) {
        let mut axial_indices = Vec::new();
        let mut equatorial_indices = Vec::new();

        for (i, &v1) in vectors_and_indices.iter().enumerate() {
            for (j, &v2) in vectors_and_indices.iter().enumerate() {
                if i < j {
                    if v1.0.is_opposite(&v2.0, 0.1) {
                        // If two vectors are nearly opposite (180 degrees apart), they form the axis
                        axial_indices.push((i, vectors_and_indices[i].1));
                        axial_indices.push((j, vectors_and_indices[j].1));
                    }
                }
            }
        }
        for (i, &v1) in vectors_and_indices.iter().enumerate() {
            if !axial_indices.contains(&(i, v1.1)) {
                equatorial_indices.push((i, v1.1));
            }
        }

        (axial_indices, equatorial_indices)
    }

    /// Determines the chirality of an octahedral structure.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the central atom in the molecule.
    /// * `neighbors` - A slice of `BondTarget` representing the neighboring atoms (ligands).
    ///
    /// # Returns
    ///
    /// A `ChiralClass` indicating the chirality of the octahedral center:
    /// - `ChiralClass::OH(n)` where n is 1-30, representing different octahedral configurations
    /// - `ChiralClass::None` if not chiral or indeterminate
    ///
    /// # Algorithm
    ///
    /// 1. Checks if there are exactly 6 neighbors (ligands).
    /// 2. Retrieves the positions of the central atom and its neighbors.
    /// 3. Identifies the shape and viewing axis based on the ligand positions.
    /// 4. Determines the winding direction (clockwise or anticlockwise) of the relevant ligands.
    /// 5. Assigns the appropriate OH number based on the shape, axis, and winding direction.
    fn determine_octahedral_chirality(&self, index: usize, neighbors: &[BondTarget]) -> ChiralClass {
        if neighbors.len() != 6 {
            return ChiralClass::None;
        }

        let central_position = match self.get_atom_position(index) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };

        let neighbor_positions: Vec<(Vector,usize)> = neighbors
            .iter()
            .map(|&neighbor| (self.get_atom_position(neighbor.target()).unwrap(),neighbor.target()))
            .collect();

        if neighbor_positions.len() != 6 {
            return ChiralClass::None;
        }

        let vectors: Vec<(Vector,usize)> = neighbor_positions.iter()
            .map(|&(pos,index)| (pos - central_position, index))
            .collect();

        let ((shape, winding), viewing_axis, relevant_ligands) = self.identify_octahedral_configuration(&vectors);
        let oh_number = self.assign_oh_number(shape, &viewing_axis, winding);

        ChiralClass::OH(oh_number)
    }

    /// Identifies the octahedral configuration (shape, viewing axis, and relevant ligands).
    fn identify_octahedral_configuration(&self, vectors: &[(Vector, usize)]) -> ((OctahedralShape, WindingDirection), (Vector, usize, usize), Vec<Vector>) {
        let pairs = self.find_opposite_pairs(vectors);
        let mut viewing_axis = (pairs[0].0.0 - pairs[0].1.0, pairs[0].0.1, pairs[0].1.1);
        if vectors[viewing_axis.1].1 > vectors[viewing_axis.2].1 {
            viewing_axis = (pairs[0].0.0, pairs[0].1.1, pairs[0].0.1);
        }

        let relevant_ligands: Vec<Vector> = vectors.iter()
            .enumerate()
            .filter(|&(index, _v)| index != viewing_axis.1 && index != viewing_axis.2)
            .map(|(_, vector)| vector.0)
            .collect();
        
        let shape = self.determine_shape(&relevant_ligands[0], &relevant_ligands[1], &relevant_ligands[2], &relevant_ligands[3], &viewing_axis.0);

        (shape, viewing_axis, relevant_ligands)
    }

    /// Determines the shape of the octahedral configuration.
    fn determine_shape(&self, v1: &Vector, v2: &Vector, v3: &Vector, v4: &Vector, axis: &Vector) -> (OctahedralShape, WindingDirection) {
        let cross_product1 = v1.cross(v2).dot(axis);
        let cross_product2 = v2.cross(v3).dot(axis);
        let cross_product3 = v3.cross(v4).dot(axis);
        match (cross_product1, cross_product2, cross_product3) {
            (f64::MIN..=-0.1, f64::MIN..=-0.1, f64::MIN..=-0.1) => (OctahedralShape::U, WindingDirection::Counterclockwise),
            (0.1..f64::MAX, 0.1..f64::MAX, 0.1..f64::MAX) => (OctahedralShape::U, WindingDirection::Clockwise),
            (f64::MIN..=-0.1, -0.1..0.1, 0.1..f64::MAX) => (OctahedralShape::Z, WindingDirection::Counterclockwise),
            (0.1..f64::MAX, -0.1..0.1, f64::MIN..=-0.1) => (OctahedralShape::Z, WindingDirection::Clockwise),
            (-0.1..0.1, f64::MIN..=-0.1, -0.1..0.1) => (OctahedralShape::Four, WindingDirection::Counterclockwise),
            (-0.1..0.1, 0.1..f64::MAX, -0.1..0.1) => (OctahedralShape::Four, WindingDirection::Clockwise),
            _ => (OctahedralShape::None, WindingDirection::None),
        }
    }

    /// Assigns the OH number based on shape, viewing axis, and winding direction.
    fn assign_oh_number(&self, shape: OctahedralShape, viewing_axis: &(Vector, usize, usize), winding: WindingDirection) -> u8 {
        let (axis_direction, axis_orientation) = self.determine_axis_direction(&viewing_axis);
        let mut winding = winding;
        if axis_orientation == AxisOrientation::Reversed {
            winding = winding.reverse();
        }
        match (shape, axis_direction, winding) {
            (OctahedralShape::U, AxisDirection::AtoF, WindingDirection::Counterclockwise) => 1,
            (OctahedralShape::U, AxisDirection::AtoF, WindingDirection::Clockwise) => 2,
            (OctahedralShape::U, AxisDirection::AtoE, WindingDirection::Counterclockwise) => 3,
            (OctahedralShape::U, AxisDirection::AtoE, WindingDirection::Clockwise) => 16,
            (OctahedralShape::U, AxisDirection::AtoD, WindingDirection::Counterclockwise) => 6,
            (OctahedralShape::U, AxisDirection::AtoD, WindingDirection::Clockwise) => 18,
            (OctahedralShape::U, AxisDirection::AtoC, WindingDirection::Counterclockwise) => 19,
            (OctahedralShape::U, AxisDirection::AtoC, WindingDirection::Clockwise) => 24,
            (OctahedralShape::U, AxisDirection::AtoB, WindingDirection::Counterclockwise) => 25,
            (OctahedralShape::U, AxisDirection::AtoB, WindingDirection::Clockwise) => 30,
            (OctahedralShape::Z, AxisDirection::AtoF, WindingDirection::Counterclockwise) => 4,
            (OctahedralShape::Z, AxisDirection::AtoF, WindingDirection::Clockwise) => 14,
            (OctahedralShape::Z, AxisDirection::AtoE, WindingDirection::Counterclockwise) => 5,
            (OctahedralShape::Z, AxisDirection::AtoE, WindingDirection::Clockwise) => 15,
            (OctahedralShape::Z, AxisDirection::AtoD, WindingDirection::Counterclockwise) => 7,
            (OctahedralShape::Z, AxisDirection::AtoD, WindingDirection::Clockwise) => 17,
            (OctahedralShape::Z, AxisDirection::AtoC, WindingDirection::Counterclockwise) => 20,
            (OctahedralShape::Z, AxisDirection::AtoC, WindingDirection::Clockwise) => 23,
            (OctahedralShape::Z, AxisDirection::AtoB, WindingDirection::Counterclockwise) => 26,
            (OctahedralShape::Z, AxisDirection::AtoB, WindingDirection::Clockwise) => 29,
            (OctahedralShape::Four, AxisDirection::AtoF, WindingDirection::Counterclockwise) => 10,
            (OctahedralShape::Four, AxisDirection::AtoF, WindingDirection::Clockwise) => 8,
            (OctahedralShape::Four, AxisDirection::AtoE, WindingDirection::Counterclockwise) => 11,
            (OctahedralShape::Four, AxisDirection::AtoE, WindingDirection::Clockwise) => 9,
            (OctahedralShape::Four, AxisDirection::AtoD, WindingDirection::Counterclockwise) => 13,
            (OctahedralShape::Four, AxisDirection::AtoD, WindingDirection::Clockwise) => 12,
            (OctahedralShape::Four, AxisDirection::AtoC, WindingDirection::Counterclockwise) => 22,
            (OctahedralShape::Four, AxisDirection::AtoC, WindingDirection::Clockwise) => 21,
            (OctahedralShape::Four, AxisDirection::AtoB, WindingDirection::Counterclockwise) => 28,
            (OctahedralShape::Four, AxisDirection::AtoB, WindingDirection::Clockwise) => 27,
            _ => 0,
        }
    }

    /// Determines the direction of the viewing axis.
    fn determine_axis_direction(&self, viewing_axis: &(Vector, usize, usize)) -> (AxisDirection, AxisOrientation) {
        match (viewing_axis.1 , viewing_axis.2) {
            (0, 5) => (AxisDirection::AtoF, AxisOrientation::Normal),
            (0, 4) => (AxisDirection::AtoE, AxisOrientation::Normal),
            (0, 3) => (AxisDirection::AtoD, AxisOrientation::Normal),
            (0, 2) => (AxisDirection::AtoC, AxisOrientation::Normal),
            (0, 1) => (AxisDirection::AtoB, AxisOrientation::Normal),
            (5, 0) => (AxisDirection::AtoF, AxisOrientation::Reversed),
            (4, 0) => (AxisDirection::AtoE, AxisOrientation::Reversed),
            (3, 0) => (AxisDirection::AtoD, AxisOrientation::Reversed),
            (2, 0) => (AxisDirection::AtoC, AxisOrientation::Reversed),
            (1, 0) => (AxisDirection::AtoB, AxisOrientation::Reversed),
            _ => (AxisDirection::None, AxisOrientation::None),
        }
    }

    /// Finds opposite pairs of ligands in an octahedral structure.
    ///
    /// # Arguments
    ///
    /// * `vectors` - A slice of `Vector` representing the vectors from the central atom to its neighbors.
    ///
    /// # Returns
    ///
    /// A vector of tuples, each containing a pair of opposite vectors.
    fn find_opposite_pairs(&self, vectors: &[(Vector,usize)]) -> Vec<((Vector,usize), (Vector,usize))> {
        let mut pairs = Vec::new();
        let tolerance = 0.1; // Tolerance for considering vectors as opposite

        for i in 0..vectors.len() {
            for j in i+1..vectors.len() {
                if vectors[i].0.is_opposite(&vectors[j].0, tolerance) {
                    pairs.push(((vectors[i].0, i), (vectors[j].0, j)));
                }
            }
        }

        pairs
    }

    /// Determines the chirality of a tetrahedral center.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the central atom.
    /// * `neighbors` - A slice of `BondTarget` representing the neighboring atoms.
    ///
    /// # Returns
    ///
    /// A `ChiralClass` indicating the chirality of the tetrahedral center.
    fn determine_tetrahedral_chirality(&self, index: usize, neighbors: &[BondTarget]) -> ChiralClass {
        if neighbors.len() != 4 {
            return ChiralClass::None;
        }

        let central_position = match self.get_atom_position(index) {
            Some(pos) => pos,
            None => return ChiralClass::None,
        };

        let neighbor_positions: Vec<Vector> = neighbors
            .iter()
            .flat_map(|&neighbor| self.get_atom_position(neighbor.target()))
            .collect();

        if neighbor_positions.len() != 4 {
            return ChiralClass::None;
        }

        // Calculate vectors from central atom to neighbors
        let vectors: Vec<Vector> = neighbor_positions.iter()
            .map(|&pos| pos - central_position)
            .collect();

        // Calculate the signed volume of the tetrahedron
        let signed_volume = vectors[0].dot(&vectors[1].cross(&vectors[2]));

        if signed_volume.abs() < f64::EPSILON {
            ChiralClass::None // Planar arrangement
        } else if signed_volume > 0.0 {
            ChiralClass::TH(2) // Right-handed (@)
        } else {
            ChiralClass::TH(1) // Left-handed (@@)
        }
    }

    /// Determines the chirality of an allene-like structure.
    ///
    /// # Arguments
    ///
    /// * `neighbors` - A slice of `BondTarget` representing the neighbors of the central atom.
    ///
    /// # Returns
    ///
    /// An `Option<ChiralClass>` indicating the chirality of the allene structure, if determinable.
    ///
    /// # Notes
    ///
    /// This function assumes that the central atom (index) is already known to be part of an
    /// allene-like structure. The `index` parameter has been removed as it was not used in the
    /// function body.
    fn determine_allene_chirality(&self, neighbors: &[BondTarget]) -> Option<ChiralClass> {
        let [left, right] = neighbors else { return None };
        
        let left_neighbors = self.get_atom_bonds(left.target())?;
        let right_neighbors = self.get_atom_bonds(right.target())?;

        if left_neighbors.len() != 2 || right_neighbors.len() != 2 {
            return None;
        }

        let left_plane = self.calculate_plane(left.target(), left_neighbors)?;
        let right_plane = self.calculate_plane(right.target(), right_neighbors)?;

        let angle = left_plane.angle_between(&right_plane);

        Some(if angle > Some(0.0) { ChiralClass::AL(1) } else { ChiralClass::AL(2) })
    }

    /// Calculates the normal vector of the plane formed by a central atom and its two neighbors.
    ///
    /// # Arguments
    ///
    /// * `center` - The index of the central atom.
    /// * `neighbors` - A slice of `BondTarget` representing the neighboring atoms.
    ///
    /// # Returns
    ///
    /// * `Option<Vector>` - The normal vector of the plane if all positions are valid and there are exactly two neighbors, or `None` otherwise.
    ///
    /// # Safety
    ///
    /// This function checks that:
    /// - The `neighbors` slice contains exactly two elements.
    /// - All atom positions (center and neighbors) are valid.
    /// - The resulting vectors are not parallel (to avoid a zero cross product).
    fn calculate_plane(&self, center: usize, neighbors: &[BondTarget]) -> Option<Vector> {
        if neighbors.len() != 2 {
            return None;
        }

        let center_pos = self.get_atom_position(center)?;
        let pos1 = self.get_atom_position(neighbors[0].target())?;
        let pos2 = self.get_atom_position(neighbors[1].target())?;

        let vec1 = pos1 - center_pos;
        let vec2 = pos2 - center_pos;

        // Check if vectors are not parallel (cross product is not zero)
        let cross_product = vec1.cross(&vec2);
        if cross_product.length() < f64::EPSILON {
            return None;
        }

        Some(cross_product)
    }

    /// Determines if an atom is a tetrahedral allene.
    ///
    /// # Arguments
    ///
    /// * `atomic_number` - The atomic number of the central atom.
    /// * `neighbors` - A slice of `BondTarget` representing the neighboring atoms.
    ///
    /// # Returns
    ///
    /// * `bool` - True if the atom is a tetrahedral allene, false otherwise.
    ///
    /// # Details
    ///
    /// A tetrahedral allene is a carbon atom (atomic number 6) with exactly two carbon neighbors,
    /// both connected by double bonds.
    fn is_tetrahedral_allene(&self, atomic_number: u8, neighbors: &[BondTarget]) -> bool {
        atomic_number == 6 // Is atom carbon 
            && neighbors.len() == 2 // Does it have two neighbors
            && neighbors
                    .iter()
                    .all(|neighbor| 
                        self.get_atomic_number(neighbor.target()) == 6 // Are all neighbors carbon 
                        && neighbor.bond_order() == BondOrder::Double  // Are all bonds double
                                                                       // bonds
                        )
    }

    /// Determines if an atom is a potential trigonal bipyramidal center.
    ///
    /// # Arguments
    ///
    /// * `neighbors` - A slice of `BondTarget` representing the neighboring atoms.
    ///
    /// # Returns
    ///
    /// * `bool` - True if the atom has exactly 5 neighbors, indicating a potential trigonal bipyramidal structure.
    fn is_potential_trigonal_bipyramidal(&self, neighbors: &[BondTarget]) -> bool {
        neighbors.len() == 5
    }

    fn is_potential_cis_trans(&self, neighbors: &[BondTarget]) -> bool {
        neighbors.len() == 3 
        && neighbors.iter().filter(|neighbor| neighbor.bond_order() == BondOrder::Double).count() == 1
        && neighbors.iter().filter(|neighbor| neighbor.bond_order() == BondOrder::Single).count() == 2
    }

    /// Determines if an atom is a potential octahedral center.
    ///
    /// # Arguments
    ///
    /// * `neighbors` - A slice of `BondTarget` representing the neighboring atoms.
    ///
    /// # Returns
    ///
    /// * `bool` - True if the atom has exactly 6 neighbors, indicating a potential octahedral structure.
    fn is_potential_octahedral(&self, neighbors: &[BondTarget]) -> bool {
        neighbors.len() == 6
    }

    fn is_potential_tetrahedral_center(&self, atomic_number: u8, neighbors: &[BondTarget]) -> bool {
        let number_of_hydrogens = neighbors
            .iter()
            .filter(|&neighbor| self.get_atomic_number(neighbor.target()) == 1)
            .count();

        neighbors.len() == 4
            && number_of_hydrogens < 2
            && atomic_number != 16 // Sulfur
            && atomic_number != 34 // Selenium
    }

    pub fn is_chiral(&self) -> bool {
        self.chiral_classes()
            .unwrap_or_default()
            .iter()
            .any(|&class| class != ChiralClass::None)
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::molecule::Molecule3D;
    use crate::chirality::ChiralClass;
    use crate::sdf::parse_sdf_file;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn test_is_chiral() {
        let mut mol = Molecule3D::from_xyz("tests/ethane.xyz");
        mol.identify_chiral_classes();
        assert!(!mol.is_chiral());
    }

    #[test]
    fn test_square_planar_chirality() {
        let mut cisplatin = Molecule3D::from_xyz("tests/cisplatin_sp2.xyz");
        cisplatin.identify_chiral_classes();
        assert_eq!(cisplatin.chiral_classes().unwrap()[0], ChiralClass::SP(3));

        let mut transplatin = Molecule3D::from_xyz("tests/transplatin_sp1.xyz");
        transplatin.identify_chiral_classes();
        assert_eq!(transplatin.chiral_classes().unwrap()[0], ChiralClass::SP(1));

        let mut transplatin = Molecule3D::from_xyz("tests/transplatin_sp2.xyz");
        transplatin.identify_chiral_classes();
        assert_eq!(transplatin.chiral_classes().unwrap()[0], ChiralClass::SP(2));

        let mut transplatin = Molecule3D::from_xyz("tests/transplatin_sp3.xyz");
        transplatin.identify_chiral_classes();
        assert_eq!(transplatin.chiral_classes().unwrap()[0], ChiralClass::SP(3));
    }

    #[test]
    fn test_trigonal_bipyramidal_chirality() {
        let file = File::open("tests/ordered_molecules.sdf").unwrap();
        let reader = BufReader::new(file);
        let sdf_entries = parse_sdf_file(reader).unwrap();
        let mut chiral_classes: Vec<ChiralClass> = Vec::new();
        for sdf_entry in sdf_entries.iter() {
            let chiral_permutation = sdf_entry.data_fields.get("ChiralPermutation").unwrap();
            let chiral_permutation = chiral_permutation.parse::<u8>().unwrap();
            let chiral_type = sdf_entry.data_fields.get("ChiralType").unwrap();
            let chiral_class = match (chiral_type.as_str(), chiral_permutation) {
                ("TH", 1..=2) => ChiralClass::TH(chiral_permutation),
                ("SP", 1..=3) => ChiralClass::SP(chiral_permutation),
                ("TB", 1..=20) => ChiralClass::TB(chiral_permutation),
                ("OH", 1..=30) => ChiralClass::OH(chiral_permutation),
                ("AL", 1..=2) => ChiralClass::AL(chiral_permutation),
                _ => ChiralClass::None,
            };
            chiral_classes.push(chiral_class);
        }
        let mut pf5: Vec<Molecule3D> = Molecule3D::from_sdf_entries(sdf_entries).unwrap();

        for (mol, chiral_class) in pf5.iter_mut().zip(chiral_classes) {
            mol.identify_chiral_classes();
            assert_eq!(mol.chiral_classes().unwrap()[0], chiral_class);
        }
    }

    #[test]
    fn test_cis_trans_chirality() {
        let molecules = Molecule3D::from_sdf("tests/cis_2_butene.sdf").unwrap();
        let mut cis_2_butene = molecules[0].clone();
        cis_2_butene.identify_chiral_classes();
        assert_eq!(cis_2_butene.chiral_classes().unwrap()[0], ChiralClass::Counterclockwise);
        assert_eq!(cis_2_butene.chiral_classes().unwrap()[1], ChiralClass::Clockwise);

        let molecules = Molecule3D::from_sdf("tests/trans_2_butene.sdf").unwrap();
        let mut trans_2_butene = molecules[0].clone();
        trans_2_butene.identify_chiral_classes();
        assert_eq!(trans_2_butene.chiral_classes().unwrap()[0], ChiralClass::Counterclockwise);
        assert_eq!(trans_2_butene.chiral_classes().unwrap()[1], ChiralClass::Counterclockwise);

        let molecules = Molecule3D::from_sdf("tests/pentene.sdf").unwrap();
        let mut pentene = molecules[0].clone();
        pentene.identify_chiral_classes();
        pentene.identify_chiral_classes();
        assert_eq!(pentene.chiral_classes().unwrap()[1], ChiralClass::Counterclockwise);
        assert_eq!(pentene.chiral_classes().unwrap()[3], ChiralClass::Counterclockwise);

    }
}