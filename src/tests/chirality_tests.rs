#[cfg(test)]
mod tests {
    use crate::chirality::{ChiralClass, ChiralClassifier};
    use crate::molecule::base::Molecule;
    use crate::molecule::molecule3d::Molecule3D;
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
        assert_eq!(
            cis_2_butene.chiral_classes().unwrap()[0],
            ChiralClass::Counterclockwise
        );
        assert_eq!(
            cis_2_butene.chiral_classes().unwrap()[1],
            ChiralClass::Clockwise
        );

        let molecules = Molecule3D::from_sdf("tests/trans_2_butene.sdf").unwrap();
        let mut trans_2_butene = molecules[0].clone();
        trans_2_butene.identify_chiral_classes();
        assert_eq!(
            trans_2_butene.chiral_classes().unwrap()[0],
            ChiralClass::Counterclockwise
        );
        assert_eq!(
            trans_2_butene.chiral_classes().unwrap()[1],
            ChiralClass::Counterclockwise
        );

        let molecules = Molecule3D::from_sdf("tests/pentene.sdf").unwrap();
        let mut pentene = molecules[0].clone();
        pentene.identify_chiral_classes();
        pentene.identify_chiral_classes();
        assert_eq!(
            pentene.chiral_classes().unwrap()[1],
            ChiralClass::Counterclockwise
        );
        assert_eq!(
            pentene.chiral_classes().unwrap()[3],
            ChiralClass::Counterclockwise
        );
    }
}
