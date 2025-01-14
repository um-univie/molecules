// Unit tests (Still needs to be improved)
#[cfg(test)]
mod kekulize_tests {
    use crate::atom::Atom;
    use crate::molecule::base::AtomLabel;
    use crate::molecule::molecular_system::extract_atom_pdb;
    use crate::prelude::*;
    use crate::vector::Vector;
    #[test]
    fn test_extract_atom() {
        let line =
            "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N  ";
        let atom = extract_atom_pdb(line).unwrap();
        let position = atom.position_vector.unwrap();
        assert_eq!(atom.atomic_number, 7);
        assert_eq!(position.x, 10.0);
        assert_eq!(position.y, 10.0);
        assert_eq!(position.z, 10.0);
    }
    #[test]
    fn test_vector_angle() {
        let v1 = Vector {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        };
        let v2 = Vector {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        };
        let angle = v1.angle_between(&v2);
        assert_eq!(angle, Some(std::f64::consts::FRAC_PI_2));
    }
    #[test]
    fn test_cross_product() {
        let v1 = Vector {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        };
        let v2 = Vector {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        };
        let v3 = v1.cross(&v2);
        assert_eq!(v3.x, 0.0);
        assert_eq!(v3.y, 0.0);
        assert_eq!(v3.z, 1.0);
    }
    #[test]
    fn test_bond_angle() {
        let atom1 = Atom::new(6).with_position((1.7, 0.0, 0.0));
        let atom2 = Atom::new(6).with_position((0.0, 0.0, 0.0));
        let atom3 = Atom::new(6).with_position((0.0, 0.0, 1.7));
        let molecule = Molecule3D::from_atoms(vec![atom1, atom2, atom3]);
        let angles = molecule.find_angles();
        println!("{:?}", angles);
        for bonds in molecule.atom_bonds.iter() {
            println!("{:?}", bonds);
        }
        assert_eq!(angles.len(), 1);
        assert_eq!(angles[0].angle(), std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_alternate_covalent_radii() {
        use chemistry_consts::*;
        let atom1 = Atom::new(6).with_position((1.7, 0.0, 0.0));
        let atom2 = Atom::new(6).with_position((0.0, 0.0, 0.0));
        let atom3 = Atom::new(6).with_position((0.0, 0.0, 1.7));
        let molecule = Molecule3D::from_atoms_alternate_covalent_radii(
            vec![atom1, atom2, atom3],
            &COVALENT_RADII,
        );
        let angles = molecule.find_angles();
        assert_eq!(angles.len(), 1);
        assert_eq!(angles[0].angle(), std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_kekulize_benzene() {
        // Test simple benzene ring
        let molecules = Molecule2D::from_smiles("c1ccccc1").unwrap();
        let mut molecule = molecules.first().unwrap().clone();
        molecule.kekulize().unwrap();

        // Check alternating single and double bonds
        let bonds = molecule.get_edges_with_type();
        let double_bonds = bonds
            .iter()
            .filter(|(_, _, bond_type)| *bond_type == BondOrder::Double)
            .count();

        assert_eq!(double_bonds, 3, "Benzene should have 3 double bonds");
        assert_eq!(molecule.to_smiles(), "C1=CC=CC=C1");
    }

    #[test]
    fn test_kekulize_pyridine() {
        // Test pyridine ring
        let molecules = Molecule2D::from_smiles("n1ccccc1").unwrap();
        let mut molecule = molecules.first().unwrap().clone();
        molecule.kekulize().unwrap();

        assert_eq!(molecule.to_smiles(), "N1=CC=CC=C1");
    }

    #[test]
    fn test_kekulize_naphthalene() {
        // Test fused rings (naphthalene)
        let molecules = Molecule2D::from_smiles("c1ccc2ccccc2c1").unwrap();
        let mut molecule = molecules.first().unwrap().clone();
        molecule.kekulize().unwrap();

        println!("{}", molecule.to_smiles());
        let bonds = molecule.get_edges_with_type();
        let double_bonds = bonds
            .iter()
            .filter(|(_, _, bond_type)| *bond_type == BondOrder::Double)
            .count();

        assert_eq!(double_bonds, 5, "Naphthalene should have 5 double bonds");
    }

    #[test]
    fn test_kekulize_non_aromatic() {
        // Test molecule with no aromatic bonds
        let molecules = Molecule2D::from_smiles("CC").unwrap();
        let mut molecule = molecules.first().unwrap().clone();
        molecule.kekulize().unwrap();

        assert_eq!(molecule.to_smiles(), "CC");
    }

    #[test]
    fn test_kekulize_furan() {
        // Test heterocyclic aromatic ring
        let molecules = Molecule2D::from_smiles("o1cccc1").unwrap();
        let mut molecule = molecules.first().unwrap().clone();
        molecule.kekulize().unwrap();

        assert_eq!(molecule.to_smiles(), "O1C=CC=C1");
    }

    #[test]
    fn test_kekulize_phenol() {
        // Test aromatic ring with substituent
        let molecules = Molecule2D::from_smiles("Oc1ccccc1").unwrap();
        let mut molecule = molecules.first().unwrap().clone();
        molecule.kekulize().unwrap();

        assert_eq!(molecule.to_smiles(), "OC1=CC=CC=C1");
    }

    #[test]
    fn test_kekulize_multiple_rings() {
        // Test molecule with multiple aromatic rings
        let molecules = Molecule2D::from_smiles("c1ccccc1Cc2ccccc2").unwrap();
        let mut molecule = molecules.first().unwrap().clone();
        molecule.kekulize().unwrap();

        let bonds = molecule.get_edges_with_type();
        let double_bonds = bonds
            .iter()
            .filter(|(_, _, bond_type)| *bond_type == BondOrder::Double)
            .count();

        assert_eq!(
            double_bonds, 6,
            "Two benzene rings should have 6 double bonds total"
        );
    }
}
#[cfg(test)]
mod aromatization_tests {
    use crate::molecule::{base::Molecule, bond::BondOrder, molecule2d::Molecule2D};
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    #[test]
    fn test_benzene_aromatization() {
        // Create benzene molecule
        let molecule = Molecule2D::from_smiles("C1=CC=CC=C1").unwrap();
        let mut molecule = molecule.first().unwrap().clone();

        molecule.simple_aromatization();

        // Check if all bonds are aromatic
        for (i, bonds) in molecule.atom_bonds().iter().enumerate() {
            for bond in bonds {
                if molecule.atomic_numbers()[bond.target()] != 1
                    && molecule.atomic_numbers()[i] != 1
                {
                    assert_eq!(bond.bond_order(), BondOrder::Aromatic);
                }
            }
        }

        assert_eq!(molecule.to_smiles(), "c1ccccc1");
    }


    #[test]
    fn test_difficult_smiles() {
        let smiles = "c1c2c(ccc1)OCC2C(=O)N(C)c3c(N)cncc3";
        let molecule = diff_smiles(smiles);
        assert_eq!(molecule, smiles);

        let smiles = "CCCCC1CCC(C(Br)C2=CC3=C(C=C2)CCC3)CC1";
        let molecule = diff_smiles(smiles);
        assert_eq!(molecule, "CCCCC1CCC(C(Br)c2cc3c(cc2)CCC3)CC1");
        
        let smiles = "C(C)(C)C1=CC(C(Cl)Cl)C2=CC=CC=C2O1";
        let molecule = diff_smiles(smiles);
        assert_eq!(molecule, "C(C)(C)C1=CC(C(Cl)Cl)c2ccccc2O1");

        let smiles = "c12c(cc(-c3c4ncn(C5CCNCC5)c4ncc3-c3cc(Cl)ccc3)cc1)cccc2";
        let molecule = diff_smiles(smiles);
        assert_eq!(molecule, "c12c(cc(c3c4ncn(C5CCNCC5)c4ncc3c6cc(Cl)ccc6)cc1)cccc2");
        
        let smiles = "c12c3cccc1C(C(c2ccc3)=O)=NNC(C(C)Nc4nc5c(cccc5)s4)=O";
        let molecule = diff_smiles(smiles);
        assert_eq!(molecule, "c12c3cccc1c(c(c2ccc3)=O)=NNC(C(C)Nc4nc5c(cccc5)s4)=O");

        let smiles = "c12c(C=CC=C1)oc3c2N(CC(Nc4cc(ccc4)C)=O)C(N(c5ccc(cc5)OCC)C3=O)=O";
        let molecule = diff_smiles(smiles);
        assert_eq!(molecule, "c12c(cccc1)oc3c2n(CC(Nc4cc(ccc4)C)=O)c(n(c5ccc(cc5)OCC)c3=O)=O");

        let smiles = "COCCOC(=O)C1=C(C)N=C2SC(=CC3=CC=C(C4=CC(Cl)=CC=C4)O3)C(=O)N2C1C1=CC=C(C(C)C)C=C1";
        let molecule = diff_smiles(smiles);
        assert_eq!(molecule, "COCCOC(=O)C1=C(C)N=c2sc(=Cc3ccc(c4cc(Cl)ccc4)o3)c(=O)n2C1c5ccc(C(C)C)cc5");

        //C(C)(C)C1=CC(C(Cl)Cl)C2=CC=CC=C2O1
    }

    fn diff_smiles(smiles: &str) -> String {
        println!("Canonizing {}", smiles);
        let molecule = parse_smiles_and_aromatize(smiles);
        molecule
    }


    #[test]
    fn test_pyridine_aromatization() {
        let molecule = Molecule2D::from_smiles("C1=CC=CN=C1").unwrap();
        let mut molecule = molecule.first().unwrap().clone();

        molecule.aromatize();

        assert_eq!(molecule.to_smiles(), "c1cccnc1");
    }

    #[test]
    fn test_naphthalene_aromatization() {
        let molecule = Molecule2D::from_smiles("C1=CC=C2C=CC=CC2=C1").unwrap();
        let mut molecule = molecule.first().unwrap().clone();

        molecule.aromatize();

        assert_eq!(molecule.to_smiles(), "c1ccc2ccccc2c1");
    }

    #[test]
    fn test_cyclopentadienyl_anion() {
        let molecule = Molecule2D::from_smiles("C1=CC=C[CH-]1").unwrap();
        let mut molecule = molecule.first().unwrap().clone();

        molecule.aromatize();

        assert_eq!(molecule.to_smiles(), "c1ccc[cH-]1");
    }

    #[test]
    fn test_non_aromatic_cycle() {
        // Cyclooctatetraene is not aromatic (8 electrons, not 4n+2)
        let molecule = Molecule2D::from_smiles("C1=CC=CC=CC=C1").unwrap();
        let mut molecule = molecule.first().unwrap().clone();

        molecule.aromatize();

        assert_eq!(molecule.to_smiles(), "C1=CC=CC=CC=C1");
    }

    #[test]
    fn test_furan_aromatization() {
        let molecule = Molecule2D::from_smiles("C1=COC=C1").unwrap();
        let mut molecule1 = molecule.first().unwrap().clone();
        molecule1.aromatize();

        assert_eq!(molecule1.to_smiles(), "c1cocc1");
    }

    #[test]
    fn test_non_cyclic_molecule() {
        // Test with ethene
        let molecule = Molecule2D::from_smiles("C=C").unwrap();
        let mut molecule = molecule.first().unwrap().clone();

        molecule.aromatize();

        // Check that the double bond remains unchanged
        assert_eq!(molecule.atom_bonds()[0][0].bond_order(), BondOrder::Double);
        assert_eq!(molecule.to_smiles(), "C=C");
    }

    #[test]
    fn test_aromatization_from_file() {
        let file = BufReader::new(File::open("tests/aromatization_tests.smi").unwrap());

        for line in file.lines() {
            let line = line.unwrap();

            // Skip comments and empty lines
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }

            // Parse the line: "kekulized, aromatized # comment"
            let parts: Vec<&str> = line
                .split('#')
                .next()
                .unwrap()
                .split(',')
                .map(|s| s.trim())
                .collect();

            if parts.len() != 2 {
                continue;
            }

            let (kekulized, expected_aromatic) = (parts[0], parts[1]);

            // Create molecule from kekulized form
            let molecule = Molecule2D::from_smiles(kekulized).unwrap();
            let mut molecule = molecule.first().unwrap().clone();

            // Try to aromatize
            molecule.simple_aromatization();

            // Compare with expected
            assert_eq!(
                molecule.to_smiles(),
                expected_aromatic,
                "Failed to correctly aromatize {}: expected {}, got {}",
                kekulized,
                expected_aromatic,
                molecule.to_smiles()
            );
        }
    }

    fn parse_smiles_and_aromatize(smiles: &str) -> String {
        let mut molecules = Molecule2D::from_smiles(smiles).unwrap();
        let mut molecule = molecules.swap_remove(0);
        molecule.simple_aromatization();
        molecule.to_smiles()
    }
}

#[cfg(test)]
mod graph_comparison_tests {
    use crate::molecule::molecule2d::Molecule2D;
    use crate::molecule::base::{Molecule, AtomLabel};
    use petgraph::graph::UnGraph;
    use petgraph::algo::isomorphism::subgraph_isomorphisms_iter;
    use petgraph::prelude::*;
    use petgraph::graph::{EdgeReference, NodeIndex};

    #[derive(Debug)]
    pub struct GraphDifferences {
        unmatched_in_first: Vec<usize>,
        unmatched_in_second: Vec<usize>,
        mapping: Option<Vec<usize>>,
    }

    /// Finds the maximum common subgraph between two molecules and returns the differences
    fn compare_molecular_graphs(mol1: &Molecule2D, mol2: &Molecule2D) -> GraphDifferences {
        let graph1 = mol1.to_ungraph_with_edge_nodes();
        let graph2 = mol2.to_ungraph_with_edge_nodes();
        let g_ref = &graph1;
        let h_ref = &graph2;
        
        let mut node_match = |node1: &AtomLabel, node2: &AtomLabel| node1 == node2;
        let mut edge_match = |edge1: &(), edge2: &()| edge1 == edge2;

        // Find all possible isomorphisms between the graphs
        let Some(iso) = subgraph_isomorphisms_iter(
            &h_ref, 
            &g_ref, 
            &mut node_match, 
            &mut edge_match
        ) else { return GraphDifferences {
                unmatched_in_first: graph1.node_indices().map(|n| n.index()).collect(),
                unmatched_in_second: graph2.node_indices().map(|n| n.index()).collect(),
                mapping: None,
            }}
        ;

        // Get the largest matching subgraph
        let best_mapping = iso.max_by_key(|mapping| mapping.len());
        
        // Identify unmatched nodes in both graphs
        if let Some(mapping) = best_mapping {
            let unmatched_g1: Vec<_> = graph1.node_indices()
                .filter(|&n| !mapping.contains(&n.index()))
                .map(|n| n.index())
                .collect();
            let unmatched_g2: Vec<_> = graph2.node_indices()
                .filter(|&n| !mapping.contains(&n.index()))
                .map(|n| n.index())
                .collect();
                
            GraphDifferences {
                unmatched_in_first: unmatched_g1,
                unmatched_in_second: unmatched_g2,
                mapping: Some(mapping),
            }
        } else {
            GraphDifferences {
                unmatched_in_first: graph1.node_indices().map(|n| n.index()).collect(),
                unmatched_in_second: graph2.node_indices().map(|n| n.index()).collect(),
                mapping: None,
            }
        }
    }

    #[test]
    fn test_graph_comparison() {
        // Test with similar but different molecules
        let mol2 = Molecule2D::from_smiles("c1ccccc1CC").unwrap();
        let mol1 = Molecule2D::from_smiles("c1ccccc1C").unwrap();
        
        let differences = compare_molecular_graphs(
            mol1.first().unwrap(),
            mol2.first().unwrap()
        );
        println!("{:?}", differences);
        
        // Should find one extra carbon in mol1
        assert_eq!(differences.unmatched_in_first.len(), 1);
        assert_eq!(differences.unmatched_in_second.len(), 0);
        
        // Visualize the differences if needed
        // visualize_graph_differences(&mol1, &mol2, &differences, "graph_diff.png");
    }
}
