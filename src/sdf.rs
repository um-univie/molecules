use std::str::FromStr;
use std::num::ParseFloatError;
use std::num::ParseIntError;
use std::collections::HashMap;
use crate::io::BondTarget;
use crate::{vector::Vector,bond::{Bond,BondOrder},molecule::Molecule3D};  
use std::io::{self, BufRead, BufReader};
use chemistry_consts::ElementProperties;
use tinyvec::ArrayVec;


/// Represents an atom in the V2000 format of an SDF file.
///
/// This struct contains all the information that can be specified for an atom
/// in the V2000 format, including its 3D coordinates, element symbol, and various
/// properties that may be used in different chemical contexts.
#[derive(Debug, Default)]
pub struct AtomV2000 {
    /// X coordinate of the atom
    pub x: f64,
    /// Y coordinate of the atom
    pub y: f64,
    /// Z coordinate of the atom
    pub z: f64,
    /// Element symbol of the atom
    pub symbol: String,
    /// Mass difference from the most abundant isotope
    pub isotope_mass_difference: Option<i8>,
    /// Formal charge of the atom
    pub formal_charge: Option<i8>,
    /// Stereo parity of the atom (0 = not stereo, 1 = odd, 2 = even, 3 = either)
    pub stereo_parity: Option<u8>,
    /// Number of hydrogen atoms attached
    pub hydrogen_count: Option<u8>,
    /// Stereo care box (0 = ignore, 1 = use in stereo perception)
    pub stereo_care_box: Option<u8>,
    /// Valence of the atom
    pub valence: Option<u8>,
    /// H0 designator for implicit hydrogen count (0 = no implicit hydrogens)
    pub implicit_hydrogen_count: Option<u8>,
    /// Atom-atom mapping number for reactions (1-based)
    pub atom_atom_mapping_number: Option<usize>,
    /// Inversion/retention flag for reactions (0 = not stereo, 1 = inversion, 2 = retention)
    pub inversion_retention_flag: Option<u8>,
    /// Exact change flag for reactions (0 = property not applied, 1 = property applied)
    pub exact_change_flag: Option<u8>,


    
}

impl FromStr for AtomV2000 {
    type Err = SDFParseError;

    /// Creates an `AtomV2000` instance from a string slice.
    ///
    /// # Arguments
    ///
    /// * `input` - A string slice containing atom data.
    ///
    /// # Errors
    ///
    /// Returns an error if the input string is too short or if parsing fails.
    fn from_str(input: &str) -> Result<Self, Self::Err> {
        if input.len() < 69 {
            return Err(SDFParseError::InvalidAtom("Atom string too short".into()));
        }

        Ok(AtomV2000 {
            x: input[0..10].trim().parse().map_err(SDFParseError::ParseFloatError)?,
            y: input[10..20].trim().parse().map_err(SDFParseError::ParseFloatError)?,
            z: input[20..30].trim().parse().map_err(SDFParseError::ParseFloatError)?,
            symbol: input[31..34].trim().to_string(),
            isotope_mass_difference: input[34..36].trim().parse().ok(),
            formal_charge: input[36..39].trim().parse().ok(),
            stereo_parity: input[39..42].trim().parse().ok(),
            hydrogen_count: input[42..45].trim().parse().ok(),
            stereo_care_box: input[45..48].trim().parse().ok(),
            valence: input[48..51].trim().parse().ok(),
            implicit_hydrogen_count: input[51..54].trim().parse().ok(),
            atom_atom_mapping_number: input[60..63].trim().parse().ok(),
            inversion_retention_flag: input[63..66].trim().parse().ok(),
            exact_change_flag: input[66..69].trim().parse().ok(),
        })
    }
}


impl Bond {
    /// Parses a bond from a MOL file line.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice that holds the line from the MOL file.
    ///
    /// # Returns
    ///
    /// * `Result<Bond, Box<dyn std::error::Error>>` - A result containing the parsed Bond or an error.
    ///
    /// # Errors
    ///
    /// * Returns an error if the line is too short or if parsing fails.
    pub fn from_mol_line(line: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 3 {
            return Err("Could not parse bond in sdf line".into());
        }

        let atom1 = parts.get(0).ok_or("Missing atom1 in bond line")?.parse::<usize>()?;
        let atom2 = parts.get(1).ok_or("Missing atom2 in bond line")?.parse::<usize>()?;

        Ok(Bond {
            atom1: atom1.checked_sub(1).ok_or("atom1 index underflow")?,
            atom2: atom2.checked_sub(1).ok_or("atom2 index underflow")?,
            bond_order: BondOrder::from_sdf_str(parts[2])?,
        })
    }
}

fn parse_counts_line(line: &str) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 2 {
        return Err("Invalid line length".into());
    }
    let num_atoms = parts[0].parse()?;
    let num_bonds = parts[1].parse()?;
    Ok((num_atoms, num_bonds))
}

impl FromStr for CountsLine {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() < 39 {
            return Err("Counts line to short".into());
        }

        let parse_u16 = |s: &str| -> Result<u16, ParseIntError> {
            s.trim().parse::<u16>()
        };

        let parse_u8 = |s: &str| -> Result<u8, ParseIntError> {
            s.trim().parse::<u8>()
        };

        Ok(CountsLine {
            num_atoms: parse_u16(&s[0..3])?,
            num_bonds: parse_u16(&s[3..6])?,
            num_atom_lists: parse_u16(&s[6..9])?,
            chiral_flag: parse_u8(&s[12..15])?,
            num_stext_entries: parse_u16(&s[15..18])?,
            num_properties: parse_u16(&s[30..33])?,
            version: s[33..39].trim().to_string(),
        })
    }
}

pub enum MolVersion {
    V2000,
    V3000
}

impl FromStr for MolVersion { 
    type Err = Box<dyn std::error::Error>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "V2000" => Ok(MolVersion::V2000),
            "V3000" => Ok(MolVersion::V3000),
            _ =>  Err("Invalid mol version {}, only V2000 and V3000 are supported".into())
        }
    }
}
/// Represents the counts line in an SDF (Structure-Data File) or MOL file.
///
/// This struct contains information about the number of atoms, bonds, and other
/// structural features of a molecule as specified in the counts line of an SDF file.
///
/// # Fields
///
/// * `num_atoms` - The number of atoms in the molecule.
/// * `num_bonds` - The number of bonds in the molecule.
/// * `num_atom_lists` - The number of atom lists (typically used for R-groups or variable attachments).
/// * `chiral_flag` - Indicates whether the molecule is chiral (1) or not (0).
/// * `num_stext_entries` - The number of Stext entries (additional text entries in the file).
/// * `num_properties` - The number of additional properties specified for the molecule.
/// * `version` - The version of the MOL file format (e.g., "V2000" or "V3000").
#[derive(Debug, Default)]
pub struct CountsLine {
    pub num_atoms: u16,
    pub num_bonds: u16,
    pub num_atom_lists: u16,
    pub chiral_flag: u8,
    pub num_stext_entries: u16,
    pub num_properties: u16,
    pub version: String,
}


/// Represents the structure of a MOL file, which is a common format for describing molecular structures.
///
/// # Fields
///
/// * `header` - A vector of strings containing the header information of the MOL file.
///              This typically includes the molecule name, user/program/date information, and comments.
/// * `counts_line` - A string representing the counts line of the MOL file.
///                   This line contains information about the number of atoms, bonds, and other structural features.
/// * `atoms` - A vector of `AtomV2000` structs, each representing an atom in the molecule.
///             This includes information such as coordinates, element type, and properties for each atom.
/// * `bonds` - A vector of `Bond` structs, each representing a bond between atoms in the molecule.
///             This includes information about which atoms are connected and the type of bond.
/// * `properties` - A vector of strings containing additional properties or attributes of the molecule.
///                  These can include various molecular properties or custom data fields.
#[derive(Debug, Default)]
pub struct MOLFile {
    pub header: Vec<String>,
    pub counts_line: String,
    pub atoms: Vec<AtomV2000>,
    pub bonds: Vec<Bond>,
    pub properties: Vec<String>,
}

#[derive(Debug)]
pub struct SDFEntry {
    pub mol_file: MOLFile,
    pub data_fields: HashMap<String, String>,
}

#[derive(Debug)]
pub enum ParseError {
    IoError(io::Error),
    Atom,
    Bond,
    UnexpectedEndOfFile,
}

#[derive(Debug)]
pub enum SDFParseError {
    MOLParseError(ParseError),
    ParseFloatError(ParseFloatError),
    ParseIntError(ParseIntError),
    InvalidAtom(String),
    IoError(io::Error),
    InvalidBond,
    UnexpectedEndOfFile,
}

impl From<ParseError> for SDFParseError {
    fn from(error: ParseError) -> Self {
        SDFParseError::MOLParseError(error)
    }
}

impl From<io::Error> for SDFParseError {
    fn from(error: io::Error) -> Self {
        SDFParseError::IoError(error)
    }
}

fn parse_mol_file<I>(lines: &mut I) -> Result<MOLFile, ParseError>
where
    I: Iterator<Item = Result<String, io::Error>>,
{
    let mut mol_file = MOLFile::default();
    let mut state = ParserState::Header(0);

    for line in lines.by_ref() {
        let line = line.map_err(ParseError::IoError)?;
        state = match state {
            ParserState::Header(count) if count < 3 => {
                mol_file.header.push(line);
                ParserState::Header(count + 1)
            }
            ParserState::Header(_) => {
                mol_file.counts_line = line;
                ParserState::Atoms
            }
            ParserState::Atoms => {
                match AtomV2000::from_str(&line) {
                    Ok(atom) => {
                        mol_file.atoms.push(atom);
                        ParserState::Atoms
                    }
                    Err(_) => {
                        match Bond::from_mol_line(&line) {
                            Ok(bond) => {
                                mol_file.bonds.push(bond);
                                ParserState::Bonds
                            }
                            Err(_) => {
                                if line.trim() == "M  END" {
                                    return Ok(mol_file);
                                } else {
                                    ParserState::Properties
                                }
                            }
                        }
                    }
                }
            }
                ParserState::Bonds => {
                    match Bond::from_mol_line(&line) {
                        Ok(bond) => {
                            mol_file.bonds.push(bond);
                            ParserState::Bonds
                        }
                        Err(_) => {
                            if line.trim() == "M  END" {
                                return Ok(mol_file);
                            } else {
                                ParserState::Properties
                            }
                        }
                    }
                }
                ParserState::Properties => {
                    println!("Property line: {}", line);
                    if line.trim() == "M  END" {
                        return Ok(mol_file);
                    } else {

                mol_file.properties.push(line);
                ParserState::Properties
            }
        }
    };
}

    Err(ParseError::UnexpectedEndOfFile)
}

pub fn parse_sdf_file<R: BufRead>(reader: R) -> Result<Vec<SDFEntry>, SDFParseError> {
    let mut entries = Vec::new();
    let mut lines = reader.lines();
    let mut current_data_fields = HashMap::new();
    let mut in_mol_block = true;
    let mut current_data_field_name = String::new();

    loop {
        if in_mol_block {
            match parse_mol_file(&mut lines) {
                Ok(mol_file) => {
                    in_mol_block = false;
                    entries.push(SDFEntry {
                        mol_file,
                        data_fields: HashMap::new(),
                    });
                }
                Err(ParseError::UnexpectedEndOfFile) => break,
                Err(e) => return Err(e.into()),
            }
        } else {
            match lines.next() {
                Some(Ok(line)) => {
                    if line.trim() == "$$$$" {
                        if let Some(entry) = entries.last_mut() {
                            entry.data_fields = std::mem::take(&mut current_data_fields);
                        }
                        in_mol_block = true;
                    } else if let Some(stripped) = line.strip_prefix("> ") {
                        // Remove < and > from the data field name
                        current_data_field_name = stripped.replace("<", "").replace(">", "").trim().to_string();
                    } else if !current_data_field_name.is_empty() {
                        current_data_fields.entry(current_data_field_name.clone())
                            .or_insert(String::new())
                            .push_str(line.trim());
                    }
                }
                Some(Err(e)) => return Err(SDFParseError::IoError(e)),
                None => break,
            }
        }
    }

    // Handle last entry if file doesn't end with $$$$
    if let Some(entry) = entries.last_mut() {
        if !current_data_fields.is_empty() {
            entry.data_fields = current_data_fields;
        }
    }

    Ok(entries)
}

enum ParserState {
    Header(usize),
    Atoms,
    Bonds,
    Properties,
}

#[derive(Debug, Default)]
pub struct BondV2000 {
    pub atom1: usize,
    pub atom2: usize,
    pub bond_type: u8,
    pub bond_stereo: u8,
    pub bond_topology: Option<u8>,
    pub reacting_center_status: Option<u8>,
}

#[derive(Debug, Default)]
pub struct SGroup {
    pub sgroup_type: String,
    pub atoms: Vec<usize>,
    pub parent: Option<usize>,
    pub component_number: Option<u8>,
    pub connectivity: Option<String>,
    pub subscript: Option<String>,
    pub expansion: Option<bool>,
    pub bracket_style: Option<u8>,
    pub class: Option<String>,
}

#[derive(Debug, Default)]
pub struct MoleculeV2000 {
    pub atoms: Vec<AtomV2000>,
    pub bonds: Vec<BondV2000>,
    pub sgroups: Vec<SGroup>,
    pub properties: Vec<String>,
    pub chiral: bool,
    pub name: String,
    pub comment: String,
    pub dimension: u8,
}


use std::path::Path;
impl Molecule3D {
    pub fn from_sdf<T: AsRef<Path>>(sdf: T) -> Result<Vec<Self>, SDFParseError>{
        let file = std::fs::File::open(sdf)?;
        let reader = BufReader::new(file);
        let sdf_entries = parse_sdf_file(reader)?; 
        Self::from_sdf_entries(sdf_entries)
    }

    pub fn from_sdf_entries(sdf_entries: Vec<SDFEntry>) -> Result<Vec<Self>, SDFParseError> {
        let mut molecules = vec![];

        for entry in sdf_entries.into_iter() {
            let number_of_atoms = entry.mol_file.atoms.len();
            let mut atomic_numbers = vec![];
            let mut charges = vec![];
            let mut positions = vec![];
            
            for atom in entry.mol_file.atoms {
                let atom_charge = atom.formal_charge.unwrap_or(0);
                atomic_numbers.push(atom.symbol.to_uppercase().as_str().atomic_number().ok_or(SDFParseError::InvalidAtom(atom.symbol.to_string()))?);
                charges.push(atom_charge);
                positions.push(Vector::new(atom.x, atom.y, atom.z));
            }

            let mut bonds = vec![ArrayVec::<[BondTarget; 10]>::new(); number_of_atoms];
            for bond in entry.mol_file.bonds {
                bonds.get_mut(bond.atom1).ok_or(SDFParseError::InvalidBond)?.push(BondTarget::new(bond.atom2, bond.bond_order));
                bonds.get_mut(bond.atom2).ok_or(SDFParseError::InvalidBond)?.push(BondTarget::new(bond.atom1, bond.bond_order));
            }
            molecules.push(
                Molecule3D {
                    charges,
                    positions,
                    atomic_numbers,
                    atom_bonds: bonds,
                    ..Default::default()
                }
            )
        }

        Ok(
            molecules
        )
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_sdf_file() {
        let sdf_content = "\
        5460033
  -OEChem-09092410012D

 11  8  0     0  0  0  0  0  0999 V2000
    3.2690    0.1614    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0
    3.9761    0.8685    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    3.9761   -0.5458    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    2.8464   -0.7450    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.5619    0.8685    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.5844   -1.3069    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.2844   -0.4829    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.4084   -1.0070    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1235    1.3069    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0003    1.3069    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.6064    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  4  6  1  0  0  0  0
  4  7  1  0  0  0  0
  4  8  1  0  0  0  0
  5  9  1  0  0  0  0
  5 10  1  0  0  0  0
  5 11  1  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
5460033

> <PUBCHEM_COMPOUND_CANONICALIZED>
1

> <PUBCHEM_CACTVS_COMPLEXITY>
7.6

> <PUBCHEM_CACTVS_HBOND_ACCEPTOR>
2

> <PUBCHEM_CACTVS_HBOND_DONOR>
2

> <PUBCHEM_CACTVS_ROTATABLE_BOND>
0

> <PUBCHEM_CACTVS_SUBSKEYS>
AAADcYADAAAGAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==

> <PUBCHEM_IUPAC_OPENEYE_NAME>
ammonia;dichloroplatinum

> <PUBCHEM_IUPAC_CAS_NAME>
ammonia;dichloroplatinum

> <PUBCHEM_IUPAC_NAME_MARKUP>
azane;dichloroplatinum

> <PUBCHEM_IUPAC_NAME>
azane;dichloroplatinum

> <PUBCHEM_IUPAC_SYSTEMATIC_NAME>
azane;bis(chloranyl)platinum

> <PUBCHEM_IUPAC_TRADITIONAL_NAME>
ammonia;dichloroplatinum

> <PUBCHEM_IUPAC_INCHI>
InChI=1S/2ClH.2H3N.Pt/h2*1H;2*1H3;/q;;;;+2/p-2

> <PUBCHEM_IUPAC_INCHIKEY>
LXZZYRPGZAFOLE-UHFFFAOYSA-L

> <PUBCHEM_EXACT_MASS>
298.955598

> <PUBCHEM_MOLECULAR_FORMULA>
Cl2H6N2Pt

> <PUBCHEM_MOLECULAR_WEIGHT>
300.05

> <PUBCHEM_REFERENCE_STANDARDIZATION>
Bypass - this structure was created from CID 5460033

> <PUBCHEM_OPENEYE_CAN_SMILES>
N.N.Cl[Pt]Cl

> <PUBCHEM_OPENEYE_ISO_SMILES>
N.N.Cl[Pt]Cl

> <PUBCHEM_CACTVS_TPSA>
2

> <PUBCHEM_MONOISOTOPIC_WEIGHT>
298.955598

> <PUBCHEM_TOTAL_CHARGE>
0

> <PUBCHEM_HEAVY_ATOM_COUNT>
5

> <PUBCHEM_ATOM_DEF_STEREO_COUNT>
0

> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_DEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_UDEF_STEREO_COUNT>
0

> <PUBCHEM_ISOTOPIC_ATOM_COUNT>
0

> <PUBCHEM_COMPONENT_COUNT>
3

> <PUBCHEM_CACTVS_TAUTO_COUNT>
-1

> <PUBCHEM_NONSTANDARDBOND>
1  4  6
1  5  6

> <PUBCHEM_COORDINATE_TYPE>
1
5
255

$$$$
";

        let reader = Cursor::new(sdf_content);
        let result = parse_sdf_file(reader);
        assert!(result.is_ok());
        
        let entries = result.unwrap();
        assert_eq!(entries.len(), 1);

        // Check first compound
        let first = &entries[0];
        assert_eq!(first.mol_file.atoms.len(), 11);
        assert_eq!(first.mol_file.bonds.len(), 8);
        assert_eq!(first.data_fields.len(), 34);
        assert_eq!(first.data_fields.get("PUBCHEM_COMPOUND_CID"), Some(&"5460033".to_string()));

    }

    #[test]
    fn test_parse_sdf_file_single_entry() {
        let sdf_content = "\
Compound 1
  -OEChem-01152411432D

  1  0  0     0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
1

$$$$
";

        let reader = Cursor::new(sdf_content);
        let result = parse_sdf_file(reader);
        assert!(result.is_ok());
        
        let entries = result.unwrap();
        assert_eq!(entries.len(), 1);

        let entry = &entries[0];
        assert_eq!(entry.mol_file.atoms.len(), 1);
        assert_eq!(entry.mol_file.bonds.len(), 0);
        assert_eq!(entry.data_fields.len(), 1);
        assert_eq!(entry.data_fields.get("PUBCHEM_COMPOUND_CID"), Some(&"1".to_string()));
    }

    #[test]
    fn test_parse_sdf_file_no_data_fields() {
        let sdf_content = "\
Compound 1
  -OEChem-01152411432D

  1  0  0     0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
";

        let reader = Cursor::new(sdf_content);
        let result = parse_sdf_file(reader);
        assert!(result.is_ok());
        
        let entries = result.unwrap();
        assert_eq!(entries.len(), 1);

        let entry = &entries[0];
        assert_eq!(entry.mol_file.atoms.len(), 1);
        assert_eq!(entry.mol_file.bonds.len(), 0);
        assert_eq!(entry.data_fields.len(), 0);
    }

    #[test]
    fn test_parse_sdf_file_missing_end_delimiter() {
        let sdf_content = "\
Compound 1
  -OEChem-01152411432D

  1  0  0     0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
1
";

        let reader = Cursor::new(sdf_content);
        let result = parse_sdf_file(reader);
        assert!(result.is_ok());
        
        let entries = result.unwrap();
        assert_eq!(entries.len(), 1);

        let entry = &entries[0];
        assert_eq!(entry.mol_file.atoms.len(), 1);
        assert_eq!(entry.mol_file.bonds.len(), 0);
        assert_eq!(entry.data_fields.len(), 1);
        assert_eq!(entry.data_fields.get("PUBCHEM_COMPOUND_CID"), Some(&"1".to_string()));
    }

    #[test]
    fn test_parse_sdf_file_empty_input() {
        let sdf_content = "";
        let reader = Cursor::new(sdf_content);
        let result = parse_sdf_file(reader);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 0);
    }

}
