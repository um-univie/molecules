from rdkit import Chem
from pathlib import Path
from typing import List, Tuple

def filter_and_write_smiles(file_path: str, smarts_pattern: str) -> str:
    """
    Read SMILES from file, filter by SMARTS pattern, and write remaining SMILES 
    to a new file with pattern name in filename.

    Args:
        file_path: Path to the input file
        smarts_pattern: SMARTS pattern to filter against

    Returns:
        Path to the output file

    Examples:
        >>> # Example content of input.smi:
        >>> # CC(=O)O,CCC(=O)O
        >>> # CN,CCO
        >>> output_file = filter_and_write_smiles('input.smi', '[OH]')
        >>> # Creates 'input_filtered_OH.smi' containing only 'CN,CCO'
        >>> print(output_file)
        'input_filtered_OH.smi'
    """
    # Compile SMARTS pattern once
    smarts_mol = Chem.MolFromSmarts(smarts_pattern)
    if not smarts_mol:
        raise ValueError(f"Invalid SMARTS pattern: {smarts_pattern}")

    # Create output filename
    input_path = Path(file_path)
    pattern_name = smarts_pattern.replace('[', '').replace(']', '')
    output_path = input_path.parent / f"{input_path.stem}_filtered_{pattern_name}{input_path.suffix}"
    count = 0
    
    with open(file_path, 'r') as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            smiles_list = [s.strip() for s in line.split('\t') if s.strip()]
            should_filter = False

            print(len(smiles_list))
            # Check each SMILES in the line
            for smiles in smiles_list:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(smarts_mol):
                    should_filter = True
                    count += 1
                    break
            
            # Write line to output file if it wasn't filtered
            if not should_filter:
                f_out.write(line)

    print(f"Filtered {count} SMILES")
    return str(output_path)

if __name__ == '__main__':
    file_path = 'smiles_tree_permutations.txt'
    smarts_pattern = 'c1(=[N])ccccc1=[C,N]'
    output_file = filter_and_write_smiles(file_path, smarts_pattern)
    #output_file = filter_and_write_smiles(file_path, smarts_pattern)

    
    #print(f"Filtered SMILES have been written to: {output_file}")