from rdkit import Chem
from typing import List

def are_smiles_equivalent(smiles_list: List[str]) -> tuple[bool, list[str]]:
    """
    Check if all SMILES strings in a given list are equivalent and identify outliers.
    
    Args:
        smiles_list (List[str]): List of SMILES strings to compare
        
    Returns:
        tuple[bool, list[str]]: (True if all equivalent, list of non-matching SMILES)
        
    Examples:
        >>> result, outliers = are_smiles_equivalent(['CC(=O)O', 'OC(C)=O', 'CC(=O)O'])
        >>> result
        True
        >>> outliers
        []
        >>> result, outliers = are_smiles_equivalent(['CC(=O)O', 'CCO', 'CC(=O)O'])
        >>> result
        False
        >>> outliers
        ['CCO']
    """
    if not smiles_list or len(smiles_list) < 2:
        return True, []
    
    outliers = []
    try:
        reference_mol = Chem.MolFromSmiles(smiles_list[0])
        if reference_mol is None:
            return False, [smiles_list[0]]
        reference_canonical = Chem.MolToSmiles(reference_mol)
        
        # Compare each SMILES with the reference
        for smiles in smiles_list[1:]:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                outliers.append(smiles)
                continue
            canonical = Chem.MolToSmiles(mol)
            if canonical != reference_canonical:
                outliers.append(smiles)
                
        return len(outliers) == 0, outliers
        
    except Exception as e:
        print(f"Error processing SMILES: {e}")
        return False, smiles_list

# Example usage
if __name__ == "__main__":
    with open("smiles_tree_permutations.txt", "r") as file:
        test_cases = [line.strip().split() for line in file]
    
    for smiles_list in test_cases:
        result, outliers = are_smiles_equivalent(smiles_list)
        if not result:
            print(f"SMILES list contains non-equivalent structures:")
            print(f"Reference SMILES: {smiles_list[0]}")
            print(f"Non-matching SMILES: {outliers}")
            print("-" * 50)
