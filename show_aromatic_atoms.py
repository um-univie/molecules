import argparse
import sys
import os
import subprocess
from typing import List

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D


def get_aromatic_atoms(smiles: str) -> List[int]:
    """
    Identify and return indices of aromatic atoms in a molecule.

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        List[int]: List of aromatic atom indices.

    Doctests:
        >>> get_aromatic_atoms("c1ccccc1")
        [0, 1, 2, 3, 4, 5]
        >>> get_aromatic_atoms("CCO")
        []
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES string provided.")

    aromatic_atoms = [
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetIsAromatic()
    ]
    return aromatic_atoms


def display_aromatic_atoms(smiles: str) -> None:
    """
    Display aromatic atoms in the molecule represented by the SMILES string.

    Args:
        smiles (str): SMILES representation of the molecule.

    Doctests:
        >>> display_aromatic_atoms("c1ccccc1")
        Aromatic atoms indices: [0, 1, 2, 3, 4, 5]
        >>> display_aromatic_atoms("CCO")
        Aromatic atoms indices: []
    """
    aromatic_indices = get_aromatic_atoms(smiles)
    print(f"Aromatic atoms indices: {aromatic_indices}")


def highlight_aromatic_atoms(smiles: str) -> None:
    """
    Generate and display a graphical representation of the molecule, highlighting aromatic atoms.

    Args:
        smiles (str): SMILES representation of the molecule.

    Doctests:
        >>> highlight_aromatic_atoms("c1ccccc1")
        Molecule image displayed with aromatic atoms highlighted.
        >>> highlight_aromatic_atoms("CCO")
        Molecule image displayed without highlighted aromatic atoms.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES string provided.")

    aromatic_indices = get_aromatic_atoms(smiles)

    # Compute 2D coordinates
    Chem.rdDepictor.Compute2DCoords(molecule)

    # Initialize drawer
    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)  # 300x300 pixels
    options = drawer.drawOptions()

    # Define highlight colors
    highlight_colors = {idx: (1.0, 0.0, 0.0) for idx in aromatic_indices}  # Red for aromatic atoms

    # Draw molecule with highlights
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        molecule,
        highlightAtoms=aromatic_indices,
        highlightAtomColors=highlight_colors
    )
    drawer.FinishDrawing()

    # Save image to a temporary file
    image_path = "temp_aromatic_molecule.png"
    with open(image_path, "wb") as img_file:
        img_file.write(drawer.GetDrawingText())
    print(f"Molecule image saved to {image_path}")

    # Display the image using the default image viewer
    open_image(image_path)


def open_image(image_path: str) -> None:
    """
    Open an image file using the default image viewer of the operating system.

    Args:
        image_path (str): Path to the image file.

    Raises:
        OSError: If the image cannot be opened.
    """
    try:
        if sys.platform.startswith('darwin'):
            subprocess.run(['open', image_path], check=True)
        elif os.name == 'nt':
            os.startfile(image_path)
        elif os.name == 'posix':
            subprocess.run(['xdg-open', image_path], check=True)
        else:
            print(f"Cannot determine how to open images on this OS: {sys.platform}")
    except Exception as e:
        print(f"Failed to open image: {e}")


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Identify and visualize aromatic atoms in a molecule given its SMILES string."
    )
    parser.add_argument(
        "smiles",
        type=str,
        help="SMILES representation of the molecule."
    )
    return parser.parse_args()


def main():
    """
    Main function to execute the aromatic atom identification and visualization.

    Example:
        >>> main()
        Aromatic atoms indices: [0, 1, 2, 3, 4, 5]
        Molecule image saved to temp_aromatic_molecule.png
        Molecule image displayed with aromatic atoms highlighted.
    """
    args = parse_arguments()
    smiles_input = args.smiles

    try:
        display_aromatic_atoms(smiles_input)
        highlight_aromatic_atoms(smiles_input)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
