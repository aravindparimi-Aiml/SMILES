"""
Cheminformatics Virtual Laboratory - 2D Molecular Depiction Engine
Generates publication-quality 2D vector SVGs and PNG bitmaps with atom/bond highlighting.
"""

import io
from typing import Optional, List, Tuple
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image


def generate_2d_svg(
    smiles: str,
    width: int = 480,
    height: int = 360,
    highlight_atoms: Optional[List[int]] = None,
    highlight_bonds: Optional[List[int]] = None,
    dark_theme: bool = True,
) -> str:
    """
    Generate clean, scalable SVG depiction with RDKit rdMolDraw2D.
    Supports stereochemical wedges/dashes, aromatic bonds, and optional highlighting.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot depict 2D structure: invalid SMILES '{smiles}'")

    # Generate 2D coordinates if not present
    if not mol.GetNumConformers():
        from rdkit.Chem import rdDepictor
        rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.2
    opts.addStereoAnnotation = True
    opts.clearBackground = False  # Transparent background

    if dark_theme:
        # High-contrast colors on dark background
        opts.setSymbolColour((0.92, 0.94, 0.96))

    h_atoms = highlight_atoms if highlight_atoms is not None else []
    h_bonds = highlight_bonds if highlight_bonds is not None else []

    drawer.DrawMolecule(mol, highlightAtoms=h_atoms, highlightBonds=h_bonds)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Clean up SVG XML header for inline embedding
    if "<?xml" in svg:
        svg = svg.split("?>", 1)[-1].strip()

    return svg


def generate_2d_png(
    smiles: str,
    width: int = 500,
    height: int = 400,
    dark_background: bool = False,
) -> bytes:
    """
    Generate PNG byte data for PDF embedding and user downloads.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot depict 2D PNG: invalid SMILES '{smiles}'")

    from rdkit.Chem import rdDepictor
    rdDepictor.Compute2DCoords(mol)

    img = Draw.MolToImage(mol, size=(width, height))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()
