"""Shared, side-effect-free validation of supported ligand graphs."""
from rdkit import Chem

def validated_molecule(smiles):
    if not isinstance(smiles, str):
        raise ValueError("SMILES must be a string")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Invalid or empty SMILES")
    if mol.GetStereoGroups():
        raise ValueError("Enhanced stereochemistry (CXSMILES stereo groups) is unsupported")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Disconnected SMILES: choose the intended ligand/protonation state explicitly")
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError("Unresolved attachment/dummy atoms")
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("Radical ligands are unsupported")
    return mol
