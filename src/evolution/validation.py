"""Shared validation for seeds, offspring and elites."""
from collections import Counter
from rdkit import Chem
from rdkit.Chem import Descriptors

class MolecularConstraints:
    def __init__(self, max_mw=500, max_logp=5, min_nitrogen=1):
        self.max_mw, self.max_logp, self.min_nitrogen = max_mw, max_logp, min_nitrogen
        self.stats = Counter()

    def validate(self, smiles):
        self.stats["total_checked"] += 1
        reason, canonical = "", None
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None or mol.GetNumAtoms() == 0: reason = "invalid_smiles"
            elif mol.GetStereoGroups(): reason = "enhanced_stereo"
            elif len(Chem.GetMolFrags(mol)) != 1: reason = "disconnected"
            elif any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()): reason = "dummy_atoms"
            elif any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()): reason = "radicals"
            elif Descriptors.MolWt(mol) > self.max_mw: reason = "mw"
            elif Descriptors.MolLogP(mol) > self.max_logp: reason = "logp"
            elif sum(a.GetAtomicNum()==7 for a in mol.GetAtoms()) < self.min_nitrogen: reason = "nitrogen"
            else:
                for a in mol.GetAtoms(): a.SetAtomMapNum(0)
                canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        except Exception:
            reason = "invalid_smiles"
        self.stats["failed_"+reason if reason else "passed"] += 1
        return canonical, reason

    def passes_constraints(self, smiles):
        canonical, reason = self.validate(smiles)
        return canonical is not None, reason

    def get_properties(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return {}
        return {"mw": Descriptors.MolWt(mol), "logp": Descriptors.MolLogP(mol),
                "nitrogen_count": sum(a.GetAtomicNum()==7 for a in mol.GetAtoms()),
                "hbd": Descriptors.NumHDonors(mol), "hba": Descriptors.NumHAcceptors(mol),
                "tpsa": Descriptors.TPSA(mol), "rotatable_bonds": Descriptors.NumRotatableBonds(mol)}

    def print_statistics(self):
        import logging
        logging.getLogger(__name__).info("Validation: %s", dict(self.stats))
