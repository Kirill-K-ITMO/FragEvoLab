"""Conservative, prevalidated graph edits for FragEvoLab.

No random selection, docking, cache or evolution state is touched here. Select an
operator externally, enumerate its plans, select a plan externally, then apply.
Final RDKit sanitization remains mandatory. See README.md for the supported scope.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
import hashlib
import json
from typing import Any

from rdkit import Chem
from rdkit.Chem import BRICS

RULES_VERSION = "1.0.0"
# Deliberate mutation policy, not a complete periodic-table valence model.
# Atom identity is (atomic number, formal charge); only closed-shell states.
VALENCE_STATES = {
    (6, 0): 4, (7, 0): 3, (8, 0): 2,
    (9, 0): 1, (17, 0): 1, (35, 0): 1, (53, 0): 1,
    (7, 1): 4, (8, -1): 1,
}
ADD_ELEMENTS = (6, 7, 8, 9)
SUBSTITUTE_ELEMENTS = (6, 7, 8)
OPERATORS = ("add_atom", "remove_atom", "substitute_atom", "change_bond")
BOND_TYPES = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}

@dataclass(frozen=True)
class MutationPlan:
    operator: str
    source_key: str
    site: tuple[int, ...]
    atomic_number: int | None = None
    bond_order: int | None = None

@dataclass(frozen=True)
class JoinPlan:
    left_key: str
    right_key: str
    left_dummy: int
    right_dummy: int
    bond_order: int

@dataclass
class PlanSet:
    plans: list[MutationPlan | JoinPlan] = field(default_factory=list)
    rejected: dict[str, int] = field(default_factory=dict)

@dataclass
class EditResult:
    mol: Any = None
    smiles: str | None = None
    reason: str = ""
    details: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.mol is not None and self.smiles is not None and not self.reason

def _copy_checked(mol, *, fragments=False):
    if mol is None or not isinstance(mol, Chem.Mol) or mol.GetNumAtoms() == 0:
        return None, "invalid_input"
    copy = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(copy)
    except (ValueError, RuntimeError):
        return None, "invalid_input"
    if len(Chem.GetMolFrags(copy)) != 1:
        return None, "disconnected_input"
    if any(a.GetNumRadicalElectrons() for a in copy.GetAtoms()):
        return None, "radical_input"
    if not fragments and any(a.GetAtomicNum() == 0 for a in copy.GetAtoms()):
        return None, "dummy_in_atom_mutation"
    if copy.GetStereoGroups():
        return None, "enhanced_stereo_unsupported"
    return copy, ""

def _key(mol):
    # Atom order is intentional: plans contain indices, not canonical atom IDs.
    atoms = [(a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope(),
              a.GetNumExplicitHs(), a.GetNumImplicitHs(), a.GetNoImplicit(),
              str(a.GetChiralTag()), a.GetIsAromatic(), a.GetAtomMapNum(),
              a.GetNumRadicalElectrons()) for a in mol.GetAtoms()]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType()),
              str(b.GetStereo()), tuple(b.GetStereoAtoms()), str(b.GetBondDir()))
             for b in mol.GetBonds()]
    return hashlib.sha256(json.dumps((atoms, bonds)).encode()).hexdigest()

def _bond_sum(atom):
    orders = [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
    if any(o not in (1, 2, 3) for o in orders):
        return None
    return int(sum(orders))

def _atom_problem(atom):
    if atom.GetIsAromatic():
        return "aromatic_site_unsupported"
    if atom.GetIsotope():
        return "isotope_site_protected"
    if any(a.GetAtomicNum() == 1 for a in atom.GetNeighbors()):
        return "explicit_h_atom_protected"
    target = VALENCE_STATES.get((atom.GetAtomicNum(), atom.GetFormalCharge()))
    if target is None:
        return "unsupported_atom_state"
    bond_sum = _bond_sum(atom)
    if bond_sum is None or bond_sum + atom.GetTotalNumHs() != target:
        return "unsupported_valence_state"
    return ""

def _stereo_bond_touched(mol, sites):
    sites = set(sites)
    for bond in mol.GetBonds():
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            involved = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), *bond.GetStereoAtoms()}
            if involved & sites:
                return True
    return False

def enumerate_mutations(mol, operator: str) -> PlanSet:
    """List supported edits without creating candidate molecules.

    Counts in rejected describe excluded proposals/sites, not failed molecules.
    Charged sites can add/remove H-replacing substituents, but element substitution
    is neutral-only. Aromatic sites and explicit hydrogen ATOMS are protected.
    """
    if operator not in OPERATORS:
        raise ValueError(f"Unknown operator: {operator}")
    mol, error = _copy_checked(mol)
    if error:
        return PlanSet(rejected={error: 1})
    plans, rejected = [], Counter()
    key = _key(mol)

    def eligible(sites):
        for i in sites:
            problem = _atom_problem(mol.GetAtomWithIdx(i))
            if problem:
                rejected[problem] += 1
                return False
        if _stereo_bond_touched(mol, sites):
            rejected["stereo_bond_protected"] += 1
            return False
        return True

    if operator == "change_bond":
        for bond in mol.GetBonds():
            if bond.GetIsAromatic() or bond.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                rejected["unsupported_bond"] += 1
                continue
            sites = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            if not eligible(sites):
                continue
            order = 2 if bond.GetBondType() == Chem.BondType.SINGLE else 1
            if order == 2 and any(mol.GetAtomWithIdx(i).GetTotalNumHs() < 1 for i in sites):
                rejected["insufficient_replaceable_h"] += 1
                continue
            plans.append(MutationPlan(operator, key, sites, bond_order=order))
    else:
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if not eligible((idx,)):
                continue
            if operator == "add_atom":
                if atom.GetTotalNumHs() < 1:
                    rejected["insufficient_replaceable_h"] += 1
                    continue
                for z in ADD_ELEMENTS:
                    plans.append(MutationPlan(operator, key, (idx,), atomic_number=z))
            elif operator == "substitute_atom":
                if atom.GetFormalCharge() != 0:
                    rejected["charged_substitution_unsupported"] += 1
                    continue
                for z in SUBSTITUTE_ELEMENTS:
                    if z == atom.GetAtomicNum():
                        continue
                    if _bond_sum(atom) > VALENCE_STATES[(z, 0)]:
                        rejected["predicted_valence_exceeded"] += 1
                        continue
                    plans.append(MutationPlan(operator, key, (idx,), atomic_number=z))
            else:
                if atom.GetDegree() != 1 or mol.GetNumHeavyAtoms() <= 2:
                    rejected["not_removable_terminal"] += 1
                    continue
                if atom.GetFormalCharge() != 0:
                    rejected["charged_deletion_unsupported"] += 1
                    continue
                neighbor = atom.GetNeighbors()[0].GetIdx()
                if eligible((neighbor,)):
                    plans.append(MutationPlan(operator, key, (idx, neighbor)))
    return PlanSet(plans, dict(sorted(rejected.items())))

def _clear_local_tetrahedral(mol, sites):
    cleared = []
    for idx in sorted(set(sites)):
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            cleared.append(idx)
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        for prop in ("_CIPCode", "_CIPRank"):
            if atom.HasProp(prop):
                atom.ClearProp(prop)
    return cleared

def _set_h(atom, count):
    # Bracket H and implicit H are pooled, then represented explicitly.
    # Hydrogen atoms/isotopes are not deleted by this helper.
    if count < 0:
        raise ValueError("Negative hydrogen count")
    atom.SetNumExplicitHs(int(count))
    atom.SetNoImplicit(True)

def _finalize(mol, details, *, allow_dummies=False):
    try:
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        if len(Chem.GetMolFrags(mol)) != 1:
            return EditResult(reason="disconnected_product", details=details)
        if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
            return EditResult(reason="radical_product", details=details)
        if not allow_dummies and any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
            return EditResult(reason="unclosed_attachments", details=details)
        mol.RemoveAllConformers()  # Old coordinates do not describe the edited graph.
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        return EditResult(mol, smiles, details=details)
    except (ValueError, RuntimeError) as error:
        return EditResult(reason="final_sanitization_failed", details={**details, "error": str(error)})

def apply_mutation(mol, plan: MutationPlan) -> EditResult:
    """Revalidate a plan, edit a copy, sanitize. Never modifies the input.

    Tetrahedral stereo on edited atoms is explicitly cleared (not guessed).
    Unaffected stereo is retained. Plans cannot be used on reordered molecules.
    """
    if not isinstance(plan, MutationPlan) or plan.operator not in OPERATORS:
        return EditResult(reason="invalid_plan")
    mol, error = _copy_checked(mol)
    if error:
        return EditResult(reason=error)
    if plan.source_key != _key(mol):
        return EditResult(reason="stale_plan")
    if plan not in enumerate_mutations(mol, plan.operator).plans:
        return EditResult(reason="disallowed_plan")
    rw = Chem.RWMol(mol)
    sites = plan.site
    details = {"rules_version": RULES_VERSION, "operator": plan.operator,
               "site": list(sites), "before": Chem.MolToSmiles(mol),
               "cleared_stereocenters": _clear_local_tetrahedral(rw, sites)}
    if plan.operator == "add_atom":
        atom = rw.GetAtomWithIdx(sites[0])
        _set_h(atom, mol.GetAtomWithIdx(sites[0]).GetTotalNumHs() - 1)
        new = Chem.Atom(plan.atomic_number)
        _set_h(new, VALENCE_STATES[(plan.atomic_number, 0)] - 1)
        idx = rw.AddAtom(new)
        rw.AddBond(sites[0], idx, Chem.BondType.SINGLE)
        details["atomic_number"] = plan.atomic_number
    elif plan.operator == "substitute_atom":
        atom = rw.GetAtomWithIdx(sites[0])
        atom.SetAtomicNum(plan.atomic_number)
        _set_h(atom, VALENCE_STATES[(plan.atomic_number, 0)] - _bond_sum(atom))
        details["atomic_number"] = plan.atomic_number
    elif plan.operator == "remove_atom":
        removed, retained = sites
        old_order = int(mol.GetBondBetweenAtoms(removed, retained).GetBondTypeAsDouble())
        _set_h(rw.GetAtomWithIdx(retained), mol.GetAtomWithIdx(retained).GetTotalNumHs() + old_order)
        rw.RemoveAtom(removed)
    else:
        bond = rw.GetBondBetweenAtoms(*sites)
        delta = plan.bond_order - int(bond.GetBondTypeAsDouble())
        for idx in sites:
            _set_h(rw.GetAtomWithIdx(idx), mol.GetAtomWithIdx(idx).GetTotalNumHs() - delta)
        bond.SetBondType(BOND_TYPES[plan.bond_order])
        details["bond_order"] = plan.bond_order
    return _finalize(rw.GetMol(), details)

def _brics_rules():
    rules = {}
    for group in BRICS.reactionDefs:
        for left, right, symbol in group:
            if symbol not in ("-", "="):
                raise RuntimeError(f"Unsupported BRICS bond symbol: {symbol}")
            a = int("".join(c for c in left if c.isdigit()))
            b = int("".join(c for c in right if c.isdigit()))
            order = 2 if symbol == "=" else 1
            rules[(a, b)] = rules[(b, a)] = order
    return rules

def _attachment_problem(atom):
    if atom.GetAtomicNum() != 0 or atom.GetDegree() != 1:
        return "invalid_attachment"
    if atom.GetFormalCharge() != 0 or atom.GetTotalNumHs() != 0:
        return "invalid_attachment"
    if atom.GetNeighbors()[0].GetAtomicNum() in (0, 1):
        return "invalid_attachment_neighbor"
    return ""

def enumerate_brics_joins(left, right) -> PlanSet:
    """Use BRICS types and replace each dummy bond with one equal-order bond.

    Matching bond order preserves each anchor's bond-order sum and H count.
    Existing valid charged/aromatic/S/P environments can thus be retained
    without inventing additional valence states. BRICS labels are trusted as
    provenance from RDKit BRICS fragmentation, not inferred from arbitrary '*'.
    """
    left, error = _copy_checked(left, fragments=True)
    if error:
        return PlanSet(rejected={"left_" + error: 1})
    right, error = _copy_checked(right, fragments=True)
    if error:
        return PlanSet(rejected={"right_" + error: 1})
    ld = [a for a in left.GetAtoms() if a.GetAtomicNum() == 0]
    rd = [a for a in right.GetAtoms() if a.GetAtomicNum() == 0]
    if not ld or not rd:
        return PlanSet(rejected={"no_attachment": 1})
    rules, plans, rejected = _brics_rules(), [], Counter()
    lk, rk = _key(left), _key(right)
    for a in ld:
        for b in rd:
            problem = _attachment_problem(a) or _attachment_problem(b)
            if problem:
                rejected[problem] += 1
                continue
            order = rules.get((a.GetIsotope(), b.GetIsotope()))
            if order is None:
                rejected["incompatible_brics_types"] += 1
                continue
            if any(x.GetBonds()[0].GetBondTypeAsDouble() != order for x in (a, b)):
                rejected["attachment_bond_order_mismatch"] += 1
                continue
            if (_stereo_bond_touched(left, (a.GetIdx(), a.GetNeighbors()[0].GetIdx())) or
                    _stereo_bond_touched(right, (b.GetIdx(), b.GetNeighbors()[0].GetIdx()))):
                rejected["stereo_bond_protected"] += 1
                continue
            plans.append(JoinPlan(lk, rk, a.GetIdx(), b.GetIdx(), order))
    return PlanSet(plans, dict(sorted(rejected.items())))

def apply_brics_join(left, right, plan: JoinPlan, *, require_complete=False) -> EditResult:
    """Connect distinct fragment copies; reject a stale or forged plan.

    require_complete=False allows intermediates with remaining BRICS points.
    Final candidates MUST be requested with require_complete=True or pass the
    generator's no-dummy validation gate before preparation/docking.
    """
    if not isinstance(plan, JoinPlan):
        return EditResult(reason="invalid_plan")
    left, error = _copy_checked(left, fragments=True)
    if error:
        return EditResult(reason="left_" + error)
    right, error = _copy_checked(right, fragments=True)
    if error:
        return EditResult(reason="right_" + error)
    if (_key(left), _key(right)) != (plan.left_key, plan.right_key):
        return EditResult(reason="stale_plan")
    if plan not in enumerate_brics_joins(left, right).plans:
        return EditResult(reason="disallowed_plan")
    n = left.GetNumAtoms()
    li, ri = plan.left_dummy, n + plan.right_dummy
    la = left.GetAtomWithIdx(plan.left_dummy).GetNeighbors()[0].GetIdx()
    ra = n + right.GetAtomWithIdx(plan.right_dummy).GetNeighbors()[0].GetIdx()
    rw = Chem.RWMol(Chem.CombineMols(left, right))
    details = {"rules_version": RULES_VERSION, "operator": "join_fragments",
               "attachments": [plan.left_dummy, plan.right_dummy], "bond_order": plan.bond_order,
               "cleared_stereocenters": _clear_local_tetrahedral(rw, (la, ra)),
               "before": [Chem.MolToSmiles(left), Chem.MolToSmiles(right)]}
    # Remove the old placeholder bonds BEFORE adding the actual chemical bond.
    rw.RemoveBond(li, la)
    rw.RemoveBond(ri, ra)
    rw.AddBond(la, ra, BOND_TYPES[plan.bond_order])
    for idx in sorted((li, ri), reverse=True):
        rw.RemoveAtom(idx)
    return _finalize(rw.GetMol(), details, allow_dummies=not require_complete)
