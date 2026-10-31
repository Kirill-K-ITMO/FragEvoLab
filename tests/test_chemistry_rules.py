from dataclasses import replace
import csv
from pathlib import Path
import random

import pytest
from rdkit import Chem
from rdkit.Chem import BRICS
from src.evolution.chemistry_rules import (
    OPERATORS, enumerate_mutations, apply_mutation,
    enumerate_brics_joins, apply_brics_join,
)

def mol(s):
    m = Chem.MolFromSmiles(s)
    assert m is not None
    return m

def smi(m):
    return Chem.MolToSmiles(m, isomericSmiles=True)

def assert_valid(result, *, fragments=False):
    assert result.ok, (result.reason, result.details)
    m = result.mol
    Chem.SanitizeMol(m)
    assert len(Chem.GetMolFrags(m)) == 1
    assert not any(a.GetNumRadicalElectrons() for a in m.GetAtoms())
    if not fragments:
        assert all(a.GetAtomicNum() > 0 for a in m.GetAtoms())
    assert Chem.MolFromSmiles(result.smiles) is not None

def test_add_replaces_explicit_h_and_clears_changed_stereo():
    m = mol('C[C@H](F)Cl')
    before = smi(m)
    p = next(p for p in enumerate_mutations(m, 'add_atom').plans if p.site == (1,) and p.atomic_number == 6)
    out = apply_mutation(m, p)
    assert_valid(out)
    center = out.mol.GetAtomWithIdx(1)
    assert center.GetTotalValence() == 4 and center.GetTotalNumHs() == 0
    assert center.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
    assert out.details['cleared_stereocenters'] == [1]
    assert smi(m) == before

def test_saturated_carbon_not_offered_for_addition_or_bond_promotion():
    m = mol('CC(C)(C)C')
    assert all(p.site != (1,) for p in enumerate_mutations(m, 'add_atom').plans)
    proposals = enumerate_mutations(m, 'change_bond')
    assert not proposals.plans
    assert proposals.rejected['insufficient_replaceable_h'] == 4

def test_substitution_excludes_impossible_oxygen_and_nitrogen():
    m = mol('CC(C)(C)C')
    proposals = enumerate_mutations(m, 'substitute_atom')
    assert not any(p.site == (1,) for p in proposals.plans)
    assert proposals.rejected['predicted_valence_exceeded'] >= 2

@pytest.mark.parametrize('s', ['C[NH3+]', 'C[NH2+]C', 'C[NH+](C)C', 'C[O-]', 'CN(C)C', 'COC', 'CCCl', 'CCBr', 'CCI'])
@pytest.mark.parametrize('op', OPERATORS)
def test_supported_charged_and_neutral_states(s, op):
    m = mol(s)
    charge = sum(a.GetFormalCharge() for a in m.GetAtoms())
    before = smi(m)
    for p in enumerate_mutations(m, op).plans:
        out = apply_mutation(m, p)
        assert_valid(out)
        assert sum(a.GetFormalCharge() for a in out.mol.GetAtoms()) == charge
    assert smi(m) == before

def test_ammonium_explicit_h_replaced_without_neutralization():
    m = mol('C[NH3+]')
    p = next(p for p in enumerate_mutations(m,'add_atom').plans if p.site == (1,) and p.atomic_number == 6)
    out = apply_mutation(m,p)
    assert_valid(out)
    n=out.mol.GetAtomWithIdx(1)
    assert n.GetFormalCharge()==1 and n.GetTotalValence()==4 and n.GetTotalNumHs()==2

def test_quaternary_ammonium_and_phenyl_not_modified_directly():
    m=mol('C[N+](C)(C)Cc1ccccc1')
    for op in OPERATORS:
        for p in enumerate_mutations(m,op).plans:
            if op != 'remove_atom':
                assert 1 not in p.site
            else:
                assert p.site[0] != 1
                out = apply_mutation(m, p)
                assert_valid(out)
                assert sum(a.GetFormalCharge() for a in out.mol.GetAtoms()) == 1
            assert not any(m.GetAtomWithIdx(i).GetIsAromatic() for i in p.site)

def test_explicit_h_atoms_and_isotopes_are_not_silently_removed():
    m=mol('[2H][C@H](F)CC')
    for op in OPERATORS:
        for p in enumerate_mutations(m,op).plans:
            assert 0 not in p.site and 1 not in p.site
            out=apply_mutation(m,p)
            assert_valid(out)
            assert sum(a.GetIsotope()==2 and a.GetAtomicNum()==1 for a in out.mol.GetAtoms())==1
    assert not enumerate_mutations(mol('[13CH4]'),'substitute_atom').plans

def test_remote_tetrahedral_stereo_retained():
    m=mol('C[C@H](F)CCCO')
    p=next(p for p in enumerate_mutations(m,'add_atom').plans if p.site==(5,) and p.atomic_number==6)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    tag=m.GetAtomWithIdx(1).GetChiralTag()
    out=apply_mutation(m,p)
    assert_valid(out)
    assert out.mol.GetAtomWithIdx(1).GetChiralTag()==tag
    assert '@' in out.smiles

def test_defined_double_stereo_protected_and_remote_stereo_kept():
    m=mol('F/C=C/CCCO')
    stereo=next(b for b in m.GetBonds() if b.GetStereo()!=Chem.BondStereo.STEREONONE)
    protected={stereo.GetBeginAtomIdx(),stereo.GetEndAtomIdx(),*stereo.GetStereoAtoms()}
    for op in OPERATORS:
        for p in enumerate_mutations(m,op).plans:
            assert not set(p.site)&protected
            out=apply_mutation(m,p)
            assert_valid(out)
            assert '/' in out.smiles or '\\' in out.smiles

def test_bond_promotion_and_reduction_balance_hydrogens():
    m=mol('CCO')
    p=next(p for p in enumerate_mutations(m,'change_bond').plans if p.site==(0,1))
    out=apply_mutation(m,p)
    assert_valid(out)
    assert out.smiles=='C=CO'
    assert sum(a.GetTotalNumHs() for a in out.mol.GetAtoms())==4
    p2=next(p for p in enumerate_mutations(out.mol,'change_bond').plans if p.site==(0,1))
    back=apply_mutation(out.mol,p2)
    assert_valid(back)
    assert back.smiles=='CCO'

def test_terminal_atom_removal_caps_neighbor_with_h():
    m=mol('CCO')
    p=next(p for p in enumerate_mutations(m,'remove_atom').plans if p.site==(2,1))
    out=apply_mutation(m,p)
    assert_valid(out)
    assert out.smiles=='CC'

def test_stale_and_forged_plans_fail_without_editing():
    m=mol('CCO')
    p=enumerate_mutations(m,'add_atom').plans[0]
    assert apply_mutation(mol('CCN'),p).reason=='stale_plan'
    assert apply_mutation(m,replace(p,atomic_number=92)).reason=='disallowed_plan'
    permuted=Chem.RenumberAtoms(m,[2,1,0])
    assert apply_mutation(permuted,p).reason=='stale_plan'
    assert smi(m)=='CCO'

@pytest.mark.parametrize('s,reason', [('CC.O','disconnected_input'),('[CH3]','radical_input'),('[3*]OCC','dummy_in_atom_mutation')])
def test_unsupported_inputs(s,reason):
    assert enumerate_mutations(mol(s),'add_atom').rejected=={reason:1}

def test_unknown_operator_and_missing_input():
    with pytest.raises(ValueError): enumerate_mutations(mol('CC'),'invent')
    assert enumerate_mutations(None,'add_atom').rejected=={'invalid_input':1}

def test_brics_connection_and_input_immutability():
    a,b=mol('[3*]OCC'),mol('[4*]CC')
    plans=enumerate_brics_joins(a,b)
    assert len(plans.plans)==1
    before=(smi(a),smi(b))
    out=apply_brics_join(a,b,plans.plans[0],require_complete=True)
    assert_valid(out)
    assert out.smiles=='CCOCC'
    assert (smi(a),smi(b))==before
    assert apply_brics_join(a,b,replace(plans.plans[0],bond_order=2)).reason=='disallowed_plan'

def test_brics_incompatible_labels():
    plans=enumerate_brics_joins(mol('[3*]OCC'),mol('[3*]OCC'))
    assert not plans.plans and plans.rejected=={'incompatible_brics_types':1}

def test_brics_double_bond_and_wrong_attachment_order():
    a,b=mol('[7*]=CC'),mol('[7*]=CC')
    plans=enumerate_brics_joins(a,b)
    assert plans.plans
    out=apply_brics_join(a,b,plans.plans[0],require_complete=True)
    assert_valid(out)
    assert out.smiles=='CC=CC'
    wrong=enumerate_brics_joins(mol('[7*]CC'),b)
    assert not wrong.plans and wrong.rejected=={'attachment_bond_order_mismatch':1}

def test_partial_assembly_cannot_be_published_as_complete():
    a,b=mol('[4*]CC[4*]'),mol('[3*]OC')
    p=enumerate_brics_joins(a,b).plans[0]
    assert_valid(apply_brics_join(a,b,p),fragments=True)
    assert apply_brics_join(a,b,p,require_complete=True).reason=='unclosed_attachments'

def test_brics_join_charged_anchor_preserves_charge_and_h():
    a,b=mol('[5*][NH2+]C'),mol('[4*]CC')
    p=enumerate_brics_joins(a,b).plans[0]
    out=apply_brics_join(a,b,p,require_complete=True)
    assert_valid(out)
    n=next(a for a in out.mol.GetAtoms() if a.GetAtomicNum()==7)
    assert n.GetFormalCharge()==1 and n.GetTotalNumHs()==2 and n.GetTotalValence()==4

def test_enhanced_stereochemistry_explicitly_unsupported():
    m=mol('C[C@H](F)Cl |&1:1|')
    assert m.GetStereoGroups()
    assert enumerate_mutations(m,'add_atom').rejected=={'enhanced_stereo_unsupported':1}

def test_all_planned_edits_on_supplied_seeds(capfd):
    path=Path(__file__).parent/'data'/'ligands.csv'
    rows=list(csv.DictReader(path.open(encoding='utf-8-sig')))
    assert len(rows)==13
    applied=0
    for row in rows:
        m=mol(row['smiles'])
        before=smi(m)
        for op in OPERATORS:
            for plan in enumerate_mutations(m,op).plans:
                assert_valid(apply_mutation(m,plan))
                applied+=1
        assert smi(m)==before
    assert applied>500
    assert "Explicit valence" not in capfd.readouterr().err

def test_brics_plans_from_real_fragmentation():
    rows=list(csv.DictReader((Path(__file__).parent/'data'/'ligands.csv').open(encoding='utf-8-sig')))
    fragments=sorted(set(f for row in rows for f in BRICS.BRICSDecompose(mol(row['smiles'])) if '*' in f))
    pairs=[(a,b) for a in fragments for b in fragments]
    random.Random(42).shuffle(pairs)
    applied=0
    for a,b in pairs[:300]:
        left,right=mol(a),mol(b)
        for p in enumerate_brics_joins(left,right).plans:
            assert_valid(apply_brics_join(left,right,p),fragments=True)
            applied+=1
    assert applied>20

def test_old_conformer_not_reused_for_new_graph():
    m = mol('CCO')
    m.AddConformer(Chem.Conformer(m.GetNumAtoms()))
    p = enumerate_mutations(m,'add_atom').plans[0]
    out = apply_mutation(m,p)
    assert_valid(out)
    assert out.mol.GetNumConformers() == 0
    assert m.GetNumConformers() == 1
