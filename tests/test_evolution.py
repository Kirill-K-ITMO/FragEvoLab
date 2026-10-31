import json
import random
import pytest
from rdkit import Chem
from src.evolution import EvolutionConfig, ProteinSpec, FitnessEvaluator, MolecularConstraints
from src.evolution.fragments import FragmentBank
from src.evolution.adaptation import Adaptation
from src.evolution.operators import Operators, join
from src.evolution.models import Candidate
from src.evolution.selection import select

def cfg(**kw):
    return EvolutionConfig(targets=[ProteinSpec("A")],min_nitrogen=0,**kw).validate()

@pytest.mark.parametrize("kwargs",[
    {"targets":[]}, {"targets":[ProteinSpec("A")],"off_targets":[ProteinSpec("A",threshold=-5)]},
    {"off_targets":[ProteinSpec("B")]}, {"mutation_steps":0},
    {"level_probabilities":{"atom":0.1,"fragment":0.1}}, {"prior_strength":0},
    {"elite_fraction":1}, {"penalty_strength":float("nan")},
    {"max_similarity":1.1}, {"targets":[ProteinSpec("../A")]},
])
def test_reject_invalid_configuration(kwargs):
    base = {"targets":[ProteinSpec("A")]}; base.update(kwargs)
    with pytest.raises(ValueError): EvolutionConfig(**base).validate()

def test_multitarget_offtarget_penalty_and_saturation():
    c = EvolutionConfig(targets=[ProteinSpec("A",threshold=-6),ProteinSpec("B",weight=3,threshold=-5)],
        off_targets=[ProteinSpec("C",threshold=-4),ProteinSpec("D",weight=2,threshold=-3)])
    ev = FitnessEvaluator(c)
    s = ev.evaluate({"A":-8,"B":-6,"C":-5,"D":-5})
    assert s.feasible and s.target_utility == 6.5 and s.off_penalty == 4 and s.fitness == 2.5
    good = ev.evaluate({"A":-8,"B":-6,"C":-3,"D":-2})
    assert good.fitness == ev.evaluate({"A":-8,"B":-6,"C":10,"D":20}).fitness == 6.5
    assert not ev.evaluate({"A":-8,"B":-6,"C":-3}).complete
    assert not ev.evaluate({"A":-8,"B":-4,"C":-3,"D":-2}).feasible
    assert not ev.evaluate({"A":-8,"B":-6,"C":float("nan"),"D":-2}).complete

def test_normalization_and_worst_target():
    c = EvolutionConfig(targets=[ProteinSpec("A",offset=-5,scale=2),ProteinSpec("B",offset=-4)],target_aggregation="worst")
    assert FitnessEvaluator(c).evaluate({"A":-9,"B":-5}).fitness == 1

@pytest.mark.parametrize("smiles,reason",[("C.C","disconnected"),("[1*]CC","dummy_atoms"),("[CH3]","radicals"),("not smiles","invalid_smiles")])
def test_shared_validation(smiles,reason):
    assert MolecularConstraints(min_nitrogen=0).validate(smiles)[1] == reason

def test_canonical_identity():
    v = MolecularConstraints(min_nitrogen=0)
    assert v.validate("OCC")[0] == v.validate("CCO")[0]
    assert v.validate("[CH3:7]CO")[0] == "CCO"
    assert v.validate("C[C@H](O)F")[0] != v.validate("C[C@@H](O)F")[0]

def test_unique_observation_and_matched_credit():
    b = FragmentBank(cfg()); parent,child = "CCOC(=O)NCC", "CCNC(=O)NCC"
    b.observe(parent,{"A":-7},0); b.observe(parent,{"A":-7},1)
    for r in b.records.values(): assert r["association"]["A"]["n"] == 1 and r["evaluated"] == 0
    added=b.add("[5*]NCC",child,1); removed=b.add("[3*]OCC",parent,0)
    events=[{"operator":"replace_fragment","added":[added],"removed":[removed]}]
    assert b.credit(parent,child,events,0.5,{"A":-7},{"A":-8})
    assert not b.credit(parent,child,events,0.5,{"A":-7},{"A":-8})
    b.update_weights()
    assert b.records[added]["evaluated"] == 1 and b.records[added]["weight"] > 1
    assert b.records[added]["energy_delta_sum"]["A"] == -1
    assert b.records[removed]["removed_evaluated"] == 1
    assert sum(b.weights().values()) == pytest.approx(1)

def test_chain_no_false_credit():
    b = FragmentBank(cfg()); key=b.add("[5*]NCC","CCN",0)
    event={"operator":"replace_fragment","added":[key]}
    b.credit("CCN","CCNC",[event,event],0.5,{"A":-7},{"A":-8})
    assert b.records[key]["evaluated"] == 0

def test_operator_zero_and_adaptation():
    a=Adaptation(cfg(adaptive_operators=True,min_operator_observations=2))
    for _ in range(10): a.record("fragment","reward",1)
    assert a.weights({"atom":1.0,"fragment":0.0}) == {"atom":1.0}
    w=a.weights({"atom":0.5,"fragment":0.5})
    assert w["fragment"]>w["atom"]>0 and sum(w.values()) == pytest.approx(1)

def test_attachment_compatibility():
    assert join(Chem.MolFromSmiles("[3*]OCC"),Chem.MolFromSmiles("[4*]CCC")) == ["CCCOCC"]
    assert join(Chem.MolFromSmiles("[3*]OCC"),Chem.MolFromSmiles("[3*]OCC")) == []

def test_local_replacement():
    c=cfg(level_probabilities={"atom":0.0,"fragment":1.0}); bank=FragmentBank(c)
    for s in ("CCOC(=O)NCC","CCNC(=O)NCCC","CCCOCC"): bank.ingest(s,0)
    op=Operators(random.Random(42),c,bank,Adaptation(c)); success=[]
    for _ in range(30):
        product,events=op.mutate("CCOC(=O)NCC",1)
        assert events[0]["operator"] == "replace_fragment"
        if product:
            assert MolecularConstraints(min_nitrogen=0).validate(product)[0]
            assert len(events[0]["added"]) == len(events[0]["removed"]) == 1
            success.append(product)
    assert success

def test_seeded_manual_operators():
    c=cfg(level_probabilities={"atom":1.0,"fragment":0.0},atom_probabilities={"add_atom":1.,"remove_atom":0.,"substitute_atom":0.,"change_bond":0.})
    runs=[]
    for _ in range(2):
        o=Operators(random.Random(17),c,FragmentBank(c),Adaptation(c))
        runs.append([o.mutate("CCN",1) for _ in range(20)])
    assert runs[0] == runs[1]
    assert all(e[0]["operator"] == "add_atom" for _,e in runs[0])

def test_selection_complete_feasible_unique():
    def candidate(s,f,feasible=True,complete=True):
        return Candidate(s,0,score={"fitness":f,"feasible":feasible,"complete":complete})
    candidates=[candidate("CCO",9),candidate("CCO",8),candidate("CCN",100,False),candidate("CCC",99,complete=False)]
    assert [p.smiles for p in select(candidates,4)] == ["CCO"]


def test_frozen_weights_until_next_generation():
    c=cfg(adaptive_operators=True,min_operator_observations=1)
    b=FragmentBank(c); key=b.add("[5*]NCC","CCN",0)
    a=Adaptation(c); op=Operators(random.Random(1),c,b,a)
    before=dict(op.levels)
    for _ in range(10): a.record("fragment","reward",1)
    assert op.levels == before
    op.freeze()
    assert op.levels["fragment"] > before["fragment"]


def test_assembly_is_bounded_and_connected():
    c=cfg(); b=FragmentBank(c)
    for s in ("CCOC(=O)NCC","CCNC(=O)NCCC","CCCOCC"): b.ingest(s,0)
    op=Operators(random.Random(7),c,b,Adaptation(c)); good=[]
    for _ in range(50):
        product,events=op.assemble()
        assert len(events[0]["added"]) <= 5
        if product:
            mol=Chem.MolFromSmiles(product)
            assert len(Chem.GetMolFrags(mol)) == 1 and all(a.GetAtomicNum()>0 for a in mol.GetAtoms())
            good.append(product)
    assert good

def test_offline_reweight_changes_roles_without_docking(tmp_path):
    from src.evolution.reweight import reweight
    checkpoint=tmp_path/"checkpoint.json"
    checkpoint.write_text(json.dumps({"archive":{
        "CCN":Candidate("CCN",0,energies={"A":-9,"B":-5}).to_dict(),
        "CCO":Candidate("CCO",0,energies={"A":-5,"B":-9}).to_dict()}}))
    c=EvolutionConfig(targets=[ProteinSpec("B")],off_targets=[ProteinSpec("A",threshold=-5)])
    r=reweight(checkpoint,c,tmp_path/"new_scores.json")
    assert r["top_molecules"][0]["smiles"] == "CCO" and r["incomplete"] == 0
    c.off_targets.append(ProteinSpec("C",threshold=-5))
    r=reweight(checkpoint,c,tmp_path/"new_scores.json")
    assert r["incomplete"] == 2 and not r["top_molecules"]
    with pytest.raises(ValueError): reweight(checkpoint,c,checkpoint)
