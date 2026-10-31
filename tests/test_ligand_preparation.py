from pathlib import Path
import pytest
from src.preparation.ligands import LigandConverter
from src.preparation import ligands as module

@pytest.fixture
def converter(tmp_path):return LigandConverter(tmp_path)

def test_uff_retry_converges_and_metadata(converter,monkeypatch):
    codes=iter([1,0]);calls=[]
    def optimize(mol,**kw):calls.append(kw['maxIters']);return next(codes)
    monkeypatch.setattr(module.AllChem,'UFFOptimizeMolecule',optimize)
    r=converter.prepare_cached('CCO')
    assert calls==[200,1000] and r['preparation_details']['converged']
    assert len(r['preparation_details']['optimization_attempts'])==2
    assert converter.prepare_cached('OCC')['cache_hit']

def test_nonconvergence_no_output_no_cache(converter,monkeypatch):
    monkeypatch.setattr(module.AllChem,'UFFOptimizeMolecule',lambda *a,**kw:1)
    p=converter.output_dir/'bad.pdbqt';p.write_text('old file')
    assert not converter.smiles_to_pdbqt_meeko('CCO','bad',p)
    assert not p.exists()
    with pytest.raises(RuntimeError,match='converge'):converter.prepare_cached('CCO')
    assert not converter.calculation_cache.statistics().get('preparation')

def test_missing_uff_parameters(converter,monkeypatch):
    monkeypatch.setattr(module.AllChem,'UFFHasAllMoleculeParams',lambda mol:False)
    with pytest.raises(RuntimeError,match='parameters'):converter.prepare_cached('CCO')

def test_embedding_failure(converter,monkeypatch):
    monkeypatch.setattr(module.AllChem,'EmbedMolecule',lambda *a,**kw:-1)
    with pytest.raises(RuntimeError,match='embedding'):converter.prepare_cached('CCO')
    assert len(converter.last_preparation_details['embedding_attempts'])==2

@pytest.mark.parametrize('smiles',['invalid','','CCN.CCO','[1*]CCN','[CH2]CN'])
def test_graph_gate(converter,smiles):
    with pytest.raises(ValueError):converter.prepare_cached(smiles)

def test_meeko_failed_writer_no_artifact(converter,monkeypatch):
    monkeypatch.setattr(module.PDBQTWriterLegacy,'write_string',lambda *a:('garbage',False,'writer failure'))
    with pytest.raises(RuntimeError,match='writer failure'):converter.prepare_cached('CCO')

@pytest.mark.parametrize('name',['../outside','CON','NUL.txt','LPT1','a:b','x.',' ','a/b','a\\b'])
def test_windows_safe_names(converter,name):
    with pytest.raises(ValueError):converter.convert_single(name,'CCO')

@pytest.mark.parametrize('header,delimiter',[('название',','),('name',';'),('name','\t')])
def test_bom_csv_and_delimiters(converter,tmp_path,header,delimiter):
    p=tmp_path/'input.csv';p.write_text(f'{header}{delimiter}smiles\nethanol{delimiter}CCO\n',encoding='utf-8-sig')
    r=converter.convert_from_csv(p)
    assert r['success']==1 and r['failed']==0 and not r['input_error']

def test_smiles_only_csv(converter,tmp_path):
    p=tmp_path/'input.csv';p.write_text('smiles\nCCO\nCCN\n')
    r=converter.convert_from_csv(p)
    assert r['success']==2 and 'ligand_000001' in converter.prepared_results

def test_duplicate_names_and_bad_rows_do_not_abort(converter):
    r=converter.convert_list([('same','CCO'),('SAME','CCN'),('bad/name','CCO'),('good','CCN')])
    assert r['total']==4 and r['failed']==3 and r['success']==1
    assert not (converter.output_dir/'same.pdbqt').exists()

def test_missing_csv_resets_previous_stats(converter,tmp_path):
    converter.convert_list([('first','CCO')])
    r=converter.convert_from_csv(tmp_path/'missing.csv')
    assert r['total']==r['success']==0 and r['input_error']

def test_empty_csv_error(converter,tmp_path):
    p=tmp_path/'empty.csv';p.touch()
    assert converter.convert_from_csv(p)['input_error']=='Empty CSV'

def test_changed_optimization_policy_invalidates_cache(converter):
    first=converter.prepare_cached('CCO');converter.optimization_iters=(300,1500)
    second=converter.prepare_cached('CCO')
    assert first['preparation_id']!=second['preparation_id']

def test_obabel_timeout_cleans_output(converter,monkeypatch):
    import subprocess
    def timeout(*a,**kw):
        assert kw['timeout']==180
        assert '-xr' not in a[0]
        raise subprocess.TimeoutExpired(a[0],180)
    monkeypatch.setattr(module.subprocess,'run',timeout)
    p=converter.output_dir/'old.pdbqt';p.write_text('old')
    assert not converter.smiles_to_pdbqt_obabel('CCO','old',p)
    assert not p.exists()

def test_removed_modules_and_menu():
    from src.ui import application as main
    assert not hasattr(main.FragEvoLabApp,'prepare_proteins')
    assert not hasattr(main.FragEvoLabApp,'view_protein_info')
    assert not hasattr(main.FragEvoLabApp,'visualize_metabolism')
    root=Path(main.__file__).resolve().parents[2]
    assert not (root/'src/protein_prep.py').exists()
    assert not (root/'src/metabolism_visualizer.py').exists()


def test_report_and_menu_share_csv_names(converter,tmp_path):
    from src.reporting.docking import ReportGenerator
    from src.core.project import ProjectUtils
    p=tmp_path/'input.csv';p.write_text('name;smiles\nethanol;CCO\n',encoding='utf-8-sig')
    assert converter.convert_from_csv(p)['success']==1
    assert ReportGenerator(tmp_path).load_smiles_data(p)=={'ethanol':'CCO'}
    assert ProjectUtils(tmp_path).load_ligands_from_csv(p)==[{'name':'ethanol','smiles':'CCO'}]

def test_seed_loader_accepts_same_csv_format(tmp_path):
    from src.evolution.engine import EvolutionaryGenerator
    from src.evolution import EvolutionConfig, MolecularConstraints
    from src.evolution.fragments import FragmentBank
    import random
    g=EvolutionaryGenerator(tmp_path);g.config=EvolutionConfig(targets=[],min_nitrogen=0)
    g.constraints=MolecularConstraints(min_nitrogen=0);g.bank=FragmentBank(g.config);g.rng=random.Random(42);g.attempts=[]
    p=tmp_path/'input.csv';p.write_text('name;smiles\nethanol;CCO\n',encoding='utf-8-sig')
    assert g._load_seeds(p)[0].smiles=='CCO'
