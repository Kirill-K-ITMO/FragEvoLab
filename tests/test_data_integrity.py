"""Data integrity, reporting, and candidate failure behavior."""
import hashlib, json, math, os, subprocess, sys, zipfile
from pathlib import Path
from types import SimpleNamespace as NS
import pytest
from src.core.calculation_cache import molecule_identity
from src.core.database import DockingDatabase
from src.core.project import ProjectUtils
from src.docking.affinity import AffinityCalculator
from src.evolution.config import EvolutionConfig
from src.evolution.validation import MolecularConstraints
from src.preparation.inputs import read_ligand_table
from src.preparation.smiles import SmilesGenerator
from src.reporting.docking import ReportGenerator
from src.ui.application import FragEvoLabApp
from src.ui.docking_menu import docking_menu


def test_history_and_provenance(tmp_path,monkeypatch):
    db=DockingDatabase(tmp_path/'history.db')
    try:
        app=NS(database=db,calculator=AffinityCalculator(),reporter=ReportGenerator(tmp_path))
        for smi,energy in [('CCO',-5),('CCN',-9),('CCN',-8)]:
            exp=db.add_experiment('same',{})
            data={'affinity':energy,'smiles':smi,'evaluation_id':'proof','artifact_sha':'sha','input_files':{'receptor':'output/inputs/hash.pdbqt'},'output_file':'pose.pdbqt','spec':{'options':{'scoring':'vinardo','seed':'123'}}}
            FragEvoLabApp._save_results_to_db(app,exp,{'same':{'A':data}},8)
        row=db.get_experiment_results(exp)[0];extra=json.loads(row['extra_json'])
        assert row['method']=='vinardo' and row['pdbqt_file']==data['input_files']['receptor']
        for k in ('evaluation_id','artifact_sha','input_files','output_file','spec'):assert extra[k]==data[k]
        captured={}
        def capture(data,smiles,*args,**kw):captured.update(data=data,smiles=smiles)
        monkeypatch.setattr(app.reporter,'generate_html_report',capture)
        app.reporter.generate_db_report(db.export_for_report())
        assert len(captured['data'])==3
        assert {(captured['smiles'][n],v['A']['affinity']) for n,v in captured['data'].items()}=={('CCO',-5),('CCN',-9),('CCN',-8)}
    finally:db.close()


def test_html_escaping_zero_and_sparse_count(tmp_path,monkeypatch):
    r=ReportGenerator(tmp_path);name='<script>window.testMarker=1</script>';protein='\" onmouseover=\"alert(1) <R>&'
    monkeypatch.setattr(r,'calculate_molecular_descriptors',lambda s:dict(formula='X',mw=0,logp=0,tpsa=0,hbd=0,hba=0,sa_score=0))
    data={name:{protein:{'affinity':-7,'smiles':'CCO','ki':0,'pki':0}},'two':{'B':{'affinity':-6}}}
    text=r.generate_html_report(data,{},title=name).read_text()
    assert name not in text and '&lt;script&gt;' in text
    assert 'data-p0-ki="0"' in text and '<option value="p0">' in text
    assert 'id="filterSa"' in text and 'data-logp="0"' in text
    assert '<div class="number">2</div><div class="label">Docking runs' in text
    assert '<td class="descriptor-value">0.00</td>' in text
    monkeypatch.setattr(r,'calculate_molecular_descriptors',lambda s:dict(formula='X',mw=None,logp=float('nan'),tpsa=None,hbd=None,hba=None,sa_score=None))
    assert 'data-logp="9999"' in r.generate_html_report(data,{}).read_text()


def test_duplicate_labels(tmp_path,monkeypatch):
    (tmp_path/'data').mkdir();(tmp_path/'data/ligands.csv').write_text('name,smiles\nx,CCO\n0003_1_x,CCC\nx,CCN\n')
    submitted={}
    def dock(**kw):submitted.update(kw['ligand_paths']);return {}
    app=NS(project_root=tmp_path,docker=NS(get_proteins=lambda:[Path('A.pdbqt')],get_config=lambda n:Path('A.txt'),run_docking=dock,_save_results_json=lambda:None),converter=NS(prepare_cached=lambda s:{'path':s}),database=NS(add_experiment=lambda **kw:1,finish_experiment=lambda e:None),_save_results_to_db=lambda *a:None)
    answers=iter(['all','all','1','1','no','yes','no']);terminal=NS(input=lambda *a:next(answers),print=lambda *a,**kw:None)
    monkeypatch.setattr('src.core.preflight.check_run',lambda *a,**k:NS(ready=True));monkeypatch.setattr('src.core.preflight.display',lambda *a:None)
    docking_menu(app,terminal)
    assert len(submitted)==3 and set(submitted.values())=={'CCO','CCC','CCN'}


@pytest.mark.parametrize('smi',['C[C@H](F)Cl |&1:1|','C[C@H](F)Cl |o1:1|'])
def test_enhanced_stereo(project,smi):
    with pytest.raises(ValueError,match='Enhanced stereo'):molecule_identity(smi)
    with pytest.raises(ValueError,match='Enhanced stereo'):project[1].prepare_cached(smi)
    assert MolecularConstraints(min_nitrogen=0).validate(smi)==(None,'enhanced_stereo')
    assert molecule_identity('CC(F)Cl')[0]


@pytest.mark.parametrize('rows',[[('same','')],[('same','CCC'),('same','CCN')]])
def test_invalid_batch_removes_alias(project,rows):
    c=project[1];assert c.convert_single('same','CCO')
    assert c.convert_list(rows)['failed']==len(rows)
    assert not (c.output_dir/'same.pdbqt').exists() and not (c.output_dir/'same.json').exists()


@pytest.mark.parametrize('contents',['name;smiles;source\nold;CCO;PubChem\n','source\tsmiles\tname\nPubChem\tCCO\told\n','smiles\nCCO\n','smiles,name,source\nCCO,old,PubChem\n'])
def test_append(tmp_path,contents):
    p=tmp_path/'ligands.csv';p.write_text(contents);g=SmilesGenerator(tmp_path);g.generated_smiles=['CCN','CCC']
    assert g.save_to_csv(output_file=p,append=True)
    fields,rows=read_ligand_table(p)
    assert [r['smiles'] for r in rows]==['CCO','CCN','CCC'] and len({r['name'] for r in rows})==3
    if 'source' in fields:assert rows[0]['source']=='PubChem'


def test_append_error_and_template_reset(tmp_path):
    p=tmp_path/'ligands.csv';p.write_text('name,smiles\na,CCO,extra\n');old=p.read_bytes()
    g=SmilesGenerator(tmp_path);g.generated_smiles=['CN']
    with pytest.raises(ValueError):g.save_to_csv(output_file=p,append=True)
    assert p.read_bytes()==old
    g.set_template('C{1}');g.add_radicals('1',['N','O']);g.add_exclusion_rule({'1':'N'});g.generate_combinations()
    g.set_template('CC{1}')
    assert not g.generated_smiles and not g.combinations_info and not g.exclusion_rules
    g.add_radicals('1',['N','O']);assert set(g.generate_combinations())=={'CCN','CCO'}


def test_export(tmp_path):
    for name in ('runs/run/results.json','inputs/receptor.pdbqt','prepared/key/ligand.pdbqt','docking_results/pose.pdbqt'):
        p=tmp_path/'output'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(name)
    db=DockingDatabase(tmp_path/'output/experiments.db')
    try:
        db.add_experiment('snapshot',{})
        with zipfile.ZipFile(ProjectUtils(tmp_path).export_results_to_zip()) as z:
            for n,m in json.loads(z.read('MANIFEST.json'))['files'].items():assert hashlib.sha256(z.read(n)).hexdigest()==m['sha256']
            for n in ('output/experiments.db','output/runs/run/results.json','output/inputs/receptor.pdbqt','output/prepared/key/ligand.pdbqt','docking_results/pose.pdbqt'):assert n in z.namelist()
            z.extract('output/experiments.db',tmp_path/'unpacked')
        copy=DockingDatabase(tmp_path/'unpacked/output/experiments.db')
        try:assert copy.get_experiment_list()[0]['name']=='snapshot'
        finally:copy.close()
    finally:db.close()


@pytest.mark.parametrize('field',['seed','preparation_seed'])
def test_seed(field):
    with pytest.raises(ValueError):EvolutionConfig.from_dict({'targets':['A'],field:2**31})


@pytest.mark.parametrize('t',[0,-1,float('nan'),float('inf')])
def test_temperature(t):
    with pytest.raises(ValueError):AffinityCalculator(t)
    with pytest.raises(ValueError):AffinityCalculator().calculate_ki(-7,t)


def test_extreme_scores():
    c=AffinityCalculator()
    for s in (-500,500):
        assert c.calculate_ki(s) is None
        p=c.calculate_pki_from_score(s);assert math.isfinite(p) and (p>0)==(s<0)
    with pytest.raises(ValueError):c.calculate_pki(0)
    data=c.process_docking_results({'a':{'A':{'affinity':-500}}});json.dumps(data,allow_nan=False)
    assert data['a']['A']['pki']>0


def test_legacy_encoding(tmp_path):
    env={**os.environ,'PYTHONUTF8':'0','PYTHONIOENCODING':'cp1251:strict'}
    p=subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'main.py'),'--project',str(tmp_path)],input=b'17\n0\n',capture_output=True,env=env,timeout=30)
    assert p.returncode==0,(p.stdout+p.stderr).decode('cp1251',errors='replace')
    assert b'Invalid choice' in p.stdout


def test_aborted_candidate_has_no_partial_checkpoint(tmp_path):
    from src.evolution.engine import EvolutionaryGenerator
    from src.evolution.fitness import FitnessEvaluator
    from src.evolution.fragments import FragmentBank
    from src.evolution.models import Candidate
    from src.evolution.reweight import reweight
    class Docker:
        stats={'cached':0,'executed':1,'total_dockings':2,'failed':1}
        def run_docking(self,**kw):return {k:{'A':{'affinity':-7}} for k in kw['ligand_paths']}
    g=EvolutionaryGenerator(tmp_path);g.config=EvolutionConfig.from_dict({'targets':['A','B'],'min_nitrogen':0})
    g.evaluator=FitnessEvaluator(g.config);g.bank=FragmentBank(g.config);g.archive={}
    c=Candidate('CCO',0);g._evaluate([c],Docker(),AffinityCalculator(),0)
    assert c.energies=={} and c.score['fitness'] is None and not c.score['complete'] and not c.score['feasible']
    assert json.loads(next(g.output_dir.glob('generation_*.json')).read_text())['candidates'][0]['missing_targets']==['B']
    cp=tmp_path/'checkpoint.json';cp.write_text(json.dumps({'archive':{'CCO':c.to_dict()}}))
    assert reweight(cp,EvolutionConfig.from_dict({'targets':['A']}),tmp_path/'reweight.json')['incomplete']==1
