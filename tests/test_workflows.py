import random
import pytest
from rdkit import Chem
from src.ui.docking_menu import parse_selection, docking_menu
from src.evolution.config import EvolutionConfig
from src.evolution.fragments import FragmentBank
from src.evolution.adaptation import Adaptation
from src.evolution.operators import Operators
from src.reporting.evolution import write_report

@pytest.mark.parametrize('text,expected', [('all',[0,1,2]),('3, 1 3',[2,0]),('2',[1])])
def test_selection(text, expected):
    assert parse_selection(text,3)==expected

@pytest.mark.parametrize('text', ['', '0', '-1', '4', 'all 1', '1 x'])
def test_bad_selection(text):
    with pytest.raises(ValueError): parse_selection(text,3)

class Terminal:
    def __init__(self, answers): self.answers=iter(answers)
    def input(self,*args): return next(self.answers)
    def print(self,*args,**kwargs): pass


def test_integrated_atom_edits_preflight(capfd):
    cfg=EvolutionConfig.from_dict({'targets':['A'],'level_probabilities':{'atom':1,'fragment':0}})
    operators=Operators(random.Random(42),cfg,FragmentBank(cfg),Adaptation(cfg))
    for smiles in ['CC(C)(C)N','CN(C)C','C[NH2+]C','C[C@H](O)CCN','NCCO']:
        for _ in range(50):
            product,events=operators.mutate(smiles,1)
            if product:
                mol=Chem.MolFromSmiles(product)
                assert mol is not None
                Chem.SanitizeMol(mol)
    assert 'Explicit valence' not in capfd.readouterr().err

def test_empty_report_and_escaped_data(tmp_path):
    report={'run_id':'<script>bad()</script>','stop_reason':'empty','unique_molecules':0,
        'total_dockings':0,'cache_hits':0,'config':{'targets':['A'],'off_targets':[]}}
    path=tmp_path/'report.html'; write_report(report,path)
    text=path.read_text()
    assert '<script>bad()' not in text
    assert '&lt;script&gt;bad()' in text
    assert '3.1.4' in text
    assert 'No data' in text

def test_ascii_banner_at_narrow_width(monkeypatch):
    import io
    from rich.console import Console
    from src.ui import console as progress_utils
    out=io.StringIO()
    monkeypatch.setattr(progress_utils,'console',Console(file=out,width=40,color_system=None))
    progress_utils.print_app_header()
    lines=out.getvalue().splitlines()
    assert 'FragEvoLab 3.1.4' in out.getvalue()
    assert all(len(line)<=40 for line in lines)

