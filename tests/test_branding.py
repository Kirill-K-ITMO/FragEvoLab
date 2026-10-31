import io
from pathlib import Path
from rich.console import Console
from src.reporting.branding import brand_html
from src.ui import console as progress_utils

def test_self_contained_brand_and_idempotence():
    html=brand_html('<html><head><title>Results</title></head><body>data</body></html>')
    assert '<title>FragEvoLab | Results</title>' in html
    assert 'data:image/png;base64,' in html
    assert 'FragEvoLab 3.1.4' in html
    assert brand_html(html)==html

def test_terminal_has_brand_without_rulers(monkeypatch):
    output=io.StringIO()
    monkeypatch.setattr(progress_utils,'console',Console(file=output,width=40,color_system=None))
    progress_utils.print_app_header()
    progress_utils.print_header('Docking','Preparation')
    progress_utils.print_section_divider('Results')
    text=output.getvalue()
    assert 'FragEvoLab 3.1.4' in text
    assert not any(token in text for token in ('---','===','══','──','╔','╚'))
    assert all(len(line)<=40 for line in text.splitlines())

def test_docking_and_database_reports_branded(project):
    from src.reporting.docking import ReportGenerator
    reporter=ReportGenerator(project[0])
    data={'one':{'A':{'affinity':-6.0,'smiles':'CCO','ki':1e-5,'pki':5.0,'ki_formatted':'10 uM'}}}
    html,_=reporter.generate_reports(data,title='Brand check')
    assert 'id="fragevolab-brand"' in Path(html).read_text()
    db=reporter.generate_db_report({'results':[{'ligand_name':'one','protein_name':'A','ligand_smiles':'CCO','delta_g':-6.,'ki':1e-5,'pki':5.,'experiment_name':'test','timestamp':'2026-10-15'}], 'total_experiments':1,'total_ligands':1,'total_proteins':1})
    assert 'id="fragevolab-brand"' in Path(db).read_text()
