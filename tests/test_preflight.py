import subprocess
import sys
from pathlib import Path
from src.core.preflight import check_run, check_project, check_dependencies
from src.core.calculation_cache import atomic_json
from src.evolution.config import EvolutionConfig


def configuration():
    return EvolutionConfig.from_dict({'targets':['A'], 'off_targets':[{'name':'B','threshold':-5}],
        'generations':0,'population_size':2,'min_nitrogen':0})


def inputs(root):
    path=root/'data/ligands.csv'
    path.write_text('name,smiles\na,CCO\nb,CCN\n')
    return path


def test_missing_dependency_explained(monkeypatch):
    from src.core import preflight
    original=preflight.importlib.import_module
    def missing(name):
        if name=='meeko': raise ImportError('test dependency missing')
        return original(name)
    monkeypatch.setattr(preflight.importlib,'import_module',missing)
    report=check_dependencies()
    assert not report.ready
    assert any('meeko' in e and 'requirements.txt' in e for e in report.errors)


def test_app_can_open_without_vina(project):
    from src.ui.application import FragEvoLabApp
    root=project[0]; inputs(root)
    atomic_json(root/'app_settings.json',{'vina_executable':str(root/'missing-vina')})
    atomic_json(root/'evolution_settings.json',configuration().to_dict())
    app=FragEvoLabApp(root)
    try:
        assert app._docker is None
        report=check_project(root)
        assert not report.ready
        assert any('Vina/CPU' in e for e in report.errors)
    finally: app.close()


def test_invalid_json_is_actionable(project):
    root=project[0]; (root/'evolution_settings.json').write_text('{broken')
    report=check_project(root)
    assert not report.ready
    assert any('Evolution settings' in e for e in report.errors)


def test_cli_check_works_without_site_packages(tmp_path):
    root=Path(__file__).resolve().parents[1]
    result=subprocess.run([sys.executable,'-S',str(root/'main.py'),'--check','--project',str(tmp_path)],
        capture_output=True,text=True,encoding='utf-8')
    assert result.returncode==1
    assert 'requirements.txt' in result.stdout
    assert 'Traceback' not in result.stderr


