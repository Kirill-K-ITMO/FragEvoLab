"""Verification CLI arguments and ligand preparation."""
import json
import subprocess
import sys
from pathlib import Path
import pytest
from scripts import verify_integration, verify_ligands
from src.evolution.config import EvolutionConfig
from src.version import APP_NAME, __version__


def test_launcher_help_and_version():
    root = Path(__file__).resolve().parents[1]
    for script in ('verify_ligands.py', 'verify_integration.py'):
        result = subprocess.run([sys.executable, str(root / script), '--help'], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
    result = subprocess.run([sys.executable, str(root / 'main.py'), '--version'], capture_output=True, text=True)
    assert result.returncode == 0
    assert result.stdout.strip() == f'{APP_NAME} {__version__}'
    assert (root / 'README.md').read_text().splitlines()[0] == f'# {APP_NAME} {__version__}'


def test_real_check_requires_explicit_inputs():
    with pytest.raises(SystemExit) as error:
        verify_integration.main([])
    assert error.value.code == 2


def test_missing_receptor_blocks_before_run(tmp_path):
    (tmp_path / 'ligands.csv').write_text('name,smiles\none,CCN\n')
    config = EvolutionConfig.from_dict({'targets': ['JAK1'], 'generations': 1})
    with pytest.raises(ValueError, match='Missing receptor'):
        verify_integration.validate_inputs(tmp_path, config)


def test_preparation_reports_failure_as_nonzero(tmp_path):
    source = tmp_path / 'input.csv'
    source.write_text('name,smiles\ngood,CCN\nbad,not-a-smiles\n')
    output = tmp_path / 'checks/preparation.json'
    assert verify_ligands.main(['--csv', str(source), '--output', str(output)]) == 1
    assert [row['success'] for row in json.loads(output.read_text())] == [True, False]


def test_preparation_refuses_to_overwrite_input(tmp_path):
    source = tmp_path / 'input.csv'
    original = 'name,smiles\none,CCN\n'
    source.write_text(original)
    with pytest.raises(SystemExit) as error:
        verify_ligands.main(['--csv', str(source), '--output', str(source)])
    assert error.value.code == 2
    assert source.read_text() == original
