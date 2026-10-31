"""Local preparation fixtures; real docking uses the separate verification CLI."""
from pathlib import Path
import pytest
from src.preparation.ligands import LigandConverter

@pytest.fixture
def project(tmp_path):
    for name in ('data/proteins', 'data/configs'):
        (tmp_path / name).mkdir(parents=True)
    converter = LigandConverter(tmp_path)
    prepared = converter.prepare_cached('CCO')
    ligand = Path(prepared['path'])
    for target in ('A', 'B'):
        (tmp_path / 'data/proteins' / (target + '.pdbqt')).write_bytes(ligand.read_bytes())
        (tmp_path / 'data/configs' / (target + '.txt')).write_text(
            'center_x=0\ncenter_y=0\ncenter_z=0\nsize_x=12\nsize_y=12\nsize_z=12\n')
    return tmp_path, converter, ligand


def pytest_collection_modifyitems(items):
    """The default suite contains local tests only; real docking has its own CLI."""
    for item in items:
        item.add_marker(pytest.mark.fast)
