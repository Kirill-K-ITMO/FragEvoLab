"""The bundled inputs are usable with the default evolution configuration."""
import csv
import hashlib
import json
from pathlib import Path
from src.evolution.config import EvolutionConfig

ROOT = Path(__file__).resolve().parents[1]


def test_default_receptors_have_matching_inputs_and_provenance():
    config = EvolutionConfig.from_dict(json.loads((ROOT / 'evolution_settings.json').read_text()))
    sources = json.loads((ROOT / 'data/receptor_sources.json').read_text())
    names = {p.name for p in config.targets + config.off_targets}
    assert names == {record['name'] for record in sources}
    for record in sources:
        for item in record['files']:
            path = ROOT / 'data' / item['path']
            assert path.is_file()
            assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256']
        assert (ROOT / 'data/proteins' / (record['name'] + '.pdbqt')).is_file()
        assert (ROOT / 'data/configs' / (record['name'] + '.txt')).is_file()


def test_default_ligands_match_source_table():
    def rows(name):
        with (ROOT / 'data' / name).open(encoding='utf-8-sig', newline='') as handle:
            return list(csv.DictReader(handle))
    ligands = rows('ligands.csv')
    sources = rows('ligand_sources.csv')
    assert ligands
    assert len({row['name'] for row in ligands}) == len(ligands)
    assert {(row['name'], row['smiles']) for row in ligands} == {
        (row['name'], row['smiles']) for row in sources
    }
    assert all(row['pubchem_cid'] and row['pubchem_url'] for row in sources)
