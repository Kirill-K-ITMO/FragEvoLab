from pathlib import Path
from src.core.cache_manager import CacheManager


def test_real_preparation_reuses_equivalent_smiles(project):
    root, converter, ligand = project
    result = converter.prepare_cached('OCC')
    assert result['cache_hit'] and Path(result['path']) == ligand
    aspirin = converter.prepare_cached('CC(=O)Oc1ccccc1C(=O)O')
    assert Path(aspirin['path']).read_text().count('ATOM') > 10
    converter.preparation_seed = 123
    assert converter.prepare_cached('CCO')['preparation_id'] != result['preparation_id']

def test_alias_collision_and_failed_preparation(project, monkeypatch):
    root, converter, _ = project
    assert converter.convert_single('gen1_mol1', 'CCO')
    old = (converter.output_dir/'gen1_mol1.pdbqt').read_bytes()
    assert converter.convert_single('gen1_mol1', 'CCN')
    assert (converter.output_dir/'gen1_mol1.pdbqt').read_bytes() != old
    monkeypatch.setattr(converter, 'smiles_to_pdbqt_meeko', lambda *a: False)
    assert not converter.convert_single('gen1_mol1', 'CCC', force=True)
    assert not (converter.output_dir/'gen1_mol1.pdbqt').exists()

def test_prepared_file_tampering_is_repaired(project):
    _, converter, path = project
    original = path.read_bytes(); path.write_bytes(b'broken')
    result = converter.prepare_cached('CCO')
    assert result['cache_hit'] and path.read_bytes() == original


def test_selection_requires_all_targets_and_accepts_zero(project):
    from src.evolution.config import EvolutionConfig
    from src.evolution.fitness import FitnessEvaluator
    from src.evolution.models import Candidate
    from src.evolution.selection import select
    evaluator=FitnessEvaluator(EvolutionConfig.from_dict({'targets':['A','B']}))
    a=Candidate('CCO',0,energies={'A':-99})
    b=Candidate('CCN',0,energies={'A':0,'B':0})
    for candidate in (a,b): candidate.score=evaluator.evaluate(candidate.energies).to_dict()
    assert [c.smiles for c in select([a,b],2)]==['CCN']


def test_restore_empty_or_corrupt_preserves_output(project,tmp_path):
    manager=CacheManager(project[0])
    marker=project[0]/'output'/'important'; marker.write_text('keep')
    empty=tmp_path/'empty'; empty.mkdir()
    assert not manager.restore_from_backup(empty,confirm=False)
    assert marker.read_text()=='keep'
    backup=manager.backup_cache()
    (backup/'output'/'important').write_text('corrupted')
    assert not manager.restore_from_backup(backup,confirm=False)
    assert marker.read_text()=='keep'

def test_restore_valid_keeps_previous_output(project):
    manager=CacheManager(project[0]); marker=project[0]/'output'/'important'
    marker.write_text('original'); backup=manager.backup_cache()
    marker.write_text('new')
    assert manager.restore_from_backup(backup,confirm=False)
    assert marker.read_text()=='original'
    assert next(project[0].glob('output_before_restore_*/important')).read_text()=='new'


def test_restore_staging_failure_leaves_current_output(project,monkeypatch):
    import shutil
    root=project[0]; manager=CacheManager(root)
    marker=root/'output'/'important'; marker.write_text('old')
    backup=manager.backup_cache(); marker.write_text('current')
    def fail(*a,**k): raise OSError('copy failed')
    monkeypatch.setattr(shutil,'copytree',fail)
    assert not manager.restore_from_backup(backup,confirm=False)
    assert marker.read_text()=='current'

def test_restore_promotion_failure_rolls_back(project,monkeypatch):
    root=project[0]; manager=CacheManager(root)
    marker=root/'output'/'important'; marker.write_text('old')
    backup=manager.backup_cache(); marker.write_text('current')
    original=Path.rename
    def fail_promotion(self,target):
        if self.name=='output' and self.parent.name.startswith('restore-'):
            raise OSError('promotion failed')
        return original(self,target)
    monkeypatch.setattr(Path,'rename',fail_promotion)
    assert not manager.restore_from_backup(backup,confirm=False)
    assert marker.read_text()=='current'
