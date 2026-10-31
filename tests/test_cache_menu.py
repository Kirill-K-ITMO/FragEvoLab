"""Cache menu selections must preserve scientific output files."""
from types import SimpleNamespace
import pytest
from src.core.cache_manager import CacheManager
from src.ui.application import FragEvoLabApp


@pytest.mark.parametrize('choice,remaining', [
    ('1', set()), ('2', {'docking'}), ('3', {'preparation'}),
])
def test_selected_cache_kind_preserves_outputs(tmp_path, monkeypatch, choice, remaining):
    manager = CacheManager(tmp_path)
    for kind in ('preparation', 'docking'):
        manager.calculation_cache.get_or_compute(kind, {}, lambda: ({}, b'artifact'))
    output = tmp_path / 'output' / 'checkpoint.json'
    output.parent.mkdir()
    output.write_text('{"saved": true}')
    from src.ui import application
    monkeypatch.setattr(application.console, 'input', lambda *args: choice)
    monkeypatch.setattr('builtins.input', lambda *args: 'yes')
    FragEvoLabApp.clear_cache(SimpleNamespace(cache_mgr=manager))
    assert set(manager.calculation_cache.statistics()) == remaining
    assert output.read_text() == '{"saved": true}'
