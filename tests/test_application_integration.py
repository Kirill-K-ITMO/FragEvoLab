import pytest
from src.ui.evolution_menu import edit_config, evolution_menu
from src.evolution.config import EvolutionConfig
from src.core.app_settings import vina_path
from src.core.calculation_cache import atomic_json

class Terminal:
    def __init__(self, answers): self.answers=iter(answers); self.messages=[]
    def input(self,*args): return next(self.answers)
    def print(self,*args,**kwargs): self.messages.append(str(args))

def test_editor_keeps_zero_probability_and_multiple_proteins():
    c=EvolutionConfig.from_dict({'targets':['A','B'],'off_targets':[{'name':'C','threshold':-5},{'name':'D','threshold':-6}]})
    keys=list(c.to_dict())
    term=Terminal([str(keys.index('level_probabilities')+1),'{"atom":0,"fragment":1}', '0'])
    updated=edit_config(term,c)
    assert updated.proteins==['A','B','C','D']
    assert updated.level_probabilities=={'atom':0,'fragment':1}

def test_editor_rejects_bad_probabilities_then_recovers():
    c=EvolutionConfig.from_dict({'targets':['A']})
    index=str(list(c.to_dict()).index('mutation_fraction')+1)
    term=Terminal([index,'2','0',index,'0.5','0'])
    assert edit_config(term,c).mutation_fraction==.5
    assert any('Parameter error' in x for x in term.messages)


def test_relative_executable_and_invalid_settings(tmp_path):
    atomic_json(tmp_path/'app_settings.json',{'vina_executable':'tools/vina'})
    assert vina_path(tmp_path)==str(tmp_path/'tools/vina')
    atomic_json(tmp_path/'app_settings.json',{'vina_executable':[]})
    with pytest.raises(ValueError): vina_path(tmp_path)


def test_zero_seed_rejected_before_docking():
    with pytest.raises(ValueError, match="seed must be positive"):
        EvolutionConfig.from_dict({'targets':['A'],'seed':0})

