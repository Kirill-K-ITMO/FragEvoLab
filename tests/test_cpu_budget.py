"""CPU budget allocation and settings."""
import json
import time
from concurrent.futures import ThreadPoolExecutor
import pytest
from src.core.cpu_budget import CpuBudget,configure_project_cpu,project_budget

@pytest.mark.parametrize('value',[0,-1,True,1.5,'2'])
def test_invalid_budget(value):
    with pytest.raises(ValueError): CpuBudget(value)

@pytest.mark.parametrize('limit,threads,workers',[(1,1,6),(4,2,8),(3,2,6),(4,1,8)])
def test_weighted_budget(limit,threads,workers,monkeypatch):
    monkeypatch.setattr('src.core.cpu_budget.available_cpus',lambda:8)
    b=CpuBudget(limit)
    def work(_):
        with b.reserve(threads):time.sleep(.01)
    with ThreadPoolExecutor(workers) as ex:list(ex.map(work,range(12)))
    assert b.peak<=limit and b.peak>0 and b.active==b.waiting==0

def test_release_on_error_and_busy_reconfigure():
    b=CpuBudget(1)
    with pytest.raises(RuntimeError):
        with b.reserve(1):b.configure(1)
    assert b.active==0
    with b.reserve(1):pass

def test_persistent_settings_and_shared_instances(tmp_path,monkeypatch):
    monkeypatch.setattr('src.core.cpu_budget.available_cpus',lambda:8)
    a=configure_project_cpu(tmp_path,4,2)
    assert project_budget(tmp_path) is a
    assert json.loads((tmp_path/'cpu_settings.json').read_text())=={'cpu_budget':4,'cpu_per_job':2}
    configure_project_cpu(tmp_path,2,1)
    assert a.total==2 and a.per_job==1
    with pytest.raises(ValueError):project_budget(tmp_path,total=4)

def test_hardware_cap(monkeypatch):
    monkeypatch.setattr('src.core.cpu_budget.available_cpus',lambda:2)
    b=CpuBudget(99,99)
    assert b.total==b.per_job==2

def test_default_reads_saved_config(tmp_path):
    (tmp_path/'cpu_settings.json').write_text('{"cpu_budget":1,"cpu_per_job":1}')
    assert project_budget(tmp_path).total==1

def test_invalid_saved_settings(tmp_path):
    (tmp_path/'cpu_settings.json').write_text('{"cpu_budget":0}')
    with pytest.raises(ValueError):project_budget(tmp_path)


