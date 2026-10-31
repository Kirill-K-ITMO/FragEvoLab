"""Shared weighted CPU budget for Vina jobs in one application process/project."""
from contextlib import contextmanager
import json
import os
from pathlib import Path
from threading import Condition, Lock
from src.core.calculation_cache import atomic_json

def positive(value, name):
    if isinstance(value,bool) or not isinstance(value,int) or value < 1:
        raise ValueError(f"{name} must be a positive integer (0/auto is not allowed)")
    return value

def available_cpus():
    try:
        return max(1,len(os.sched_getaffinity(0)))
    except (AttributeError,OSError):
        return max(1,os.cpu_count() or 1)

class CpuBudget:
    def __init__(self, total, per_job=1):
        self.condition = Condition()
        self.active = self.waiting = self.peak = 0
        self.total = min(positive(total,'cpu_budget'),available_cpus())
        self.per_job = min(positive(per_job,'cpu_per_job'),self.total)

    def configure(self, total, per_job=1):
        total=min(positive(total,'cpu_budget'),available_cpus())
        per_job=min(positive(per_job,'cpu_per_job'),total)
        with self.condition:
            if self.active or self.waiting:
                raise RuntimeError('Cannot change CPU settings while docking jobs are running or waiting')
            self.total,self.per_job = total,per_job
            self.peak=0
        return self

    @contextmanager
    def reserve(self, count):
        positive(count,'task CPU')
        with self.condition:
            if count>self.total:
                raise ValueError('CPU budget changed; restart this docking request')
            self.waiting+=1
            try:
                self.condition.wait_for(lambda:self.active+count<=self.total)
                self.active+=count
                self.peak=max(self.peak,self.active)
            finally:
                self.waiting-=1
        try:
            yield
        finally:
            with self.condition:
                self.active-=count
                self.condition.notify_all()

    def snapshot(self):
        with self.condition:
            return dict(cpu_budget=self.total,cpu_per_job=self.per_job,
                        active_cpu=self.active,waiting_jobs=self.waiting,peak_cpu=self.peak)

_registry = {}
_registry_lock = Lock()

def project_budget(root, *, total=None, per_job=None):
    root=Path(root).resolve()
    with _registry_lock:
        if root not in _registry:
            path=root/'cpu_settings.json'
            values=json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
            unknown=set(values)-{'cpu_budget','cpu_per_job'}
            if unknown: raise ValueError(f'Unknown CPU settings: {sorted(unknown)}')
            _registry[root]=CpuBudget(total if total is not None else values.get('cpu_budget',min(4,available_cpus())),
                                      per_job if per_job is not None else values.get('cpu_per_job',1))
        else:
            pool=_registry[root]
            if total is not None or per_job is not None:
                target_total=pool.total if total is None else min(positive(total,'cpu_budget'),available_cpus())
                target_per=pool.per_job if per_job is None else min(positive(per_job,'cpu_per_job'),target_total)
                if (target_total,target_per)!=(pool.total,pool.per_job):
                    # Explicit reconfiguration belongs to the settings action, not
                    # construction of another utility's private docker instance.
                    raise ValueError('Conflicting CPU settings in one project; use configure_project_cpu')
        return _registry[root]

def configure_project_cpu(root, total, per_job=1, *, save=True):
    root=Path(root).resolve()
    with _registry_lock:
        pool=_registry.get(root)
        if pool is None:
            pool=CpuBudget(total,per_job)
        else:
            pool.configure(total,per_job)
        if save:
            atomic_json(root/'cpu_settings.json',{'cpu_budget':pool.total,'cpu_per_job':pool.per_job})
        _registry[root]=pool
        return pool

def child_environment(cpu):
    env=os.environ.copy()
    # Prevent inherited OMP/BLAS settings from bypassing the task's requested cap.
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[key]=str(cpu)
    env['OMP_DYNAMIC']='FALSE'
    return env
