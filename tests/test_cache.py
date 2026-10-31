import multiprocessing as mp
import os
from pathlib import Path
import sqlite3
import time
from concurrent.futures import ThreadPoolExecutor
import pytest
from src.core.calculation_cache import CalculationCache, molecule_identity, atomic_write

def test_persistent_cache_and_order_independent_keys(tmp_path):
    cache = CalculationCache(tmp_path)
    entry = cache.get_or_compute('x', {'a': 1, 'b': 2}, lambda: ({'score': -7}, b'pose'))
    assert not entry.hit
    second = CalculationCache(tmp_path).get_or_compute('x', {'b': 2, 'a': 1},
                                                      lambda: pytest.fail('recomputed'))
    assert second.hit and second.key == entry.key
    assert cache.get('x', {'a': 2, 'b': 2}) is None
    assert cache.get('y', {'a': 1, 'b': 2}) is None

@pytest.mark.parametrize('column,value', [('artifact', b'damaged'), ('payload', '{'), ('spec','{}')])
def test_corruption_is_miss_and_repaired(tmp_path, column, value):
    cache = CalculationCache(tmp_path)
    entry = cache.get_or_compute('x', {'a': 1}, lambda: ({'n': 1}, b'old'))
    with sqlite3.connect(cache.database) as db:
        db.execute(f'UPDATE entries SET {column}=? WHERE key=?', (value, entry.key))
    assert cache.get('x', {'a': 1}) is None
    result = cache.get_or_compute('x', {'a': 1}, lambda: ({'n': 2}, b'new'))
    assert not result.hit and result.artifact == b'new'

def test_failure_retried_and_force_does_not_return_stale(tmp_path):
    cache = CalculationCache(tmp_path)
    def fail():
        raise RuntimeError('failure')
    with pytest.raises(RuntimeError):
        cache.get_or_compute('x', {}, fail)
    assert cache.get('x', {}) is None
    cache.get_or_compute('x', {}, lambda: ({}, b'valid'))
    with pytest.raises(RuntimeError):
        cache.get_or_compute('x', {}, fail, force=True)
    # A forced failure does not destroy an earlier independently valid result.
    assert cache.get('x', {}).artifact == b'valid'

def test_same_key_threads_compute_once(tmp_path):
    cache = CalculationCache(tmp_path)
    calls = []
    def compute():
        calls.append(1)
        time.sleep(.05)
        return {}, b'pose'
    with ThreadPoolExecutor(8) as pool:
        results = list(pool.map(lambda _: cache.get_or_compute('x', {}, compute), range(16)))
    assert len(calls) == 1
    assert sum(not r.hit for r in results) == 1

def worker(root, calls):
    cache = CalculationCache(Path(root))
    def compute():
        with open(calls, 'ab') as stream:
            stream.write(b'compute\n')
        time.sleep(.05)
        return {}, b'pose'
    cache.get_or_compute('x', {}, compute)

def crash_worker(root):
    cache = CalculationCache(Path(root))
    cache.get_or_compute('crash', {}, lambda: os._exit(9))

def test_process_lock_and_crash_recovery(tmp_path):
    ctx = mp.get_context('spawn')
    cache = CalculationCache(tmp_path / 'cache')
    calls = tmp_path / 'calls'
    processes = [ctx.Process(target=worker, args=(cache.root, calls)) for _ in range(4)]
    for p in processes: p.start()
    for p in processes:
        p.join(20)
        assert p.exitcode == 0
    assert calls.read_text() == 'compute\n'
    p = ctx.Process(target=crash_worker, args=(cache.root,))
    p.start(); p.join(20)
    assert p.exitcode == 9
    assert not cache.get_or_compute('crash', {}, lambda: ({}, b'retry')).hit
    with sqlite3.connect(cache.database) as db:
        assert db.execute("SELECT count(*) FROM attempts WHERE status='interrupted'").fetchone()[0] == 1

def test_distinct_keys_compute_in_parallel(tmp_path):
    import threading
    cache = CalculationCache(tmp_path)
    barrier = threading.Barrier(2)
    def compute():
        barrier.wait(timeout=3)
        return {}, b'pose'
    with ThreadPoolExecutor(2) as pool:
        results = list(pool.map(lambda i: cache.get_or_compute('x', {'i': i}, compute), range(2)))
    assert len(results) == 2

def test_nonfinite_and_empty_not_cached(tmp_path):
    cache = CalculationCache(tmp_path)
    for data, artifact in [({'score': float('nan')}, b'pose'), ({}, b'')]:
        with pytest.raises(ValueError):
            cache.get_or_compute('x', {}, lambda: (data, artifact))
        assert cache.get('x', {}) is None

def test_chemical_identity():
    assert molecule_identity('CCO') == molecule_identity('OCC')
    assert molecule_identity('[CH3:1][CH2:2][OH:3]') == molecule_identity('CCO')
    assert molecule_identity('C[C@H](O)F') != molecule_identity('C[C@@H](O)F')
    assert molecule_identity('[13CH3]CO') != molecule_identity('CCO')
    assert molecule_identity('CC[NH3+]') != molecule_identity('CCN')
    with pytest.raises(ValueError): molecule_identity('invalid')

def test_future_schema_rejected(tmp_path):
    cache = CalculationCache(tmp_path)
    with sqlite3.connect(cache.database) as db: db.execute('PRAGMA user_version=99')
    with pytest.raises(RuntimeError): CalculationCache(tmp_path)

def test_atomic_write_preserves_old_on_failed_replace(tmp_path, monkeypatch):
    path = tmp_path / 'file'; path.write_bytes(b'old')
    def fail(*args): raise OSError('disk failure')
    monkeypatch.setattr(os, 'replace', fail)
    with pytest.raises(OSError): atomic_write(path, b'new')
    assert path.read_bytes() == b'old'
    assert list(tmp_path.iterdir()) == [path]
