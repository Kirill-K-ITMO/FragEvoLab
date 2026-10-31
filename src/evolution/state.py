"""JSON checkpoints only, never pickle. Resume after a committed generation."""
import hashlib
import json
import platform
import sys
from importlib.metadata import version
from pathlib import Path
from rdkit import rdBase
from src.core.calculation_cache import atomic_json

SCHEMA = 1

def fingerprint(config, docker):
    values = config.to_dict()
    values.pop("generations")  # extending a run is allowed
    inputs = {}
    for protein in config.proteins:
        for folder,ext in (("proteins_dir","pdbqt"),("configs_dir","txt")):
            path = (Path(docker.proteins_dir)/f"{protein}.pdbqt" if folder == "proteins_dir"
                    else docker.get_config(protein))
            if path is None: raise ValueError(f"Missing configuration for {protein}")
            inputs[f"{folder}/{protein}.{ext}"] = hashlib.sha256(path.read_bytes()).hexdigest()
    engine = getattr(docker,"vina_executable","")
    import shutil
    engine_path = shutil.which(engine) or engine
    if Path(engine_path).is_file(): inputs["engine"] = hashlib.sha256(Path(engine_path).read_bytes()).hexdigest()
    source_root = Path(__file__).parent.parent
    source_files = sorted(Path(__file__).parent.glob("*.py")) + [source_root/n for n in (
        "preparation/ligands.py", "preparation/validation.py", "docking/vina.py", "core/calculation_cache.py",
        "docking/affinity.py", "core/cpu_budget.py", "preparation/inputs.py")]
    values["cpu_budget"] = docker.cpu_pool.total
    values.update(inputs=inputs, rdkit=rdBase.rdkitVersion,
                  runtime={"python":sys.version, "platform":platform.platform(), **{p:version(p) for p in ("meeko","numpy","scipy","gemmi")}},
                  engine_options={k:getattr(docker,k,None) for k in ("seed","cpu","timeout")},
                  implementation={str(p.relative_to(source_root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files})
    return hashlib.sha256(json.dumps(values,sort_keys=True).encode()).hexdigest()

def tuples(value):
    return tuple(tuples(v) for v in value) if isinstance(value,list) else value

def save(path, payload): atomic_json(Path(path), {"schema":SCHEMA, **payload})

def load(path, expected):
    data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if data.get("schema") != SCHEMA or data.get("fingerprint") != expected:
        raise ValueError("Checkpoint settings, implementation or docking inputs changed; start a new run")
    return data
