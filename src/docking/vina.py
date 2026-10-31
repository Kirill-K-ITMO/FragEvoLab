"""vina - FragEvoLab module."""

import math
import shutil
import tempfile
import platform
import uuid
from threading import Lock
from src.core.cpu_budget import project_budget, positive, child_environment
from src.core.calculation_cache import (CalculationCache, digest, atomic_write, atomic_json,
                                   validate_pdbqt)
import subprocess
import json
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
from datetime import datetime
import time

class AutoDockVina:
    """
    Perform molecular docking with AutoDock Vina.
    
    Supports:
    - Multi-target docking (multiple receptors)
    - Parallel processing
    - Configurable exhaustiveness
    - Saving results
    - Progress bars and timers
    """
    
    def __init__(self, project_root: Optional[Path] = None, *, vina_executable=None,
                 seed: int = 42, cpu: int | None = None, timeout: float = 3600,
                 cpu_budget: int | None = None):
        """
        Initialize the docking module.
        
        Args:
            project_root: Project root directory (default: current directory)
        """
        self.project_root = Path(project_root or Path.cwd()).resolve()
        from src.core.app_settings import vina_path
        self.vina_executable = str(vina_executable) if vina_executable is not None else vina_path(self.project_root)
        self.cpu_pool = project_budget(self.project_root, total=cpu_budget)
        self._cpu_override = positive(cpu,"cpu") if cpu is not None else None
        self.seed, self.timeout = seed, timeout
        self._stats_lock = Lock()
        self.calculation_cache = CalculationCache(self.project_root / ".calculation_cache" / "v1")
        self.last_batch_dir = None
        
        # Directory paths
        self.data_dir = self.project_root / "data"
        self.proteins_dir = self.data_dir / "proteins"
        self.configs_dir = self.data_dir / "configs"
        self.ligands_dir = self.project_root / "output" / "ligands_pdbqt"
        self.results_dir = self.project_root / "output" / "docking_results"
        self.cache_dir = self.project_root / "output" / "cache"
        
        # Create the required directories
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Docking statistics
        self.stats = {
            'total_dockings': 0,
            'successful': 0,
            'failed': 0,
            'cached': 0,
            'errors': [],
            'elapsed_time': 0
        }
        
        # Docking results
        self.docking_results = {}  # {ligand_name: {protein_name: {affinity, energy, ...}}}
        
        # Configure logging
        self._setup_logging()
        
        # Check the Vina installation
        self._check_vina()
    
    @property
    def cpu(self):
        return min(self._cpu_override or self.cpu_pool.per_job,self.cpu_pool.total)

    def _setup_logging(self):
        """Configure logging for docking."""
        logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
        self.logger = logging.getLogger(__name__)
    
    def _check_vina(self):
        executable = shutil.which(self.vina_executable)
        if executable is None:
            raise RuntimeError("AutoDock Vina executable not found")
        self.vina_executable = str(Path(executable).resolve())
        result = subprocess.run([self.vina_executable, "--version"],
                                capture_output=True, text=True, check=True, timeout=30)
        self.vina_version = (result.stdout or result.stderr).strip()
        if not self.vina_version:
            raise RuntimeError("Vina did not report a version")
        self.logger.info("Vina: %s", self.vina_version)

    def get_proteins(self) -> List[Path]:
        """
        Return all receptor PDBQT files.
        
        Returns:
            List of receptor file paths
        """
        proteins = sorted(self.proteins_dir.glob("*.pdbqt"))
        
        if not proteins:
            self.logger.warning(f"No receptors found in: {self.proteins_dir}")
        
        return proteins
    
    def get_ligands(self) -> List[Path]:
        """
        Return all ligand PDBQT files.
        
        Returns:
            List of ligand file paths
        """
        ligands = sorted(self.ligands_dir.glob("*.pdbqt"))
        
        if not ligands:
            self.logger.warning(f"No ligands found in: {self.ligands_dir}")
        
        return ligands
    
    def get_config(self, protein_name: str) -> Optional[Path]:
        """
        Return the receptor configuration file.
        
        Args:
            protein_name: Receptor name (without extension)
            
        Returns:
            Configuration file path or None
        """
        # Try supported extensions
        for ext in ['.txt', '.conf', '.config']:
            config_file = self.configs_dir / f"{protein_name}{ext}"
            if config_file.exists():
                return config_file
        
        self.logger.warning(f"Configuration not found for receptor: {protein_name}")
        return None
    
    def parse_config(self, config_path: Path) -> Dict[str, str]:
        params = {}
        for raw in Path(config_path).read_text(encoding="utf-8-sig").splitlines():
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            if "=" not in line:
                raise ValueError(f"Invalid Vina config line: {raw}")
            key, value = (part.strip() for part in line.split("=", 1))
            if not key or not value or key in params:
                raise ValueError(f"Empty or duplicate Vina option: {key}")
            params[key] = value
        return params

    def _effective_options(self, config_path, exhaustiveness):
        params = self.parse_config(config_path)
        # The supplied receptor and exhaustiveness override the config values.
        params.pop("receptor", None)
        params["exhaustiveness"] = str(exhaustiveness)
        params.setdefault("seed", str(self.seed))
        params.setdefault("cpu", str(self.cpu))
        params.setdefault("scoring", "vina")
        integers = {"seed", "cpu", "exhaustiveness", "num_modes", "max_evals", "verbosity"}
        floats = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
                  "energy_range", "min_rmsd", "spacing"}
        # Unknown options fail closed: e.g. --maps/--flex would add unhashed inputs.
        unknown = set(params) - integers - floats - {"scoring"}
        if unknown:
            raise ValueError(f"Unsupported config options: {sorted(unknown)}")
        for key in integers & params.keys():
            params[key] = str(int(params[key]))
        for key in floats & params.keys():
            number = float(params[key])
            if not math.isfinite(number):
                raise ValueError(f"Non-finite option: {key}")
            params[key] = format(number, ".17g")
        required = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z"}
        if not required <= params.keys():
            raise ValueError("Docking box requires center_x/y/z and size_x/y/z")
        if any(float(params[k]) <= 0 for k in ("size_x", "size_y", "size_z")):
            raise ValueError("Docking box sizes must be positive")
        if not -2147483648 <= int(params["seed"]) <= 2147483647:
            raise ValueError("Vina seed must fit a signed 32-bit integer")
        if int(params["seed"]) == 0:
            raise ValueError("Use a nonzero explicit seed; seed=0 requests an uncontrolled replicate")
        if int(params["cpu"]) < 1 or int(params["exhaustiveness"]) < 1:
            raise ValueError("cpu and exhaustiveness must be positive")
        if params["scoring"] not in ("vina", "vinardo"):
            raise ValueError("This cache integration supports rigid receptor Vina/Vinardo docking")
        requested_cpu = int(params["cpu"])
        effective_cpu = min(requested_cpu, self.cpu, self.cpu_pool.total)
        if effective_cpu != requested_cpu:
            self.logger.warning("Config cpu=%s limited to %s by application CPU settings",requested_cpu,effective_cpu)
        params["cpu"] = str(effective_cpu)
        return params

    def _count(self, field, error=None):
        with self._stats_lock:
            self.stats[field] = self.stats.get(field, 0) + 1
            if error is not None:
                self.stats['errors'].append(error)

    def run_single_docking(self, ligand_path: Path, protein_path: Path,
                           config_path: Path, exhaustiveness: int = 8,
                           force: bool = False) -> Optional[Dict]:
        ligand_path, protein_path = Path(ligand_path), Path(protein_path)
        ligand_name, protein_name = ligand_path.stem, protein_path.stem
        try:
            # Read once and execute on these snapshots, so the key describes the
            # exact input bytes even if someone edits an original during docking.
            ligand_bytes, protein_bytes = ligand_path.read_bytes(), protein_path.read_bytes()
            validate_pdbqt(ligand_bytes)
            validate_pdbqt(protein_bytes)
            options = self._effective_options(config_path, exhaustiveness)
            executable_sha = digest(Path(self.vina_executable).read_bytes())
            spec = {"ligand_sha": digest(ligand_bytes), "receptor_sha": digest(protein_bytes),
                    "options": options, "engine_sha": executable_sha,
                    "engine_version": self.vina_version, "platform": platform.platform(),
                    "implementation_sha": digest(Path(__file__).read_bytes())}
            def compute():
                with tempfile.TemporaryDirectory(prefix="vina-") as directory:
                    temp = Path(directory)
                    ligand, receptor, pose = (temp / n for n in ("ligand.pdbqt", "receptor.pdbqt", "pose.pdbqt"))
                    ligand.write_bytes(ligand_bytes)
                    receptor.write_bytes(protein_bytes)
                    command = [self.vina_executable, "--ligand", str(ligand),
                               "--receptor", str(receptor), "--out", str(pose)]
                    for key, value in sorted(options.items()):
                        command.extend(["--" + key, value])
                    if digest(Path(self.vina_executable).read_bytes()) != executable_sha:
                        raise RuntimeError("Vina executable changed while preparing the task")
                    with self.cpu_pool.reserve(int(options['cpu'])):
                        self._count('executed')
                        process = subprocess.run(command, capture_output=True, text=True,
                                                 check=True, timeout=self.timeout,
                                                 env=child_environment(int(options['cpu'])))
                    artifact = pose.read_bytes()
                    validate_pdbqt(artifact)
                    parsed = self._parse_docking_output(pose, ligand_name, protein_name)
                    if parsed is None or not math.isfinite(parsed['affinity']):
                        raise ValueError("No finite Vina score in the produced PDBQT")
                    return {"affinity": parsed['affinity'], "mode": parsed['mode'],
                            "timestamp": parsed['timestamp'], "spec": spec,
                            "stdout": process.stdout, "stderr": process.stderr}, artifact
            entry = self.calculation_cache.get_or_compute("docking", spec, compute, force=force)
            output_file = self.results_dir / "by_id" / entry.key / (digest(entry.artifact) + ".pdbqt")
            if not output_file.exists() or digest(output_file.read_bytes()) != digest(entry.artifact):
                atomic_write(output_file, entry.artifact)
            inputs = {}
            for role, content in (("ligand", ligand_bytes), ("receptor", protein_bytes)):
                saved_input = self.project_root / "output" / "inputs" / (digest(content) + ".pdbqt")
                if not saved_input.exists() or digest(saved_input.read_bytes()) != digest(content):
                    atomic_write(saved_input, content)
                inputs[role] = str(saved_input)
            self._count('cached' if entry.hit else 'successful')
            return {**entry.payload, 'input_files': inputs, 'ligand': ligand_name, 'protein': protein_name,
                    'evaluation_id': entry.key, 'artifact_sha': digest(entry.artifact),
                    'output_file': str(output_file), 'cache_hit': entry.hit}
        except Exception as error:
            self._count('failed', (ligand_name, protein_name, str(error)))
            self.logger.error("Docking %s / %s: %s", ligand_name, protein_name, error)
            return None

    def _parse_vina_output(self, stdout: str, stderr: str, ligand_name: str, protein_name: str) -> Optional[Dict]:
        """
        Parse AutoDock Vina output to extract results.
        
        Args:
            stdout: Vina standard output
            stderr: Vina standard error
            ligand_name: Ligand name
            protein_name: Receptor name
            
        Returns:
            Result dictionary or None
        """
        # Combine stdout and stderr for parsing
        full_output = stdout + "\n" + stderr
        lines = full_output.split('\n')
        
        # Find the results table
        for i, line in enumerate(lines):
            # Format 1: standard Vina table
            if 'mode |   affinity | dist from best mode' in line.lower():
                if i + 2 < len(lines):
                    result_line = lines[i + 2].strip()
                    parts = result_line.split()
                    
                    if len(parts) >= 2:
                        try:
                            mode = int(parts[0])
                            affinity = float(parts[1])
                            
                            return {
                                'ligand': ligand_name,
                                'protein': protein_name,
                                'mode': mode,
                                'affinity': affinity,
                                'timestamp': datetime.now().isoformat()
                            }
                        except (ValueError, IndexError) as e:
                            self.logger.debug(f"Failed to parse format 1: {e}")
            
            # Format 2: line containing "Affinity:"
            if 'affinity:' in line.lower():
                match = re.search(r'affinity:\s*([-+]?\d+\.?\d*)', line, re.IGNORECASE)
                if match:
                    try:
                        affinity = float(match.group(1))
                        return {
                            'ligand': ligand_name,
                            'protein': protein_name,
                            'mode': 1,
                            'affinity': affinity,
                            'timestamp': datetime.now().isoformat()
                        }
                    except ValueError as e:
                        self.logger.debug(f"Failed to parse format 2: {e}")
            
            # Format 3: line beginning with "1   "
            if re.match(r'^\s*1\s+[-+]?\d+\.?\d*', line):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mode = int(parts[0])
                        affinity = float(parts[1])
                        
                        return {
                            'ligand': ligand_name,
                            'protein': protein_name,
                            'mode': mode,
                            'affinity': affinity,
                            'timestamp': datetime.now().isoformat()
                        }
                    except (ValueError, IndexError) as e:
                        self.logger.debug(f"Failed to parse format 3: {e}")
        
        self.logger.warning(f"Failed to parse Vina output for {ligand_name} × {protein_name}")
        return None
    
    def _parse_docking_output(self, output_file: Path, ligand_name: str, protein_name: str) -> Optional[Dict]:
        """
        Parse the output PDBQT file to extract the docking score.
        
        Args:
            output_file: Output PDBQT file path
            ligand_name: Ligand name
            protein_name: Receptor name
            
        Returns:
            Result dictionary or None
        """
        try:
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('REMARK VINA RESULT:'):
                        parts = line.split()
                        if len(parts) >= 4:
                            try:
                                affinity = float(parts[3])
                                
                                return {
                                    'ligand': ligand_name,
                                    'protein': protein_name,
                                    'mode': 1,
                                    'affinity': affinity,
                                    'timestamp': datetime.now().isoformat()
                                }
                            except (ValueError, IndexError) as e:
                                self.logger.error(f"Failed to parse PDBQT: {e}")
        except Exception as e:
            self.logger.error(f"Failed to read {output_file}: {e}")
        
        return None
    
    def run_docking(self, ligands: Optional[List[str]] = None,
                    proteins: Optional[List[str]] = None, exhaustiveness: int = 8,
                    num_workers: int = 4, force: bool = False, *,
                    ligand_paths: Optional[Dict[str, Path]] = None) -> Dict:
        """Batch API preserved; explicit paths avoid mutable generation-name aliases.

        None means all available inputs; an empty list means no inputs.
        One instance runs one batch at a time. Task-level work is parallel.
        """
        positive(num_workers,"num_workers")
        effective_workers=min(num_workers,self.cpu_pool.total)
        self.logger.info("CPU budget=%s; requested jobs=%s; worker cap=%s; max CPU/job=%s",
                         self.cpu_pool.total,num_workers,effective_workers,self.cpu)
        started = time.monotonic()
        self.stats = {'total_dockings': 0, 'successful': 0, 'failed': 0,
                      'cached': 0, 'executed': 0, 'errors': [], 'elapsed_time': 0,
                      'cpu_budget': self.cpu_pool.total, 'cpu_per_job_cap': self.cpu,
                      'requested_workers': num_workers, 'effective_workers': effective_workers}
        self.docking_results = {}
        available_proteins = {p.stem: p for p in self.get_proteins()}
        targets = list(available_proteins) if proteins is None else list(dict.fromkeys(proteins))
        if ligand_paths is None:
            available_ligands = {p.stem: p for p in self.get_ligands()}
        else:
            available_ligands = {name: Path(path) for name, path in ligand_paths.items()}
        labels = list(available_ligands) if ligands is None else list(dict.fromkeys(ligands))
        self.stats['total_dockings'] = len(targets) * len(labels)
        tasks = []
        for protein in targets:
            config = self.get_config(protein) if protein in available_proteins else None
            for label in labels:
                if config is None or label not in available_ligands:
                    self._count('failed', (label, protein, 'Missing input/config'))
                    continue
                tasks.append((label, protein, available_ligands[label], available_proteins[protein], config))
        from src.ui.progress import WorkProgress
        with WorkProgress("Docking", self.stats['total_dockings']) as progress:
            progress.update(completed=self.stats['failed'])
            with ThreadPoolExecutor(max_workers=effective_workers) as executor:
                futures = {executor.submit(self.run_single_docking, lig, prot, conf, exhaustiveness, force):
                           (label, protein) for label, protein, lig, prot, conf in tasks}
                for future in as_completed(futures):
                    label, protein = futures[future]
                    try:
                        result = future.result()
                        if result is not None:
                            result = {**result, 'ligand': label, 'protein': protein}
                            self.docking_results.setdefault(label, {})[protein] = result
                    except Exception as error:
                        self._count('failed', (label, protein, str(error)))
                    finally:
                        progress.advance()
                        progress.update(description=f"Docking: cached {self.stats['cached']}; failed {self.stats['failed']}")
        self.stats['elapsed_time'] = time.monotonic() - started
        self._save_results_json()
        return self.docking_results

    def _save_results_json(self):
        # Preserve an immutable batch record; docking_results.json is a latest-result
        # export for existing callers, never an input to the calculation cache.
        self.last_batch_dir = self.project_root / "output" / "runs" / ("docking-" + uuid.uuid4().hex)
        atomic_json(self.last_batch_dir / "results.json", self.docking_results)
        atomic_json(self.last_batch_dir / "statistics.json", self.stats)
        atomic_json(self.cache_dir / "docking_results.json", self.docking_results)

    def load_results_json(self) -> Dict:
        """
        Load docking results from a JSON file.
        
        Returns:
            Result dictionary
        """
        json_file = self.cache_dir / "docking_results.json"
        
        if not json_file.exists():
            self.logger.warning("Results file not found")
            return {}
        
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                results = json.load(f)
            
            self.docking_results = results
            self.logger.info(f"Loaded results for {len(results)} ligands")
            return results
        
        except Exception as e:
            self.logger.error(f"Failed to load JSON: {e}")
            return {}
    
    def get_top_ligands(self, protein: str, n: int = 10) -> List[Tuple[str, float]]:
        """
        Return the top N ligands for a receptor ranked by docking score.
        
        Args:
            protein: Receptor name
            n: Number of ligands
            
        Returns:
            List of (ligand_name, affinity) tuples
        """
        ligands_affinities = []
        
        for ligand_name, proteins_data in self.docking_results.items():
            if protein in proteins_data:
                affinity = proteins_data[protein]['affinity']
                ligands_affinities.append((ligand_name, affinity))
        
        # Sort by docking score (lower is better)
        ligands_affinities.sort(key=lambda x: x[1])
        
        return ligands_affinities[:n]

