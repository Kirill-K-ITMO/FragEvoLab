"""ligands - FragEvoLab module."""

from __future__ import annotations

import csv
import math
from threading import local
import importlib.metadata
import platform
import tempfile
import shutil
from src.core.calculation_cache import (CalculationCache, molecule_identity, digest,
                                   atomic_write, atomic_json, validate_pdbqt)
from src.preparation.inputs import read_ligand_rows
from filelock import FileLock
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit is not installed. Functionality is limited.")

try:
    from meeko import MoleculePreparation
    from meeko import PDBQTWriterLegacy
    MEEKO_AVAILABLE = True
except ImportError:
    MEEKO_AVAILABLE = False
    print("Warning: Meeko is not installed. OpenBabel will be used as a less reliable fallback.")

class LigandConverter:
    """
    Convert SMILES strings to PDBQT for AutoDock Vina.
    
    Supports two conversion methods:
    1. RDKit + Meeko (recommended, more accurate)
    2. OpenBabel (fallback)
    """
    
    def __init__(self, project_root: Optional[Path] = None):
        """
        Initialize the converter.
        
        Args:
            project_root: Project root directory (default: current directory)
        """
        self.project_root = Path(project_root or Path.cwd()).resolve()
        self.calculation_cache = CalculationCache(self.project_root / ".calculation_cache" / "v1")
        self.preparation_seed = 42
        self.optimization_iters = (200, 1000)
        self.preparation_timeout = 180
        self._preparation_local = local()
        self.prepared_results = {}
        self.data_dir = self.project_root / "data"
        self.output_dir = self.project_root / "output" / "ligands_pdbqt"
        self.cache_dir = self.project_root / "output" / "cache"
        
        # Create the required directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Conversion statistics
        self.stats = {
            'total': 0,
            'success': 0,
            'failed': 0,
            'cached': 0,
            'errors': []
        }
        
        # Configure logging
        self._setup_logging()
        
        # Check tool availability
        self._check_tools()
    
    def _setup_logging(self):
        """Configure logging to track conversion progress."""
        logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
        self.logger = logging.getLogger(__name__)
    
    def _check_tools(self):
        """Check the availability of conversion tools."""
        self.use_meeko = RDKIT_AVAILABLE and MEEKO_AVAILABLE
        self.use_obabel = False
        
        if self.use_meeko:
            self.logger.info("Using RDKit + Meeko (recommended)")
        else:
            # Check OpenBabel availability
            try:
                subprocess.run(
                    ['obabel', '-V'],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    check=True
                )
                self.use_obabel = True
                self.logger.info("Using OpenBabel (fallback)")
            except (subprocess.CalledProcessError, FileNotFoundError):
                self.logger.error(
                    "Neither RDKit+Meeko nor OpenBabel is available. "
                    "Install one of the conversion backends."
                )
    
    @property
    def last_preparation_details(self):
        return getattr(self._preparation_local, "details", {})

    def _validated_molecule(self, smiles):
        if not RDKIT_AVAILABLE:
            raise RuntimeError("RDKit is required for validated ligand preparation")
        from src.preparation.validation import validated_molecule
        return validated_molecule(smiles)

    def smiles_to_pdbqt_meeko(self, smiles, name, output_path):
        output_path=Path(output_path)
        self._preparation_local.details={"backend":"meeko", "status":"started"}
        details=self.last_preparation_details
        try:
            if not MEEKO_AVAILABLE: raise RuntimeError("Meeko unavailable")
            if not self.optimization_iters or any(type(n) is not int or n<1 for n in self.optimization_iters):
                raise ValueError("Optimization iteration limits must be positive integers")
            mol=Chem.AddHs(self._validated_molecule(smiles))
            details["embedding_attempts"]=[]
            for random_coords in (False,True):
                params=AllChem.ETKDGv3()
                params.randomSeed=self.preparation_seed
                params.useRandomCoords=random_coords
                params.numThreads=1
                params.maxIterations=1000
                code=AllChem.EmbedMolecule(mol,params)
                details["embedding_attempts"].append({"random_coords":random_coords,"return_code":code})
                if code==0: break
            if code!=0: raise ValueError("3D embedding failed after two ETKDGv3 attempts")
            if not AllChem.UFFHasAllMoleculeParams(mol):
                raise ValueError("UFF parameters are unavailable for this ligand")
            details["optimization_attempts"]=[]
            for iterations in self.optimization_iters:
                code=AllChem.UFFOptimizeMolecule(mol,maxIters=iterations)
                details["optimization_attempts"].append({"max_iters":iterations,"return_code":code})
                if code==0: break
                if code!=1: raise ValueError(f"UFF failed with code {code}")
                self.logger.info("%s: UFF not converged after %s iterations; continuing within configured retry budget",name,iterations)
            if code!=0: raise ValueError("UFF did not converge within the permitted optimization attempts")
            if not all(math.isfinite(v) for xyz in mol.GetConformer().GetPositions() for v in xyz):
                raise ValueError("Non-finite coordinates")
            setups=MoleculePreparation().prepare(mol)
            if not setups or len(setups)!=1: raise ValueError("Expected exactly one Meeko setup")
            text,success,error=PDBQTWriterLegacy.write_string(setups[0])
            if not success: raise ValueError(f"Meeko writer: {error}")
            artifact=text.encode("utf-8")
            validate_pdbqt(artifact)
            atomic_write(output_path,artifact)
            details.update(status="success",optimization="UFF",converged=True)
            return True
        except Exception as error:
            output_path.unlink(missing_ok=True)
            details.update(status="failed",error=str(error))
            self.logger.error("%s: %s",name,error)
            return False

    def smiles_to_pdbqt_obabel(self, smiles, name, output_path):
        output_path=Path(output_path)
        self._preparation_local.details={"backend":"obabel","status":"started"}
        try:
            self._validated_molecule(smiles)
            with tempfile.TemporaryDirectory(prefix="obabel-ligand-") as directory:
                source=Path(directory)/"ligand.smi"; target=Path(directory)/"ligand.pdbqt"
                canonical,_=molecule_identity(smiles)
                source.write_text(canonical+" ligand\n",encoding="utf-8")
                command=["obabel",str(source),"-O",str(target),"--gen3d","-h"]
                subprocess.run(command,capture_output=True,text=True,check=True,timeout=self.preparation_timeout)
                artifact=target.read_bytes(); validate_pdbqt(artifact)
                atomic_write(output_path,artifact)
            self.last_preparation_details.update(status="success",converged=None,
                optimization="OpenBabel backend; UFF convergence not verified by this adapter")
            return True
        except Exception as error:
            output_path.unlink(missing_ok=True)
            self.last_preparation_details.update(status="failed",error=str(error))
            self.logger.error("%s: %s",name,error)
            return False

    def preparation_spec(self, smiles: str) -> dict:
        if type(self.preparation_seed) is not int or not 0 <= self.preparation_seed <= 2147483647:
            raise ValueError("preparation_seed must be in 0..2147483647")
        self._validated_molecule(smiles)
        canonical, molecule_id = molecule_identity(smiles)
        backend = "meeko" if self.use_meeko else "obabel" if self.use_obabel else None
        if backend is None:
            raise RuntimeError("No ligand preparation backend available")
        versions = {}
        for package in ("rdkit", "meeko", "numpy", "gemmi", "scipy"):
            try:
                versions[package] = importlib.metadata.version(package)
            except importlib.metadata.PackageNotFoundError:
                versions[package] = None
        spec = {"molecule_id": molecule_id, "canonical_smiles": canonical,
                "backend": backend, "versions": versions,
                "platform": platform.platform(), "python": platform.python_version(),
                "implementation_sha": digest(Path(__file__).read_bytes()),
                "graph_policy_sha": digest(Path(__file__).with_name("validation.py").read_bytes()),
                "seed": self.preparation_seed, "hydrogens": "add",
                "embedding": "ETKDGv3_then_random_ETKDGv3", "max_attempts": 1000,
                "optimization": "UFF_strict" if backend == "meeko" else "backend_defined", "max_iters": list(self.optimization_iters),
                "preparation_timeout": self.preparation_timeout,
                "obabel_rigid": False}
        if backend == "obabel":
            executable = shutil.which("obabel")
            if not executable:
                raise RuntimeError("Open Babel executable not found")
            spec["obabel_binary_sha"] = digest(Path(executable).read_bytes())
        return spec

    def prepare_cached(self, smiles: str, force: bool = False) -> dict:
        """Return an exact, validated artifact path; never infer success from an alias."""
        spec = self.preparation_spec(smiles)
        def compute():
            with tempfile.TemporaryDirectory(prefix="prepare-") as directory:
                path = Path(directory) / "ligand.pdbqt"
                method = (self.smiles_to_pdbqt_meeko if spec["backend"] == "meeko"
                          else self.smiles_to_pdbqt_obabel)
                if not method(spec["canonical_smiles"], spec["molecule_id"], path):
                    raise RuntimeError(self.last_preparation_details.get("error", "Ligand preparation failed"))
                artifact = path.read_bytes()
                validate_pdbqt(artifact)
                return {"spec": spec, "molecule_id": spec["molecule_id"],
                        "canonical_smiles": spec["canonical_smiles"],
                        "preparation_details": dict(self.last_preparation_details)}, artifact
        entry = self.calculation_cache.get_or_compute("preparation", spec, compute, force=force)
        # Content hash also distinguishes forced, potentially stochastic preparations.
        path = (self.project_root / "output" / "prepared" / entry.key /
                (digest(entry.artifact) + ".pdbqt"))
        if not path.exists() or digest(path.read_bytes()) != digest(entry.artifact):
            atomic_write(path, entry.artifact)
        return {**entry.payload, "path": str(path), "preparation_id": entry.key,
                "artifact_sha": digest(entry.artifact), "cache_hit": entry.hit}

    def convert_single(self, name: str, smiles: str, force: bool = False) -> bool:
        """Compatibility API: publish a verified alias only after preparation succeeds."""
        self.validate_name(name)
        output_file = self.output_dir / f"{name}.pdbqt"
        sidecar = output_file.with_suffix(".json")
        alias_key = digest(str(output_file).encode())
        with FileLock(str(self.calculation_cache.locks / f"alias-{alias_key}.lock")):
            self.prepared_results.pop(name, None)
            try:
                prepared = self.prepare_cached(smiles, force=force)
                atomic_write(output_file, Path(prepared["path"]).read_bytes())
                atomic_json(sidecar, prepared)
                self.prepared_results[name] = prepared
                self.stats['cached' if prepared['cache_hit'] else 'success'] += 1
                return True
            except Exception as error:
                # An old alias must not survive a failed conversion under the same label.
                output_file.unlink(missing_ok=True)
                sidecar.unlink(missing_ok=True)
                self.stats['failed'] += 1
                self.stats['errors'].append((name, str(error)))
                self.logger.error("%s: %s", name, error)
                return False

    @staticmethod
    def validate_name(name):
        if not isinstance(name,str) or not name or name!=name.strip() or name.endswith(("."," ")):
            raise ValueError("Empty or unsafe ligand filename")
        if any(ord(c)<32 or c in '<>:"/\\|?*' for c in name):
            raise ValueError("Ligand name contains forbidden filename characters")
        if name.split('.')[0].upper() in {"CON","PRN","AUX","NUL",*[f"COM{i}" for i in range(1,10)],*[f"LPT{i}" for i in range(1,10)]}:
            raise ValueError("Reserved Windows filename")
        return name

    def _reset_stats(self):
        self.stats={"total":0,"success":0,"failed":0,"cached":0,"errors":[],"input_error":None}
        self.prepared_results={}

    def convert_from_csv(self, csv_path=None, force=False):
        self._reset_stats()
        path=Path(csv_path) if csv_path is not None else self.data_dir/"ligands.csv"
        try:
            rows=read_ligand_rows(path)
            return self.convert_list(rows,force=force)
        except (OSError,UnicodeError,ValueError,csv.Error) as error:
            self.stats["input_error"]=str(error)
            self.logger.error("CSV input error: %s",error)
            return self.stats

    def convert_list(self, ligands, force=False):
        from collections import Counter
        self._reset_stats()
        rows=list(ligands);self.stats["total"]=len(rows)
        counts=Counter(str(name).strip().casefold() for name,smiles in rows)
        for name,smiles in rows:
            try:
                name=self.validate_name(name)
                if counts[name.casefold()]>1: raise ValueError("Duplicate ligand name in this batch")
                if not isinstance(smiles,str) or not smiles.strip(): raise ValueError("Empty SMILES")
                self.convert_single(name,smiles.strip(),force=force)
            except (ValueError,TypeError) as error:
                try:
                    safe_name = self.validate_name(name)
                except (ValueError, TypeError):
                    safe_name = None
                if safe_name is not None:
                    output_file = self.output_dir / f"{safe_name}.pdbqt"
                    alias_key = digest(str(output_file).encode())
                    with FileLock(str(self.calculation_cache.locks / f"alias-{alias_key}.lock")):
                        output_file.unlink(missing_ok=True)
                        output_file.with_suffix(".json").unlink(missing_ok=True)
                    self.prepared_results.pop(safe_name, None)
                self.stats["failed"]+=1;self.stats["errors"].append((str(name),str(error)))
                self.logger.error("%s: %s",name,error)
        return self.stats

    def get_converted_ligands(self) -> List[Path]:
        """
        Return all converted ligand PDBQT files.
        
        Returns:
            List of PDBQT file paths
        """
        pdbqt_files = sorted(self.output_dir.glob("*.pdbqt"))
        return pdbqt_files
    
    def clear_cache(self, confirm: bool = False) -> bool:
        """
        Clear the cache (delete all ligand PDBQT files).
        
        Args:
            confirm: If True, delete without requesting confirmation
            
        Returns:
            True if the cleanup was completed
        """
        pdbqt_files = self.get_converted_ligands()
        
        if not pdbqt_files and not self.calculation_cache.statistics().get("preparation"):
            self.logger.info("No files to delete")
            return True
        
        if not confirm:
            response = input(f"\nDelete {len(pdbqt_files)} ligand PDBQT files? (yes/no): ")
            if response.lower() not in ['да', 'yes', 'y', 'д']:
                self.logger.info("Cleanup cancelled")
                return False
        
        self.calculation_cache.clear("preparation")
        # Delete files
        deleted = 0
        for file in pdbqt_files:
            try:
                file.unlink()
                deleted += 1
            except Exception as e:
                self.logger.error(f"Failed to delete {file.name}: {e}")
        
        self.logger.info(f"✓ Files deleted: {deleted}/{len(pdbqt_files)}")
        return True

def convert_smiles_to_pdbqt(smiles_list: List[Tuple[str, str]], 
                            project_root: Optional[Path] = None,
                            force: bool = False) -> Dict:
    """
    Convenience function to convert a list of SMILES.
    
    Args:
        smiles_list: List of (name, SMILES) tuples
        project_root: Project root directory
        force: Overwrite existing files
        
    Returns:
        Dictionary of conversion statistics
        
    Example:
        >>> ligands = [('aspirin', 'CC(=O)Oc1ccccc1C(=O)O'), ('caffeine', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C')]
        >>> stats = convert_smiles_to_pdbqt(ligands)
    """
    converter = LigandConverter(project_root=project_root)
    return converter.convert_list(ligands=smiles_list, force=force)

