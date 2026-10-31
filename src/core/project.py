"""project - FragEvoLab module."""

import zipfile
from pathlib import Path
from typing import List, Dict, Optional
from datetime import datetime
import logging

class ProjectUtils:
    """
    Project utility functions.
    """
    
    def __init__(self, project_root: Optional[Path] = None):
        """
        Initialize project utilities.
        
        Args:
            project_root: Project root directory (default: current directory)
        """
        self.project_root = project_root or Path.cwd()
        
        # Directory paths
        self.data_dir = self.project_root / "data"
        self.ligands_csv = self.data_dir / "ligands.csv"
        self.ligands_pdbqt_dir = self.project_root / "output" / "ligands_pdbqt"
        self.docking_results_dir = self.project_root / "output" / "docking_results"
        self.reports_dir = self.project_root / "output" / "reports"
        self.cache_dir = self.project_root / "output" / "cache"
        
        # Configure logging
        self.logger = logging.getLogger(__name__)
    
    def export_results_to_zip(self,
                              output_file: Optional[Path] = None,
                              include_ligands_csv: bool = True,
                              include_ligands_pdbqt: bool = False,
                              include_docking_pdbqt: bool = True,
                              include_reports: bool = True,
                              include_logs: bool = False) -> Path:
        """
        Create a ZIP archive containing project results.
        
        Archive contents:
        - Ligand CSV file
        - PDBQT files containing docked ligand poses
        - HTML and CSV reports
        - Run metadata, immutable inputs/preparations and a SQLite snapshot
        - Optional: logs and legacy ligand PDBQT aliases

        Export after calculations finish. MANIFEST.json maps original paths to
        archive paths. This is a result collection, not a runnable project.
        
        Args:
            output_file: ZIP file path (project root by default)
            include_ligands_csv: Include the ligand CSV file
            include_ligands_pdbqt: Include input ligand PDBQT files
            include_docking_pdbqt: Include docking PDBQT results
            include_reports: Include reports
            include_logs: Include logs
            
        Returns:
            Path to the created ZIP file
        """
        import hashlib
        import json
        import sqlite3
        import tempfile
        from contextlib import closing
        from src.core.calculation_cache import atomic_write

        if output_file is None:
            output_file = self.project_root / f"fragevolab_results_{datetime.now():%Y%m%d_%H%M%S}.zip"
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        selected = {}
        def collect(folder, prefix, pattern="*"):
            if folder.exists():
                for path in sorted(folder.rglob(pattern)):
                    if path.is_file() and path.resolve() != output_file.resolve():
                        selected[prefix + path.relative_to(folder).as_posix()] = path

        # Export after calculations finish. Immutable inputs are provenance, not
        # optional legacy aliases. This is a result collection, not a cache backup.
        for name in ("runs", "inputs", "prepared"):
            collect(self.project_root / "output" / name, f"output/{name}/")
        if include_ligands_csv and self.ligands_csv.exists():
            selected["data/ligands.csv"] = self.ligands_csv
        if include_ligands_pdbqt:
            collect(self.ligands_pdbqt_dir, "ligands_pdbqt/", "*.pdbqt")
        if include_docking_pdbqt:
            collect(self.docking_results_dir, "docking_results/", "*.pdbqt")
        if include_reports:
            collect(self.reports_dir, "reports/")
            collect(self.project_root / "output/evolutionary", "evolutionary/")
        if include_logs:
            collect(self.cache_dir, "logs/", "*.log")
            collect(self.cache_dir, "logs/", "*.json")
        for name in ("evolution_settings.json", "cpu_settings.json", "app_settings.json"):
            path = self.project_root / name
            if path.exists(): selected[name] = path
        manifest = {"format": 1, "kind": "result_export", "source_root": str(self.project_root.resolve()),
                    "created_at": datetime.now().isoformat(),
                    "note": "Export after calculations finish. Absolute source paths are mapped by the files table. This is not a runnable project or a calculation-cache backup.",
                    "options": {"include_ligands_csv": include_ligands_csv,
                                "include_ligands_pdbqt": include_ligands_pdbqt,
                                "include_docking_pdbqt": include_docking_pdbqt,
                                "include_reports": include_reports, "include_logs": include_logs},
                    "files": {}}
        with tempfile.TemporaryDirectory(prefix="fragevolab-export-") as temp:
            temp = Path(temp)
            db = self.project_root / "output/experiments.db"
            if db.exists():
                snapshot = temp / "experiments.db"
                with closing(sqlite3.connect(db.resolve().as_uri() + "?mode=ro", uri=True)) as source:
                    with closing(sqlite3.connect(snapshot)) as destination:
                        source.backup(destination)
                selected["output/experiments.db"] = snapshot
            archive = temp / "results.zip"
            with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zipf:
                for name, path in selected.items():
                    data = path.read_bytes()
                    zipf.writestr(name, data)
                    source_path = db if name == "output/experiments.db" else path
                    manifest["files"][name] = {"sha256": hashlib.sha256(data).hexdigest(),
                                               "size": len(data), "source_path": str(source_path.resolve())}
                zipf.writestr("MANIFEST.json", json.dumps(manifest, indent=2))
                zipf.writestr("CONTENTS.txt", "FragEvoLab result export\n\n" + "\n".join(zipf.namelist()) + "\n")
            atomic_write(output_file, archive.read_bytes())
        self.logger.info("ZIP archive created: %s", output_file)
        return output_file

    def load_ligands_from_csv(self, csv_path: Optional[Path] = None) -> List[Dict[str, str]]:
        """
        Load ligands from a CSV file.
        
        Args:
            csv_path: CSV file path
            
        Returns:
            List of ligand data dictionaries
        """
        if csv_path is None:
            csv_path = self.ligands_csv
        
        if not csv_path.exists():
            self.logger.error(f"CSV file not found: {csv_path}")
            return []
        
        ligands = []
        
        try:
            from src.preparation.inputs import read_ligand_rows, unambiguous_rows
            ligands=[{"name":name,"smiles":smiles} for name,smiles in unambiguous_rows(read_ligand_rows(csv_path))]

            self.logger.info(f"Ligands loaded: {len(ligands)}")
        
        except Exception as e:
            self.logger.error(f"Failed to read CSV: {e}")
        
        return ligands

    def validate_project_structure(self) -> Dict[str, bool]:
        """
        Check the project structure and required files.
        
        Returns:
            Dictionary of check results
        """
        checks = {
            'data_dir': self.data_dir.exists(),
            'proteins_dir': (self.data_dir / "proteins").exists(),
            'configs_dir': (self.data_dir / "configs").exists(),
            'ligands_csv': self.ligands_csv.exists(),
            'output_dir': (self.project_root / "output").exists(),
            'src_dir': (self.project_root / "src").exists()
        }
        
        # Additional checks
        if checks['proteins_dir']:
            proteins = list((self.data_dir / "proteins").glob("*.pdbqt"))
            checks['proteins_count'] = len(proteins)
        else:
            checks['proteins_count'] = 0
        
        if checks['configs_dir']:
            configs = list((self.data_dir / "configs").glob("*.txt"))
            checks['configs_count'] = len(configs)
        else:
            checks['configs_count'] = 0
        
        return checks
    
    def print_project_status(self):
        """Display project status."""
        checks = self.validate_project_structure()
        
        print("PROJECT STATUS")
        
        print("\nStructure:")
        print(f"  {'✓' if checks['data_dir'] else '✗'} Directory data/")
        print(f"  {'✓' if checks['proteins_dir'] else '✗'} Directory data/proteins/ (receptors: {checks['proteins_count']})")
        print(f"  {'✓' if checks['configs_dir'] else '✗'} Directory data/configs/ (configurations: {checks['configs_count']})")
        print(f"  {'✓' if checks['ligands_csv'] else '✗'} File data/ligands.csv")
        print(f"  {'✓' if checks['output_dir'] else '✗'} Directory output/")
        print(f"  {'✓' if checks['src_dir'] else '✗'} Directory src/")
        
        # Statistics
        ligands = self.load_ligands_from_csv()
        ligands_pdbqt = list(self.ligands_pdbqt_dir.glob("*.pdbqt")) if self.ligands_pdbqt_dir.exists() else []
        docking_results = list(self.docking_results_dir.rglob("*.pdbqt")) if self.docking_results_dir.exists() else []
        reports = list(self.reports_dir.glob("*.html")) if self.reports_dir.exists() else []
        
        print("\nStatistics:")
        print(f"  Ligands in CSV: {len(ligands)}")
        print(f"  PDBQT files generated: {len(ligands_pdbqt)}")
        print(f"  Docking results: {len(docking_results)}")
        print(f"  Reports: {len(reports)}")

    def get_available_proteins(self) -> List[str]:
        """
        Return available receptors.
        
        Returns:
            List of receptor names (without extensions)
        """
        proteins_dir = self.data_dir / "proteins"
        
        if not proteins_dir.exists():
            return []
        
        proteins = [p.stem for p in proteins_dir.glob("*.pdbqt")]
        return sorted(proteins)

def export_results(output_file: Optional[Path] = None,
                  include_all: bool = True,
                  project_root: Optional[Path] = None) -> Path:
    """
    Convenience function to export results to ZIP.
    
    Args:
        output_file: ZIP file path
        include_all: Include all available data
        project_root: Project root directory
        
    Returns:
        Path to the created ZIP file
    """
    utils = ProjectUtils(project_root=project_root)
    
    return utils.export_results_to_zip(
        output_file=output_file,
        include_ligands_csv=True,
        include_ligands_pdbqt=False,
        include_docking_pdbqt=True,
        include_reports=True,
        include_logs=include_all
    )

