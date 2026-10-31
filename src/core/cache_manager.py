"""cache manager - FragEvoLab module."""

import os
import tempfile
import uuid
from src.core.calculation_cache import CalculationCache, digest, atomic_json
import shutil
import json
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime, timedelta
import logging

class CacheManager:
    """
    Manage the project cache.
    
    Features:
    - Cache usage analysis
    - Selective calculation-cache cleanup by type
    - Storage size statistics
    - Backup
    """
    
    def __init__(self, project_root: Optional[Path] = None):
        """
        Initialize the cache manager.
        
        Args:
            project_root: Project root directory (default: current directory)
        """
        self.project_root = Path(project_root or Path.cwd()).resolve()
        self.calculation_cache = CalculationCache(self.project_root / ".calculation_cache" / "v1")
        
        # Paths to cache directories
        self.cache_locations = {
            'ligands_pdbqt': self.project_root / "output" / "ligands_pdbqt",
            'docking_results': self.project_root / "output" / "docking_results",
            'cache': self.project_root / "output" / "cache",
            'reports': self.project_root / "output" / "reports"
        }
        
        # Configure logging
        self.logger = logging.getLogger(__name__)
    
    def get_directory_size(self, directory: Path) -> int:
        """
        Calculate directory size in bytes.
        
        Args:
            directory: Directory path
            
        Returns:
            Size in bytes
        """
        total_size = 0
        
        try:
            for dirpath, dirnames, filenames in os.walk(directory):
                for filename in filenames:
                    filepath = os.path.join(dirpath, filename)
                    if os.path.exists(filepath):
                        total_size += os.path.getsize(filepath)
        except Exception as e:
            self.logger.error(f"Failed to calculate the size of {directory}: {e}")
        
        return total_size
    
    def format_size(self, size_bytes: int) -> str:
        """
        Format a size as a human-readable string.
        
        Args:
            size_bytes: Size in bytes
            
        Returns:
            Formatted string
        """
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if size_bytes < 1024.0:
                return f"{size_bytes:.2f} {unit}"
            size_bytes /= 1024.0
        
        return f"{size_bytes:.2f} PB"
    
    def get_file_count(self, directory: Path, pattern: str = "*") -> int:
        """
        Count files in a directory.
        
        Args:
            directory: Directory path
            pattern: Filename pattern (for example, "*.pdbqt")
            
        Returns:
            Number of files
        """
        if not directory.exists():
            return 0
        
        return len(list(directory.glob(pattern)))
    
    def get_cache_statistics(self) -> Dict:
        """
        Collect statistics for all cache directories.
        
        Returns:
            Statistics dictionary
        """
        stats = {
            'total_size': 0,
            'total_files': 0,
            'locations': {}
        }
        
        for name, path in self.cache_locations.items():
            if path.exists():
                size = self.get_directory_size(path)
                files = self.get_file_count(path)
                
                stats['locations'][name] = {
                    'path': str(path),
                    'size': size,
                    'size_formatted': self.format_size(size),
                    'files': files
                }
                
                stats['total_size'] += size
                stats['total_files'] += files
            else:
                stats['locations'][name] = {
                    'path': str(path),
                    'size': 0,
                    'size_formatted': '0 B',
                    'files': 0
                }
        
        stats['calculation_entries'] = self.calculation_cache.statistics()
        stats['calculation_cache_bytes'] = self.get_directory_size(self.calculation_cache.root)
        stats['total_size_formatted'] = self.format_size(stats['total_size'])
        
        return stats
    
    def print_statistics(self):
        """Display cache statistics."""
        stats = self.get_cache_statistics()
        
        print("CACHE STATISTICS")
        print(f"\nTotal size: {stats['total_size_formatted']}")
        print(f"Total files: {stats['total_files']}")
        print(f"Reusable calculation entries: {stats['calculation_entries']}")
        print(f"SQLite cache size: {self.format_size(stats['calculation_cache_bytes'])}")
        print("\nBreakdown by directory:")
        
        for name, data in stats['locations'].items():
            print(f"\n{name}:")
            print(f"  Path: {data['path']}")
            print(f"  Size: {data['size_formatted']}")
            print(f"  Files: {data['files']}")

    def get_old_files(self, directory: Path, days: int = 7) -> List[Path]:
        """
        Find files older than the specified number of days.
        
        Args:
            directory: Directory to search
            days: File age in days
            
        Returns:
            List of paths to old files
        """
        if not directory.exists():
            return []
        
        cutoff_date = datetime.now() - timedelta(days=days)
        old_files = []
        
        for file in directory.rglob("*"):
            if file.is_file():
                file_time = datetime.fromtimestamp(file.stat().st_mtime)
                if file_time < cutoff_date:
                    old_files.append(file)
        
        return old_files
    
    def clear_cache(self, location: Optional[str] = None,
                    older_than_days: Optional[int] = None, confirm: bool = True) -> Dict:
        """Delete reusable records, never scientific outputs. Stop workers first."""
        kinds = {None: None, 'cache': None, 'ligands_pdbqt': 'preparation',
                 'docking_results': 'docking'}
        if location not in kinds:
            raise ValueError("Reports are experiment history, not disposable cache")
        if older_than_days is not None:
            raise ValueError("Age-based eviction is not implemented; use an explicit cache clear")
        if confirm and input("Clear reusable calculation cache? (yes/no): ").strip().lower() not in (
                'yes', 'y', 'да', 'д'):
            return {'deleted_files': 0, 'deleted_entries': 0, 'freed_space': 0}
        count = self.calculation_cache.clear(kinds[location])
        # SQLite can reuse free pages. File size need not shrink after DELETE.
        return {'deleted_files': 0, 'deleted_entries': count, 'freed_space': 0}

    def clear_ligands_cache(self, confirm: bool = True) -> Dict:
        """
        Clear the prepared ligand cache.
        
        Args:
            confirm: Request confirmation
            
        Returns:
            Deletion statistics
        """
        return self.clear_cache(location='ligands_pdbqt', confirm=confirm)
    
    def clear_docking_cache(self, confirm: bool = True) -> Dict:
        """
        Clear the docking results cache.
        
        Args:
            confirm: Request confirmation
            
        Returns:
            Deletion statistics
        """
        return self.clear_cache(location='docking_results', confirm=confirm)
    
    @staticmethod
    def _inventory(directory: Path) -> Dict[str, str]:
        if not directory.is_dir():
            raise ValueError("Backup output directory missing")
        result = {}
        for path in sorted(directory.rglob("*")):
            if path.is_symlink():
                raise ValueError("Symlinks are not supported in output backups")
            if path.is_file():
                result[path.relative_to(directory).as_posix()] = digest(path.read_bytes())
        if not result:
            raise ValueError("Refusing an empty output backup")
        return result

    def backup_cache(self, backup_dir: Optional[Path] = None) -> Path:
        """Offline snapshot of scientific outputs. Reusable cache is not copied."""
        output = self.project_root / 'output'
        before = self._inventory(output)
        parent = Path(backup_dir or self.project_root / 'backups').resolve()
        if parent == output or output in parent.parents:
            raise ValueError("Backup destination must be outside output")
        parent.mkdir(parents=True, exist_ok=True)
        destination = parent / ('cache_backup_' + uuid.uuid4().hex)
        destination.mkdir()
        try:
            shutil.copytree(output, destination / 'output')
            if self._inventory(destination / 'output') != before or self._inventory(output) != before:
                raise ValueError("Output changed during backup; stop workers first")
            atomic_json(destination / 'metadata.json', {'format': 1, 'sha256': before,
                        'datetime': datetime.now().isoformat()})
        except Exception:
            shutil.rmtree(destination)
            raise
        return destination

    def restore_from_backup(self, backup_path: Path, confirm: bool = True) -> bool:
        """Validate and stage first; preserve prior output and roll back on failure.

        Call only after all calculations stop. Legacy backups without hashes are
        rejected; inspect/copy them manually instead of silently trusting them.
        """
        staging = None
        output = self.project_root / 'output'
        previous = self.project_root / ('output_before_restore_' + uuid.uuid4().hex)
        try:
            backup_path = Path(backup_path).resolve()
            if backup_path == output or output in backup_path.parents:
                raise ValueError("Backup must be outside the output being replaced")
            manifest = json.loads((backup_path / 'metadata.json').read_text(encoding='utf-8'))
            if manifest.get('format') != 1 or not manifest.get('sha256'):
                raise ValueError("Backup has no verifiable manifest")
            if self._inventory(backup_path / 'output') != manifest['sha256']:
                raise ValueError("Backup checksum mismatch")
            staging = Path(tempfile.mkdtemp(prefix='restore-', dir=self.project_root))
            shutil.copytree(backup_path / 'output', staging / 'output')
            if self._inventory(staging / 'output') != manifest['sha256']:
                raise ValueError("Staged backup checksum mismatch")
            if confirm and input("Restore verified backup? (yes/no): ").strip().lower() not in (
                    'yes', 'y', 'да', 'д'):
                return False
            if output.exists():
                output.rename(previous)
            try:
                (staging / 'output').rename(output)
            except Exception:
                if previous.exists():
                    previous.rename(output)
                raise
            self.logger.info("Restored output; previous files retained at %s", previous)
            return True
        except Exception as error:
            self.logger.error("Restore failed without deleting current output: %s", error)
            return False
        finally:
            if staging is not None:
                shutil.rmtree(staging, ignore_errors=True)

    def optimize_cache(self) -> Dict:
        """Offline SQLite compaction; never remove files from experiment outputs."""
        from contextlib import closing
        with closing(self.calculation_cache._connect()) as db:
            db.execute('VACUUM')
        return {'removed_files': 0, 'removed_temp_files': 0, 'freed_space': 0,
                'freed_space_formatted': self.format_size(0), 'vacuum_completed': True}

    def get_cache_info_json(self) -> str:
        """
        Return cache information as JSON.
        
        Returns:
            JSON string containing statistics
        """
        stats = self.get_cache_statistics()
        return json.dumps(stats, indent=2, ensure_ascii=False)

