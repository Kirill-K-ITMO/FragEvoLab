from src.version import APP_NAME, __version__
"""application - FragEvoLab module."""

import logging
from pathlib import Path
from typing import Optional, Dict

from src.core.cpu_budget import configure_project_cpu, project_budget, available_cpus
from src.preparation.smiles import SmilesGenerator
from src.preparation.ligands import LigandConverter
from src.docking.vina import AutoDockVina
from src.docking.affinity import AffinityCalculator
from src.reporting.docking import ReportGenerator
from src.core.cache_manager import CacheManager
from src.core.project import ProjectUtils
from src.core.database import DockingDatabase
from src.ui.console import (
    print_app_header, print_header, print_success,
    print_error, print_info, console
)

class FragEvoLabApp:
    """
    Main FragEvoLab application class.
    Integrates the modules and provides the terminal interface.
    """
    
    def __init__(self, project_root: Optional[Path] = None):
        """
        Initialize the application.
        
        Args:
            project_root: Project root directory
        """
        self.project_root = Path(project_root or Path(__file__).resolve().parents[2]).resolve()
        
        project_budget(self.project_root)

        # Initialize modules
        self.smiles_gen = SmilesGenerator(project_root=self.project_root)
        self.converter = LigandConverter(project_root=self.project_root)
        self._docker = None
        self.calculator = AffinityCalculator()
        self.reporter = ReportGenerator(project_root=self.project_root)
        self.cache_mgr = CacheManager(project_root=self.project_root)
        self.utils = ProjectUtils(project_root=self.project_root)
        self.database = DockingDatabase(
            db_file=self.project_root / "output" / "experiments.db"
        )
        
        # Configure logging
        self._setup_logging()
        
        self.logger = logging.getLogger(__name__)
    
    def _setup_logging(self):
        """Configure application logging."""
        log_file = self.project_root / "output" / "cache" / "fragevolab.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        
        self._log_handler = logging.FileHandler(log_file, encoding="utf-8")
        self._log_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logging.getLogger().addHandler(self._log_handler)

    @property
    def docker(self):
        if self._docker is None:
            self._docker = AutoDockVina(project_root=self.project_root)
        return self._docker

    def close(self):
        if getattr(self, "_closed", False):
            return
        self._closed = True
        self.database.close()
        logging.getLogger().removeHandler(self._log_handler)
        self._log_handler.close()

    def show_menu(self):
        """Display the main menu."""
        pool = project_budget(self.project_root)
        console.print("\n[bold cyan]FragEvoLab - main menu[/bold cyan]")
        console.print(f"[dim]CPU: budget {pool.total}, maximum per Vina job {pool.per_job}[/dim]")
        groups = [
            ("Calculations", [(21, "Evolutionary generation"), (4, "Docking")]),
            ("Ligands", [(1, "SMILES generation"), (2, "SMILES → PDBQT"), (3, "Ligand list")]),
            ("Results", [(9, "HTML/CSV reports"), (10, "Statistics"), (11, "Export ZIP"),
                             (12, "Experiment history"), (13, "Database report")]),
            ("Settings and tools", [(26, "Readiness check"), (25, "CPU budget"), (14, "Cache statistics"),
                 (15, "Clear cache"), (16, "Project status")])]
        for title, items in groups:
            console.print(f"[bold]{title}[/bold]")
            for number, label in items:
                console.print(f"  {number:>2}  {label}")
        console.print("   0  Exit")

    def check_readiness(self):
        from src.core.preflight import check_project, display
        display(check_project(self.project_root), console)

    def configure_cpu(self):
        pool=project_budget(self.project_root)
        try:
            total=int(console.input(f"Total logical CPU budget (1-{available_cpus()}, currently {pool.total}): ").strip() or pool.total)
            per_job=int(console.input(f"Maximum threads per Vina job (currently {pool.per_job}): ").strip() or pool.per_job)
            pool=configure_project_cpu(self.project_root,total,per_job)
            print_success(f"Saved: total CPU budget {pool.total}, maximum per Vina job {pool.per_job}")
        except (ValueError,RuntimeError) as error:
            print_error(str(error))

    def generate_smiles(self):
        """Start template-based SMILES generation."""
        print_header("📝 SMILES library generation", "Generate molecules from templates")
        self.smiles_gen.interactive_generation()
    
    def convert_smiles(self):
        """Start SMILES-to-PDBQT conversion."""
        print_header("🔄 SMILES → PDBQT conversion", "Prepare ligands for docking")
        
        ligands_csv = self.project_root / "data" / "ligands.csv"
        if not ligands_csv.exists():
            print_error("File data/ligands.csv not found!")
            console.print("[dim]  Create a ligand file or generate SMILES first (option 1)[/dim]")
            return
        
        force_option = console.input("\n[cyan]Overwrite existing PDBQT files? (yes/no, default no):[/cyan] ").strip().lower()
        force = force_option in ['да', 'yes', 'y', 'д']
        
        stats=self.converter.convert_from_csv(force=force)
        if stats.get("input_error"):
            print_error(f"CSV error: {stats['input_error']}")
            return
        console.print(f"Prepared: {stats['success']}; cache hits: {stats['cached']}; errors: {stats['failed']}")
        for name,reason in stats['errors']:
            print_error(f"{name}: {reason}")
        if stats['failed']==0:
            print_success("Conversion complete!")
    
    def view_ligands(self):
        """Display the ligand list."""
        print_header("📋 Ligand list", "All molecules in the project")
        
        ligands = self.utils.load_ligands_from_csv()
        
        if not ligands:
            print_error("No ligands found in data/ligands.csv")
            return
        
        console.print(f"\n[bold]Total ligands:[/bold] [cyan]{len(ligands)}[/cyan]")
        console.print("\n[dim]First 20 ligands:[/dim]")
        
        from rich.table import Table
        
        table = Table(show_header=True, header_style="bold cyan")
        table.add_column("№", justify="right", style="dim")
        table.add_column("Name", style="cyan")
        table.add_column("SMILES", style="white")
        
        for i, lig in enumerate(ligands[:20], 1):
            smiles_short = lig['smiles'][:50] + "..." if len(lig['smiles']) > 50 else lig['smiles']
            table.add_row(str(i), lig['name'], smiles_short)
        
        console.print(table)
        
        if len(ligands) > 20:
            console.print(f"\n[dim]... and another {len(ligands) - 20} ligands[/dim]")

    def docking(self):
        from src.ui.docking_menu import docking_menu
        docking_menu(self, console)

    def generate_reports(self, exhaustiveness: int = None):
        """Generate HTML and CSV reports."""
        print_header("📄 Report generation", "HTML and CSV results")
        
        results_json = self.project_root / "output" / "cache" / "docking_results.json"
        
        if not results_json.exists():
            print_error("Docking results not found!")
            console.print("[dim]  Run docking first (option 4)[/dim]")
            return
        
        import json
        with open(results_json, 'r', encoding='utf-8') as f:
            docking_results = json.load(f)
        
        if not docking_results:
            print_error("Docking results are empty")
            return
        
        console.print(f"\n[dim]Results found for {len(docking_results)} ligands[/dim]")
        console.print("[dim]Estimating Ki/pKi from docking scores at 298.15 K (not experimental values)...[/dim]")
        
        enriched_results = self.calculator.process_docking_results(docking_results)
        
        title = console.input("\n[cyan]Report title (Enter = default):[/cyan] ").strip()
        if not title:
            title = "Docking results"
        
        html_path, csv_path = self.reporter.generate_reports(
            enriched_results, title=title, exhaustiveness=exhaustiveness
        )
        
        print_success("Reports created:")
        console.print(f"  📄 [cyan]HTML:[/cyan] {html_path}")
        console.print(f"  📊 [cyan]CSV:[/cyan] {csv_path}")
    
    def view_statistics(self):
        """Display result statistics."""
        print_header("📊 Result statistics", "Docking data analysis")
        
        results_json = self.project_root / "output" / "cache" / "docking_results.json"
        
        if not results_json.exists():
            print_error("Docking results not found")
            return
        
        import json
        with open(results_json, 'r', encoding='utf-8') as f:
            docking_results = json.load(f)
        
        enriched_results = self.calculator.process_docking_results(docking_results)
        
        proteins = set()
        for ligand_data in enriched_results.values():
            proteins.update(ligand_data.keys())
        
        proteins = sorted(proteins)
        
        console.print(f"\n[bold]Ligands:[/bold] [cyan]{len(enriched_results)}[/cyan]")
        console.print(f"[bold]Receptors:[/bold] [cyan]{len(proteins)}[/cyan]")
        console.print(f"[bold]Total docking calculations:[/bold] [yellow]{sum(len(values) for values in enriched_results.values())}[/yellow]")
        
        console.print(f"Estimated Ki/pKi from docking score at {self.calculator.temperature:g} K; not experimental affinity.")
        from rich.table import Table
        
        for protein in proteins:
            stats = self.calculator.calculate_statistics(enriched_results, protein)
            
            table = Table(title=f"[bold cyan]{protein}[/bold cyan]", show_header=True)
            table.add_column("Parameter", style="cyan")
            table.add_column("Value", style="white")
            
            if stats['delta_g']:
                table.add_row("ΔG (min)", f"{stats['delta_g']['min']:.2f} kcal/mol")
                table.add_row("ΔG (max)", f"{stats['delta_g']['max']:.2f} kcal/mol")
                table.add_row("ΔG (mean)", f"{stats['delta_g']['mean']:.2f} kcal/mol")
            
            if stats['ki']:
                table.add_row("Ki (min)", stats['ki']['min_formatted'])
                table.add_row("Ki (max)", stats['ki']['max_formatted'])
            
            if stats['pki']:
                table.add_row("pKi (min)", f"{stats['pki']['min']:.2f}")
                table.add_row("pKi (max)", f"{stats['pki']['max']:.2f}")
                table.add_row("pKi (mean)", f"{stats['pki']['mean']:.2f}")
            
            console.print(table)
            console.print()
            
            top = self.calculator.get_best_ligands_by_pki(enriched_results, protein, n=5)
            console.print(f"[bold]Top 5 ligands for {protein}:[/bold]")
            for i, (name, dg, pki) in enumerate(top, 1):
                console.print(f"  {i}. [cyan]{name}[/cyan]: ΔG={dg:.2f}, pKi={pki:.2f}")
            console.print()

    def export_to_zip(self):
        """Export results to ZIP."""
        print_header("📦 Export to ZIP", "Archive results")
        
        include_logs = console.input("\n[cyan]Include logs and JSON? (yes/no):[/cyan] ").strip().lower() in ['да', 'yes', 'y', 'д']
        
        output_name = console.input("[cyan]ZIP filename (Enter = automatic):[/cyan] ").strip()
        output_file = Path(output_name) if output_name else None
        
        try:
            zip_path = self.utils.export_results_to_zip(
                output_file=output_file,
                include_logs=include_logs
            )
            
            print_success("ZIP archive created:")
            console.print(f"  📦 [cyan]{zip_path}[/cyan]")
            console.print(f"  💾 Size: [yellow]{zip_path.stat().st_size / (1024*1024):.2f} MB[/yellow]")
        except Exception as e:
            print_error(f"Failed to create archive: {e}")
    
    def view_db_experiments(self):
        """Display experiment history from the database."""
        print_header("💾 Experiment history", "All runs stored in the database")
        
        experiments = self.database.get_experiment_list()
        
        if not experiments:
            print_info("No saved experiments")
            return
        
        from rich.table import Table
        
        table = Table(show_header=True, header_style="bold cyan")
        table.add_column("ID", justify="right")
        table.add_column("Name")
        table.add_column("Start")
        table.add_column("End")
        table.add_column("Status")
        
        for exp in experiments:
            status = "[green]✓[/green]" if exp['finished_at'] else "[yellow]⏳[/yellow]"
            table.add_row(
                str(exp['id']),
                exp['name'],
                exp['started_at'][:16] if exp['started_at'] else 'N/A',
                exp['finished_at'][:16] if exp['finished_at'] else 'N/A',
                status
            )
        
        console.print(table)
        console.print(f"\n[bold]Total experiments:[/bold] [cyan]{len(experiments)}[/cyan]")
    
    def generate_db_report(self):
        """Generate db-report.html from the database."""
        print_header("📊 Database report", "History of all docking calculations")
        
        experiments = self.database.get_experiment_list()
        
        if not experiments:
            print_error("No data in the database")
            return
        
        console.print(f"\n[dim]Experiments in the database: {len(experiments)}[/dim]")
        console.print("[dim]Generating report...[/dim]")
        
        db_data = self.database.export_for_report()
        
        report_path = self.reporter.generate_db_report(db_data)
        
        print_success("Database report created:")
        console.print(f"  📄 [cyan]{report_path}[/cyan]")
    
    def cache_statistics(self):
        """Cache statistics."""
        print_header("📊 Cache statistics", "Disk usage")
        self.cache_mgr.print_statistics()
    
    def clear_cache(self):
        """Clear cache."""
        print_header("🗑️  Clear cache", "Remove reusable calculation entries")
        
        console.print("\n[bold cyan]Select what to clear:[/bold cyan]")
        console.print("  1. All preparation and docking cache entries")
        console.print("  2. Preparation cache entries only")
        console.print("  3. Docking cache entries only")
        console.print("Saved poses, reports, and checkpoints are retained.")
        console.print("  0. Cancel")
        
        choice = console.input("\n[cyan]Your choice:[/cyan] ").strip()
        
        if choice == '1':
            self.cache_mgr.clear_cache(confirm=True)
        elif choice == '2':
            self.cache_mgr.clear_ligands_cache(confirm=True)
        elif choice == '3':
            self.cache_mgr.clear_docking_cache(confirm=True)
        elif choice == '0':
            print_info("Cancel")
        else:
            print_error("Invalid choice")
    
    def project_status(self):
        """Project status."""
        print_header("ℹ️  Project status", "File and data overview")
        self.utils.print_project_status()
    
    def evolutionary_denovo(self):
        from src.ui.evolution_menu import evolution_menu
        evolution_menu(self, console)

    def _save_results_to_db(self, exp_id: int, results: Dict, exhaustiveness: int):
        """Save docking results to the database."""
        smiles_dict = self.reporter.load_smiles_data()
        
        for ligand_name, proteins_data in results.items():
            smiles = next((row['smiles'] for row in proteins_data.values() if row.get('smiles')), smiles_dict.get(ligand_name, ''))
            
            descriptors = self.reporter.calculate_molecular_descriptors(smiles)
            
            ligand_dict = {
                'name': ligand_name,
                'smiles': smiles,
                'formula': descriptors.get('formula'),
                'mw': descriptors.get('mw'),
                'logp': descriptors.get('logp'),
                'tpsa': descriptors.get('tpsa'),
                'hbd': descriptors.get('hbd'),
                'hba': descriptors.get('hba'),
                'sa_score': descriptors.get('sa_score')
            }
            ligand_id = self.database.add_or_get_ligand(ligand_dict)
            
            for protein_name, data in proteins_data.items():
                protein_id = self.database.add_or_get_protein({
                    'name': protein_name,
                    'pdbqt_file': data.get('input_files', {}).get('receptor', f"{protein_name}.pdbqt")
                })
                
                delta_g = data.get('affinity')
                ki = self.calculator.calculate_ki(delta_g) if delta_g is not None else None
                pki = self.calculator.calculate_pki_from_score(delta_g) if delta_g is not None else None
                
                result_dict = {
                    **data,
                    'delta_g': delta_g,
                    'ki': ki,
                    'pki': pki,
                    'rank': data.get('mode', 1),
                    'method': data.get('spec', {}).get('options', {}).get('scoring', 'vina'),
                    'exhaustiveness': exhaustiveness
                }
                
                self.database.add_docking_result(exp_id, ligand_id, protein_id, result_dict)
    

    def run(self):
        """Main application loop."""
        print_app_header()
        
        while True:
            self.show_menu()
            
            choice = console.input("\n[bold yellow]Your choice:[/bold yellow] ").strip()
            
            try:
                if choice == '1':
                    self.generate_smiles()
                elif choice == '2':
                    self.convert_smiles()
                elif choice == '3':
                    self.view_ligands()
                elif choice == '4':
                    self.docking()
                elif choice == '9':
                    self.generate_reports()
                elif choice == '10':
                    self.view_statistics()
                elif choice == '11':
                    self.export_to_zip()
                elif choice == '12':
                    self.view_db_experiments()
                elif choice == '13':
                    self.generate_db_report()
                elif choice == '14':
                    self.cache_statistics()
                elif choice == '15':
                    self.clear_cache()
                elif choice == '16':
                    self.project_status()
                elif choice == '21':
                    self.evolutionary_denovo()
                elif choice == '26':
                    self.check_readiness()
                elif choice == '25':
                    self.configure_cpu()
                elif choice == '0':
                    console.print(f"[bold cyan]Thank you for using {APP_NAME} {__version__}![/bold cyan]")
                    self.close()
                    break
                else:
                    print_error("Invalid choice. Please try again.")
            
            except KeyboardInterrupt:
                console.print("\n\n[yellow]⚠ Operation interrupted by the user[/yellow]")
                continue
            
            except Exception as e:
                self.logger.error(f"Error: {e}", exc_info=True)
                print_error(f"An error occurred: {e}")
                console.print("[dim]  See the log for details: output/cache/fragevolab.log[/dim]")

