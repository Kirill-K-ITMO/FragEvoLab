"""Readiness checks. No preparation, docking or evolutionary state is created."""
from dataclasses import dataclass, field
from pathlib import Path
import importlib
import json
import sys
import tempfile


@dataclass
class Readiness:
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    summary: dict = field(default_factory=dict)

    @property
    def ready(self):
        return not self.errors

    def require(self):
        if self.errors:
            raise ValueError("Calculation is not ready:\n" + "\n".join(self.errors))

    def to_dict(self):
        return dict(ready=self.ready, errors=self.errors, warnings=self.warnings, summary=self.summary)


def check_dependencies():
    report = Readiness()
    report.summary['Python'] = sys.version.split()[0]
    if sys.version_info[:2] != (3, 12):
        report.warnings.append('This build is tested with Python 3.12; the current version differs.')
    for name in ('rdkit.Chem', 'rdkit.Chem.AllChem', 'meeko', 'gemmi', 'numpy', 'scipy', 'rich', 'filelock', 'jinja2'):
        try:
            importlib.import_module(name)
        except Exception as error:
            report.errors.append(f'{name}: {error}. Install requirements.txt in the active environment.')
    return report


def _probe_directory(path):
    path.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryFile(dir=path) as stream:
        stream.write(b'fragevolab-readiness')
        stream.flush()


def check_run(root, *, config=None, seeds=None, resume=None, proteins=None,
              rows=None, exhaustiveness=8, num_workers=4, docker=None,
              dependencies=True):
    """Check selected inputs and effective Vina options, using production validators.

    Invalid evolutionary seeds are reported as excluded, as in the engine. If all
    seeds fail, readiness fails. For ordinary docking every selected SMILES must
    be valid. A checkpoint is checked against the exact same fingerprint as resume.
    """
    root = Path(root).resolve()
    report = check_dependencies() if dependencies else Readiness()
    report.summary['Project'] = str(root)
    if report.errors:
        return report
    if config is not None:
        try:
            config.validate()
            proteins = config.proteins
            exhaustiveness, num_workers = config.exhaustiveness, config.num_workers
            report.summary.update({'Targets': ', '.join(p.name for p in config.targets),
                'Off-targets': ', '.join(p.name for p in config.off_targets) or 'none',
                'Population': config.population_size, 'Generations after G0': config.generations,
                'Mutation levels': config.level_probabilities, 'Atom-level operators': config.atom_probabilities,
                'Fragment / operator adaptation': f'{config.adaptive_fragments} / {config.adaptive_operators}',
                'Force preparation / docking': f'{config.force_preparation} / {config.force_docking}'})
        except Exception as error:
            report.errors.append(f'Evolution settings: {error}')
            return report
    from src.core.cpu_budget import positive
    try:
        positive(num_workers, 'num_workers'); positive(exhaustiveness, 'exhaustiveness')
    except ValueError as error:
        report.errors.append(str(error))
    for path in (root/'output', root/'.calculation_cache'/'v1'):
        try:
            _probe_directory(path)
        except OSError as error:
            report.errors.append(f'Cannot write to {path}: {error}')
    if docker is None:
        try:
            from src.docking.vina import AutoDockVina
            docker = AutoDockVina(root, seed=config.seed if config else 42)
        except Exception as error:
            report.errors.append(f'Vina/CPU: {error}. Check app_settings.json, cpu_settings.json, and the Vina executable.')
    if docker is not None:
        report.summary.update({'Vina': docker.vina_version, 'Exhaustiveness': exhaustiveness,
            'CPU budget': docker.cpu_pool.total, 'Maximum CPU threads per Vina job': docker.cpu,
            'Concurrent jobs (upper limit)': min(num_workers, docker.cpu_pool.total)})
        from src.core.calculation_cache import validate_pdbqt
        if not proteins:
            report.errors.append('No receptors selected.')
        for name in proteins or []:
            try:
                if not isinstance(name,str) or Path(name).name != name or any(c in name for c in '/\\'):
                    raise ValueError('A receptor name must be a filename stem without a path or extension')
                receptor = root/'data'/'proteins'/f'{name}.pdbqt'
                validate_pdbqt(receptor.read_bytes())
                conf = docker.get_config(name)
                if conf is None:
                    raise ValueError('No .txt, .conf, or .config configuration file in data/configs')
                docker._effective_options(conf, exhaustiveness)
            except Exception as error:
                report.errors.append(f'Receptor {name}: {error}')
    if resume is not None:
        try:
            if config is None or docker is None:
                raise ValueError('Valid settings and Vina are required to validate a checkpoint')
            from src.evolution import state
            data = state.load(resume, state.fingerprint(config, docker))
            if config.generations < data['generation']:
                raise ValueError('The generation limit is below the saved generation')
            if not data['population'] or any(s not in data['archive'] for s in data['population']):
                raise ValueError('The checkpoint population list is corrupted')
            report.summary['Resume after generation'] = data['generation']
            report.summary['Results'] = str(Path(resume).resolve().parent)
            _probe_directory(Path(resume).resolve().parent)
        except Exception as error:
            report.errors.append(f'Checkpoint: {error}')
    else:
        try:
            if rows is None:
                from src.preparation.inputs import read_ligand_rows
                rows = read_ligand_rows(Path(seeds) if seeds else root/'data'/'ligands.csv')
            if not rows:
                raise ValueError('The CSV or molecule selection is empty')
            from src.core.calculation_cache import molecule_identity
            from src.evolution.validation import MolecularConstraints
            constraints = MolecularConstraints(config.max_mw,config.max_logp,config.min_nitrogen) if config else None
            seen, excluded = set(), []
            for name, smiles in rows:
                try:
                    if constraints:
                        canonical, reason = constraints.validate(smiles)
                        if not canonical: raise ValueError(reason)
                    else:
                        from src.preparation.validation import validated_molecule
                        validated_molecule(smiles)
                        canonical = molecule_identity(smiles)[0]
                    if canonical in seen:
                        excluded.append(f'{name}: duplicate')
                    seen.add(canonical)
                except Exception as error:
                    message=f'{name}: {error}'
                    if config: excluded.append(message)
                    else: report.errors.append(message)
            report.summary['Input rows / unique valid molecules'] = f'{len(rows)} / {len(seen)}'
            if not seen:
                report.errors.append('No valid seeds satisfy the molecular constraints: no suitable molecules.')
            if excluded:
                report.warnings.append(f'Excluded/duplicate seed structures: {len(excluded)}. ' + '; '.join(excluded[:5]))
            if config and len(seen) < config.population_size:
                report.warnings.append('G0 will be smaller than population_size: insufficient unique seed structures.')
            if config and len(seen) > config.population_size:
                report.warnings.append('A subset of population_size seed molecules will be selected using the specified seed.')
            if not config:
                report.summary['Molecule-receptor pairs'] = len(rows)*len(proteins or [])
        except Exception as error:
            report.errors.append(f'Seed molecules: {error}')
        report.summary['Results'] = str(root/'output'/('evolutionary' if config else 'runs'))
    return report


def check_project(root, *, config_path=None, seeds=None, resume=None):
    root = Path(root).resolve()
    report = check_dependencies()
    if report.errors:
        return report
    try:
        from src.evolution.config import EvolutionConfig
        path = Path(config_path) if config_path else root/'evolution_settings.json'
        config = EvolutionConfig.from_dict(json.loads(path.read_text(encoding='utf-8-sig')))
    except Exception as error:
        report.errors.append(f'Evolution settings: {error}. Check evolution_settings.json or --config.')
        return report
    return check_run(root, config=config, seeds=seeds, resume=resume, dependencies=False)


def display(report, console=None):
    def emit(message):
        if console is None: print(message)
        else: console.print(message, markup=False)
    emit('FragEvoLab readiness check')
    for key, value in report.summary.items(): emit(f'{key}: {value}')
    for message in report.warnings: emit(f'Note: {message}')
    for message in report.errors: emit(f'Error: {message}')
    emit('Ready to run.' if report.ready else 'Calculation blocked. Resolve the errors and repeat the check.')
