"""Opt-in real-Vina verification using explicitly supplied data and configuration.

Uses a temporary project. Does not modify input data or the application's cache.
This checks software behavior, not biological activity or scientific accuracy.
"""
import argparse
import json
import logging
import platform
import shutil
import tempfile
from contextlib import ExitStack
from pathlib import Path
from src.core.calculation_cache import atomic_json
from src.core.cpu_budget import configure_project_cpu
from src.core.preflight import check_run
from src.docking.affinity import AffinityCalculator
from src.docking.vina import AutoDockVina
from src.evolution.config import EvolutionConfig
from src.evolution.engine import EvolutionaryGenerator
from src.evolution.reweight import reweight
from src.version import __version__


def validate_inputs(data_dir, config):
    """Check required files without starting Vina or preparing any molecules."""
    if config.generations < 1:
        raise ValueError('Integration verification requires generations >= 1')
    if config.force_preparation or config.force_docking:
        raise ValueError('Set force_preparation and force_docking to false for cache verification')
    if not (data_dir / 'ligands.csv').is_file():
        raise ValueError(f'Missing ligand CSV: {data_dir / "ligands.csv"}')
    files = [Path('ligands.csv')]
    for receptor in [*config.targets, *config.off_targets]:
        name = receptor.name
        if not name or Path(name).name != name or '/' in name or '\\' in name or name in ('.', '..'):
            raise ValueError('Receptor names must be filename stems, not paths')
        protein = Path('proteins') / (name + '.pdbqt')
        if not (data_dir / protein).is_file():
            raise ValueError(f'Missing receptor: {data_dir / protein}')
        box = next((Path('configs') / (name + suffix) for suffix in ('.txt', '.conf', '.config')
                    if (data_dir / 'configs' / (name + suffix)).is_file()), None)
        if box is None:
            raise ValueError(f'Missing box configuration for {name}')
        files.extend([protein, box])
    return files


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vina', required=True, help='Vina executable path or command on PATH')
    parser.add_argument('--data-dir', required=True, type=Path,
                        help='Directory containing ligands.csv, proteins/, and configs/')
    parser.add_argument('--config', required=True, type=Path, help='Evolution JSON with matching receptor names')
    parser.add_argument('--output', type=Path, default=Path('integration_results.json'))
    parser.add_argument('--cpu-budget', type=int, default=2)
    parser.add_argument('--timeout', type=float, default=3600)
    args = parser.parse_args(argv)
    try:
        if args.cpu_budget < 1 or args.timeout <= 0:
            raise ValueError('CPU budget and timeout must be positive')
        data_dir = args.data_dir.resolve()
        config = EvolutionConfig.from_dict(json.loads(args.config.read_text(encoding='utf-8-sig')))
        files = validate_inputs(data_dir, config)
        executable = shutil.which(args.vina)
        if executable is None:
            raise ValueError('AutoDock Vina executable not found; supply --vina')
        executable = str(Path(executable).resolve())
        output = args.output.resolve()
        protected = [args.config.resolve(), *((data_dir / file).resolve() for file in files)]
        destinations = [output, output.with_suffix('.html'), output.with_name(output.stem + '_report.json')]
        if output.suffix.lower() != '.json' or any(path in protected for path in destinations):
            raise ValueError('Use a separate .json output path that does not overwrite an input file')
    except (ValueError, TypeError, OSError) as error:
        parser.error(str(error))

    result = {'fragevolab_version': __version__, 'platform': platform.platform(),
              'config': config.to_dict(), 'status': 'failed', 'checks': {},
              'note': 'Technical verification only; not a validation of biological activity or selectivity.'}
    try:
        with tempfile.TemporaryDirectory(prefix='fragevolab-integration-') as folder, ExitStack() as cleanup:
            cleanup.callback(logging.shutdown)
            root = Path(folder)
            for file in files:
                target = root / 'data' / file
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(data_dir / file, target)
            configure_project_cpu(root, args.cpu_budget, 1)
            docker = AutoDockVina(root, vina_executable=executable, seed=config.seed, timeout=args.timeout)
            result['vina'] = docker.vina_version
            check_run(root, config=config, seeds=root / 'data/ligands.csv', docker=docker).require()
            calculator = AffinityCalculator()

            def run(**kwargs):
                return EvolutionaryGenerator(root).run_evolution(
                    root / 'data/ligands.csv', docker, calculator, config=config, **kwargs)

            paused = run(stop_after_generation=0)
            first = run(resume_from=paused['checkpoint'])
            repeat = run()
            first_state = json.loads(Path(first['checkpoint']).read_text(encoding='utf-8'))
            repeat_state = json.loads(Path(repeat['checkpoint']).read_text(encoding='utf-8'))
            rescored = reweight(first['checkpoint'], config, root / 'rescored.json')
            checks = {
                'resume_after_G0': first['total_generations'] >= 2,
                'identical_repeat': first_state['archive'] == repeat_state['archive'],
                'repeat_uses_cache': repeat['total_dockings'] == 0 and repeat['cache_hits'] > 0,
                'new_docking_performed': first['total_dockings'] > 0,
                'cpu_budget_respected': docker.cpu_pool.peak <= args.cpu_budget,
                'all_candidates_complete': bool(first_state['archive']) and all(
                    candidate['score']['complete'] for candidate in first_state['archive'].values()),
                'reweight': len(rescored['candidates']) == len(first_state['archive']),
            }
            result.update(checks=checks, status='passed' if all(checks.values()) else 'failed',
                          repeat_new_dockings=repeat['total_dockings'], repeat_cache_hits=repeat['cache_hits'],
                          first_new_dockings=first['total_dockings'], cpu_peak_allocated=docker.cpu_pool.peak,
                          candidates=list(first_state['archive'].values()), generations=first['generations'],
                          fragments=first['fragments'])
            output.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(first['html_report'], output.with_suffix('.html'))
            shutil.copyfile(first['json_report'], output.with_name(output.stem + '_report.json'))
    except Exception as error:
        result['error'] = str(error)
    atomic_json(output, result)
    print(json.dumps({'status': result['status'], 'checks': result['checks'],
                      'error': result.get('error')}, indent=2))
    return 0 if result['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
