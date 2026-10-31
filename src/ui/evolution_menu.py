"""Terminal settings editor and the same workflow used by the CLI."""
import json
from pathlib import Path
from src.core.calculation_cache import atomic_json
from src.evolution.config import EvolutionConfig
from src.evolution.reweight import reweight

HELP = {
 'targets': 'Targets: a list of names or objects {name, weight, threshold, offset, scale}',
 'off_targets': 'Off-targets: a list of objects; each requires a threshold in kcal/mol',
 'penalty_strength': 'Off-target penalty multiplier',
 'target_aggregation': 'mean = weighted mean; worst = worst target',
 'generations': 'Number of generations after the initial G0',
 'mutation_fraction': 'Fraction of non-elite slots assigned to mutation, from 0 to 1',
 'level_probabilities': 'Probabilities of atom and fragment; sum = 1',
 'atom_probabilities': 'Probabilities of add_atom, remove_atom, substitute_atom, change_bond; sum = 1',
 'adaptive_fragments': 'Adaptive fragment weights: true / false',
 'adaptive_operators': 'Adaptive operator probabilities: true / false; zero probabilities remain zero',
 'mutation_steps': 'Number of sequential mutation steps',
 'exploration': 'Uniform exploration component for fragments, from 0 to 1',
}

def resolve_input(root, value):
    path = Path(value.strip().strip('"')).expanduser()
    return path if path.is_absolute() else Path(root)/path

def edit_config(console, config):
    values = config.to_dict()
    keys = list(values)
    console.print('Enter JSON values: numbers, true/false, lists, and objects. Enter keeps the current value.')
    while True:
        groups = {
            "Receptors and selection": ["targets", "off_targets", "penalty_strength", "target_aggregation"],
            "Population and generations": ["population_size", "generations", "parent_count", "elite_fraction", "max_attempts_per_slot"],
            "Mutations and adaptation": ["mutation_fraction", "level_probabilities", "atom_probabilities", "mutation_steps", "adaptive_fragments", "adaptive_operators", "exploration", "prior_strength", "reward_scale", "min_operator_observations"],
            "Molecular constraints": ["max_similarity", "max_mw", "max_logp", "min_nitrogen"]}
        listed = {key for group in groups.values() for key in group}
        groups["Docking and reproducibility"] = [key for key in keys if key not in listed]
        for title, group in groups.items():
            console.print(title, markup=False)
            for key in group:
                if key in values:
                    console.print(f'{keys.index(key)+1}. {key}: {json.dumps(values[key], ensure_ascii=False)}', markup=False)
        choice = console.input('Parameter number; 0 = validate and save: ').strip()
        if choice == '0':
            try: return EvolutionConfig.from_dict(values)
            except (ValueError, TypeError, AttributeError) as error:
                console.print(f'Parameter error: {error}', markup=False)
                continue
        try:
            index = int(choice)-1
            if not 0 <= index < len(keys): raise ValueError('Invalid number')
            key = keys[index]
            console.print(HELP.get(key, key), markup=False)
            value = console.input('New value (JSON): ').strip()
            if value: values[key] = json.loads(value)
        except (ValueError, TypeError) as error:
            console.print(f'Input error: {error}', markup=False)

def evolution_menu(app, console):
    root = app.project_root
    settings = root/'evolution_settings.json'
    while True:
        console.print('Evolution: 1 - settings; 2 - import JSON; 3 - new run; 4 - resume; 5 - rescore; 6 - check readiness; 0 - back')
        choice = console.input('Choice: ').strip()
        if choice == '0': return
        try:
            if choice == '2':
                source = resolve_input(root, console.input('JSON file: '))
                config = EvolutionConfig.from_dict(json.loads(source.read_text(encoding='utf-8-sig')))
                atomic_json(settings, config.to_dict())
                continue
            if settings.exists():
                config = EvolutionConfig.from_dict(json.loads(settings.read_text(encoding='utf-8-sig')))
            else:
                proteins = sorted(p.stem for p in app.docker.proteins_dir.glob('*.pdbqt'))
                console.print('Available receptors: '+', '.join(proteins), markup=False)
                names = [n.strip() for n in console.input('Target receptor names, separated by commas: ').split(',') if n.strip()]
                config = EvolutionConfig.from_dict({'targets': names})
            if choice == '6':
                from src.core.preflight import check_run, display
                display(check_run(root, config=config), console)
            elif choice == '1':
                atomic_json(settings, edit_config(console, config).to_dict())
            elif choice in ('3', '4'):
                checkpoint = resolve_input(root, console.input('checkpoint.json: ')) if choice == '4' else None
                seeds = None if checkpoint else resolve_input(root, console.input('Ligand CSV [data/ligands.csv]: ') or 'data/ligands.csv')
                from src.evolution.engine import EvolutionaryGenerator
                # Same reproducible seed as the CLI.
                old_seed = app.docker.seed
                try:
                    app.docker.seed = config.seed
                    from src.core.preflight import check_run, display
                    readiness = check_run(root, config=config, seeds=seeds, resume=checkpoint, docker=app.docker)
                    display(readiness, console)
                    if not readiness.ready:
                        continue
                    if console.input('Start calculation? (yes/no): ').strip().lower() not in ('да', 'д', 'yes', 'y'):
                        continue
                    report = EvolutionaryGenerator(root).run_evolution(seeds, app.docker, app.calculator, config=config, resume_from=checkpoint)
                finally:
                    app.docker.seed = old_seed
                for key in ('stop_reason', 'html_report', 'json_report', 'checkpoint'):
                    console.print(f'{key}: {report[key]}', markup=False)
            elif choice == '5':
                checkpoint = resolve_input(root, console.input('checkpoint.json: '))
                output = resolve_input(root, console.input('New JSON report file: '))
                reweight(checkpoint, config, output)
                console.print(f'Report: {output}', markup=False)
            else: console.print('Invalid choice')
        except (ValueError, TypeError, AttributeError, OSError, RuntimeError) as error:
            console.print(f'Error: {error}', markup=False)
