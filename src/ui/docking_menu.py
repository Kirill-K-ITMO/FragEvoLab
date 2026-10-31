"""One ordinary docking workflow; selections never modify the user's CSV."""
from src.preparation.inputs import read_ligand_rows
from src.core.calculation_cache import molecule_identity

def parse_selection(text, size):
    text = text.strip().lower()
    if text == "all":
        return list(range(size))
    tokens = text.replace(",", " ").split()
    if not tokens:
        raise ValueError("Enter numbers or all")
    indices = []
    for token in tokens:
        value = int(token) - 1
        if not 0 <= value < size:
            raise ValueError(f"The number must be between 1 and {size}")
        if value not in indices:
            indices.append(value)
    return indices

def docking_menu(app, console):
    proteins = sorted(p.stem for p in app.docker.get_proteins())
    if not proteins:
        raise ValueError("No receptors found in data/proteins")
    for i, protein in enumerate(proteins, 1):
        console.print(f"{i}. {protein}", markup=False)
    selected = [proteins[i] for i in parse_selection(console.input("Receptors: numbers or all: "), len(proteins))]
    for name in selected:
        if app.docker.get_config(name) is None:
            raise ValueError(f"No configuration file for {name}")
    path = app.project_root / "data" / "ligands.csv"
    rows = read_ligand_rows(path) if path.exists() else []
    for i, (name, smiles) in enumerate(rows, 1):
        console.print(f"{i}. {name}: {smiles}", markup=False, soft_wrap=True)
    choice = console.input("SMILES: numbers, all, or write smiles: ").strip()
    if choice.lower() == "write smiles":
        rows = []
        console.print("Enter one SMILES per line; an empty line finishes input. 0 = cancel.")
        while True:
            value = console.input(f"SMILES {len(rows)+1}: ").strip()
            if value == "0":
                return
            if not value:
                break
            # Validate immediately, before any expensive calculation.
            canonical, _ = molecule_identity(value)
            rows.append((f"manual_{len(rows)+1}", canonical))
    else:
        rows = [rows[i] for i in parse_selection(choice, len(rows))]
    if not rows:
        raise ValueError("No molecules selected")
    checked = [(name, molecule_identity(smi)[0]) for name, smi in rows]
    exhaustiveness = int(console.input("Exhaustiveness [8]: ").strip() or 8)
    workers = int(console.input("Concurrent jobs [4], within the CPU budget: ").strip() or 4)
    if min(exhaustiveness, workers) < 1:
        raise ValueError("Exhaustiveness and the job count must be positive")
    force = console.input("Overwrite docking results? (yes/no): ").strip().lower() in ("да", "yes", "y", "д")
    from src.core.preflight import check_run, display
    readiness = check_run(app.project_root, proteins=selected, rows=checked,
                          exhaustiveness=exhaustiveness, num_workers=workers, docker=app.docker)
    display(readiness, console)
    if not readiness.ready:
        return
    console.print(f"Molecules: {len(checked)}; receptors: {len(selected)}; pairs: {len(checked)*len(selected)}")
    if console.input("Start docking? (yes/no): ").strip().lower() not in ("да", "yes", "y", "д"):
        return
    from src.ui.progress import WorkProgress
    prepared = {}
    selected_smiles = {}
    selected_metadata = {}
    failures = []
    with WorkProgress("Ligand preparation", len(checked)) as progress:
        for i, (name, smiles) in enumerate(checked, 1):
            try:
                record = app.converter.prepare_cached(smiles)
                # Row prefix keeps duplicate names independent and deterministic.
                label = name
                suffix = 0
                while label in prepared:
                    suffix += 1
                    label = f"{i:04d}_{suffix}_{name}"
                prepared[label] = record["path"]
                selected_smiles[label] = smiles
                selected_metadata[label] = {"input_name": name, "selection_row_id": i}
            except Exception as error:
                failures.append((name, str(error)))
            finally:
                progress.advance()
    for name, error in failures:
        console.print(f"Preparation failed for {name}: {error}", markup=False)
    if not prepared:
        raise ValueError("No prepared ligands")
    exp_id = app.database.add_experiment(name=f"Docking {len(prepared)}x{len(selected)}",
        user_params={"exhaustiveness":exhaustiveness,"workers":workers,"proteins":selected,
                     "molecules":checked}, description="Unified docking workflow")
    results = app.docker.run_docking(proteins=selected, ligand_paths=prepared,
        exhaustiveness=exhaustiveness, num_workers=workers, force=force)
    for name, proteins_data in results.items():
        for row in proteins_data.values():
            row["smiles"] = selected_smiles[name]
            row.update(selected_metadata[name])
    app.docker._save_results_json()
    app._save_results_to_db(exp_id, results, exhaustiveness)
    app.database.finish_experiment(exp_id)
    if console.input("Generate reports? (yes/no): ").strip().lower() in ("да", "yes", "y", "д"):
        app.generate_reports(exhaustiveness=exhaustiveness)
