# Testing FragEvoLab

## 1. Developer local suite (no installed Vina)

End users do not need pytest. Follow the README installation instructions, supply prepared receptors and box configurations, and run a small docking/evolution example through the menu. The checks below are for development and troubleshooting.

Run from the project root after installing Python 3.12 and the dependencies:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt -r requirements-test.txt
.\.venv\Scripts\python.exe -m pytest -q
```

Explicit equivalent:

```powershell
.\.venv\Scripts\python.exe -m pytest -q -m fast
```

On Linux use `./.venv/bin/python` instead of `.\.venv\Scripts\python.exe`.

The local suite does not execute docking and requires no installed Vina. RDKit/Meeko and the other Python dependencies are required. Tests cover molecular validation and mutation rules, real ligand preparation, cache utilities, CPU allocation, configuration, reports and bundled-data consistency. Small local fixtures and isolated unit-test mocks are used; no Vina executable is generated or simulated. These checks do not establish docking accuracy or end-to-end behavior with Vina.

`pytest.ini` restricts default collection to `tests/`. The real docking entry point is the separate `verify_integration.py` command below; ordinary pytest collection does not start it. Local tests are marked `fast` in `tests/conftest.py`.

## 2. Preparation of a supplied ligand dataset (no Vina)

This opt-in command prepares every row in a supplied CSV and verifies reuse of its preparation cache:

```powershell
.\.venv\Scripts\python.exe verify_ligands.py --csv data/ligands.csv --output checks/ligand_preparation.json
```

Exit status is 0 only when all rows succeed; 1 indicates preparation/cache failures, and 2 indicates invalid command arguments or input. Results include per-ligand errors. The cache is temporary and does not alter the application's cache. Preparation statistics are not docking results.

## 3. Real-Vina integration verification (opt-in)

Required inputs:

- An installed executable for the operating system, verified using `--version`.
- A data directory with `ligands.csv`, prepared `proteins/<name>.pdbqt` files, and matching `configs/<name>.txt` (or `.conf`/`.config`) files.
- An explicit evolution JSON configuration whose target/off-target names match these filenames.
- `generations >= 1`, `force_preparation=false`, and `force_docking=false`.

Use a small population and low exhaustiveness for an initial technical smoke run. Ensure constraints and thresholds permit feasible seeds and offspring; a run that cannot progress beyond G0 fails the resumption check. The real integration command does not silently reduce your configured workload.

Windows:

```powershell
.\.venv\Scripts\python.exe verify_integration.py --vina tools/vina.exe --data-dir data --config evolution_settings.json --cpu-budget 2 --timeout 3600 --output checks/real_integration.json
```

Linux:

```bash
./.venv/bin/python verify_integration.py --vina tools/vina --data-dir data --config evolution_settings.json --cpu-budget 2 --timeout 3600 --output checks/real_integration.json
```

Relative paths are resolved from the repository root. All three of `--vina`, `--data-dir`, and `--config` are required. The bundled example contains prepared JAK1 3EYG and JAK2 6BBV receptors. See [test_data.md](test_data.md) for provenance.

The script copies only selected inputs into a temporary project, then checks:

1. Preflight readiness.
2. G0 checkpoint creation and continuation into at least one subsequent generation.
3. Repeating the same run in the same environment with no new docking calculations and an identical archive.
4. Complete docking results for archived candidates.
5. The allocated Vina CPU budget.
6. Offline reweighting without docking.

Archive identity here tests reuse within the same process/environment; it does not claim bitwise reproducibility between operating systems or independent Vina runs. Reweighting currently checks the reweighting path and candidate count using the supplied configuration; scientific sensitivity to alternative weights requires a separate experiment.

Outputs are a summary JSON with Boolean checks, candidate records, scores, and version information, plus adjacent HTML and full JSON reports. Exit code 0 means all checks passed, 1 means execution or verification failed, and 2 means invalid arguments/inputs. Runtime failures are written to the summary JSON. Argument validation errors are printed before calculation starts.

The temporary project is deleted after execution. Report artifact paths referring to that temporary project are not persistent pose/checkpoint exports. For a scientific run whose poses and checkpoint must be retained, use the regular evolution CLI described in the README.

## Scope and versions

The application version is defined in `src/version.py`. Use `python main.py --version` to display it without loading scientific dependencies. Dependency and Vina versions are separate from the application version.

The example contains 40 ligands. Its configured population size is six, so the integration check does not dock all 40 inputs. Passing technical tests does not establish receptor preparation quality, docking accuracy, biological activity or selectivity.
