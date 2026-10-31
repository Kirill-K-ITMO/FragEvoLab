# JAK1/JAK2 input dataset

The repository includes 40 ligand SMILES in `data/ligands.csv`, prepared JAK1 (3EYG) and JAK2 (6BBV) receptors in `data/proteins/`, and matching docking boxes in `data/configs/`. The root `evolution_settings.json` selects JAK1 as target and JAK2 as off-target. No copying or renaming of the supplied inputs is necessary.

Install Python dependencies and the Vina executable as described in the [main README](../README.md). Run these commands from the repository root. Vina is installed separately.

## Run evolution and retain outputs

Windows PowerShell:

```powershell
.\.venv\Scripts\python.exe -m src.evolution.cli --project . --config evolution_settings.json --seeds data/ligands.csv --vina tools/vina.exe --cpu-budget 2 --cpu 1
```

Linux:

```bash
./.venv/bin/python -m src.evolution.cli --project . --config evolution_settings.json --seeds data/ligands.csv --vina tools/vina --cpu-budget 2 --cpu 1
```

The CLI performs readiness checks and stores outputs under `output/`. To resume, repeat the command with `--resume <checkpoint-path>`, using the checkpoint path printed by the program. Keep the same implementation, settings and inputs for resumption.

For menu-based use, set `vina_executable` in `app_settings.json` to `tools/vina.exe` on Windows or `tools/vina` on Linux, then run `python main.py` with your environment's Python interpreter. Menu item 4 starts docking; item 21 opens evolution. The prepared receptors and ligand table are available directly from `data/`.

The bundled evolution configuration uses population size 6, one offspring generation and exhaustiveness 1. It is a technical example, not a scientific study protocol. It does not dock all 40 input ligands. Use the docking menu and select all ligands to run the entire table.

## Check preparation and integration

Windows PowerShell:

```powershell
.\.venv\Scripts\python.exe verify_ligands.py --csv data/ligands.csv --output output/verification/preparation.json
.\.venv\Scripts\python.exe verify_integration.py --vina tools/vina.exe --data-dir data --config evolution_settings.json --cpu-budget 2 --output output/verification/integration.json
```

On Linux use `./.venv/bin/python` and `--vina tools/vina`.

The preparation command checks all 40 rows and preparation-cache reuse. The integration command tests docking, resumption, cache reuse and offline reweighting using the configured population. It retains summary JSON and HTML reports, but deletes its temporary calculation project, including poses and checkpoints. Use the regular CLI above to retain calculation artifacts. See [testing instructions](../docs/testing.md).

## Sources and customization

- [ligand_sources.csv](ligand_sources.csv) links each input to its PubChem record.
- [receptor_sources.json](receptor_sources.json) identifies source structures, preparation metadata and file hashes.
- [Data provenance](../docs/test_data.md) records preparation limitations and unspecified stereochemistry.
- [manual_mutations.json](../examples/manual_mutations.json) illustrates fixed mutation probabilities with these receptor names.
- [multiple_targets.json](../examples/multiple_targets.json) uses placeholder names and requires your own matching receptor files.

Keep an unchanged copy of the supplied dataset before replacing inputs. Update provenance when changing structures; the supplied hashes apply only to the supplied receptors/configs. Target/off-target roles and thresholds are demonstration settings, not experimental activity labels. Estimated Ki/pKi values, when reported, use 298.15 K.
