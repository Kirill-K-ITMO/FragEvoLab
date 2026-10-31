"""Generation adapter around the preparation/docking cache."""
import math
import uuid
from typing import List, Dict
from src.core.calculation_cache import molecule_identity, atomic_json

class DockingAdapter:
    def _dock_generation(self, smiles_list: List[str], docker, calculator,
                         proteins: List[str], generation: int, num_workers: int = 3) -> Dict:
        """Every request is checked against the preparation and per-target caches.

        All-or-nothing candidate policy: any failed required receptor evaluation
        excludes the candidate from selection. Partial scores remain only in
        the per-calculation cache/diagnostic JSON, not the candidate checkpoint.
        Offline reweighting cannot recover a rejected candidate's partial scores.
        Derived Ki/pKi are recomputed with the current calculator on every call.
        """
        from src.preparation.ligands import LigandConverter
        if not proteins:
            raise ValueError("An explicit nonempty target list is required")
        proteins = list(dict.fromkeys(proteins))
        converter = LigandConverter(project_root=self.project_root)
        converter.preparation_seed = self.preparation_seed
        prepared, rows, aliases, by_identity = {}, [], {}, {}
        for index, smiles in enumerate(smiles_list, 1):
            label = f"gen{generation}_mol{index}"
            row = {'label': label, 'input_smiles': smiles, 'status': 'preparation_failed'}
            try:
                canonical, molecule_id = molecule_identity(smiles)
                row.update(canonical_smiles=canonical, molecule_id=molecule_id)
                if molecule_id not in by_identity:
                    # One preparation per unique structure in this generation.
                    try:
                        record = converter.prepare_cached(canonical, force=self.force_preparation)
                    except Exception as error:
                        by_identity[molecule_id] = (None, str(error))
                    else:
                        by_identity[molecule_id] = (record, None)
                        prepared[molecule_id] = record['path']
                record, error = by_identity[molecule_id]
                if error is not None:
                    raise RuntimeError(error)
                aliases[label] = molecule_id
                row.update(preparation=record, status='prepared')
            except Exception as error:
                row['error'] = str(error)
                self.logger.error("%s: %s", label, error)
            rows.append(row)
        # Empty lists must never accidentally trigger docking of the shared folder.
        if prepared:
            raw = docker.run_docking(proteins=proteins, exhaustiveness=self.exhaustiveness,
                                     num_workers=num_workers, force=self.force_docking,
                                     ligand_paths=prepared)
            stats = dict(docker.stats)
        else:
            raw, stats = {}, {'cached': 0, 'executed': 0, 'total_dockings': 0, 'failed': 0}
        enriched = calculator.process_docking_results(raw)
        complete = {}
        for row in rows:
            label = row['label']
            if label not in aliases:
                continue
            data = enriched.get(aliases[label], {})
            row['results'] = {target: {**value, 'ligand': label}
                              for target, value in data.items() if target in proteins}
            row['missing_targets'] = [target for target in proteins if target not in data]
            valid = all(target in data and data[target].get('pki') is not None
                        and math.isfinite(data[target]['pki']) for target in proteins)
            row['status'] = 'complete' if valid else 'incomplete'
            if valid:
                complete[label] = row['results']
        hits, executed = stats.get('cached', 0), stats.get('executed', 0)
        misses = stats.get('total_dockings', 0) - hits
        self.last_generation_stats = {'hits': hits, 'misses': misses, 'executed': executed,
                                      'failed': stats.get('failed', 0),
                                      'unique_prepared': len(prepared),
                                      'complete_candidates': len(complete)}
        for field in ('hits', 'misses', 'executed'):
            self.cache_stats[field] += self.last_generation_stats[field]
        # Unique filename also preserves repeated manual evaluations of one generation.
        atomic_json(self.output_dir / f"generation_{generation}_{uuid.uuid4().hex}.json",
                    {'generation': generation, 'targets': proteins, 'candidates': rows,
                     'statistics': self.last_generation_stats,
                     'temperature': getattr(calculator, 'temperature', None)})
        return complete

