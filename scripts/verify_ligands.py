"""Real RDKit/Meeko preparation and cache checks; no Vina is required."""
import argparse
import tempfile
from pathlib import Path
from src.preparation.ligands import LigandConverter
from src.preparation.inputs import read_ligand_rows
from src.core.calculation_cache import atomic_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.csv.resolve() == args.output.resolve():
            raise ValueError('Output must not overwrite the input CSV')
        rows = read_ligand_rows(args.csv)
        if not rows:
            raise ValueError('The ligand CSV is empty')
    except (OSError, ValueError) as error:
        parser.error(str(error))
    results = []
    with tempfile.TemporaryDirectory(prefix='fragevolab-preparation-') as directory:
        converter = LigandConverter(Path(directory))
        for name, smiles in rows:
            try:
                first = converter.prepare_cached(smiles)
                again = converter.prepare_cached(smiles)
                if not again['cache_hit'] or first['artifact_sha'] != again['artifact_sha']:
                    raise RuntimeError('Repeated preparation did not reuse the same cached artifact')
                results.append({'name': name, 'smiles': smiles, 'success': True,
                                'details': first['preparation_details'], 'repeat_cache_hit': True})
            except Exception as error:
                results.append({'name': name, 'smiles': smiles, 'success': False,
                                'error': str(error), 'details': converter.last_preparation_details})
    atomic_json(args.output, results)
    count = sum(row['success'] for row in results)
    print('Prepared', count, 'of', len(results))
    return 0 if count == len(results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
