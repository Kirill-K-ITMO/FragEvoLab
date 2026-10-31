"""Check translated exports and scientific unit formatting."""
import csv
import re

from src.docking.affinity import AffinityCalculator
from src.reporting.docking import ReportGenerator
from src.preparation.smiles import SmilesGenerator


def test_english_docking_exports_preserve_numeric_values(project):
    reporter = ReportGenerator(project[0])
    data = {'example': {'A': {'affinity': -6.0, 'smiles': 'CCO',
                            'ki': 1e-5, 'pki': 5.0, 'ki_formatted': '10 µM'}}}
    html, table = reporter.generate_reports(data)
    text = html.read_text(encoding='utf-8')
    assert 'lang="en"' in text
    assert not re.search('[\u0400-\u04ff]', text)
    assert 'Docking score (kcal/mol)' in text
    assert '298.15 K' in text
    with table.open(encoding='utf-8', newline='') as stream:
        row = next(csv.DictReader(stream))
    assert row['name'] == 'example'
    assert float(row['A_delta_g']) == -6.0
    assert float(row['A_ki']) == 1e-5
    assert float(row['A_pki']) == 5.0


def test_english_smiles_export_header(tmp_path):
    generator = SmilesGenerator(tmp_path)
    generator.set_template('C{1}')
    generator.add_radicals('1', ['N'])
    generator.generate_combinations()
    output = tmp_path / 'molecules.csv'
    assert generator.save_to_csv(output_file=output)
    with output.open(encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == ['name', 'smiles']
        assert next(reader)['smiles'] == 'CN'


def test_ki_unit_prefixes_are_international():
    calculator = AffinityCalculator()
    for value, unit in [(1.0, 'M'), (1e-3, 'mM'), (1e-6, 'µM'),
                        (1e-9, 'nM'), (1e-12, 'pM')]:
        assert calculator.format_ki(value).split()[-1] == unit
