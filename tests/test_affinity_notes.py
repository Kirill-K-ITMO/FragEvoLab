import csv
import math
from pathlib import Path
import pytest
from src.docking.affinity import AffinityCalculator
from src.reporting.docking import ReportGenerator


def test_standard_temperature_and_units():
    calculator = AffinityCalculator()
    assert calculator.temperature == 298.15
    # A score corresponding to a 1 micromolar estimate at the standard temperature.
    score = -6 * math.log(10) * 8.314463 * 298.15 / 4184
    ki = calculator.calculate_ki(score)
    assert ki == pytest.approx(1e-6)
    assert calculator.calculate_pki(ki) == pytest.approx(6)


def test_html_csv_notes_keep_values_and_schema(project):
    reporter = ReportGenerator(project[0])
    data = AffinityCalculator().process_docking_results({'sample': {'A': {'affinity': -7.0}}})
    html, table = reporter.generate_reports(data, smiles_dict={'sample': 'CCO'})
    text = Path(html).read_text()
    assert '298.15 K' in text and 'not experimental affinity' in text
    assert 'Imported results or custom API temperatures' in text
    note = Path(str(table) + '.notes.txt').read_text()
    assert '298.15 K' in note and 'kcal/mol' in note
    with Path(table).open() as stream:
        row = next(csv.DictReader(stream))
    assert float(row['A_ki']) == data['sample']['A']['ki']
    assert float(row['A_pki']) == data['sample']['A']['pki']
