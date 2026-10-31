"""SA-score colors follow the reported numerical category boundaries."""
import pytest
from src.reporting.docking import ReportGenerator


@pytest.mark.parametrize('score,category', [
    (2.0, 'sa-good'), (3.0, 'sa-moderate'),
    (5.9, 'sa-moderate'), (6.0, 'sa-hard'), (8.0, 'sa-hard'),
])
def test_sa_category_in_report(tmp_path, monkeypatch, score, category):
    reporter = ReportGenerator(tmp_path)
    descriptors = reporter.calculate_molecular_descriptors('CCO')
    descriptors['sa_score'] = score
    monkeypatch.setattr(reporter, 'calculate_molecular_descriptors', lambda _: descriptors)
    path = reporter.generate_html_report({'ethanol': {'target': {'affinity': -5.0}}},
                                         {'ethanol': 'CCO'})
    assert f'<span class="{category}">{score:.2f}</span>' in path.read_text(encoding='utf-8')
