import subprocess
import sys
from pathlib import Path


def test_removed_calculator_choice_is_rejected_and_menu_exits(tmp_path):
    launcher = Path(__file__).resolve().parents[1] / 'main.py'
    result = subprocess.run(
        [sys.executable, str(launcher), '--project', str(tmp_path)],
        input='17\n0\n', text=True, capture_output=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'Ki calculator' not in result.stdout
    assert 'Invalid choice' in result.stdout
    assert 'Thank you for using FragEvoLab' in result.stdout
