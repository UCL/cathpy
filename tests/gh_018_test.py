import subprocess
from pathlib import Path

SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "cath-funfhmmer-api"


def test_script():
    result = subprocess.run([SCRIPT_PATH, '-h'])
    assert result.returncode == 0
