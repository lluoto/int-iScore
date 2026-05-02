"""Legacy input-data generation wrapper.

This exposes the historical `input_data_generate.py` workflow from inside the package
 tree while preserving its standalone nature and external template/module assumptions.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[3]
    script = repo_root / "input_data_generate.py"
    if not script.exists():
        raise FileNotFoundError(f"Legacy input-data generator not found: {script}")
    result = subprocess.run([sys.executable, str(script)], cwd=str(repo_root))
    return int(result.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
