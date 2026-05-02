"""Legacy frustratometer runner wrapper.

This keeps the historical `frustratometer2-master/run.py` workflow reachable from the
package tree without claiming full integration into the packaged mode/API layer.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[3]
    script = repo_root / "frustratometer2-master" / "run.py"
    if not script.exists():
        raise FileNotFoundError(f"Legacy frustratometer runner not found: {script}")
    result = subprocess.run([sys.executable, str(script)], cwd=str(repo_root))
    return int(result.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
