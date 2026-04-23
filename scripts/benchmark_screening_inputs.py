from __future__ import annotations

import time
from pathlib import Path
from typing import Any


def benchmark_paths(input_paths: list[str]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for input_path in input_paths:
        path = Path(input_path).expanduser().resolve()
        start = time.time()
        size_bytes = path.stat().st_size if path.is_file() else 0
        rows.append(
            {
                "input_path": str(path),
                "elapsed_seconds": time.time() - start,
                "size_bytes": size_bytes,
            }
        )
    return {"benchmarks": rows}
