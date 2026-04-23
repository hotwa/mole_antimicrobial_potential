import re
from pathlib import Path

content = Path("src/batch_screening.py").read_text()

replacement = """from dataclasses import dataclass
from typing import List, Dict, Any

@dataclass
class ScreeningSummary:
    normalized_rows: int
    predicted_rows: int
    normalized_input_path: Path
    predictions_all_path: Path
    grouped_outputs: List[Dict[str, str]]
    grouping_mode: str = "auto"
    cpu_workers_selected: int = 1
    prefetch_queue_size_selected: int = 4
    work_unit_count: int = 0
    target_rows_per_group: int = 0
    target_bytes_per_group: int = 0
"""

content = re.sub(
    r"from dataclasses import dataclass\nfrom typing import List, Dict\n\n@dataclass\nclass ScreeningSummary:.*?grouped_outputs: List\[Dict\[str, str\]\]",
    replacement,
    content,
    flags=re.DOTALL
)

Path("src/batch_screening.py").write_text(content)
