import unittest
from pathlib import Path
from src.batch_screening import _append_tsv
import pandas as pd
import tempfile

class TestAppend(unittest.TestCase):
    def test_append(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.tsv"
            df1 = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
            _append_tsv(df1, path, True)
            df2 = pd.DataFrame({"a": [5], "b": [6]})
            _append_tsv(df2, path, False)
            res = pd.read_csv(path, sep='\t')
            self.assertEqual(len(res), 3)
            self.assertEqual(list(res["a"]), [1, 2, 5])

if __name__ == "__main__":
    unittest.main()
