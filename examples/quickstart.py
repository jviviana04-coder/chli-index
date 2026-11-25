"""Build the CHLI on the synthetic example cohort and print the report.

Run from the repository root:  python examples/quickstart.py
"""
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import chli  # noqa: E402

cfg = chli.config.load(ROOT / "configs" / "example_synthetic.json")
data = chli.read_table(ROOT / "examples" / "synthetic_cohort.csv")
res = chli.build(data, cfg, id_column="id")

print(chli.to_markdown(res.report))
print(res.scores[["id", "D1", "D2", "D3", "D4", "CHLI", "phenotype_tertile"]].head(10).round(3).to_string(index=False))
