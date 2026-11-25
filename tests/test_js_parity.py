"""The web calculator (web/chli-core.js) must give the same results as the Python package."""
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import chli

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
TOL = 1e-9


def _run_js(csv_path: Path, cfg_path: Path) -> dict:
    out = subprocess.run(["node", str(ROOT / "web" / "run_node.js"), str(csv_path), str(cfg_path)],
                         capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


def _compare(df: pd.DataFrame, cfg_path: Path, csv_path: Path) -> None:
    py = chli.build(df, chli.config.load(cfg_path))
    js = _run_js(csv_path, cfg_path)
    s_js = pd.DataFrame(js["scores"])
    for c in ["D1", "D2", "D3", "D4", "CHLI_equal", "CHLI", "CHLI_pca2"]:
        assert np.max(np.abs(py.scores[c].values - s_js[c].values)) < TOL, c
    for c in ["phenotype_tertile", "phenotype_median"]:
        assert (py.scores[c].values == s_js[c].values).all(), c
    for k, d in py.report["domains"].items():
        j = js["report"]["domains"][k]
        assert d["items_retained"] == j["items_retained"]
        assert d["excluded_screening"] == j["excluded_screening"]
        assert d["n_complete_cases"] == j["n_complete_cases"]
        assert abs(d["alpha"] - j["alpha"]) < TOL
        assert abs(d["alpha_before_screening"] - j["alpha_before_screening"]) < TOL
        for nm, r in d["item_total"].items():
            assert abs(r - j["item_total"][nm]) < TOL, nm
        if d["pc1_variance"] is not None:
            assert abs(d["pc1_variance"] - j["pc1_variance"]) < TOL
        for nm, w in d["weights"].items():
            assert abs(w - j["weights"][nm]) < TOL, nm
    for key in ["inter_domain_spearman"]:
        for k, v in py.report[key].items():
            assert abs(v - js["report"][key][k]) < TOL, k
    s = py.report["second_level_pca"]
    assert abs(s["pc1_variance"] - js["report"]["second_level_pca"]["pc1_variance"]) < TOL
    p, q = py.report["phenotypes"], js["report"]["phenotypes"]
    assert abs(p["tertile_cutpoints"]["P33"] - q["tertile_cutpoints"]["P33"]) < TOL
    assert p["tertile_counts"] == q["tertile_counts"]


def test_parity_synthetic():
    csv = ROOT / "examples" / "synthetic_cohort.csv"
    _compare(pd.read_csv(csv), ROOT / "configs" / "example_synthetic.json", csv)


@pytest.mark.skipif(not os.environ.get("CHLI_PRIVATE_DATA"), reason="private cohort data not available")
def test_parity_colombian_cohort():
    df = pd.read_excel(os.environ["CHLI_PRIVATE_DATA"])
    with tempfile.TemporaryDirectory() as tmp:
        csv = Path(tmp) / "cohort.csv"
        df.to_csv(csv, index=False)
        _compare(pd.read_csv(csv), ROOT / "configs" / "colombian_centenarians_manuscript.json", csv)
