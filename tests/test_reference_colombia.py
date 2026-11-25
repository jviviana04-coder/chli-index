"""Reproduction of the published figures from the Colombian cohort.

The individual data are not public, so this test only runs when the path to the cohort
master dataset is given in the environment variable ``CHLI_PRIVATE_DATA``. The expected
values below are the aggregate figures reported in the manuscript; they contain no
individual information.

    CHLI_PRIVATE_DATA=/path/to/base_nueva_centenario_con_fenotipo_funcional.xlsx pytest
"""
import os
from pathlib import Path

import pandas as pd
import pytest

import chli

DATA = os.environ.get("CHLI_PRIVATE_DATA")
pytestmark = pytest.mark.skipif(not DATA, reason="private cohort data not available")
ROOT = Path(__file__).resolve().parents[1]

EXPECTED = {
    "alpha": {"D1": 0.8742377812910753, "D2": 0.1369010735274067, "D3": 0.7675525723968123, "D4": 0.624545733520408},
    "alpha_ci95": {"D1": (0.8388486857384151, 0.8996449673963304), "D2": (-0.24843497475306614, 0.3975585379078175),
                   "D3": (0.6959685264039356, 0.8155491459292569), "D4": (0.4916364707836132, 0.7126668371436997)},
    # omega from a converged one-factor model (tol=1e-8). Manuscript v14 reported 0.876, 0.267,
    # 0.773 and 0.644, obtained with the default tolerance of FactorAnalysis, which stops
    # before convergence (audit 2026-11-24, PAQ-10); v15 reports the values below.
    "omega": {"D1": 0.876357063515342, "D2": 0.2692761135930623, "D3": 0.7717641975377156, "D4": 0.6427351013954125},
    "pc1": {"D1": 0.6164106696420484, "D3": 0.45261892263571063, "D4": 0.3772002673762249},
    "n_complete": {"D1": 146, "D2": 83, "D3": 129, "D4": 115},
    "removed": {"D1": [], "D2": [], "D3": ["GAI-SF"], "D4": ["Vitamin D", "Vitamin B12", "Serum creatinine"]},
    "pca2_var": 0.5103938884222563,
    "pca2": {"D1": 0.6016586737198126, "D2": 0.14306289715142884, "D3": 0.5589314828148958, "D4": 0.552390663674449},
    # With n = 163 both tertile cut-points fall exactly on observed values. np.percentile
    # with q = 100/3 and 200/3 computes the virtual indices 54.00000000000001 and
    # 108.00000000000001, so P33 and P67 come out about 1e-16 above the observed values at
    # ranks 55 and 109. Compared exactly, the participant at rank 109 would fall just below
    # P67 and stay Intermediate (54 / 54 / 55, as in manuscript v13-v14). Comparing values
    # rounded to 10 decimals applies the rule as written (Healthy if CHLI >= P67,
    # Non-healthy if CHLI <= P33): 55 / 53 / 55. See docs/REPRODUCIBILITY.md.
    "tertiles": {"Healthy": 55, "Intermediate": 53, "Non-healthy": 55},
    "median_split": {"Healthy": 82, "Non-healthy": 81},
}


@pytest.fixture(scope="module")
def res():
    cfg = chli.config.load(ROOT / "configs" / "colombian_centenarians_manuscript.json")
    return chli.build(pd.read_excel(DATA), cfg)


def test_reliability(res):
    for k, d in res.report["domains"].items():
        assert d["alpha"] == pytest.approx(EXPECTED["alpha"][k], abs=1e-12)
        assert d["alpha_ci95"] == pytest.approx(EXPECTED["alpha_ci95"][k], abs=1e-12)
        assert d["omega"] == pytest.approx(EXPECTED["omega"][k], abs=1e-6)
        assert d["n_complete_cases"] == EXPECTED["n_complete"][k]
        assert d["excluded_screening"] == EXPECTED["removed"][k]


def test_aggregation(res):
    for k, v in EXPECTED["pc1"].items():
        assert res.report["domains"][k]["pc1_variance"] == pytest.approx(v, abs=1e-12)
    s = res.report["second_level_pca"]
    assert s["pc1_variance"] == pytest.approx(EXPECTED["pca2_var"], abs=1e-12)
    for k, v in EXPECTED["pca2"].items():
        assert s["loadings"][k] == pytest.approx(v, abs=1e-12)


def test_phenotypes(res):
    assert res.report["phenotypes"]["tertile_counts"] == EXPECTED["tertiles"]
    assert res.report["phenotypes"]["median_counts"] == EXPECTED["median_split"]
