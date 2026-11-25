"""Behaviour of the pipeline on the synthetic example cohort."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import chli

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def run():
    df = pd.read_csv(ROOT / "examples" / "synthetic_cohort.csv")
    cfg = chli.config.load(ROOT / "configs" / "example_synthetic.json")
    return df, chli.build(df, cfg, id_column="id")


def test_screening_removes_only_the_weak_items(run):
    _, r = run
    assert r.report["domains"]["D1"]["excluded_screening"] == []
    assert r.report["domains"]["D3"]["excluded_screening"] == ["Short anxiety scale"]
    assert r.report["domains"]["D4"]["excluded_screening"] == ["Vitamin B12"]
    assert r.report["domains"]["D2"]["excluded_screening"] == []       # formative: never screened


def test_screening_raises_alpha(run):
    _, r = run
    for k in ("D3", "D4"):
        d = r.report["domains"][k]
        assert d["alpha"] > d["alpha_before_screening"]


def test_scores_are_standardized(run):
    _, r = run
    for k in ("D1", "D2", "D3", "D4", "CHLI"):
        assert abs(r.scores[k].mean()) < 1e-12
        assert abs(r.scores[k].std(ddof=1) - 1) < 1e-12
    assert np.isclose(r.scores["CHLI_equal"].values, r.scores[["D1", "D2", "D3", "D4"]].mean(axis=1).values).all()


def test_polarity_is_respected(run):
    df, r = run
    rho = lambda col: r.scores["D1"].corr(df[col], method="spearman")
    assert rho("barthel") > 0.5
    assert rho("frailty") < -0.3


def test_pca_weights_are_unit_vectors_and_formative_weights_equal(run):
    _, r = run
    for k in ("D1", "D3", "D4"):
        w = np.array(list(r.report["domains"][k]["weights"].values()))
        assert abs((w ** 2).sum() - 1) < 1e-12
        assert (w > 0).all()
    w2 = list(r.report["domains"]["D2"]["weights"].values())
    assert np.allclose(w2, 1 / 3)


def test_phenotype_labels(run):
    _, r = run
    c = r.report["phenotypes"]["tertile_counts"]
    assert sum(c.values()) == 180 and max(c.values()) - min(c.values()) <= 2
    p = r.report["phenotypes"]["tertile_cutpoints"]
    s = r.scores
    assert (s.loc[s.phenotype_tertile == "Healthy", "CHLI_equal"] >= p["P67"]).all()
    assert (s.loc[s.phenotype_tertile == "Non-healthy", "CHLI_equal"] <= p["P33"]).all()


def test_reproducible(run):
    df, r = run
    cfg = chli.config.load(ROOT / "configs" / "example_synthetic.json")
    r2 = chli.build(df, cfg, id_column="id")
    pd.testing.assert_frame_equal(r.scores, r2.scores)
    assert r.report["domains"]["D1"]["alpha_ci95"] == r2.report["domains"]["D1"]["alpha_ci95"]


def test_missing_column_is_reported():
    df = pd.read_csv(ROOT / "examples" / "synthetic_cohort.csv").drop(columns=["sppb"])
    cfg = chli.config.load(ROOT / "configs" / "example_synthetic.json")
    with pytest.raises(KeyError, match="sppb"):
        chli.build(df, cfg)
