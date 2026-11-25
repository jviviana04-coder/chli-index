"""Small, dependency-light statistics used by the pipeline.

Every function here has a line-by-line counterpart in ``web/chli.js``; the two are
checked against each other by ``tests/test_js_parity.py``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def zscore(s: pd.Series) -> pd.Series:
    """Standardize with the sample SD (ddof = 1), ignoring missing values."""
    return (s - s.mean()) / s.std(ddof=1)


def cronbach_alpha(X: pd.DataFrame) -> float:
    """Cronbach's alpha on complete cases (NaN when it cannot be estimated)."""
    Xc = X.dropna()
    k = Xc.shape[1]
    if k < 2 or len(Xc) < 3:
        return float("nan")
    tot = Xc.sum(axis=1).var(ddof=1)
    if not tot > 0:
        return float("nan")
    return float((k / (k - 1)) * (1 - Xc.var(ddof=1).sum() / tot))


def alpha_bootstrap_ci(X: pd.DataFrame, rng: np.random.Generator, n_boot: int = 2000,
                       level: float = 0.95) -> tuple[float, float]:
    """Percentile bootstrap CI for alpha, resampling complete cases (None if not estimable)."""
    Xc = X.dropna().values
    n, k = Xc.shape
    if n < 3 or k < 2:
        return None
    out = []
    for _ in range(n_boot):
        s = Xc[rng.integers(0, n, n)]
        t = s.sum(axis=1).var(ddof=1)
        if t > 0:
            out.append((k / (k - 1)) * (1 - s.var(axis=0, ddof=1).sum() / t))
    if not out:
        return None
    a = (1 - level) / 2 * 100
    lo, hi = np.percentile(out, [a, 100 - a])
    return float(lo), float(hi)


def mcdonald_omega(X: pd.DataFrame) -> float:
    """McDonald's omega from a one-factor maximum-likelihood model on complete cases.

    Needs scikit-learn. Returns NaN with fewer than three items or too few complete cases.
    """
    Xc = X.dropna()
    if Xc.shape[1] < 3 or len(Xc) <= Xc.shape[1]:
        return float("nan")
    from sklearn.decomposition import FactorAnalysis
    try:
        # a tight tolerance: with the default (tol=1e-2) the EM algorithm stops before the
        # maximum-likelihood solution and omega can differ in the third decimal
        fa = FactorAnalysis(n_components=1, tol=1e-8, max_iter=10000, random_state=0).fit(Xc)
    except (ValueError, np.linalg.LinAlgError):
        return float("nan")
    L = fa.components_[0]
    u = fa.noise_variance_
    return float((L.sum() ** 2) / ((L.sum() ** 2) + u.sum()))


TIE_DECIMALS = 10


def _stable(s: pd.Series) -> pd.Series:
    """Round before ranking so that rows with identical inputs stay tied.

    Different numerical paths (NumPy, a browser) can leave identical rows a few units
    of 1e-16 apart, which would otherwise break their ties in a rank correlation.
    """
    return s.round(TIE_DECIMALS)


def item_total(X: pd.DataFrame, method: str = "spearman") -> dict[str, float]:
    """Corrected item-total correlation: each item against the sum of the others.

    Computed on complete cases of the full item set, as in the manuscript.
    """
    Xc = X.dropna()
    return {c: float(_stable(Xc[c]).corr(_stable(Xc.drop(columns=[c]).sum(axis=1)), method=method))
            for c in Xc.columns}


def pca(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """PCA by eigendecomposition of the covariance matrix.

    Returns (explained variance ratio, components as rows, scores). Equivalent to
    ``sklearn.decomposition.PCA`` up to the sign of each component.
    """
    Xc = X - X.mean(axis=0)
    C = np.cov(Xc, rowvar=False, ddof=1)
    w, v = np.linalg.eigh(C)
    order = np.argsort(w)[::-1]
    w, v = w[order], v[:, order]
    return w / w.sum(), v.T, Xc @ v


def spearman(x: pd.Series, y: pd.Series) -> float:
    o = x.notna() & y.notna()
    return float(_stable(x[o]).corr(_stable(y[o]), method="spearman"))
