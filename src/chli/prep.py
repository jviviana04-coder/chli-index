"""Helpers to prepare inputs that the CHLI expects but that cohorts rarely store as such."""
from __future__ import annotations

import numpy as np
import pandas as pd


def phenoage_residual(phenoage: pd.Series, age: pd.Series, sex: pd.Series) -> pd.Series:
    """Age- and sex-adjusted phenotypic age, rank-based inverse-normal transformed.

    Regresses phenotypic age (Levine et al., 2018) on chronological age and sex by
    ordinary least squares, then applies Blom's rank-based inverse-normal transform to
    the residuals. Higher values mean an older phenotype than expected for age and sex,
    so the item enters the CHLI with polarity -1.

    In the Colombian centenarian cohort this function correlates r = 0.9998 with the
    variable used in the manuscript; the small difference comes from the sample on
    which the original residual was fitted. Computing phenotypic age itself requires the
    nine Levine biomarkers and is outside the scope of this package.
    """
    from scipy.stats import norm, rankdata

    d = pd.DataFrame({"y": pd.to_numeric(phenoage, errors="coerce"),
                      "age": pd.to_numeric(age, errors="coerce"),
                      "sex": pd.to_numeric(sex, errors="coerce")}).dropna()
    X = np.column_stack([np.ones(len(d)), d["age"].values, d["sex"].values])
    beta, *_ = np.linalg.lstsq(X, d["y"].values, rcond=None)
    res = d["y"].values - X @ beta
    c = 3 / 8
    out = pd.Series(np.nan, index=phenoage.index, dtype=float)
    out.loc[d.index] = norm.ppf((rankdata(res) - c) / (len(res) - 2 * c + 1))
    return out
