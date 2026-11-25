"""The CHLI construction pipeline.

The steps follow the manuscript (Sections 2.5 to 2.6.5) in order:

1. Prepare every item: numeric conversion, optional transform, polarity, z-score.
2. Exclude reflective items with too much missing data (default > 30%).
3. Reliability screening of reflective domains: drop items whose corrected
   item-total correlation is below the threshold (default Spearman < 0.20), in one pass.
4. Aggregate each domain: PC1 for reflective domains, mean of z-scores for formative
   ones, after per-variable median imputation. Domain scores are standardized.
5. Build the index: CHLI_equal (mean of the domain scores) and CHLI_pca2 (second-level PCA).
6. Label phenotypes: tertiles and median split of CHLI_equal.

Every decision is recorded in ``Result.report`` so that a run can be audited.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import dataio, stats
from .config import Config, Item


@dataclass
class Result:
    scores: pd.DataFrame
    report: dict = field(default_factory=dict)


MIN_COMPLETE = 3      # below this, reliability cannot be estimated at all
FEW_COMPLETE = 30     # below this, estimates are reported with a warning
OUTPUT_PREFIX = "n_items_observed_"


def _raw_item(df: pd.DataFrame, it: Item) -> tuple[pd.Series, list[str]]:
    if it.column not in df.columns:
        raise KeyError(f"Column '{it.column}' (item '{it.name}') is not in the data")
    x, bad = dataio.to_numeric(df[it.column])
    if it.transform == "neg_abs_distance":
        x = -(x - it.target).abs()
    return it.polarity * x, bad


def build(df: pd.DataFrame, cfg: Config, id_column: str | None = None) -> Result:
    """Run the full pipeline on ``df`` with configuration ``cfg``.

    Rows with no usable value in any configured column are kept in the output with empty
    scores and are not used to build the index (they are counted in the report).
    """
    st = cfg.settings
    rng = np.random.default_rng(st.seed)
    dupcols = sorted({c for c in df.columns[df.columns.duplicated()] if c in set(cfg.columns)})
    if dupcols:
        raise ValueError(f"Repeated column names in the data: {', '.join(map(str, dupcols))}")
    if id_column is not None and id_column not in df.columns:
        raise KeyError(f"Identifier column '{id_column}' is not in the data")
    reserved = {"CHLI", "CHLI_equal", "CHLI_pca2", "phenotype_tertile", "phenotype_median"}
    if id_column is not None and (id_column in reserved or id_column in {d.id for d in cfg.domains}
                                  or str(id_column).startswith(OUTPUT_PREFIX)):
        raise ValueError(f"The identifier column cannot be called '{id_column}'; it needs another name")

    raw: dict[str, dict[str, pd.Series]] = {}
    unreadable: dict[str, dict] = {}
    for dom in cfg.domains:
        raw[dom.id] = {}
        for it in dom.items:
            x, bad = _raw_item(df, it)
            raw[dom.id][it.name] = x
            if bad:
                unreadable[it.column] = {"count": len(bad), "examples": bad[:3]}
    format_warnings = []
    for col in dict.fromkeys(cfg.columns):
        w = dataio.number_format_warning(col, df[col].tolist())
        if w:
            format_warnings.append(w)
    has_data = pd.concat([x.notna() for d in raw.values() for x in d.values()], axis=1).any(axis=1)
    keep = has_data.values
    n_in = int(len(df))
    sub_index = df.index[keep]
    if keep.sum() < MIN_COMPLETE:
        raise ValueError(f"Only {int(keep.sum())} rows have any usable value; the index needs more participants")

    rep: dict = {"config": cfg.name, "n": int(keep.sum()), "n_rows_in_file": n_in,
                 "n_rows_without_data": int(n_in - keep.sum()), "settings": vars(st).copy(),
                 "unreadable_cells": unreadable, "warnings": [], "domains": {}}
    for col, u in unreadable.items():
        rep["warnings"].append(f"Column '{col}': {u['count']} cell(s) could not be read as numbers and were "
                               f"treated as missing (e.g. {', '.join(repr(e) for e in u['examples'])})")
    rep["warnings"] += format_warnings
    if rep["n_rows_without_data"]:
        rep["warnings"].append(f"{rep['n_rows_without_data']} row(s) have no usable value and were left out")

    final_X: dict[str, pd.DataFrame] = {}
    for dom in cfg.domains:
        X = pd.DataFrame({k: v[keep] for k, v in raw[dom.id].items()}, index=sub_index)
        d = {"label": dom.label, "type": dom.type, "items_considered": list(X.columns),
             "missing_pct": {c: round(100 * float(X[c].isna().mean()), 1) for c in X.columns},
             "excluded_unusable": [], "excluded_missing": [], "excluded_screening": [],
             "item_total": {}, "warnings": []}

        # 1 · items without variation (or almost no data) carry no information
        # (no variation is judged on rounded values, so that floating-point noise in a constant
        # column, e.g. 0.1 + 0.2, is not taken for variation; same rule in the web calculator)
        unusable = [c for c in X.columns if X[c].notna().sum() < MIN_COMPLETE
                    or X[c].dropna().round(stats.TIE_DECIMALS).nunique() < 2]
        d["excluded_unusable"] = unusable
        X = X.drop(columns=unusable)
        X = X.apply(stats.zscore)

        # 2 · missingness rule (reflective domains only)
        if dom.type == "reflective":
            too = [c for c in X.columns if X[c].isna().mean() > st.max_missing_reflective]
            d["excluded_missing"] = too
            X = X.drop(columns=too)
        else:
            d["warnings"] += [f"{c}: {d['missing_pct'][c]}% missing (formative domain, kept)"
                              for c in X.columns if X[c].isna().mean() > st.max_missing_reflective]
        if X.shape[1] < 2:
            raise ValueError(f"Domain {dom.id}: fewer than two usable items "
                             f"(without variation: {', '.join(unusable) or 'none'}; "
                             f"too much missing: {', '.join(d['excluded_missing']) or 'none'})")

        # 3 · reliability screening (reflective domains only, one pass)
        n_cc = int(len(X.dropna()))
        if dom.type == "reflective":
            if n_cc < MIN_COMPLETE:
                raise ValueError(f"Domain {dom.id}: only {n_cc} participant(s) have every item, so the "
                                 "reliability screening cannot be done. The column mapping should be checked, "
                                 "or the item with the most missing values removed.")
            d["alpha_before_screening"] = stats.cronbach_alpha(X)
            it_tot = stats.item_total(X, method=st.screening_method)
            d["item_total"] = it_tot
            drop = [c for c, r in it_tot.items() if not r >= st.screening_threshold]
            d["excluded_screening"] = drop
            X = X.drop(columns=drop)
            if X.shape[1] < 2:
                raise ValueError(f"Domain {dom.id}: fewer than two items left after the reliability screening")
        else:
            d["alpha_before_screening"] = stats.cronbach_alpha(X)
        d["items_retained"] = list(X.columns)
        d["n_complete_cases"] = int(len(X.dropna()))
        if d["n_complete_cases"] < FEW_COMPLETE:
            d["warnings"].append(f"only {d['n_complete_cases']} complete cases: reliability estimates are unstable")
        final_X[dom.id] = X
        rep["domains"][dom.id] = d

    # reliability of the final domains; the bootstrap stream runs over domains in config order
    for dom in cfg.domains:
        X = final_X[dom.id]
        d = rep["domains"][dom.id]
        d["alpha"] = stats.cronbach_alpha(X)
        ci = stats.alpha_bootstrap_ci(X, rng, st.bootstrap_n) if st.bootstrap_n else None
        d["alpha_ci95"] = list(ci) if ci else None
        try:
            d["omega"] = stats.mcdonald_omega(X)
        except ImportError:
            d["omega"] = None

    # 4 · aggregation
    S: dict[str, pd.Series] = {}
    nobs: dict[str, pd.Series] = {}
    for dom in cfg.domains:
        X = final_X[dom.id]
        d = rep["domains"][dom.id]
        nobs[dom.id] = X.notna().sum(axis=1)
        Xi = X.fillna(X.median())
        if dom.type == "reflective":
            ratio, comps, sc = stats.pca(Xi.values)
            s, comp = sc[:, 0], comps[0]
            if np.corrcoef(s, Xi.mean(axis=1))[0, 1] < 0:
                s, comp = -s, -comp
            S[dom.id] = stats.zscore(pd.Series(s, index=sub_index))
            d["aggregation"] = "PC1"
            d["pc1_variance"] = float(ratio[0])
            d["pc1_meets_threshold"] = bool(ratio[0] >= st.pc1_min_variance)
            d["weights"] = {c: float(w) for c, w in sorted(zip(X.columns, comp), key=lambda t: -abs(t[1]))}
        else:
            S[dom.id] = stats.zscore(Xi.mean(axis=1))
            d["aggregation"] = "mean of z-scores (formative)"
            d["pc1_variance"] = None
            d["weights"] = {c: 1 / X.shape[1] for c in X.columns}
            d["descriptive_pca_variance"] = [float(v) for v in stats.pca(Xi.values)[0]]
        d["participants_without_any_item"] = int((nobs[dom.id] == 0).sum())

    CH = pd.DataFrame(S)

    # 5 · index
    chli_equal = CH.mean(axis=1)
    ratio2, comps2, sc2 = stats.pca(CH.values)
    L, pca2 = comps2[0], sc2[:, 0]
    if L.sum() < 0:
        L, pca2 = -L, -pca2

    # 6 · phenotypes
    p33, p67 = np.percentile(chli_equal, [100 / 3, 200 / 3])
    med = float(chli_equal.median())
    # compare on rounded values so that a participant sitting exactly on a cut-point is
    # classified the same way regardless of floating-point noise
    r = stats.TIE_DECIMALS
    ce = chli_equal.round(r)
    tert = pd.Series("Intermediate", index=sub_index)
    tert[ce <= round(p33, r)] = "Non-healthy"
    tert[ce >= round(p67, r)] = "Healthy"   # checked last: Healthy wins if P33 == P67
    binl = pd.Series(np.where(ce >= round(med, r), "Healthy", "Non-healthy"), index=sub_index)
    if round(p33, r) == round(p67, r):
        rep["warnings"].append("The two tertile cut-points are equal (many tied index values); "
                               "tertile labels are not informative in this cohort")

    out = pd.DataFrame(index=sub_index)
    for k in CH.columns:
        out[k] = CH[k]
    out["CHLI_equal"] = chli_equal
    out["CHLI"] = stats.zscore(chli_equal)
    out["CHLI_pca2"] = pca2
    out["phenotype_tertile"] = tert
    out["phenotype_median"] = binl
    for k in CH.columns:
        out[f"{OUTPUT_PREFIX}{k}"] = nobs[k]
    out = out.reindex(df.index)
    for k in CH.columns:
        out[f"{OUTPUT_PREFIX}{k}"] = out[f"{OUTPUT_PREFIX}{k}"].fillna(0).astype(int)
    if id_column is not None:
        out.insert(0, id_column, df[id_column].values)

    ids = list(CH.columns)
    rep["second_level_pca"] = {"pc1_variance": float(ratio2[0]), "loadings": dict(zip(ids, map(float, L))),
                               "spearman_with_CHLI_equal": stats.spearman(chli_equal, pd.Series(pca2, index=sub_index))}
    rep["inter_domain_spearman"] = {f"{a}-{b}": stats.spearman(CH[a], CH[b])
                                    for i, a in enumerate(ids) for b in ids[i + 1:]}
    loo = {}
    for k in ids:
        rest = stats.zscore(CH.drop(columns=[k]).mean(axis=1))
        loo[k] = {"rho_with_full_index": stats.spearman(stats.zscore(chli_equal), rest),
                  "rho_with_excluded_domain": stats.spearman(rest, CH[k])}
    rep["leave_one_domain_out"] = loo
    rep["CHLI_equal"] = {"mean": float(chli_equal.mean()), "sd": float(chli_equal.std(ddof=1)),
                         "median": med, "min": float(chli_equal.min()), "max": float(chli_equal.max())}
    rep["phenotypes"] = {"tertile_cutpoints": {"P33": float(p33), "P67": float(p67)},
                         "tertile_counts": {k: int(v) for k, v in tert.value_counts().items()},
                         "median_cutpoint": med,
                         "median_counts": {k: int(v) for k, v in binl.value_counts().items()}}
    return Result(scores=out, report=rep)
