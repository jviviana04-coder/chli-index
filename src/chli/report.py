"""Human-readable (Markdown) and machine-readable (JSON) reports of a CHLI run."""
from __future__ import annotations

import json
import math


def _clean(o):
    """Replace NaN and infinities by None so that the JSON is valid."""
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if hasattr(o, "item"):  # numpy scalars
        return _clean(o.item())
    return o


def to_json(rep: dict) -> str:
    return json.dumps(_clean(rep), indent=1, ensure_ascii=False, allow_nan=False)


def _f(x, nd=3):
    return "n/a" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{nd}f}"


def to_markdown(rep: dict) -> str:
    L = [f"# CHLI report · {rep['config']}", "", f"Participants: **{rep['n']}**", ""]
    if rep.get("warnings"):
        L += ["## Please check", ""] + [f"- {w}" for w in rep["warnings"]] + [""]
    L += ["## Domains", "", "| Domain | Type | Items retained | Complete cases | Alpha (95% CI) | Omega | PC1 variance |",
          "|---|---|---|---|---|---|---|"]
    for k, d in rep["domains"].items():
        ci = d.get("alpha_ci95")
        a = _f(d["alpha"]) + (f" ({_f(ci[0])} to {_f(ci[1])})" if ci else "")
        pc = "n/a (formative)" if d["pc1_variance"] is None else f"{100 * d['pc1_variance']:.1f}%"
        if d["pc1_variance"] is not None and not d["pc1_meets_threshold"]:
            pc += " ⚠ below threshold"
        L.append(f"| {k} · {d['label']} | {d['type']} | {len(d['items_retained'])}/{len(d['items_considered'])} "
                 f"| {d['n_complete_cases']} | {a} | {_f(d['omega'])} | {pc} |")
    L += ["", "## Decisions taken by the pipeline", ""]
    for k, d in rep["domains"].items():
        if d.get("excluded_unusable"):
            L.append(f"- **{k}**: excluded because they have no variation or almost no data: {', '.join(d['excluded_unusable'])}")
        if d["excluded_missing"]:
            L.append(f"- **{k}**: excluded for missing data: {', '.join(d['excluded_missing'])}")
        for c in d["excluded_screening"]:
            L.append(f"- **{k}**: removed by reliability screening: {c} (item-total r = {_f(d['item_total'][c])})")
        for w in d.get("warnings", []):
            L.append(f"- **{k}**: {w}")
        if d["participants_without_any_item"]:
            L.append(f"- **{k}**: {d['participants_without_any_item']} participants had no item observed; "
                     f"their score comes entirely from median imputation")
        if d["excluded_screening"]:
            L.append(f"- **{k}**: alpha {_f(d['alpha_before_screening'])} before screening, {_f(d['alpha'])} after")
    L += ["", "## Weights", ""]
    for k, d in rep["domains"].items():
        L.append(f"- **{k}** ({d['aggregation']}): " + ", ".join(f"{c} {w:+.3f}" for c, w in d["weights"].items()))
    s = rep["second_level_pca"]
    L += ["", "## Index", "",
          f"- CHLI_equal: mean {rep['CHLI_equal']['mean']:.3f}, SD {rep['CHLI_equal']['sd']:.3f}, "
          f"median {rep['CHLI_equal']['median']:.3f}, range {rep['CHLI_equal']['min']:.3f} to {rep['CHLI_equal']['max']:.3f}",
          f"- Second-level PCA: PC1 explains {100 * s['pc1_variance']:.1f}% · loadings "
          + ", ".join(f"{k} {v:+.3f}" for k, v in s["loadings"].items()),
          f"- Spearman CHLI_equal vs CHLI_pca2: {s['spearman_with_CHLI_equal']:.3f}", "",
          "## Inter-domain Spearman correlations", ""]
    L += [f"- {k}: {v:+.3f}" for k, v in rep["inter_domain_spearman"].items()]
    L += ["", "## Leave-one-domain-out", ""]
    L += [f"- without {k}: ρ with full index {v['rho_with_full_index']:.3f}; ρ with {k} {v['rho_with_excluded_domain']:+.3f}"
          for k, v in rep["leave_one_domain_out"].items()]
    p = rep["phenotypes"]
    L += ["", "## Phenotypes", "",
          f"- Tertiles: P33 = {p['tertile_cutpoints']['P33']:.3f}, P67 = {p['tertile_cutpoints']['P67']:.3f} · "
          + ", ".join(f"{k} {v}" for k, v in p["tertile_counts"].items()),
          f"- Median split at {p['median_cutpoint']:.3f} · " + ", ".join(f"{k} {v}" for k, v in p["median_counts"].items()),
          "", "---", "",
          "The CHLI characterizes phenotypes within the cohort on which it is built. It is not a "
          "diagnostic or prognostic instrument, and scores are not comparable across cohorts."]
    return "\n".join(L) + "\n"
