"""Generate the synthetic example cohort shipped with the repository.

THE DATA PRODUCED HERE ARE ENTIRELY SIMULATED. No value is copied from, sampled from or
fitted to any real participant. The generative model below was written by hand, with
scale ranges taken from the published instruments, so that the example behaves like a
plausible centenarian cohort: three correlated reflective domains, one formative domain
whose markers barely covary, a few weak items that the reliability screening should
remove, and realistic missingness.

Run:  python examples/make_synthetic.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

rng = np.random.default_rng(2026)
N = 180

# latent structure: a general health factor shared by D1, D3 and D4
g = rng.normal(size=N)
f1 = 0.80 * g + 0.60 * rng.normal(size=N)   # physical function
f3 = 0.70 * g + 0.71 * rng.normal(size=N)   # cognitive-psychological
f4 = 0.55 * g + 0.84 * rng.normal(size=N)   # metabolic-nutritional


def item(f, loading, noise=None):
    noise = np.sqrt(1 - loading ** 2) if noise is None else noise
    return loading * f + noise * rng.normal(size=N)


def scale(x, lo, hi, step=1.0):
    """Map a standard-normal-like signal onto an instrument range, rounding to its step."""
    u = (x - x.min()) / (x.max() - x.min())
    v = lo + u * (hi - lo)
    return np.clip(np.round(v / step) * step, lo, hi)


df = pd.DataFrame({"id": [f"S{i:03d}" for i in range(1, N + 1)]})
df["age"] = rng.integers(100, 112, N)
df["sex"] = rng.binomial(1, 0.28, N)

# D1 · physical function (higher = better, except frailty and sarcopenia)
df["barthel"] = scale(item(f1, 0.85), 0, 100, 5)
df["lawton_brody"] = scale(item(f1, 0.70), 0, 8)
df["sppb"] = scale(item(f1, 0.80), 0, 12)
df["handgrip_kg"] = scale(item(f1, 0.65), 4, 30, 0.5)
df["frailty"] = scale(-item(f1, 0.65), 0, 2)
df["sarcopenia"] = scale(-item(f1, 0.60), 1, 2)

# D3 · cognitive-psychological (CDR, GDS, ISI and anxiety: higher = worse)
df["mmse"] = scale(item(f3, 0.70), 0, 30)
df["cdr_sb"] = scale(-item(f3, 0.75), 0, 18, 0.5)
df["gds15"] = scale(-item(f3, 0.70), 0, 15)
df["swls"] = scale(item(f3, 0.60), 5, 35)
df["isi"] = scale(-item(f3, 0.40), 0, 28)
df["anxiety_short"] = scale(-item(f3, 0.05), 0, 5)          # weak item: screening should drop it

# D4 · metabolic-nutritional
df["mna"] = scale(item(f4, 0.75), 8, 30, 0.5)
bmi_dist = 1.5 * np.log1p(np.exp(-1.6 * item(f4, 0.65)))           # farther from 24 when f4 is low
df["bmi"] = np.round(24 + bmi_dist * np.where(rng.random(N) < 0.75, -1, 1), 1)
df["hemoglobin"] = np.round(13 + 1.4 * item(f4, 0.60), 1)
df["albumin"] = np.round(3.8 + 0.4 * item(f4, 0.60), 2)
df["crp_mg_l"] = np.round(np.exp(1.0 - 0.9 * item(f4, 0.45)), 2)
df["vitamin_b12"] = np.round(np.exp(6.0 + 0.4 * rng.normal(size=N)))  # unrelated: screening should drop it

# D2 · molecular markers: nearly independent of each other and of the other domains
df["mtdna_cn"] = np.round(np.exp(5.0 + 0.35 * (0.10 * g + rng.normal(size=N))), 1)
df["telomere_ts"] = np.round(0.9 + 0.18 * (0.25 * np.log(df["mtdna_cn"]) + rng.normal(size=N)), 3)
df["phenoage"] = np.round(df["age"] - 8 + 6 * (-0.08 * g + rng.normal(size=N)), 1)

# missingness, completely at random, at rates similar to a field cohort
for col, rate in {"handgrip_kg": .06, "frailty": .05, "mmse": .15, "cdr_sb": .08, "gds15": .07,
                  "swls": .09, "isi": .03, "anxiety_short": .06, "mna": .01, "bmi": .04,
                  "hemoglobin": .08, "albumin": .25, "crp_mg_l": .11, "vitamin_b12": .12,
                  "mtdna_cn": .10, "telomere_ts": .28, "phenoage": .30}.items():
    df.loc[rng.random(N) < rate, col] = np.nan

# the CHLI uses phenotypic age adjusted for age and sex (see chli.phenoage_residual)
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chli import phenoage_residual  # noqa: E402

df["phenoage_residual"] = np.round(phenoage_residual(df["phenoage"], df["age"], df["sex"]), 4)

out = Path(__file__).with_name("synthetic_cohort.csv")
df.to_csv(out, index=False)
print(f"{N} synthetic participants written to {out}")
