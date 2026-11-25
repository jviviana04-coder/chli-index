"""CHLI: the Centenarians Healthy Longevity Index.

Build the index in your own cohort::

    import pandas as pd, chli
    cfg = chli.config.load("configs/example_synthetic.json")
    res = chli.build(pd.read_csv("examples/synthetic_cohort.csv"), cfg, id_column="id")
    res.scores.head(); res.report["domains"]["D1"]["alpha"]
"""
from . import config, dataio, stats
from .dataio import read_table
from .pipeline import Result, build
from .prep import phenoage_residual
from .report import to_json, to_markdown

__all__ = ["build", "Result", "config", "dataio", "stats", "read_table", "to_json", "to_markdown",
           "phenoage_residual"]
__version__ = "1.0.0"
