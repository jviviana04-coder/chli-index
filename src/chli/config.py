"""Domain configuration for the CHLI pipeline.

A configuration says which columns of the user's table feed each domain, in which
direction each variable points, and how each domain is aggregated. It is plain JSON so
that the same file drives the Python package and the web calculator; the browser applies
the same checks (web/chli-core.js, ``validateConfig``).
"""
from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any

VALID_TYPES = {"reflective", "formative"}
VALID_TRANSFORMS = {None, "neg_abs_distance"}
VALID_METHODS = {"spearman", "pearson"}
RESERVED_IDS = {"CHLI", "CHLI_equal", "CHLI_pca2", "phenotype_tertile", "phenotype_median"}
_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,31}$")


def _polarity(v: Any, name: str) -> int:
    try:
        p = float(v)
    except (TypeError, ValueError):
        p = math.nan
    if p not in (1.0, -1.0) or isinstance(v, bool):
        raise ValueError(f"Item '{name}': direction (polarity) must be 1 (higher is better) "
                         f"or -1 (higher is worse), got {v!r}")
    return int(p)


@dataclass
class Item:
    name: str
    column: str
    polarity: int = 1
    transform: str | None = None
    target: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("Every item needs a non-empty 'name'")
        if not isinstance(self.column, str) or not self.column.strip():
            raise ValueError(f"Item '{self.name}': 'column' must be the name of a column in the data")
        self.polarity = _polarity(self.polarity, self.name)
        if self.transform not in VALID_TRANSFORMS:
            raise ValueError(f"Item '{self.name}': unknown transform {self.transform!r} "
                             "(use none or 'neg_abs_distance')")
        if self.transform == "neg_abs_distance":
            try:
                self.target = float(self.target)
            except (TypeError, ValueError):
                self.target = math.nan
            if not math.isfinite(self.target):
                raise ValueError(f"Item '{self.name}': transform 'neg_abs_distance' needs a numeric 'target'")


@dataclass
class Domain:
    id: str
    label: str
    type: str
    items: list[Item]

    def __post_init__(self) -> None:
        if not isinstance(self.id, str) or not _ID.match(self.id):
            raise ValueError(f"Domain id {self.id!r}: use a short name with letters, digits or '_' "
                             "that starts with a letter (for example D1)")
        if self.id in RESERVED_IDS or self.id.startswith("n_items_observed_"):
            raise ValueError(f"Domain id '{self.id}' is reserved for an output column; another id is needed")
        t = str(self.type).strip().lower()
        if t not in VALID_TYPES:
            raise ValueError(f"Domain '{self.id}': type must be 'reflective' or 'formative', got {self.type!r}")
        self.type = t
        if len(self.items) < 2:
            raise ValueError(f"Domain '{self.id}': needs at least two items")
        names = [it.name for it in self.items]
        dup = sorted({n for n in names if names.count(n) > 1})
        if dup:
            raise ValueError(f"Domain '{self.id}': repeated item names: {', '.join(dup)}")


@dataclass
class Settings:
    screening_method: str = "spearman"
    screening_threshold: float = 0.20
    max_missing_reflective: float = 0.30
    pc1_min_variance: float = 0.20
    bootstrap_n: int = 2000
    seed: int = 20261020

    def __post_init__(self) -> None:
        if self.screening_method not in VALID_METHODS:
            raise ValueError("screening_method must be 'spearman' or 'pearson'")
        for k, lo, hi in (("screening_threshold", -1, 1), ("max_missing_reflective", 0, 1), ("pc1_min_variance", 0, 1)):
            v = getattr(self, k)
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not lo <= v <= hi:
                raise ValueError(f"{k} must be a number between {lo} and {hi}, got {v!r}")
        for k, lo in (("bootstrap_n", 0), ("seed", 0)):
            v = getattr(self, k)
            if isinstance(v, bool) or not isinstance(v, int) or v < lo:
                raise ValueError(f"{k} must be a whole number of at least {lo}, got {v!r}")


@dataclass
class Config:
    name: str
    domains: list[Domain]
    settings: Settings = field(default_factory=Settings)
    description: str = ""

    @property
    def columns(self) -> list[str]:
        return [it.column for d in self.domains for it in d.items]


def _item(d: dict[str, Any]) -> Item:
    return Item(name=d.get("name"), column=d.get("column"), polarity=d.get("polarity", 1),
                transform=d.get("transform"), target=d.get("target"))


def from_dict(d: dict[str, Any]) -> Config:
    if not isinstance(d, dict) or not isinstance(d.get("domains"), list):
        raise ValueError("The configuration needs a list of 'domains'")
    doms = [Domain(id=x.get("id"), label=x.get("label") or x.get("id"), type=x.get("type"),
                   items=[_item(i) for i in x.get("items", [])]) for x in d["domains"]]
    if len(doms) < 2:
        raise ValueError("The index needs at least two domains")
    ids = [x.id for x in doms]
    if len(set(ids)) != len(ids):
        raise ValueError("Domain ids must be unique")
    raw = d.get("settings", {}) or {}
    allowed = {f.name for f in fields(Settings)}
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError(f"Unknown settings: {', '.join(unknown)} (allowed: {', '.join(sorted(allowed))})")
    return Config(name=d.get("name", "CHLI"), domains=doms, settings=Settings(**raw),
                  description=d.get("description", ""))


def load(path: str | Path) -> Config:
    with open(path, encoding="utf-8") as fh:
        return from_dict(json.load(fh))
