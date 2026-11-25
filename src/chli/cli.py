"""Command-line interface: ``chli build --data cohort.csv --config config.json --out results/``."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import config as cfgmod
from .dataio import read_table
from .pipeline import build
from .report import to_json, to_markdown


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="chli", description="Build the Centenarians Healthy Longevity Index")
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="build the index in a cohort")
    b.add_argument("--data", required=True, help="CSV or Excel file, one row per participant")
    b.add_argument("--config", required=True, help="JSON configuration (see configs/)")
    b.add_argument("--out", default="chli_results", help="output folder")
    b.add_argument("--id", default=None, help="name of the participant identifier column")
    b.add_argument("--encoding", default="utf-8", choices=["utf-8", "windows-1252"],
                   help="text encoding of a CSV file (default utf-8; Excel 'CSV UTF-8' is utf-8)")
    a = ap.parse_args(argv)

    try:
        cfg = cfgmod.load(a.config)
        res = build(read_table(a.data, encoding=a.encoding), cfg, id_column=a.id)
    except (ValueError, KeyError, OSError) as e:
        msg = e.args[0] if isinstance(e, KeyError) and e.args else e
        sys.exit(f"chli: {msg}")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res.scores.to_csv(out / "chli_scores.csv", index=False)
    (out / "chli_report.json").write_text(to_json(res.report), encoding="utf-8")
    (out / "chli_report.md").write_text(to_markdown(res.report), encoding="utf-8")
    for w in res.report.get("warnings", []):
        print(f"warning: {w}", file=sys.stderr)
    print(f"{res.report['n']} participants · scores, report.json and report.md written to {out}/")


if __name__ == "__main__":
    main()
