"""Messy files, odd configurations and edge cases.

Real cohorts arrive as spreadsheets exported in different locales, with blank rows,
units typed into cells and repeated headers. These tests make sure that the package and
the web calculator either agree on what to do or stop with a clear message.
"""
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import chli
from chli.dataio import parse_number, read_csv_text

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "configs" / "example_synthetic.json").read_text(encoding="utf-8"))
DATA = pd.read_csv(ROOT / "examples" / "synthetic_cohort.csv")
needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")


def cfg(**changes):
    d = json.loads(json.dumps(CFG))
    for path, value in changes.items():
        obj, keys = d, path.split("__")
        for k in keys[:-1]:
            obj = obj[int(k)] if k.isdigit() else obj[k]
        last = keys[-1]
        obj[int(last) if last.isdigit() else last] = value
    return d


# ------------------------------------------------------------------ numbers
CELLS = ["12", " 12 ", "-0.5", "+3", ".5", "5.", "1e-3", "12,5", "-0,25", "", "NA", "n/a", "NULL", ".", "-",
         "1,250.5", "1.250,5", "12 mg", "abc", "inf", "1e999", "0x18", "TRUE"]
EXPECTED = [12, 12, -0.5, 3, 0.5, 5, 0.001, 12.5, -0.25, None, None, None, None, None, None,
            "bad", "bad", "bad", "bad", "bad", "bad", "bad", "bad"]


@pytest.mark.parametrize("cell,expected", list(zip(CELLS, EXPECTED)))
def test_parse_number(cell, expected):
    v, bad = parse_number(cell)
    if expected == "bad":
        assert bad and math.isnan(v)
    elif expected is None:
        assert not bad and math.isnan(v)
    else:
        assert not bad and v == expected


@needs_node
def test_parse_number_same_in_browser():
    js = ("const C=require(process.argv[1]);const cells=JSON.parse(process.argv[2]);"
          "console.log(JSON.stringify(cells.map(c=>{const [v,b]=C.parseNumber(c);return [Number.isNaN(v)?null:v,b];})))")
    out = subprocess.run(["node", "-e", js, str(ROOT / "web" / "chli-core.js"), json.dumps(CELLS)],
                         capture_output=True, text=True, check=True)
    for cell, (v, b) in zip(CELLS, json.loads(out.stdout)):
        pv, pb = parse_number(cell)
        assert b == pb, cell
        assert (v is None and math.isnan(pv)) or v == pv, cell


# ------------------------------------------------------------------ CSV reading
def test_semicolon_file_with_trailing_separators_and_blank_rows():
    text = "id;a;b;\r\n1;1,5;2;\r\n;;;\r\n2;2,5;3;\r\n\r\n"
    df = read_csv_text(text)
    assert len(df) == 2 and list(df.columns[:3]) == ["id", "a", "b"]
    assert df.loc[1, "a"] == "2,5"


def test_repeated_header_is_an_error():
    with pytest.raises(ValueError, match="Repeated column names"):
        read_csv_text("a,b,a\n1,2,3\n")


def test_row_longer_than_header_is_an_error():
    with pytest.raises(ValueError, match="Row 3"):
        read_csv_text("a,b\n1,2\n3,4,5\n")


def test_bom_and_quotes():
    df = read_csv_text('﻿"a","b c"\n"1","x, y"\n')
    assert list(df.columns) == ["a", "b c"] and df.loc[0, "b c"] == "x, y"


# ------------------------------------------------------------------ configuration
@pytest.mark.parametrize("change,ok", [
    ({"domains__0__items__4__polarity": "-1"}, True),
    ({"domains__0__type": "Reflective"}, True),
    ({"domains__0__items__4__polarity": 2}, False),
    ({"domains__0__type": "mixed"}, False),
    ({"domains__0__id": "CHLI"}, False),
    ({"domains__0__id": "1st"}, False),
    ({"domains__3__items__1__target": "<img src=x onerror=alert(1)>"}, False),
    ({"settings": {"threshold": 0.2}}, False),
    ({"settings": {"bootstrap_n": -5}}, False),
])
def test_config_validation(change, ok):
    d = cfg(**change)
    if ok:
        chli.config.from_dict(d)
    else:
        with pytest.raises(ValueError):
            chli.config.from_dict(d)


def test_repeated_item_names_are_an_error():
    d = cfg()
    d["domains"][0]["items"][1]["name"] = d["domains"][0]["items"][0]["name"]
    with pytest.raises(ValueError, match="repeated item names"):
        chli.config.from_dict(d)


@needs_node
@pytest.mark.parametrize("change", [
    {"domains__0__items__4__polarity": 2}, {"domains__0__type": "mixed"}, {"domains__0__id": "CHLI"},
    {"domains__3__items__1__target": "<img src=x onerror=alert(1)>"}, {"settings": {"threshold": 0.2}},
])
def test_browser_rejects_the_same_configurations(change):
    js = "const C=require(process.argv[1]);try{C.validateConfig(JSON.parse(process.argv[2]));console.log('ok')}catch(e){console.log('error')}"
    out = subprocess.run(["node", "-e", js, str(ROOT / "web" / "chli-core.js"), json.dumps(cfg(**change))],
                         capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "error"


# ------------------------------------------------------------------ degenerate data
def test_constant_item_is_set_aside_and_reported():
    df = DATA.copy()
    df["mtdna_cn"] = 1.0
    r = chli.build(df, chli.config.from_dict(CFG))
    assert r.report["domains"]["D2"]["excluded_unusable"] == ["mtDNA-CN"]
    assert r.scores["CHLI"].notna().all()


def test_too_few_complete_cases_gives_a_clear_error():
    df = DATA.copy()
    for k, col in enumerate(["barthel", "lawton_brody", "sppb", "handgrip_kg"]):   # 25% each, never together
        df.loc[k * 45:k * 45 + 44, col] = np.nan
    with pytest.raises(ValueError, match="Domain D1"):
        chli.build(df, chli.config.from_dict(CFG))


def test_rows_without_any_value_are_kept_but_not_scored():
    df = pd.concat([DATA, pd.DataFrame([{"id": "empty"}])], ignore_index=True)
    r = chli.build(df, chli.config.from_dict(CFG), id_column="id")
    assert r.report["n"] == len(DATA) and r.report["n_rows_without_data"] == 1
    assert len(r.scores) == len(df) and math.isnan(r.scores["CHLI"].iloc[-1])
    assert r.scores["CHLI"].iloc[:-1].tolist() == pytest.approx(chli.build(DATA, chli.config.from_dict(CFG))
                                                                  .scores["CHLI"].tolist(), abs=1e-12)


def test_unreadable_cells_are_reported():
    df = DATA.copy().astype({"handgrip_kg": object})
    df.loc[:4, "handgrip_kg"] = "20 kg"
    r = chli.build(df, chli.config.from_dict(CFG))
    assert r.report["unreadable_cells"]["handgrip_kg"]["count"] == 5
    assert any("handgrip_kg" in w for w in r.report["warnings"])


def test_report_json_is_valid_even_with_missing_estimates():
    r = chli.build(DATA, chli.config.from_dict(CFG))
    r.report["domains"]["D1"]["omega"] = float("nan")
    json.loads(chli.to_json(r.report))   # would fail on NaN


# ------------------------------------------------------------------ browser parity on a messy file
@needs_node
def test_browser_and_package_agree_on_a_messy_file(tmp_path):
    df = DATA.copy()
    lines = [";".join(df.columns) + ";"]
    for i, row in df.iterrows():
        cells = [str(v).replace(".", ",") if isinstance(v, float) else str(v) for v in row]
        if i == 3:
            cells[list(df.columns).index("handgrip_kg")] = "18 kg"
        lines.append(";".join(cells) + ";")
        if i == 10:
            lines.append(";" * len(df.columns))
    csv = tmp_path / "messy.csv"
    csv.write_text("﻿" + "\r\n".join(lines) + "\r\n", encoding="utf-8")
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps(CFG), encoding="utf-8")
    py = chli.build(chli.read_table(csv), chli.config.from_dict(CFG))
    js = json.loads(subprocess.run(["node", str(ROOT / "web" / "run_node.js"), str(csv), str(cfg_path)],
                                   capture_output=True, text=True, check=True).stdout)
    s = pd.DataFrame(js["scores"])
    assert py.report["n"] == js["report"]["n"] == len(DATA)
    assert np.max(np.abs(py.scores["CHLI"].values - s["CHLI"].astype(float).values)) < 1e-9
    assert (py.scores["phenotype_tertile"].values == s["phenotype_tertile"].values).all()
    assert py.report["unreadable_cells"] == js["report"]["unreadable_cells"]


# ------------------------------------------------------------------ audit 2026-11-24
# Findings of the independent audit (PAR-01 to PAR-07): both tools must read the same file
# in the same way, or both must stop with the same message.
def _node(fn_body: str, *args):
    js = "const C=require(process.argv[1]);const a=JSON.parse(process.argv[2]);" + fn_body
    out = subprocess.run(["node", "-e", js, str(ROOT / "web" / "chli-core.js"), json.dumps(list(args))],
                         capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


QUOTE_CASES = [
    'a,b\n12",3\n4,5\n',           # stray quote inside a value: an ordinary character
    'a,b\n"q""q",3\n"",4\n',       # escaped quote and empty quoted value
    'a,b\n"open,3\n4,5\n',         # quote never closed: error
    'a,b\n"x"y,3\n',               # text after a closing quote: error
]


def _py_csv(text):
    try:
        return {"rows": read_csv_text(text).values.tolist()}
    except ValueError as e:
        return {"error": str(e)}


@pytest.mark.parametrize("text", QUOTE_CASES)
def test_quotes_python(text):
    r = _py_csv(text)
    if "open" in text or '"x"y' in text:
        assert "error" in r
    else:
        assert "rows" in r and len(r["rows"]) == (2 if '""' in text else 2)


@needs_node
def test_quotes_same_in_browser():
    js = ("console.log(JSON.stringify(a.map(t=>{try{const p=C.parseCSV(t);"
          "return {rows:p.rows.map(r=>p.columns.map(c=>r[c]))}}catch(e){return {error:e.message}}})))")
    for text, j in zip(QUOTE_CASES, _node(js, *QUOTE_CASES)):
        assert _py_csv(text) == j, text


def test_stray_quote_does_not_swallow_the_file():
    rows = "\n".join(f"{i},{i}" for i in range(1, 51))
    df = read_csv_text('a,b\n12",3\n' + rows + "\n")
    assert len(df) == 51


ENC_TEXT = "id;LongTelomérica;x\n1;1,5;2\n2;2,5;3\n"


def test_encoding_is_not_guessed():
    from chli.dataio import decode_bytes
    raw = ENC_TEXT.encode("cp1252")
    with pytest.raises(ValueError, match=r"not UTF-8 text \(byte 0xE9 on line 1\)"):
        decode_bytes(raw)
    assert decode_bytes(raw, "windows-1252") == ENC_TEXT
    assert decode_bytes(b"\xef\xbb\xbf" + ENC_TEXT.encode("utf-8")) == ENC_TEXT


@needs_node
def test_encoding_same_in_browser():
    from chli.dataio import decode_bytes
    raw = list(ENC_TEXT.encode("cp1252"))
    js = ("const b=Uint8Array.from(a[0]);const out=[];"
          "try{C.decodeBytes(b,'utf-8');out.push(null)}catch(e){out.push(e.message)}"
          "out.push(C.decodeBytes(b,'windows-1252'));console.log(JSON.stringify(out))")
    err, txt = _node(js, raw)
    with pytest.raises(ValueError) as e:
        decode_bytes(bytes(raw))
    assert err == str(e.value)
    assert txt == decode_bytes(bytes(raw), "windows-1252")


def test_cli_encoding_option(tmp_path):
    from chli.cli import main
    p = tmp_path / "cohort.csv"
    p.write_bytes(DATA.rename(columns={"id": "código"}).to_csv(index=False, sep=";").encode("cp1252"))
    with pytest.raises(SystemExit) as e:
        main(["build", "--data", str(p), "--config", str(ROOT / "configs" / "example_synthetic.json"),
              "--out", str(tmp_path / "o1"), "--id", "código"])
    assert "not UTF-8 text" in str(e.value.code)
    main(["build", "--data", str(p), "--config", str(ROOT / "configs" / "example_synthetic.json"),
          "--out", str(tmp_path / "o2"), "--id", "código", "--encoding", "windows-1252"])
    assert (tmp_path / "o2" / "chli_scores.csv").exists()


def _build_both(df, config, id_column=None):
    """Python and browser results for the same table (written as CSV)."""
    import tempfile
    try:
        py = chli.build(df, chli.config.from_dict(config), id_column=id_column)
        py = {"scores": py.scores, "report": py.report}
    except (ValueError, KeyError) as e:
        py = {"error": e.args[0] if e.args else str(e)}
    with tempfile.TemporaryDirectory() as tmp:
        csv, cj = Path(tmp) / "d.csv", Path(tmp) / "c.json"
        df.to_csv(csv, index=False)
        cj.write_text(json.dumps(config), encoding="utf-8")
        args = ["node", str(ROOT / "web" / "run_node.js"), str(csv), str(cj), "utf-8"] + ([id_column] if id_column else [])
        js = json.loads(subprocess.run(args, capture_output=True, text=True, check=True).stdout)
    return py, js


@needs_node
@pytest.mark.parametrize("idcol", ["no_such_column", "CHLI", "n_items_observed_D1", "D1"])
def test_identifier_column_checks(idcol):
    df = DATA.copy()
    if idcol != "no_such_column":
        df[idcol] = range(len(df))
    py, js = _build_both(df, cfg(), idcol)
    assert "error" in py and "error" in js
    assert py["error"] == js["error"]


@needs_node
def test_constant_column_with_float_noise_is_unusable():
    df = DATA.copy()
    col = CFG["domains"][0]["items"][0]["column"]
    df[col] = [0.1 + 0.2 if i % 2 else 0.30000000000000004 - 1e-17 for i in range(len(df))]
    py, js = _build_both(df, cfg())
    d = CFG["domains"][0]["id"]
    name = CFG["domains"][0]["items"][0]["name"]
    assert name in py["report"]["domains"][d]["excluded_unusable"]
    assert js["report"]["domains"][d]["excluded_unusable"] == py["report"]["domains"][d]["excluded_unusable"]


@needs_node
def test_missing_rule_uses_the_exact_fraction():
    # 30.005% missing is above the 30% rule even though it rounds to 30.0%
    df = pd.concat([DATA] * 20, ignore_index=True).iloc[:19994].copy()   # n chosen so that k/n = 0.30005...
    n = len(df)
    k = int(np.ceil(0.30005 * n))
    col = CFG["domains"][0]["items"][1]["column"]
    df.loc[: k - 1, col] = np.nan
    py, js = _build_both(df, cfg(**{"settings": {"bootstrap_n": 0}}))
    d = CFG["domains"][0]["id"]
    assert py["report"]["domains"][d]["excluded_missing"] == js["report"]["domains"][d]["excluded_missing"]
    assert py["report"]["domains"][d]["excluded_missing"] == [CFG["domains"][0]["items"][1]["name"]]


@needs_node
def test_number_format_warnings_same_in_both():
    df = DATA.astype(object).copy()
    c1 = CFG["domains"][0]["items"][0]["column"]
    c2 = CFG["domains"][0]["items"][2]["column"]
    df.loc[0, c1] = "1,234"
    df.loc[1, c2] = "2,5"
    df.loc[2, c2] = "2.75"
    py, js = _build_both(df, cfg())
    pw = [w for w in py["report"]["warnings"] if "decimal" in w]
    jw = [w for w in js["report"]["warnings"] if "decimal" in w]
    assert len(pw) == 2 and len(jw) == 2
    assert [w.split(":")[0] for w in pw] == [w.split(":")[0] for w in jw]


@needs_node
@pytest.mark.parametrize("values,warns", [
    (["1,234", "2,500", "3,750"], True),            # every comma value could be a thousands separator
    (["1,147", "0,982", "1,282"], False),           # 0,982 shows that the column uses a decimal comma
    (["1,147", "1,15", "1,282"], False),            # so does 1,15
    (["1.5", "2,5", "3.25"], True),                 # decimal comma and decimal point mixed
])
def test_number_format_rule(values, warns):
    from chli.dataio import number_format_warning
    py = number_format_warning("x", values)
    js = _node("console.log(JSON.stringify(C.numberFormatWarning('x', a)))", *values)
    assert (py is not None) == warns
    assert py == js
