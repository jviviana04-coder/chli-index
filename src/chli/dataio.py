"""Reading tables and turning cells into numbers, with the same rules as the web calculator.

The browser (web/chli-core.js) and this module must read a file in exactly the same way,
otherwise the two tools would give different indices for the same cohort. The rules are
deliberately strict and visible:

* A cell is a number only if it looks like one: ``12``, ``-0.5``, ``1e-3``, or a decimal
  comma with nothing else (``12,5``). Thousands separators, units (``12 mg``) or text are
  not guessed: the cell is treated as missing and counted in the report, so that the user
  can see it and fix the file.
* Empty cells and the usual missing-value codes (``NA``, ``N/A``, ``NaN``, ``NULL``,
  ``None``, ``#N/A``, ``.``, ``-``) are missing without a warning.
* Rows that are completely empty (for example ``;;;;`` exported by a spreadsheet) are
  ignored. A row with more filled cells than the header is an error.
* Repeated column names are an error: it would be ambiguous which one to use.
* A quotation mark opens a quoted value only as the first character of a cell; anywhere
  else it is an ordinary character. A quoted value that is never closed is an error.
* Text files are read as UTF-8 (with or without BOM). Bytes that are not valid UTF-8 stop
  the reading with a message; the encoding is never guessed. ``windows-1252`` (the
  default of older Spanish-language spreadsheet programs) can be declared explicitly.
* Values that look like a thousands separator (``1,234``) or columns that mix a decimal
  comma and a decimal point are read by the rules above and reported as a warning.
"""
from __future__ import annotations

import codecs
import csv
import io
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

MISSING_TOKENS = {"", "na", "n/a", "nan", "null", "none", "#n/a", "#na", "<na>", "-nan", ".", "-"}
_NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
_DEC_COMMA = re.compile(r"^[+-]?\d+,\d+$")
_THOUSANDS_COMMA = re.compile(r"^[+-]?[1-9]\d{0,2},\d{3}$")
_DEC_POINT = re.compile(r"^[+-]?(\d+\.\d*|\.\d+)$")
ENCODINGS = ("utf-8", "windows-1252")

# WHATWG windows-1252 (what a browser uses) maps the five bytes that Python's cp1252 leaves
# undefined to the C1 control characters; do the same so that both tools read the same text.
codecs.register_error("chli_c1", lambda e: (chr(e.object[e.start]), e.start + 1))


def parse_number(v) -> tuple[float, bool]:
    """Return (value, unreadable). Missing cells give (nan, False)."""
    if v is None or isinstance(v, bool):
        return math.nan, v is not None
    if isinstance(v, (int, float, np.integer, np.floating)):
        f = float(v)
        if math.isnan(f):
            return math.nan, False
        return (f, False) if math.isfinite(f) else (math.nan, True)
    if v is pd.NA or v is pd.NaT:
        return math.nan, False
    t = str(v).strip()
    if t.lower() in MISSING_TOKENS:
        return math.nan, False
    if _NUM.match(t):
        f = float(t)
    elif _DEC_COMMA.match(t):
        f = float(t.replace(",", "."))
    else:
        return math.nan, True
    return (f, False) if math.isfinite(f) else (math.nan, True)


def to_numeric(s: pd.Series) -> tuple[pd.Series, list[str]]:
    """Convert a column with :func:`parse_number`; also return the unreadable cells."""
    vals, bad = [], []
    for v in s.tolist():
        f, unreadable = parse_number(v)
        vals.append(f)
        if unreadable:
            bad.append(str(v))
    return pd.Series(vals, index=s.index, dtype=float), bad


def number_format_warning(column: str, values) -> str | None:
    """Warning for cells whose number format is ambiguous (same text as the web calculator)."""
    texts = [str(v).strip() for v in values if isinstance(v, str)]
    thousands = [t for t in texts if _THOUSANDS_COMMA.match(t)]
    comma_cells = [t for t in texts if _DEC_COMMA.match(t)]
    comma = bool(comma_cells)
    point = any(_DEC_POINT.match(t) for t in texts)
    # the comma is only ambiguous if every comma value in the column could be a thousands
    # separator; a value such as 0,982 or 1,15 shows that the column uses a decimal comma
    if thousands and len(thousands) == len(comma_cells):
        ex = ", ".join(json.dumps(t, ensure_ascii=False) for t in thousands[:3])
        return (f"Column '{column}': {len(thousands)} cell(s) such as {ex} were read as decimal numbers; "
                "if the comma is a thousands separator, the file needs to be corrected")
    if comma and point:
        return (f"Column '{column}': some cells use a decimal comma and others a decimal point; "
                "the values should be checked")
    return None


def decode_bytes(data: bytes, encoding: str = "utf-8") -> str:
    """Decode a text file without guessing its encoding."""
    enc = encoding.lower().replace("_", "-")
    if enc in ("utf8", "utf-8", "utf-8-sig"):
        if data.startswith(b"\xef\xbb\xbf"):
            data = data[3:]
        try:
            return data.decode("utf-8")
        except UnicodeDecodeError as e:
            line = data[:e.start].count(b"\n") + 1
            raise ValueError(f"The file is not UTF-8 text (byte 0x{data[e.start]:02X} on line {line}). "
                             "It can be saved again as 'CSV UTF-8', or the encoding can be set to windows-1252.") from None
    if enc in ("windows-1252", "cp1252"):
        return data.decode("cp1252", errors="chli_c1")
    raise ValueError(f"Unknown encoding '{encoding}' (use utf-8 or windows-1252)")


def detect_delimiter(first_line: str) -> str:
    counts = [(d, first_line.count(d)) for d in (",", ";", "\t")]
    return max(counts, key=lambda t: t[1])[0]  # first one wins a tie, as in the browser


def read_csv_text(text: str) -> pd.DataFrame:
    text = text.lstrip("﻿")
    if not text.strip():
        raise ValueError("The file is empty")
    delim = detect_delimiter(re.split(r"\r\n|\r|\n", text, maxsplit=1)[0])
    reader = csv.reader(io.StringIO(text, newline=""), delimiter=delim, strict=True)
    rows = []
    try:
        for r in reader:
            rows.append(r)
    except csv.Error as e:
        if "unexpected end of data" in str(e):
            raise ValueError("The file ends inside a value that opens with a quotation mark (\") and is never "
                             "closed. A quotation mark may be missing or there may be one too many.") from None
        raise ValueError(f"Line {reader.line_num}: after a closing quotation mark (\") a separator or the end "
                         "of the line is expected.") from None
    head = [h.strip() for h in rows[0]]
    named = [h for h in head if h]
    dup = sorted({h for h in named if named.count(h) > 1})
    if dup:
        raise ValueError(f"Repeated column names: {', '.join(dup)}. Each column needs a unique name.")
    body = []
    for i, r in enumerate(rows[1:], start=2):
        if all(not f.strip() for f in r):
            continue
        if len(r) > len(head) and any(f.strip() for f in r[len(head):]):
            raise ValueError(f"Row {i} has {len(r)} cells but the header has {len(head)}. "
                             "The row may have an extra separator or an unquoted comma inside a value.")
        body.append([(r[j] if j < len(r) else "") for j in range(len(head))])
    return pd.DataFrame(body, columns=head, dtype=object)


def read_table(path: str | Path, encoding: str = "utf-8") -> pd.DataFrame:
    """Read a CSV (comma, semicolon or tab) or an Excel file.

    ``encoding`` applies to text files only: ``utf-8`` (default, with or without BOM) or
    ``windows-1252``.
    """
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(p)
        df.columns = [str(c).strip() for c in df.columns]
        return df.dropna(how="all").reset_index(drop=True)
    return read_csv_text(decode_bytes(p.read_bytes(), encoding))
