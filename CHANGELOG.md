# Changelog

## 1.0.0 · unreleased

First public version: Python package, bilingual web calculator (English and Spanish) and guided
tour, with the method of the article (manuscript v15).

Changes after the independent audit of 24 November 2026, before release:

- **Reading files.** A quotation mark only opens a quoted value at the start of a cell; an
  unclosed quoted value stops the run with a message, in both tools. Text files are read as
  strict UTF-8 with a clear message on invalid bytes; Windows-1252 can be chosen explicitly
  (`--encoding` and a selector in the calculator). Values such as `1,234` and columns that mix
  decimal comma and point are reported as warnings.
- **Same rules in both tools.** The 30% missing-data rule uses the exact fraction in the
  calculator; "no variation" is judged on values rounded to 10 decimals; the identifier column
  is checked in both.
- **Omega** is computed from a converged factor model (tolerance 10⁻⁸). D2, D3 and D4 change in
  the third decimal with respect to manuscript v14.
- **Calculator.** Complete HTML document with viewport; visible focus on the language buttons;
  Spanish text for every message of the shared core; the status line and the example follow a
  change of language; the same Spanish group labels in the calculator and the tour; fonts bundled
  in `web/fonts/` so that the pages make no request to other servers.
- **File preparation.** Step 1 of the calculator explains how the file should be prepared, shows the first
  rows of the example file and offers the example file and its configuration for download. A file
  whose columns do not match the example asks for the column assignment instead of failing.
  Step 2 opens with a short guide to the column assignment: a saved configuration can be loaded,
  the configuration of the article applied with one button (it is also applied automatically to a
  file with the columns of the Colombian cohort), and a new assignment downloaded as `.json`.
- **Report.** "Download report" saves a self-contained HTML file with the results for the cohort,
  the configuration used (variable, column, direction) and the notes of the page; it opens in
  any browser and can be printed or saved as PDF.
- **Repository.** GitHub Actions workflows for the tests and for GitHub Pages; `CONTRIBUTING.md`;
  `openpyxl>=3.0.7`; SPDX license expression in `pyproject.toml`.

Reproduction of the article and the figures of v14 that were corrected in v15:
`docs/REPRODUCIBILITY.md`.
