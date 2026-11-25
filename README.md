# CHLI · Centenarians Healthy Longevity Index

The CHLI summarizes how a very old person is aging across four domains: **physical function**,
**molecular biological aging**, **cognitive-psychological status** and **metabolic-nutritional
status**. It was developed in a cohort of 163 Colombian centenarians.

This repository lets any group **build the CHLI in its own cohort** with the same method as the
article, either in a browser or in Python. It does not contain the Colombian participant data.

- **Web calculator:** <https://jviviana04-coder.github.io/chli-index/>
- **Guided tour:** <https://jviviana04-coder.github.io/chli-index/tour.html>
- **Code and documentation:** <https://github.com/jviviana04-coder/chli-index>

> **What the CHLI is not.** It is not a diagnostic test or a prognostic score. It ranks people
> *within the cohort on which it is built*, so scores from two cohorts are not on the same scale.
> In the original cohort the physical function domain (D1) alone discriminated mortality as well as
> the full index.

---

## Two ways to use it

| | Who it is for | What you need |
|---|---|---|
| **Guided tour** · `web/tour.html` | Anyone who wants to understand the index first | A browser |
| **Web calculator** · `web/index.html` | Anyone. No installation, no programming | A CSV file with your cohort |
| **Python package** · `src/chli` | Analysts who want to script, audit or extend the method | Python ≥ 3.9 |

Both give the same numbers: `tests/test_js_parity.py` checks the browser code against the Python
code on the example cohort and, when the data are available, on the Colombian cohort.

### Web calculator

Open `web/index.html` in a browser, or publish the `web/` folder with GitHub Pages: the workflow
`.github/workflows/pages.yml` does it on every push to `main` (*Settings → Pages → Source: GitHub
Actions*). The calculation runs entirely in the browser: **the data you load are processed on
your computer and are not sent anywhere**. The pages make no request to other servers: the
fonts are included in `web/fonts/` (SIL Open Font License). The page opens with a simulated
cohort already loaded so you can see what it does. Both the calculator and the guided tour are
in **English and Spanish** (button at the top; `?lang=es` or `?lang=en` in the address also works).

1. Load a CSV with one row per participant. From Excel, save it as **CSV UTF-8**; a file saved
   as plain "CSV" by an older Spanish-language Excel is Windows-1252, which can be chosen in the
   *File encoding* list.
2. Map your columns to the variables of each domain and set their direction.
3. Press **Build the index**, then copy or download the scores and the report.

**How files are read.** The calculator and the package apply the same rules
(`src/chli/dataio.py`, `parseNumber` in `web/chli-core.js`): a cell is a number only if it
looks like one (`12`, `-0.5`, `1e-3`, or a decimal comma such as `12,5`); thousands separators,
units (`12 kg`) or text are treated as missing **and listed in the report**, never guessed.
Empty rows are ignored, repeated column names stop the run with a message, and variables
without variation are set aside and reported. A quotation mark only opens a quoted value at the
start of a cell, and a quoted value that is never closed stops the run. Files are read as UTF-8
unless Windows-1252 is chosen; the encoding is never guessed. Values such as `1,234`, which could
carry a thousands separator, are read as decimals and flagged in the report.

The guided tour (`web/tour.html`) explains the method in nine short chapters, running the
real pipeline on the simulated cohort; it can be used to present the index.

### Python

```bash
pip install .            # from a clone of this repository
pip install ".[test]"    # with the test dependencies (pytest)
```

The browser–Python comparison tests need [Node.js](https://nodejs.org) on the PATH; without it they
are skipped. On Windows PowerShell, set the variable for the reference test with
`$env:CHLI_PRIVATE_DATA = "C:\path\to\data.xlsx"; pytest`.

```python
import pandas as pd
import chli

cfg = chli.config.load("configs/example_synthetic.json")
res = chli.build(chli.read_table("examples/synthetic_cohort.csv"), cfg, id_column="id")

res.scores.head()                       # one row per participant
res.report["domains"]["D3"]             # alpha, omega, screening, weights...
print(chli.to_markdown(res.report))     # a readable report
res.report["warnings"]                  # anything in the file that needs checking
```

Or from the command line:

```bash
chli build --data examples/synthetic_cohort.csv --config configs/example_synthetic.json --out results/ --id id
```

This writes `chli_scores.csv`, `chli_report.json` and `chli_report.md`. With your own file, `--id`
is the name of its identifier column; add `--encoding windows-1252` for a CSV saved in that
encoding.

---

## Your data

- One row per participant, one column per variable. CSV or Excel (Excel only in Python).
- Missing values may be left empty. They are imputed by the variable's median, as in the article,
  and the output says how many variables each participant actually had in each domain.
- **Phenotypic age** must enter as an age- and sex-adjusted residual computed in your own cohort:
  `chli.phenoage_residual(phenoage, age, sex)`. Computing phenotypic age itself needs the nine
  biomarkers of Levine et al. (2018) and is outside this package.
- **Mitochondrial DNA copy number and telomere length** measured by qPCR are relative values
  specific to the assay and laboratory. That is one reason the index is rebuilt in each cohort
  rather than scored against reference values.

## The configuration

A JSON file says which column feeds each variable and in which direction it points. Copy
`configs/example_synthetic.json` and edit the `column` fields:

```json
{
  "name": "My cohort",
  "domains": [
    {"id": "D1", "label": "Physical function", "type": "reflective", "items": [
      {"name": "Barthel index", "column": "barthel", "polarity": 1},
      {"name": "Frailty score", "column": "frailty", "polarity": -1}
    ]},
    {"id": "D4", "label": "Metabolic-nutritional status", "type": "reflective", "items": [
      {"name": "BMI, distance to optimum", "column": "bmi", "polarity": 1,
       "transform": "neg_abs_distance", "target": 24}
    ]}
  ],
  "settings": {"screening_threshold": 0.20, "max_missing_reflective": 0.30,
               "pc1_min_variance": 0.20, "bootstrap_n": 2000, "seed": 20261020}
}
```

- `polarity`: `1` if higher values mean better health, `-1` if they mean worse.
- `type`: `reflective` for a scale whose variables are expected to move together (summarized by
  its first principal component), `formative` for a set of distinct markers that together define
  the domain without having to correlate (summarized by their mean).
- `transform: neg_abs_distance` turns a variable whose optimum lies in the middle, such as body
  mass index, into its negative distance to `target`.

`configs/colombian_centenarians_manuscript.json` is the exact configuration of the article, with
the column names of the original dataset.

## What the pipeline does

The steps follow the article (Sections 2.5 to 2.6.5):

1. **Prepare.** Each variable is converted to numbers, transformed if requested, oriented so that
   higher means healthier, and standardized (z-score, sample SD).
2. **Missing data.** In reflective domains, variables with more than 30% missing values are
   excluded. Formative domains keep them and report a warning.
3. **Reliability screening.** In reflective domains, every variable whose corrected item-total
   correlation (Spearman, complete cases) is below 0.20 is removed, in a single pass. Cronbach's
   alpha with a 2000-replicate bootstrap CI and McDonald's omega are reported for every domain.
4. **Aggregation.** After per-variable median imputation, reflective domains are summarized by
   their first principal component, oriented to agree with the mean of the variables; the 20%
   explained-variance quality threshold is checked. Formative domains are the mean of their
   standardized markers. Every domain score is standardized.
5. **Index.** `CHLI_equal` is the mean of the four domain scores; `CHLI` is its standardized
   version; `CHLI_pca2` is the first component of a second-level PCA on the domain scores.
6. **Phenotypes.** Tertiles of `CHLI_equal` (Healthy ≥ P67, Non-healthy ≤ P33, Intermediate
   otherwise) and a median split (Healthy ≥ median).

It also reports the Spearman correlations between domains and a leave-one-domain-out check of
how much the index depends on each domain. `docs/METHOD.md` gives the formulas.

## Output columns

| Column | Meaning |
|---|---|
| `D1` … `D4` | Standardized domain scores |
| `CHLI_equal` | Mean of the domain scores |
| `CHLI` | `CHLI_equal` standardized (mean 0, SD 1) |
| `CHLI_pca2` | Second-level PCA score |
| `phenotype_tertile` | Healthy / Intermediate / Non-healthy |
| `phenotype_median` | Healthy / Non-healthy |
| `n_items_observed_D1` … | Variables actually observed for that participant in each domain |

---

## Reproducing the article

The participant data are not public. With access to the original dataset:

```bash
CHLI_PRIVATE_DATA=/path/to/the_original_dataset.xlsx pytest
```

The package reproduces the published reliability figures, bootstrap intervals, omega, removed
variables, principal-component weights and variances, second-level PCA and domain scores to
machine precision. `docs/REPRODUCIBILITY.md` lists the figures of manuscript v14 that were
corrected in v15 and their cause.

Without the private data, `pytest` runs the tests on the synthetic cohort and the browser–Python
comparison.

## Repository

```
src/chli/        Python package (config, dataio, stats, pipeline, report, cli, prep)
web/             published pages: calculator (index.html), guided tour (tour.html), shared
                 pipeline (chli-core.js), bundled fonts (fonts/)
web-src/         source of the calculator (app.html); web/build_web.py writes web/index.html
configs/         the article's configuration and the example configuration
examples/        synthetic cohort and the script that generates it
tests/           synthetic, reference, browser-parity and messy-file tests
docs/            method and reproducibility notes
```

The example cohort in `examples/synthetic_cohort.csv` is simulated from a model written by hand
(`examples/make_synthetic.py`). No value comes from a real participant.

## Contributing and support

See `CONTRIBUTING.md` for how to report a problem, ask a question or propose a change.
Changes between versions are listed in `CHANGELOG.md`.

## Citation

If you use the CHLI, please cite the article (see `CITATION.cff`).

## License

See `LICENSE`.
