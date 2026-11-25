# Contributing to the CHLI

Thank you for using the Centenarians Healthy Longevity Index.

## Reporting a problem

Open an issue in this repository with:

- what you did (calculator or Python package, and the version shown in `CHANGELOG.md`);
- what you expected and what happened, including the exact message;
- a small file that reproduces the problem. **Never attach real participant data**: make a
  synthetic file with the same columns (see `examples/make_synthetic.py`) or replace the values.

## Asking a question

Questions about how to use the index in a cohort, or about the method, can also be opened as
issues. Questions about the article can be sent to the corresponding authors listed in it.

## Proposing a change

1. Open an issue first, so that the change can be discussed.
2. The package and the web calculator must stay identical: a change to the method in
   `src/chli/` needs the same change in `web/chli-core.js`, and a test in
   `tests/test_js_parity.py` or `tests/test_robustness.py` that checks both.
3. Run the tests before proposing the change: `pip install ".[test]"` and `pytest`
   (Node.js is needed for the browser–Python comparison).
4. After changing `web-src/app.html`, `examples/` or `configs/`, run `python web/build_web.py`.
5. Changes that alter published results must be documented in `docs/REPRODUCIBILITY.md` and
   in `CHANGELOG.md`.
