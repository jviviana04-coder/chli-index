# Method

Notation: *x<sub>ij</sub>* is variable *j* for participant *i*; *n* participants; a domain has *k*
variables after the missing-data rule.

## 0 · Reading the data

Text files are read as UTF-8 (with or without BOM). The encoding is never guessed: a byte that is
not valid UTF-8 stops the reading with a message giving the byte and the line, and
`windows-1252` can be declared explicitly (`--encoding` in the command line, a selector in the
calculator). A quotation mark opens a quoted value only as the first character of a cell; a
quoted value that is never closed is an error. Cells are converted to numbers with strict rules
shared by the package and the calculator (see the README): unreadable cells become missing and
are listed in `report["unreadable_cells"]` and `report["warnings"]`. Cells that look like a
thousands separator (`1,234`) and columns that mix a decimal comma and a decimal point are read
by the same rules and reported as a warning. Rows with no usable value in any configured column are kept in the
output with empty scores and are not used to build the index.

## 1 · Preparation

A variable with fewer than three observed values, or with no variation (fewer than two distinct
values after rounding to 10 decimals), carries no information and is set aside
(`excluded_unusable`). For each remaining variable:

- optional transform: `neg_abs_distance` gives −|x − target| (body mass index, target 24 kg/m²);
- orientation: multiply by the polarity (+1 or −1) so that higher always means healthier;
- standardization: z = (x − mean) / SD, with the sample SD (n − 1), over the observed values.

## 2 · Missing data

In a reflective domain, a variable with more than `max_missing_reflective` (default 0.30) of its
values missing (the exact fraction, not a rounded percentage) is excluded before screening. A formative domain keeps all its markers.
A reflective domain needs at least three participants with every variable observed; otherwise
the run stops with a message naming the domain.

## 3 · Reliability

On complete cases of the domain:

- **Cronbach's alpha** = k / (k − 1) · (1 − Σ var(z<sub>j</sub>) / var(Σ z<sub>j</sub>)).
- **95% CI**: percentile bootstrap of alpha, resampling complete cases, `bootstrap_n` replicates
  (default 2000). One random stream (seed 20261020) runs over the domains in configuration order.
- **McDonald's omega** = (Σλ)² / ((Σλ)² + Σψ), from a one-factor maximum-likelihood model
  (scikit-learn `FactorAnalysis`, `tol=1e-8`, `max_iter=10000`, `random_state=0`; with the
  default tolerance the algorithm stops before convergence). Python only.
- **Corrected item-total correlation**: Spearman correlation between each variable and the sum of
  the other variables of the domain, on complete cases of the full variable set.

**Screening** (reflective domains only): variables with a corrected item-total correlation below
`screening_threshold` (default 0.20) are removed, all at once. Alpha is reported before and after.

## 4 · Aggregation

Missing values of each standardized variable are replaced by that variable's median.

- **Reflective domain**: principal component analysis by eigendecomposition of the covariance
  matrix of the imputed variables. The score is the projection on the first component, with its
  sign chosen so that it correlates positively with the row mean of the variables. Its
  explained-variance share is compared with `pc1_min_variance` (default 0.20).
- **Formative domain**: the row mean of the imputed standardized markers.

Each domain score is then standardized (mean 0, SD 1).

Why a formative domain: principal components extract variance the variables share. When the
markers of a domain index distinct processes and do not covary (in the Colombian cohort, alpha =
0.137 for mitochondrial DNA copy number, telomere length and the phenotypic age residual), the
first component has no substantive meaning, and equal weights are the rule that follows from the
construct.

## 5 · Index

- `CHLI_equal` = mean of the domain scores (equal weights).
- `CHLI` = standardized `CHLI_equal`.
- `CHLI_pca2`: projection on the first component of a PCA of the domain scores, with the
  loadings' sign chosen so that they sum to a positive number.

## 6 · Phenotypes

- Tertiles of `CHLI_equal` with NumPy's default linear interpolation: Healthy if ≥ P67, Non-healthy
  if ≤ P33, Intermediate otherwise.
- Median split: Healthy if ≥ median, Non-healthy otherwise.

Comparisons with cut-points and rank correlations are made on values rounded to 10 decimals, so
that participants with identical inputs stay tied whatever the numerical path.

## 7 · Diagnostics

- Spearman correlation between every pair of domain scores.
- Leave-one-domain-out: for each domain, the index rebuilt without it is correlated with the full
  index and with the excluded domain.
