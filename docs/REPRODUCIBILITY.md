# Reproducibility against the article

Run on the original dataset with `configs/colombian_centenarians_manuscript.json`, the package
reproduces to machine precision (differences below 10⁻¹²):

- the four domain scores and the CHLI of all 163 participants;
- Cronbach's alpha and its bootstrap 95% CI for the four domains;
- the number of complete cases per domain;
- the four variables removed by the reliability screening (GAI-SF, vitamin D, vitamin B12,
  serum creatinine) and their item-total correlations;
- the principal-component weights and explained variance of D1, D3 and D4 (Table S2);
- the second-level PCA (51.0% of variance; loadings 0.602, 0.143, 0.559, 0.552);
- the leave-D2-out correlations (0.921 and 0.079);
- the median split (82 Healthy, 81 Non-healthy).

An independent audit (24 November 2026) reimplemented the method from the text of the article
and reproduced the package scores to 2.7 × 10⁻¹⁵. Supplementary Table S5 of the article (v15)
compares every reported estimate with the value returned by the package.

## Differences with manuscript v14, corrected in v15

**Tertile counts.** With n = 163, both tertile cut-points fall exactly on an observed value. The
article's rule (Healthy if CHLI ≥ P67, Non-healthy if CHLI ≤ P33) therefore puts both boundary
participants in the extreme groups: **55 Healthy, 53 Intermediate, 55 Non-healthy**. Version v14
reported 54 / 54 / 55. The noise is in the cut-points, not in the scores: `np.percentile` with
q = 100/3 and 200/3 computes the virtual indices 54.00000000000001 and 108.00000000000001, which
places P33 and P67 about 10⁻¹⁶ above the observed values at ranks 55 and 109. Compared exactly,
the participant at rank 109 stays just below P67. The package compares values rounded to 10
decimals and so applies the rule as written. (`np.quantile(x, [1/3, 2/3])` gives exact virtual
indices and the same 55 / 53 / 55.)

**Spearman correlation D1–D2.** D1 has groups of participants with identical physical-function
inputs after imputation, and therefore tied D1 scores. With ties kept, ρ = 0.0894 (0.089);
v14 reported 0.090. Breaking one of these ties by a difference of the order of 10⁻¹⁶ gives
0.0896 in some orderings, which is compatible with the reported value, but the original
ordering cannot be recovered. The correlations of D1 with D3 and D4 shift only in the fourth
decimal.

**McDonald's omega.** Omega is computed from a one-factor maximum-likelihood model. The values
in v14 (0.876, 0.267, 0.773, 0.644) came from `FactorAnalysis` with its default tolerance
(10⁻²), which stops before convergence. The package uses a tolerance of 10⁻⁸: 0.876, 0.269,
0.772 and 0.643.

**Other figures.** The mean missingness of the retained variables of D3 and D4 (8.6% and
10.4%), the variance explained by PC1 before the reliability screening (37.6% and 24.3%), the
agreement of the median split with and without the screening (157 of 163, 96.3%) and the sample
sizes of Table 1 were corrected in v15 to the values returned by the package.

None of these differences changes a conclusion of the article.

## Phenotypic age residual

The article uses an age- and sex-adjusted, rank-inverse-normal phenotypic age residual computed
before this package existed. `chli.phenoage_residual` follows the same recipe (ordinary least
squares on age and sex, Blom transform of the residuals). In the Colombian cohort it correlates
r = 0.9998 with the original variable; the difference comes from the sample on which the
original residual was fitted. The reference test uses the original variable.
