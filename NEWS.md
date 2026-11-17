# GIR 0.2.1

* Updated both quantreg references in the guide to use the canonical CRAN
  package URL. Replaced the README hyperlink to excluded maintainer instructions
  with a plain-text reference to their location in the source repository.
  These changes address CRAN's review of the first submission and do not alter
  estimation or inference.

# GIR 0.2.0

* Prepared the first CRAN release. The detailed guide is now a PDF vignette
  with an executable introductory example. Revised documentation wording and
  added `citation("GIR")` entries for the software and its underlying working
  paper. The paper has not yet been published on arXiv.
* Added full release checks that build the vignette and PDF function reference,
  plus Linux checks on GitHub Actions. These changes do not alter estimation
  or inference calculations.
* Made the polynomial design-grouping test portable across numerical libraries.
  It now checks exact reconstruction of the design and its actual distinct
  rows, in addition to agreement with the independent coefficient calculation.

* Empty estimation samples after `na.omit`, `na.exclude` or a custom action
  now give an informative error before matrix construction or the routine
  removal warning. Positive weights are checked in the retained sample,
  and entirely missing weights follow the requested missing-value policy.
* First-stage support checks now use retained observations. If omission leaves
  an undeclared numeric covariate with at most two distinct values, an explicit
  `discrete` or `continuous` declaration is required. Custom actions report
  missing versus complete removals accurately, preserve their own error
  context and reject invalid row identifiers. Estimator formulas are unchanged.
* Added regression tests for sample validation, row alignment and empty
  prediction/interval shapes. Replaced all broad warning suppression in the
  point-estimator and bootstrap/prediction tests with checks of the complete
  expected warning set or an explicit absence of warnings. User-facing
  numerical and inference warnings remain enabled.
* Added deterministic mathematical validation of both estimation stages and
  bootstrap refits against independent QR-based threshold regressions for
  full-rank Legendre and spline tensor bases. Added pairwise group-Gini,
  positive outcome-scale and row-order invariance, and GIR-R-squared identity
  tests. Existing expected bootstrap warnings are now asserted explicitly.
  The new `GIR-validation` help topic explains the tests and their limits;
  estimator calculations, defaults and returned results are unchanged.
* First-stage ridge penalties now scale with the average diagonal of the
  weighted Gram matrix, without the former absolute floor of 1. The initial
  and retry penalties both respect common weight rescaling. Unusable scales,
  penalties and adjusted matrices give informative errors, and successful
  ridge retries are identified in the diagnostics. The `ridge = 1e-10`
  default is unchanged and documented as a numerical heuristic. Ridge-based
  point estimates and inference can change where the old floor was active;
  affected analyses should refit and rerun their bootstrap.
* Added `cluster_order = "auto"`, `"original"`, or `"estimation"` to
  `bboot.gir()`, also available through automatic bootstraps in `predict()`.
  Automatic ordering now rejects cluster vectors when the original data and
  estimation sample have equal lengths but different row mappings. Explicit
  ordering prevents silent cluster misassignment after rows are reordered or
  repeated. Retained-column formulas are already aligned. Affected analyses
  should rerun their bootstrap; original point estimates are unchanged.
* Fixed multi-column `cluster` formulas accidentally merging distinct groups
  when their combined text labels collided (for example, `("a.b", "c")` and
  `("a", "b.c")`). Grouping now uses exact value combinations, including
  numeric identifiers that print alike, and retains first-occurrence order.
  Affected analyses should rerun their bootstrap: draws and inference can
  change; original point estimates are unaffected. Added regression tests
  and clarified that multi-column formulas specify intersection clusters.
* Added bboot.gir(): common normalized exponential multipliers in both stages,
  basic/type-7 intervals, bootstrap covariance and SEs, centered marginal tests,
  optional joint bootstrap Wald tests, smooth coefficient transformations,
  supplied multipliers and a one-draw mode for Monte Carlo studies.
* Added predict.gir(): training/newdata linear Gini indices, pointwise bootstrap
  confidence intervals and standard errors, reusable coefficient draws,
  missing-data handling and estimability checks for aliased point predictions.
* Fits retain a resolved estimation snapshot; bootstrap bases and ordering
  are cached without changing the sums used in estimation.
* Added help files and edge-case tests.
* bboot.gir() gains a `cluster` argument: one Exp(1) multiplier per cluster,
  shared by all its observations in both stages (a Bayesian cluster
  bootstrap), given as a one-sided formula over retained columns or as a
  vector aligned to the estimation rows or to the original data. The result
  reports the number of clusters and whether the first-stage cells are nested
  in the clusters, in which case only the second stage varies across draws.
* The second stage accumulates its sums over the distinct regressor rows, in
  gir() and in every bootstrap draw. Estimates are unchanged. Fits and
  bootstrap draws with cell-level regressors (industry and year effects in
  firm data, for example) are stored and refitted at the size of the distinct
  rows instead of the sample; `x = TRUE` still returns the full design.
* The automatic sieve rule now targets the total tensor-product dimension
  (the paper's K) for splines as well as Legendre polynomials; a
  knots-per-margin alternative remains available as
  sieve_size_type = "knots". Every Legendre margin has at least
  two columns, and the diagnostics report when that floor, or the spline
  minimum, exceeded the requested size.

# GIR 0.0.0.9000

* First-stage diagnostics retain named marginal dimension and interior-knot
  vectors. Basis choices are legendre and spline only.
