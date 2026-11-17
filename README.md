# GIR: Gini Index Regression

Linear Gini Index Regression (GIR) in R. `gir()` fits the two-stage
estimator of Bhardwaj and Diaz-Villarejo (2026): a least-squares series
first stage for the conditional mid-distribution function and a weighted
linear projection for the conditional Gini index. `bboot.gir()` provides
Bayesian-bootstrap inference (basic confidence intervals, standard errors,
covariance, tests of smooth coefficient transformations) and `predict.gir()`
evaluates the fitted conditional Gini index with pointwise uncertainty.

## Installation

Before the first CRAN release, installation from GitHub requires access to
the private repository:

```r
# install.packages("remotes")
remotes::install_github("uptoaparameter/GIR", build_vignettes = TRUE)
```

Building the guide from GitHub requires LaTeX. A source package produced by
`R CMD build` already contains the guide PDF and can be installed without
LaTeX or roxygen2. To install from GitHub without building the guide, omit
`build_vignettes = TRUE`.

Use `?gir`, `?bboot.gir` and `?predict.gir` for function help. Open the detailed
guide with `vignette("GIR-guide", package = "GIR")` and retrieve the software
and paper references with `citation("GIR")`.

## Example

```r
library(GIR)
set.seed(150)
n <- 2000L
dat <- data.frame(x1 = runif(n, 0.1, 1),
                  x2 = rbinom(n, 1, 0.5), x3 = rbinom(n, 1, 0.5))
m <- with(dat, 0.075*x1 + 0.100*x2 + 0.125*x3)
dat$y <- rbeta(n, (1-m)/(2*m), 1)

fit <- gir(y ~ 0 + x1 + x2 + x3, data = dat,
           continuous = "x1", discrete = c("x2", "x3"),
           basis = "spline", basis_bounds = list(x1 = c(0.1, 1)))
summary(fit)

boot <- bboot.gir(fit, R = 999, seed = 42)
summary(boot)
confint(boot)
vcov(boot)
summary(boot, transform = c(x2 = 1, x3 = -1), null = 0)

nd <- data.frame(x1 = c(0.2, 0.8), x2 = 1, x3 = 0)
predict(fit, nd, interval = "confidence", se.fit = TRUE, bootstrap = boot)
```

The first stage is estimated separately in every cell formed by the
`discrete` variables, with a tensor-product Legendre or spline basis in the
`continuous` ones; the user declares which variables are which. The bootstrap
refits both stages with the same normalized exponential multipliers and keeps
the estimation sample, cells, basis sizes and bounds fixed. `t0` and `t`
contain the untransformed coefficients and draws, which tests and predictions
reuse.

For clustered multipliers, `cluster = ~ industry + region` creates one cluster
for each distinct industry-region combination. This is one-way clustering at
the intersections, not two-way clustering by the two variables separately.
Clusters are identified by the values of their defining variables and numbered
in order of first occurrence. The paper's bootstrap theorem assumes
independent observations; clustered inference requires separate justification.

For a cluster vector, use `cluster_order = "original"` when its labels follow
the original dataset's order at fitting time. The stored `fit$row.index`
aligns it to the estimation sample. Use `cluster_order = "estimation"` for
labels already in final estimation-row order. The default `"auto"` infers
the order from length, but requires an explicit choice when both lengths
match and the row mapping has changed (for example, a reordered `subset`).
Formulas such as `cluster = ~ industry` already use aligned retained columns;
leave `cluster_order = "auto"` for formulas or when no clustering is requested.
These controls also work in `predict(..., se.fit = TRUE)` when it generates
the bootstrap; the cluster labels describe the fitting sample, not `newdata`.

`predict()` estimates the conditional **Gini index**, not individual outcomes.
Its intervals are pointwise basic intervals, not simultaneous bands. For an
externally supplied `midcdf`, a refitting callback or an explicit fixed
first-stage choice is required for bootstrap inference.

## Missing data and diagnostics

`na.omit` and `na.exclude` remove rows with missing values in either stage,
the case weights or supplied mid-CDF values. Fitting requires at least one
retained observation with positive weight. An empty sample or loss of all
positive weights produces a specific error before the routine removal
warning. Missing weights follow the selected `na.action` policy.

Support checks use the retained sample: if omission leaves an undeclared
numeric first-stage variable with at most two distinct values, declare it
explicitly in `discrete` or `continuous`. Custom missing-value actions report
how many removed rows had missing values and how many were complete; their
own failures are identified without assuming missing data caused them.

Prediction preserves empty outputs: no prediction rows gives an empty vector
or a zero-row confidence-interval matrix. `na.exclude` restores missing
positions as `NA`. Numerical and inference warnings identify issues that need
attention before interpreting a fit.

## Numerical stabilization

The default `first_stage_singular = "error"` stops on a singular or
ill-conditioned first stage. If you explicitly select `"ridge"`, the default
`ridge = 1e-10` is a numerical stabilization heuristic. The initial diagonal
penalty is `ridge * mean(diag(G))`, where `G` is the weighted cell Gram matrix;
it scales with the weights. A failed solve receives an additional 100 times
that base penalty. Inspect `fit$first_stage$diagnostics$ridge_lambda` for the
total amount used. This parameter is not selected by cross-validation or
prescribed by the paper's theory, and its effect can be substantial in poorly
identified directions. Assess sensitivity when using ridge. Common weight
scaling preserves internal fits up to floating-point error within usable
numerical ranges; external first-stage callbacks need their own guarantees.

## Mathematical validation

The test suite checks both stages against independent threshold regressions,
saturated group estimates against the pairwise empirical Gini, outcome-unit
and row-order invariance, and the GIR-R-squared loss decomposition. It also
checks paired bootstrap draws and the arithmetic of basic confidence intervals.
The mathematical tests use full-rank designs without ridge or clipping;
cluster and ridge behavior have separate tests.

Run `testthat::test_local(filter = "math-")` from the source directory for
these checks, or `testthat::test_local()` for the complete suite. The tests
also run during `R CMD check` and need neither the draft paper nor network
access. See `help("GIR-validation", package = "GIR")` for the correspondence
with the paper's lemmas and algorithms, numerical tolerances and limitations.
These finite-sample checks do not establish asymptotic validity or nominal
confidence-interval coverage; those require separate theoretical and
simulation assessment.

## Development

The help files are generated from the roxygen comments in `R/`; after
editing them run `roxygen2::roxygenize()` in the package directory.
Run `testthat::test_local()` for the tests. Release checks also build and check
the PDF guide and function reference. Maintainer build instructions are in
`tools/README.md` in the source repository.
The guide's editable source is `vignettes/GIR-guide.Rnw`.

## Reference

Bhardwaj, V. and Diaz-Villarejo, D. (2026).
[Gini Index Regression](https://ddiazvil.com/assets/papers/JMP_GIR_DiazVillarejo.pdf). Working
paper, Universidad Carlos III de Madrid.
