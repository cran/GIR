ridge_scaling_data <- function() {
  data.frame(y = rep(1:10, 2), x = rep(c(0, 1), each = 10))
}

ridge_scaling_fit <- function(d, basis, ridge, a = rep(1, nrow(d))) {
  expect_warning(expect_warning(
    fit <- gir(y ~ x, data = d, weights = a, continuous = "x", basis = basis,
                sieve_size = if (basis == "legendre") 3L else 4L,
                first_stage_singular = "ridge", ridge = ridge),
    "A marginal basis exceeds", fixed = TRUE),
    "First-stage singularity/ill-conditioning", fixed = TRUE)
  fit
}

ridge_scaling_basis <- function(d, basis) {
  if (basis == "spline") {
    return(splines::bs(d$x, degree = 3L, intercept = TRUE,
                       Boundary.knots = c(0, 1)))
  }
  z <- 2 * d$x - 1
  cbind(1, sqrt(3) * z, sqrt(5) * (3 * z^2 - 1) / 2)
}

ridge_scaling_reference <- function(y, P, a, ridge) {
  # Direct threshold regressions, independent of the cumulative-score code.
  H <- outer(y, y, function(u, v) (u < v) + 0.5 * (u == v))
  G <- crossprod(P, P * a)
  lambda <- ridge * mean(diag(G))
  M <- diag(P %*% solve(G + diag(lambda, ncol(P)), crossprod(P, H * a)))
  list(midcdf = M, lambda = lambda)
}

test_that("ridge bootstrap refits are invariant to common weight scales", {
  d <- ridge_scaling_data()
  X <- model.matrix(y ~ x, d)
  scales <- c(1e-12, 1e-6, 1, 1e6, 1e12)
  W <- outer(1:20, scales)
  for (basis in c("legendre", "spline")) {
    P <- ridge_scaling_basis(d, basis)
    for (ridge in c(1e-8, 1e-10, 1e-12)) {
      fit <- ridge_scaling_fit(d, basis, ridge)
      expect_warning(bb <- bboot.gir(fit, bootstrap_weights = W,
        normalize_weights = FALSE, inference = FALSE, keep_midcdf = TRUE,
        keep_weights = TRUE), "Bootstrap diagnostics:", fixed = TRUE)
      expect_warning(normalized <- bboot.gir(fit, bootstrap_weights = W,
        normalize_weights = TRUE, inference = FALSE, keep_midcdf = TRUE),
        "Bootstrap diagnostics:", fixed = TRUE)
      ref <- ridge_scaling_reference(d$y, P, 1:20, ridge)
      theta <- solve(crossprod(X, X * ((1:20) * d$y)),
                     colSums(X * ((1:20) * d$y * (2 * ref$midcdf - 1))))
      expect_identical(unname(bb$bootstrap_weights), W)
      for (b in seq_along(scales)) {
        expect_equal(unname(bb$midcdf[, b]), ref$midcdf, tolerance = 1e-10)
        expect_equal(unname(bb$t[b, ]), unname(theta), tolerance = 1e-10)
      }
      expect_equal(bb$midcdf, normalized$midcdf, tolerance = 1e-10)
      expect_equal(bb$t, normalized$t, tolerance = 1e-10)
    }
  }
})

test_that("ridge point fits and their penalties scale with observation weights", {
  d <- ridge_scaling_data()
  X <- model.matrix(y ~ x, d)
  for (basis in c("legendre", "spline")) {
    P <- ridge_scaling_basis(d, basis)
    fit <- ridge_scaling_fit(d, basis, 1e-10, 1:20)
    for (scale in c(1e-12, 1e12)) {
      a <- (1:20) * scale
      scaled <- ridge_scaling_fit(d, basis, 1e-10, a)
      ref <- ridge_scaling_reference(d$y, P, a, 1e-10)
      theta <- solve(crossprod(X, X * (a * d$y)),
                     colSums(X * (a * d$y * (2 * ref$midcdf - 1))))
      expect_equal(unname(midcdf(scaled)), ref$midcdf, tolerance = 1e-10)
      expect_equal(unname(coef(scaled)), unname(theta), tolerance = 1e-10)
      expect_equal(coef(scaled), coef(fit), tolerance = 1e-10)
      # Compare ratios so an absolute tolerance cannot hide a tiny penalty.
      expect_equal(scaled$first_stage$diagnostics$ridge_lambda / ref$lambda, 1,
                   tolerance = 1e-12)
    }
  }
})

test_that("full-rank ill-conditioned first stages also use relative ridge", {
  P <- cbind(1, seq(-1, 1, length.out = 12) * 1e-4)
  reference <- ridge_scaling_reference(1:12, P, 1:12, 1e-8)
  for (scale in c(1e-12, 1, 1e12)) {
    fit <- GIR:::.gir_fit_cell(1:12, P, (1:12) * scale, "ill-conditioned",
                               "ridge", 1e-7, 1e-8)
    expect_true(fit$ok)
    expect_equal(fit$rank, 2L)
    expect_identical(fit$action, "ridge")
    expect_equal(fit$fitted, reference$midcdf, tolerance = 1e-10)
    expect_equal(fit$ridge_lambda / (reference$lambda * scale), 1,
                 tolerance = 1e-12)
  }
})

test_that("ridge retries scale the increment and report the actual total", {
  # The initial increment rounds away; the larger retry allows a solve.
  P <- cbind(rep(1, 8), rep(1, 8))
  for (scale in c(1e-12, 1, 1e12)) {
    a <- (1:8) * scale
    fit <- GIR:::.gir_fit_cell(1:8, P, a, "retry", "ridge", 1e-12, 1e-17)
    expected <- ridge_scaling_reference(1:8, P, a, 101e-17)
    expect_true(fit$ok)
    expect_identical(fit$action, "ridge")
    expect_equal(fit$ridge_lambda / expected$lambda, 1, tolerance = 1e-12)
    expect_equal(fit$fitted, expected$midcdf, tolerance = 1e-10)
  }

  # A stricter solve threshold can trigger a retry without an initial penalty.
  P <- cbind(1, seq(-1, 1, length.out = 12) * 1e-9)
  for (scale in c(1e-12, 1, 1e12)) {
    a <- (1:12) * scale
    fit <- GIR:::.gir_fit_cell(1:12, P, a, "late retry", "ridge", 1e-25, 1e-10)
    expected <- ridge_scaling_reference(1:12, P, a, 100e-10)
    expect_false(fit$singular)
    expect_identical(fit$action, "ridge")
    expect_equal(fit$ridge_lambda / expected$lambda, 1, tolerance = 1e-12)
    expect_equal(fit$fitted, expected$midcdf, tolerance = 1e-10)
  }
})

test_that("unusable ridge scales, penalties and adjusted matrices fail clearly", {
  adjust <- GIR:::.gir_ridge_adjustment
  for (scale in c(0, Inf, NA_real_)) {
    G <- diag(scale, 2)
    expect_error(adjust(G, G, 1e-10, "invalid scale"),
                 "average Gram diagonal must be finite and positive", fixed = TRUE)
  }
  tiny <- diag(1e-320, 2)
  expect_error(adjust(tiny, tiny, 1e-10, "underflow"),
               "penalty underflowed or overflowed", fixed = TRUE)
  G <- diag(2, 2)
  expect_error(adjust(G, G, .Machine$double.xmax, "overflow"),
               "penalty underflowed or overflowed", fixed = TRUE)
  expect_error(adjust(G, G, .Machine$double.xmax / 10, "retry overflow", 100),
               "penalty underflowed or overflowed", fixed = TRUE)
  huge <- diag(1e308, 2)
  expect_error(adjust(huge, huge, 1, "matrix overflow"),
               "Ridge-adjusted Gram matrix contains non-finite values", fixed = TRUE)
})

test_that("ridge remains opt-in and does not change well-conditioned fits", {
  expect_identical(formals(gir)$ridge, 1e-10)
  d <- ridge_scaling_data()
  expect_error(gir(y ~ x, d, continuous = "x", sieve_size = 3),
               "has 2 distinct value(s)", fixed = TRUE)
  P <- cbind(1, seq(-1, 1, length.out = 20))
  plain <- GIR:::.gir_fit_cell(1:20, P, 1:20, "ordinary", "none", 1e-12, 1e-10)
  for (method in c("error", "ridge")) {
    fit <- GIR:::.gir_fit_cell(1:20, P, 1:20, "ordinary", method, 1e-12, 1e-10)
    expect_identical(fit, plain)
  }
  expect_identical(plain$action, "none")
  expect_identical(plain$ridge_lambda, 0)
})
