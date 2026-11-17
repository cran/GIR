test_that("Algorithm 3 uses the same conditional ranks in both losses", {
  d <- math_data()
  X <- model.matrix(~ x + z + g, d)
  for (basis in c("legendre", "spline")) {
    fit <- math_fit(d, basis)
    reference <- math_reference(d, math_basis(d, basis))
    z <- 2 * reference$midcdf - 1
    restricted <- stats::lm.wfit(matrix(1, nrow(d), 1), z, w = d$y)$coefficients[1]
    residual <- z - reference$fitted
    total <- sum(d$y * (z - restricted)^2)
    unexplained <- sum(d$y * residual^2)
    explained <- sum(d$y * (reference$fitted - restricted)^2)
    r2 <- fit$GIR.R2
    expect_gt(total, 0)
    expect_equal(r2$restricted_coefficient, unname(restricted), tolerance = 1e-10)
    expect_equal(r2$restricted_loss, total, tolerance = 1e-10)
    expect_equal(r2$unrestricted_loss, unexplained, tolerance = 1e-10)
    expect_equal(r2$value, 1 - unexplained / total, tolerance = 1e-10)
    expect_equal(r2$value, explained / total, tolerance = 1e-10)
    expect_equal(total, unexplained + explained, tolerance = 1e-10)
    expect_equal(as.numeric(crossprod(X, d$y * residual) / sum(d$y)),
                 rep(0, ncol(X)), tolerance = 1e-10)
    expect_gte(r2$value, -1e-12)
    expect_lte(r2$value, 1 + 1e-12)

    # The conditional benchmark generally differs from an unconditional Gini.
    unconditional <- sum(abs(outer(d$y, d$y, "-"))) /
      (2 * nrow(d) * sum(d$y))
    expect_gt(abs(restricted - unconditional), 1e-5)

    # Hold the nuisance estimates and sample fixed as the second stage expands.
    m <- reference$midcdf
    fits <- lapply(list(y ~ 1, y ~ x, y ~ x + z, y ~ x + z + g), function(f) {
      gir(f, data = d, midcdf = m, GIR = TRUE, singular.ok = FALSE)
    })
    values <- vapply(fits, function(f) f$GIR.R2$value, 0)
    benchmarks <- vapply(fits, function(f) f$GIR.R2$restricted_loss, 0)
    expect_equal(values[1], 0, tolerance = 1e-12)
    expect_true(all(diff(values) >= -1e-12))
    expect_equal(benchmarks, rep(total, length(fits)), tolerance = 1e-10)
    expect_equal(values[length(values)], r2$value, tolerance = 1e-10)
  }
})

test_that("zero benchmark variation gives an explicitly undefined GIR R-squared", {
  d <- data.frame(y = rep(2, 12), x = seq_len(12))
  expect_warning(fit <- gir(y ~ x, data = d, first_stage = ~ 1, GIR = TRUE),
                 "GIR-R^2 is undefined because the restricted loss is zero.",
                 fixed = TRUE)
  expect_true(is.na(fit$GIR.R2$value))
  expect_equal(fit$GIR.R2$restricted_loss, 0, tolerance = 1e-12)
  expect_equal(fit$GIR.R2$unrestricted_loss, 0, tolerance = 1e-12)
})

test_that("an unnested no-intercept model can have negative GIR R-squared", {
  d <- data.frame(y = rep(1, 12), x = seq_len(12),
                  m = 0.9 + 0.02 * sin(seq_len(12)))
  expect_warning(fit <- gir(y ~ 0 + x, data = d, midcdf = m, GIR = TRUE),
                 "constant-only benchmark is not nested", fixed = TRUE)
  z <- 2 * d$m - 1
  reference <- stats::lm.wfit(matrix(d$x, ncol = 1), z, w = d$y)
  expected <- 1 - sum(reference$residuals^2) / sum((z - mean(z))^2)
  expect_equal(fit$GIR.R2$value, expected, tolerance = 1e-10)
  expect_lt(fit$GIR.R2$value, 0)
})
