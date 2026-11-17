test_that("Algorithms 1 and 2 agree with independent tensor threshold regressions", {
  d <- math_data()
  W <- math_multipliers(nrow(d))
  for (basis in c("legendre", "spline")) {
    P <- math_basis(d, basis)
    # These tests exercise unregularized, full-rank fits in every cell.
    for (rows in split(seq_len(nrow(d)), d$g)) {
      expect_equal(qr(P[rows, ])$rank, ncol(P))
      expect_gt(rcond(crossprod(P[rows, ])), 1e-7)
    }
    # Unit weights are Algorithm 1; unequal weights check the package extension.
    for (a in list(rep(1, nrow(d)), 0.7 + (seq_len(nrow(d)) %% 7) / 10)) {
      fit <- math_fit(d, basis, a)
      reference <- math_reference(d, P, a)
      expect_equal(unname(midcdf(fit)), reference$midcdf, tolerance = 1e-10)
      expect_equal(coef(fit), reference$coefficients, tolerance = 1e-10)
      expect_equal(unname(fitted(fit)), unname(reference$fitted), tolerance = 1e-10)
      expect_true(all(fit$first_stage$diagnostics$ridge_lambda == 0))

      if (length(unique(a)) == 1L) {
        bb <- math_boot(fit, W)
      } else {
        expect_warning(bb <- math_boot(fit, W),
                       "Observation weights are multiplied", fixed = TRUE)
      }
      expected <- t(vapply(seq_len(ncol(W)), function(b) {
        draw <- math_reference(d, P, a * W[, b])
        expect_equal(unname(bb$midcdf[, b]), draw$midcdf, tolerance = 1e-10)
        draw$coefficients
      }, coef(fit)))
      expect_equal(unname(bb$t), unname(expected), tolerance = 1e-10)
      expect_equal(unname(bb$std.error), unname(apply(expected, 2, sd)),
                   tolerance = 1e-10)
      expect_equal(unname(vcov(bb)), unname(cov(expected)), tolerance = 1e-10)
      roots <- sweep(expected, 2, reference$coefficients, "-")
      q <- apply(roots, 2, quantile, probs = c(0.25, 0.75), type = 7)
      ci <- cbind(reference$coefficients - q[2, ], reference$coefficients - q[1, ])
      expect_equal(unname(bb$conf.int), unname(ci), tolerance = 1e-10)
      # Ensure a fixed-first-stage implementation could not pass this fixture.
      expect_gt(max(abs(bb$midcdf[, 1] - reference$midcdf)), 1e-4)
      X <- model.matrix(~ x + z + g, d)
      fixed <- stats::lm.wfit(X, 2 * reference$midcdf - 1,
                             w = a * W[, 1] * d$y)$coefficients
      expect_gt(max(abs(expected[1, ] - fixed)), 1e-5)
      # Basic and percentile intervals must be distinguishable on this fixture.
      percentile <- t(apply(expected, 2, quantile, probs = c(0.25, 0.75)))
      expect_gt(max(abs(ci - percentile)), 1e-5)
    }
  }
})
