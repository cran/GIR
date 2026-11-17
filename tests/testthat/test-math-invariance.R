math_invariance_fit <- function(d, basis) {
  if (basis != "discrete") return(math_fit(d, basis))
  gir(y ~ x + z + g, data = d, first_stage = ~ g, discrete = "g",
      first_stage_singular = "error", clip = FALSE, singular.ok = FALSE,
      GIR = TRUE)
}

test_that("positive outcome rescaling preserves estimates and paired inference", {
  d <- math_data()
  W <- math_multipliers(nrow(d))
  for (basis in c("discrete", "legendre", "spline")) {
    fit <- math_invariance_fit(d, basis)
    bb <- math_boot(fit, W)
    for (scale in c(0.01, 100)) {
      scaled_data <- d
      scaled_data$y <- scale * d$y
      scaled <- math_invariance_fit(scaled_data, basis)
      scaled_bb <- math_boot(scaled, W)
      expect_equal(midcdf(scaled), midcdf(fit), tolerance = 1e-10)
      expect_equal(coef(scaled), coef(fit), tolerance = 1e-10)
      expect_equal(fitted(scaled), fitted(fit), tolerance = 1e-10)
      expect_equal(predict(scaled, newdata = d), predict(fit, newdata = d),
                   tolerance = 1e-10)
      expect_equal(scaled$GIR.R2$value, fit$GIR.R2$value, tolerance = 1e-10)
      expect_equal(scaled_bb$midcdf, bb$midcdf, tolerance = 1e-10)
      expect_equal(scaled_bb$t, bb$t, tolerance = 1e-10)
      expect_equal(scaled_bb$std.error, bb$std.error, tolerance = 1e-10)
      expect_equal(vcov(scaled_bb), vcov(bb), tolerance = 1e-10)
      expect_equal(scaled_bb$conf.int, bb$conf.int, tolerance = 1e-10)
    }
  }
})

test_that("row permutations preserve fits and bootstrap draws paired by observation", {
  d <- math_data()
  W <- math_multipliers(nrow(d))
  permutation <- order((seq_len(nrow(d)) * 137) %% nrow(d))
  undo <- order(permutation)
  for (basis in c("discrete", "legendre", "spline")) {
    fit <- math_invariance_fit(d, basis)
    permuted <- math_invariance_fit(d[permutation, ], basis)
    expect_equal(coef(permuted), coef(fit), tolerance = 1e-10)
    expect_equal(unname(midcdf(permuted)[undo]), unname(midcdf(fit)),
                 tolerance = 1e-10)
    expect_equal(unname(fitted(permuted)[undo]), unname(fitted(fit)),
                 tolerance = 1e-10)
    expect_equal(predict(permuted, newdata = d), predict(fit, newdata = d),
                 tolerance = 1e-10)
    expect_equal(permuted$GIR.R2$value, fit$GIR.R2$value, tolerance = 1e-10)
    # The random seed alone would not pair observation-level weights after a
    # permutation. Move the supplied multipliers with their observations.
    bb <- math_boot(fit, W)
    permuted_bb <- math_boot(permuted, W[permutation, ])
    expect_equal(unname(permuted_bb$midcdf[undo, ]), unname(bb$midcdf),
                 tolerance = 1e-10)
    expect_equal(permuted_bb$t, bb$t, tolerance = 1e-10)
    expect_equal(permuted_bb$std.error, bb$std.error, tolerance = 1e-10)
    expect_equal(vcov(permuted_bb), vcov(bb), tolerance = 1e-10)
    expect_equal(permuted_bb$conf.int, bb$conf.int, tolerance = 1e-10)
  }
})
