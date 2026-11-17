math_group_data <- function() {
  outcomes <- list(tied = c(0, 1, 1, 2), constant = rep(3, 5),
                   unequal = c(0, 0, 1, 2, 5, 9), concentrated = c(rep(0, 6), 7))
  data.frame(y = unlist(outcomes, use.names = FALSE),
             g = factor(rep(names(outcomes), lengths(outcomes)),
                        levels = names(outcomes)))
}

math_pairwise_gini <- function(y, a = rep(1, length(y))) {
  sum(outer(a, a) * abs(outer(y, y, "-"))) / (2 * sum(a) * sum(a * y))
}

test_that("saturated GIR equals the uncorrected pairwise group Gini", {
  d <- math_group_data()
  groups <- split(seq_len(nrow(d)), d$g)
  expected <- vapply(groups, function(rows) math_pairwise_gini(d$y[rows]), 0)
  expect_equal(unname(expected["tied"]), 0.375, tolerance = 1e-12)
  expect_equal(unname(expected["constant"]), 0, tolerance = 1e-12)
  expect_equal(unname(expected["concentrated"]), 6 / 7, tolerance = 1e-12)

  fit <- gir(y ~ 0 + g, data = d, first_stage = ~ g, discrete = "g",
             minimum_cell_size = 1L, singular.ok = FALSE)
  expect_equal(unname(coef(fit)), unname(expected), tolerance = 1e-12)
  expected_rows <- unname(expected[as.character(d$g)])
  expect_equal(unname(fitted(fit)), expected_rows, tolerance = 1e-12)
  expect_equal(unname(predict(fit, newdata = d)), expected_rows, tolerance = 1e-12)
  expect_true(all(expected >= 0 & expected <= 1 - 1 / lengths(groups)))
  corrected <- expected * lengths(groups) / (lengths(groups) - 1)
  expect_gt(max(abs(coef(fit) - corrected)), 0.1)

  # Treatment coding changes coefficients, but not the fitted group Ginis.
  treatment <- gir(y ~ g, data = d, first_stage = ~ g, discrete = "g",
                   minimum_cell_size = 1L, singular.ok = FALSE)
  expect_equal(unname(coef(treatment)),
               unname(c(expected[1], expected[-1] - expected[1])), tolerance = 1e-12)
  expect_equal(unname(fitted(treatment)), expected_rows, tolerance = 1e-12)
})

test_that("saturated bootstrap draws equal weighted pairwise group Ginis", {
  d <- math_group_data()
  fit <- gir(y ~ 0 + g, data = d, first_stage = ~ g, discrete = "g",
             minimum_cell_size = 1L, singular.ok = FALSE)
  W <- math_multipliers(nrow(d), R = 7L)
  bb <- bboot.gir(fit, bootstrap_weights = W, normalize_weights = FALSE,
                  inference = FALSE)
  groups <- split(seq_len(nrow(d)), d$g)
  for (b in seq_len(ncol(W))) {
    expected <- vapply(groups, function(rows) {
      math_pairwise_gini(d$y[rows], W[rows, b])
    }, 0)
    expect_equal(unname(bb$t[b, ]), unname(expected), tolerance = 1e-12)
  }
})
