empty_sample_data <- function() {
  data.frame(y = rep(c(1, 2, 3, 5), 6), x = seq_len(24),
             z = seq_len(24), a = rep(1, 24), m = rep(c(.1, .3, .6, .9), 6))
}

test_that("empty input and subset errors remain explicit and quiet", {
  d <- empty_sample_data()
  expect_gir_warnings(expect_error(gir(y ~ x, data = d[FALSE, ]),
                                   "data has no rows.", fixed = TRUE))
  for (rows in list(integer(), rep(FALSE, nrow(d)), -seq_len(nrow(d)))) {
    expect_gir_warnings(expect_error(gir(y ~ x, data = d, subset = rows),
                                     "subset selects no observations.", fixed = TRUE))
  }
})

test_that("removing every row gives one sample error before support or fitting", {
  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)
  for (action in list(na.omit, na.exclude)) {
    action_name <- if (identical(action, na.omit)) "na.omit" else "na.exclude"
    pattern <- paste0("No observations remain after ", action_name,
                       ": all 24 selected observations were removed.")
    for (column in c("y", "x", "z", "a", "m")) {
      d <- empty_sample_data()
      d[[column]][] <- NA_real_
      args <- list(formula = y ~ x, data = d, weights = d$a, na.action = action)
      if (column == "z") args$first_stage <- ~ z
      if (column == "m") args$midcdf <- d$m
      # No calling handler: under warn=2 any warning would become the error
      # and fail this match, even if it preceded the intended sample error.
      expect_error(do.call(gir, args), pattern, fixed = TRUE)
    }
    d <- empty_sample_data()
    d$y[] <- NA_real_
    expect_error(gir(y ~ x, data = d, subset = 1:7, na.action = action),
                 paste0("No observations remain after ", action_name,
                        ": all 7 selected observations were removed."), fixed = TRUE)
    expect_error(bboot.gir(y ~ x, data = d, na.action = action,
                           R = 2, inference = FALSE), pattern, fixed = TRUE)
    # The result is the same when the all-missing covariate is declared.
    d$x[] <- NA_real_
    expect_error(gir(y ~ x, data = d, continuous = "x", na.action = action),
                 pattern, fixed = TRUE)
  }
})

test_that("missing and zero weights are checked on the retained sample", {
  d <- empty_sample_data()
  d$a[] <- NA_real_
  expect_gir_warnings(expect_error(gir(y ~ x, data = d, weights = a),
                                   "Missing values are present", fixed = TRUE))
  expect_gir_warnings(expect_error(gir(y ~ x, data = d, weights = a,
                                       na.action = na.pass),
                                   "Missing values remain after na.action", fixed = TRUE))
  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)
  d$a[] <- 0
  expect_error(gir(y ~ x, data = d, weights = a),
               "positive in the retained estimation sample", fixed = TRUE)
  d$a[1] <- 1
  d$y[1] <- NA_real_
  for (action in list(na.omit, na.exclude)) {
    expect_error(gir(y ~ x, data = d, weights = a, na.action = action),
                 "positive in the retained estimation sample", fixed = TRUE)
    expect_error(gir(y ~ x, data = d, weights = a, midcdf = m, na.action = action),
                 "positive in the retained estimation sample", fixed = TRUE)
  }
  # Nonempty, positive-weight, zero-outcome data are a different problem.
  d <- empty_sample_data()
  d$y[] <- 0
  expect_error(gir(y ~ x, data = d, first_stage = ~ 1),
               "second stage has zero numerical rank", fixed = TRUE)
})

test_that("partial omissions preserve coefficients, metadata and zero case weights", {
  d <- empty_sample_data()
  d$a[2] <- 0
  d$y[3] <- NA_real_
  d$x[5] <- NA_real_
  d$z[7] <- NA_real_
  d$a[9] <- NA_real_
  selected <- c(24:13, 1:12)
  keep <- selected[complete.cases(d[selected, ])]
  for (action in list(na.omit, na.exclude)) {
    fit <- expect_gir_warnings(gir(y ~ x, data = d, subset = selected, weights = a,
      first_stage = ~ z, continuous = "z", sieve_size = 2, na.action = action),
      "4 observation(s) with missing values were removed")
    ref <- expect_gir_warnings(gir(y ~ x, data = d[keep, ], weights = a,
      first_stage = ~ z, continuous = "z", sieve_size = 2))
    expect_equal(coef(fit), coef(ref), tolerance = 1e-12)
    expect_equal(midcdf(fit), midcdf(ref), tolerance = 1e-12)
    expect_identical(fit$row.index, keep)
    expect_equal(nobs(fit), length(keep))
    expect_equal(weights(fit), d$a[keep])
    expect_true(any(weights(fit) == 0))
    expect_s3_class(fit$na.action, if (identical(action, na.exclude)) "exclude" else "omit")
  }
})

test_that("custom actions report removals and failures without inventing missingness", {
  d <- empty_sample_data()
  drop_all <- function(x) x[FALSE, , drop = FALSE]
  expect_gir_warnings(expect_error(gir(y ~ x, data = d, na.action = drop_all),
    "No observations remain after custom na.action: all 24 selected observations were removed.",
    fixed = TRUE))
  expect_gir_warnings(expect_error(gir(y ~ x, data = d,
    na.action = function(x) stop("policy failed")),
    "custom na.action failed: policy failed", fixed = TRUE))
  drop_first <- function(x) x[-1, , drop = FALSE]
  fit <- expect_gir_warnings(gir(y ~ x, data = d, first_stage = ~ 1,
    na.action = drop_first),
    "1 observation(s) were removed by custom na.action (0 with missing values, 1 complete).")
  ref <- expect_gir_warnings(gir(y ~ x, data = d[-1, ], first_stage = ~ 1))
  expect_equal(coef(fit), coef(ref), tolerance = 1e-12)
  d$y[2] <- NA_real_
  mixed <- function(x) x[-c(1, 2), , drop = FALSE]
  expect_gir_warnings(gir(y ~ x, data = d, first_stage = ~ 1, na.action = mixed),
    "2 observation(s) were removed by custom na.action (1 with missing values, 1 complete).")
  for (index in list(NA_real_, 1.5, 0, 25, Inf, "1", rep(1, 24))) {
    bad <- function(x) { x$.GIR_row <- rep(index, length.out = nrow(x)); x }
    expect_gir_warnings(expect_error(gir(y ~ x, data = d, na.action = bad),
      "na.action returned invalid row indices.", fixed = TRUE))
  }
})

test_that("support diagnostics describe the retained estimation sample", {
  d <- empty_sample_data()
  d$x <- rep(1:3, 8)
  d$y[d$x == 3] <- NA_real_
  expect_gir_warnings(expect_error(gir(y ~ x, data = d, na.action = na.omit),
    "has only 2 distinct value(s) and was not declared", fixed = TRUE))
  # A deliberate discrete declaration permits this nonempty two-point support.
  fit <- expect_gir_warnings(gir(y ~ x, data = d, discrete = "x", na.action = na.omit),
    "8 observation(s) with missing values were removed")
  ref <- expect_gir_warnings(gir(y ~ x, data = d[!is.na(d$y), ], discrete = "x"))
  expect_equal(coef(fit), coef(ref), tolerance = 1e-12)
})

test_that("empty prediction samples preserve vector and interval shapes", {
  d <- empty_sample_data()
  d$g <- factor(rep(c("a", "b"), each = 12))
  fit <- gir(y ~ g, data = d, first_stage = ~ g, discrete = "g")
  bb <- bboot.gir(fit, bootstrap_weights = math_multipliers(nrow(d)), level = 0.5)
  empty <- d[FALSE, ]
  missing <- data.frame(g = factor(rep(NA_character_, 3), levels = levels(d$g)))
  expect_gir_warnings(expect_length(predict(fit, empty), 0L))
  for (action in list(na.omit, na.exclude, na.pass)) {
    n <- if (identical(action, na.omit)) 0L else 3L
    p <- expect_gir_warnings(predict(fit, missing, na.action = action))
    expect_length(p, n)
    expect_true(all(is.na(p)))
    for (nd in list(empty, missing)) {
      expected_n <- if (nrow(nd) == 0) 0L else n
      pr <- expect_gir_warnings(predict(fit, nd, na.action = action, bootstrap = bb,
        interval = "confidence", se.fit = TRUE, level = 0.5))
      expect_equal(dim(pr$fit), c(expected_n, 3L))
      expect_identical(colnames(pr$fit), c("fit", "lwr", "upr"))
      expect_length(pr$se.fit, expected_n)
      expect_true(all(is.na(pr$fit)))
      expect_true(all(is.na(pr$se.fit)))
    }
  }
  expect_gir_warnings(expect_error(predict(fit, missing, na.action = na.fail),
                                   "missing values", fixed = TRUE))
})
