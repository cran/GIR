test_that("supplied mid-CDF reproduces the direct second-stage formula", {
  set.seed(11)
  n <- 100
  dat <- data.frame(
    y = runif(n, 0.1, 2),
    x = rnorm(n),
    d = rbinom(n, 1, 0.5),
    m = runif(n)
  )
  X <- model.matrix(y ~ x + d, dat)
  expected <- solve(
    crossprod(X, X * dat$y) / n,
    colMeans(X * (dat$y * (2 * dat$m - 1)))
  )
  fit <- gir(y ~ x + d, data = dat, midcdf = m, singular.ok = FALSE)
  expect_equal(unname(coef(fit)), unname(expected), tolerance = 1e-12)
  expect_equal(unname(midcdf(fit)), dat$m)
})

test_that("ties use the mid-distribution response", {
  dat <- data.frame(y = c(1, 1, 2, 2), group = factor(c("a", "a", "a", "a")))
  fit <- gir(
    y ~ 1,
    data = dat,
    first_stage = ~ group,
    discrete = "group",
    minimum_cell_size = 1,
    singular.ok = FALSE
  )
  expect_equal(unname(midcdf(fit)), c(0.25, 0.25, 0.75, 0.75),
               tolerance = 1e-12)
})

test_that("first-stage diagnostics retain marginal dimensions and knot counts", {
  set.seed(19)
  n <- 100
  dat <- data.frame(
    y = runif(n, 0.1, 2),
    x1 = runif(n),
    x2 = runif(n)
  )

  fit <- expect_gir_warnings(gir(
    y ~ x1 + x2,
    data = dat,
    first_stage = ~ x1 + x2,
    continuous = c("x1", "x2"),
    basis = "spline",
    sieve_size_type = "knots",
    sieve_size = c(2L, 1L),
    first_stage_singular = "ridge",
    singular.ok = FALSE
  ),
    "2 mid-CDF value(s) from series estimation are outside [0,1]")

  diagnostics <- fit$first_stage$diagnostics
  expect_true(all(c(
    "marginal_dimensions",
    "marginal_knot_counts",
    "basis_dimension"
  ) %in% names(diagnostics)))
  expect_equal(diagnostics$marginal_dimensions[[1L]], c(x1 = 6L, x2 = 5L))
  expect_equal(diagnostics$marginal_knot_counts[[1L]], c(x1 = 2L, x2 = 1L))
  expect_equal(diagnostics$basis_dimension[1L], 30L)
})

test_that("missing weights follow na.action like any other variable", {
  set.seed(31)
  d <- data.frame(y = runif(60, 0.1, 2), x = runif(60), a = runif(60, 1, 2))
  d$a[5] <- NA_real_
  expect_error(gir(y ~ x, data = d, weights = a, continuous = "x",
                   sieve_size = 2), "Missing values")
  fit <- expect_gir_warnings(gir(y ~ x, data = d, weights = a,
                              na.action = na.omit, continuous = "x",
                              sieve_size = 2),
    "1 observation(s) with missing values were removed by na.omit.")
  expect_equal(nobs(fit), 59L)
  expect_false(5L %in% fit$row.index)
  ref <- expect_gir_warnings(gir(y ~ x, data = d[-5, ], weights = a,
                              continuous = "x", sieve_size = 2))
  expect_equal(coef(fit), coef(ref), tolerance = 1e-12)
  expect_equal(unname(midcdf(fit)), unname(midcdf(ref)), tolerance = 1e-12)
  excl <- expect_gir_warnings(gir(y ~ x, data = d, weights = a,
                               na.action = na.exclude, continuous = "x",
                               sieve_size = 2),
    "1 observation(s) with missing values were removed by na.exclude.")
  expect_s3_class(excl$na.action, "exclude")
  expect_equal(nobs(excl), 59L)
  d$a[5] <- Inf
  expect_error(gir(y ~ x, data = d, weights = a, na.action = na.omit,
                   continuous = "x", sieve_size = 2), "infinite")
})

test_that("defaulted numeric first-stage columns are reported by a message", {
  set.seed(41)
  d <- data.frame(y = runif(80, 0.1, 2), x = runif(80), z = runif(80))
  expect_message(expect_gir_warnings(gir(y ~ x + z, data = d, sieve_size = 4)),
                 "x, z were not declared")
  expect_no_message(expect_gir_warnings(gir(y ~ x + z, data = d,
    continuous = c("x", "z"), sieve_size = 4)))
  expect_no_message(expect_gir_warnings(gir(y ~ x, data = d, first_stage = ~ z,
    continuous = "z", sieve_size = 3)))
})

test_that("an undeclared two-valued numeric column is refused, a declared one accepted", {
  set.seed(42)
  d <- data.frame(y = runif(80, 0.1, 2), x = runif(80), b = rbinom(80, 1, 0.5))
  expect_error(gir(y ~ x + b, data = d, sieve_size = 4), "distinct value")
  fit <- expect_gir_warnings(gir(y ~ x + b, data = d, continuous = "x",
                              discrete = "b", sieve_size = 3),
    "1 mid-CDF value(s) from series estimation are outside [0,1]")
  expect_equal(nrow(fit$first_stage$diagnostics), 2L)
  fit2 <- expect_gir_warnings(gir(y ~ x + b, data = d,
                               continuous = c("x", "b"), sieve_size = 4))
  expect_equal(unname(fit2$first_stage$diagnostics$marginal_dimensions[[1L]]),
               c(2L, 2L))
})

test_that("a marginal basis larger than the variable's support is diagnosed", {
  set.seed(44)
  d <- data.frame(y = runif(120, 0.1, 2), s = sample(1:3, 120, TRUE))
  expect_error(suppressMessages(gir(y ~ s, data = d, continuous = "s",
                                    sieve_size = 5)),
               "3 distinct value")
  fit <- expect_gir_warnings(gir(y ~ s, data = d, continuous = "s",
    sieve_size = 5, first_stage_singular = "drop"),
    c("singular by construction and was handled with 'drop': s has 3 distinct",
      "First-stage singularity/ill-conditioning was handled with 'drop'"))
  dg <- fit$first_stage$diagnostics
  expect_equal(dg$basis_dimension, 5L)
  expect_lt(dg$effective_dimension, dg$basis_dimension)
  expect_true(all(is.finite(midcdf(fit))))
})

test_that("the size rule targets the total tensor dimension for both bases", {
  set.seed(43)
  n <- 500
  d <- data.frame(y = runif(n, 0.1, 2), x1 = runif(n), x2 = runif(n),
                  x3 = runif(n))
  # One continuous variable: floor(500^0.3) = 6 columns, so 2 interior knots.
  sp <- expect_gir_warnings(gir(y ~ x1, data = d, continuous = "x1",
                             basis = "spline", first_stage_singular = "ridge"),
    "2 mid-CDF value(s) from series estimation are outside [0,1]")
  dg <- sp$first_stage$diagnostics
  expect_equal(dg$size_type, "dimension")
  expect_equal(unname(dg$marginal_dimensions[[1L]]), 6L)
  expect_equal(unname(dg$marginal_knot_counts[[1L]]), 2L)
  expect_false(dg$floored)
  # The knots-per-margin alternative: 6 knots, 10 columns.
  kn <- expect_gir_warnings(gir(y ~ x1, data = d, continuous = "x1",
                             basis = "spline", sieve_size_type = "knots",
                             first_stage_singular = "ridge"),
    "2 mid-CDF value(s) from series estimation are outside [0,1]")
  expect_equal(unname(kn$first_stage$diagnostics$marginal_dimensions[[1L]]), 10L)
  expect_equal(unname(kn$first_stage$diagnostics$marginal_knot_counts[[1L]]), 6L)
  # Two continuous variables: target 6 is below the cubic minimum 4 x 4 = 16.
  sp2 <- expect_gir_warnings(gir(y ~ x1 + x2, data = d,
    continuous = c("x1", "x2"), basis = "spline", first_stage_singular = "ridge"),
    "2 mid-CDF value(s) from series estimation are outside [0,1]")
  dg2 <- sp2$first_stage$diagnostics
  expect_equal(unname(dg2$marginal_dimensions[[1L]]), c(4L, 4L))
  expect_equal(unname(dg2$marginal_knot_counts[[1L]]), c(0L, 0L))
  expect_true(dg2$floored)
  expect_false(dg2$capped)
  # Legendre never drops a variable: target 6 over three margins gives (2, 2, 2).
  lg <- expect_gir_warnings(gir(y ~ x1 + x2 + x3, data = d,
    continuous = c("x1", "x2", "x3"), sieve_size = 6),
    "1 mid-CDF value(s) from series estimation are outside [0,1]")
  dg3 <- lg$first_stage$diagnostics
  expect_equal(unname(dg3$marginal_dimensions[[1L]]), c(2L, 2L, 2L))
  expect_equal(dg3$basis_dimension, 8L)
  expect_true(dg3$floored)
  lg12 <- expect_gir_warnings(gir(y ~ x1 + x2 + x3, data = d,
    continuous = c("x1", "x2", "x3"), sieve_size = 12))
  expect_equal(unname(lg12$first_stage$diagnostics$marginal_dimensions[[1L]]),
               c(3L, 2L, 2L))
  expect_false(lg12$first_stage$diagnostics$floored)
  # A cell too small for a linear term in every variable is an error.
  expect_gir_warnings(expect_error(gir(y ~ x1 + x2 + x3, data = d[1:9, ],
    continuous = c("x1", "x2", "x3"), minimum_cell_size = 1L), "too small"))
})

# ---- second-stage sums over distinct regressor rows ------------------------
cell_gir_data <- function(seed = 21) {
  set.seed(seed)
  cells <- expand.grid(ind = factor(c("a", "b", "c", "d")),
                       year = factor(c("t1", "t2")))
  cells$dose <- c(0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4) *
    (cells$year == "t2")
  size <- rep(c(20, 30, 40, 50), 2)
  rows <- rep(seq_len(nrow(cells)), size)
  dat <- cells[rows, ]
  dat$y <- rgamma(nrow(dat), shape = 2, rate = 1)
  rownames(dat) <- NULL
  dat
}
midrank_gini <- function(y) {
  y <- sort(y); n <- length(y)
  m <- (findInterval(y, y, left.open = TRUE) + findInterval(y, y)) / (2 * n)
  sum(y * (2 * m - 1)) / sum(y)
}

test_that("second-stage sums over distinct rows reproduce the observation-level formulas", {
  dat <- cell_gir_data()
  fit <- gir(y ~ ind + year + dose, data = dat, first_stage = ~ ind + year,
             discrete = c("ind", "year"), x = TRUE, GIR = TRUE)
  X <- model.matrix(y ~ ind + year + dose, dat)
  m <- midcdf(fit)
  expected <- solve(crossprod(X, X * dat$y) / nrow(X),
                    colMeans(X * (dat$y * (2 * m - 1))))
  expect_equal(unname(coef(fit)), unname(expected), tolerance = 1e-12)
  expect_equal(nrow(fit$estimation$X), 8L)
  expect_equal(length(fit$estimation$group), nrow(dat))
  expect_equal(fit$estimation$X[fit$estimation$group, ], unname(X)[, ],
               ignore_attr = TRUE)
  expect_equal(unname(fit$x), unname(X))
  expect_equal(attr(fit$x, "assign"), attr(X, "assign"))
  expect_equal(unname(fitted(fit)), unname(drop(X %*% expected)), tolerance = 1e-12)
  expect_equal(unname(residuals(fit)), unname(2 * m - 1 - drop(X %*% expected)),
               tolerance = 1e-12)
  expect_equal(unname(model.matrix(fit)), unname(X))
  # GIR-R^2 from observation-level losses
  z <- 2 * m - 1; f <- drop(X %*% expected)
  expect_equal(fit$GIR.R2$value,
               1 - sum(dat$y * (z - f)^2) /
                 sum(dat$y * (z - sum(dat$y * z) / sum(dat$y))^2),
               tolerance = 1e-12)
})

test_that("a design with all rows distinct keeps one row per observation", {
  set.seed(4)
  dat <- data.frame(y = runif(60, 0.5, 2), x = rnorm(60), m = runif(60))
  fit <- gir(y ~ x, data = dat, midcdf = m)
  expect_equal(nrow(fit$estimation$X), 60L)
  expect_equal(fit$estimation$group, seq_len(60))
  expect_equal(fit$estimation$version, 2L)
})

test_that("cell-level regressors give the weighted regression of cell Ginis", {
  dat <- cell_gir_data()
  fit <- gir(y ~ ind + year + dose, data = dat, first_stage = ~ ind + year,
             discrete = c("ind", "year"))
  cells <- aggregate(y ~ ind + year + dose, dat,
                     function(v) c(G = midrank_gini(v), W = sum(v)))
  cells <- data.frame(cells[c("ind", "year", "dose")],
                      G = cells$y[, "G"], W = cells$y[, "W"])
  wls <- lm(G ~ ind + year + dose, data = cells, weights = W)
  expect_equal(unname(coef(fit)), unname(coef(wls)), tolerance = 1e-10)
})

test_that("an aliased second stage with duplicated rows keeps the weighted fit", {
  dat <- cell_gir_data()
  dat$dose2 <- 2 * dat$dose
  fit <- expect_gir_warnings(gir(y ~ ind + year + dose + dose2, data = dat,
                              first_stage = ~ ind + year,
                              discrete = c("ind", "year")),
    "The second-stage design is rank deficient; 1 coefficient(s) are aliased")
  expect_true(any(fit$aliased))
  X <- model.matrix(y ~ ind + year + dose + dose2, dat)
  z <- 2 * midcdf(fit) - 1
  reference <- lm(z ~ X - 1, weights = dat$y)
  expect_equal(unname(fitted(fit)), unname(fitted(reference)), tolerance = 1e-10)
})

test_that("predictions on the estimation sample expand the distinct rows", {
  dat <- cell_gir_data()
  fit <- gir(y ~ ind + year + dose, data = dat, first_stage = ~ ind + year,
             discrete = c("ind", "year"))
  expect_equal(predict(fit), fitted(fit))
  expect_warning(bb <- bboot.gir(fit, R = 50L, seed = 3L),
                 "Fewer than ten expected draws per confidence-interval tail.",
                 fixed = TRUE)
  from_sample <- expect_gir_warnings(predict(fit, se.fit = TRUE,
                                          interval = "confidence", bootstrap = bb),
    "Few bootstrap draws per prediction-interval tail")
  from_newdata <- expect_gir_warnings(predict(fit, newdata = dat, se.fit = TRUE,
                                           interval = "confidence", bootstrap = bb),
    "Few bootstrap draws per prediction-interval tail")
  expect_equal(unname(from_sample$fit), unname(from_newdata$fit), tolerance = 1e-12)
  expect_equal(unname(from_sample$se.fit), unname(from_newdata$se.fit),
               tolerance = 1e-12)
  expect_equal(names(from_sample$se.fit), names(fitted(fit)))
})

test_that("a formula with data-dependent terms groups exactly", {
  dat <- cell_gir_data()
  fit <- gir(y ~ ind + poly(dose, 2), data = dat, first_stage = ~ ind + year,
             discrete = c("ind", "year"))
  X <- model.matrix(y ~ ind + poly(dose, 2), dat)
  m <- midcdf(fit)
  expected <- solve(crossprod(X, X * dat$y), colSums(X * (dat$y * (2 * m - 1))))
  expect_equal(unname(coef(fit)), unname(expected), tolerance = 1e-10)
  # QR rounding in poly() can differ across platforms, so the number of
  # bitwise-distinct rows is not fixed. Group the actual numeric rows and
  # verify that expanding the stored design recovers every value exactly.
  rows <- lapply(seq_len(nrow(X)), function(i) unname(X[i, ]))
  expect_equal(nrow(fit$estimation$X), length(unique(rows)))
  expanded <- fit$estimation$X[fit$estimation$group, , drop = FALSE]
  expect_identical(as.vector(expanded), as.vector(X))
  expect_lt(nrow(fit$estimation$X), nrow(dat))
})
