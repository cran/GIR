# Independent small-matrix tests; no helper here calls the package first stage.
small_gir_data <- function() {
  data.frame(y = c(0, 1, 1, 3, 1, 2, 2, 4),
             d = factor(rep(c("a", "b"), each = 4)))
}
small_weights <- function() {
  W <- cbind(1:8, c(2, 1, 4, 3, 5, 8, 6, 7), c(8, 2, 6, 4, 7, 1, 5, 3))
  sweep(W, 2, colMeans(W), "/")
}
naive_mid <- function(y, group, a) {
  vapply(seq_along(y), function(i) {
    rows <- group == group[i]
    sum(a[rows] * ((y[rows] < y[i]) + 0.5 * (y[rows] == y[i]))) / sum(a[rows])
  }, numeric(1))
}
naive_theta <- function(X, y, mid, a) {
  solve(crossprod(X, X * (a * y)), colSums(X * (a * y * (2 * mid - 1))))
}
small_fit <- function() {
  gir(y ~ d, data = small_gir_data(), first_stage = ~ d,
      discrete = "d", minimum_cell_size = 1L, singular.ok = FALSE)
}

test_that("both bootstrap stages equal independent weighted midrank moments", {
  d <- small_gir_data(); fit <- small_fit(); W <- small_weights()
  bb <- expect_gir_warnings(bboot.gir(fit, bootstrap_weights = W,
    normalize_weights = FALSE, keep_midcdf = TRUE, keep_weights = TRUE),
    "Fewer than ten expected draws per confidence-interval tail.")
  X <- model.matrix(y ~ d, d)
  expected <- matrix(NA_real_, ncol(W), ncol(X))
  for (b in seq_len(ncol(W))) {
    m <- naive_mid(d$y, d$d, W[, b])
    expected[b, ] <- naive_theta(X, d$y, m, W[, b])
    expect_equal(unname(bb$midcdf[, b]), m, tolerance = 1e-12)
  }
  expect_equal(unname(bb$t), expected, tolerance = 1e-12)
  expect_equal(unname(bb$std.error), apply(expected, 2, sd), tolerance = 1e-12)
  expect_equal(unname(vcov(bb)), cov(expected), tolerance = 1e-12)
  expect_equal(unname(bb$bootstrap_weights), unname(W))
  expect_true(any(abs(bb$midcdf[, 1] - midcdf(fit)) > 1e-8))
  delta <- sweep(expected, 2, coef(fit), "-")
  lo <- apply(delta, 2, quantile, probs = 0.025, names = FALSE, type = 7)
  hi <- apply(delta, 2, quantile, probs = 0.975, names = FALSE, type = 7)
  expect_equal(unname(bb$conf.int), unname(cbind(coef(fit) - hi, coef(fit) - lo)),
               tolerance = 1e-12)
})

test_that("existing observation weights enter both stages exactly once", {
  d <- small_gir_data(); d$a <- c(0, 1, 2, 1, 2, 1, 3, 1)
  fit <- gir(y ~ d, data = d, weights = a, discrete = "d",
              minimum_cell_size = 1L, singular.ok = FALSE)
  W <- small_weights()
  bb <- expect_gir_warnings(bboot.gir(fit, bootstrap_weights = W,
                  normalize_weights = FALSE, keep_midcdf = TRUE),
    c("Observation weights are multiplied by iid bootstrap multipliers.",
      "Fewer than ten expected draws per confidence-interval tail."))
  X <- model.matrix(y ~ d, d)
  for (b in seq_len(ncol(W))) {
    a <- d$a * W[, b]; m <- naive_mid(d$y, d$d, a)
    expect_equal(unname(bb$midcdf[, b]), m, tolerance = 1e-12)
    expect_equal(unname(bb$t[b, ]), as.numeric(naive_theta(X, d$y, m, a)),
                 tolerance = 1e-12)
  }
})

test_that("named contrasts and centered hypothesis tests reuse raw draws", {
  fit <- small_fit()
  bb <- expect_gir_warnings(bboot.gir(fit, bootstrap_weights = small_weights(),
                                    normalize_weights = FALSE),
    "Fewer than ten expected draws per confidence-interval tail.")
  s <- expect_gir_warnings(summary(bb, transform = c(db = 1), null = 0.1),
    "Fewer than ten expected draws per confidence-interval tail.")
  roots <- bb$t[, "db"] - coef(fit)["db"]
  observed <- coef(fit)["db"] - 0.1
  p <- min(1, 2 * min((1 + sum(roots <= observed)) / 4,
                      (1 + sum(roots >= observed)) / 4))
  expect_equal(unname(s$test$p.value), p)
  expect_equal(unname(s$test$std.error), sd(bb$t[, "db"]), tolerance = 1e-12)
  f <- function(b) c(slope = unname(exp(b["db"])))
  sf <- expect_gir_warnings(summary(bb, transform = f, null = 1),
    "Fewer than ten expected draws per confidence-interval tail.")
  expect_equal(unname(sf$test$std.error), sd(exp(bb$t[, "db"])), tolerance = 1e-12)
  expect_equal(coef(bb), coef(fit))
  expect_warning(
    expect_error(summary(bb, transform = function(b) NA_real_), "finite"),
    "Fewer than ten expected draws per confidence-interval tail.", fixed = TRUE)
})

test_that("R=1 is explicitly draw-only and supplied multipliers consume no RNG", {
  fit <- small_fit(); w <- small_weights()[, 1]
  expect_error(bboot.gir(fit, R = 1, bootstrap_weights = w), "R=1")
  set.seed(10); before <- .Random.seed
  b <- bboot.gir(fit, R = 1, inference = FALSE,
                 bootstrap_weights = w, normalize_weights = FALSE)
  expect_identical(.Random.seed, before)
  expect_null(b$std.error)
  expect_null(b$conf.int)
  expect_error(summary(b), "single draw")
  expect_error(confint(b), "two successful draws")
  # A common rescaling of the multipliers leaves the draw unchanged.
  scaled <- bboot.gir(fit, R = 1, inference = FALSE,
                      bootstrap_weights = w * 2, normalize_weights = FALSE)
  expect_equal(scaled$t, b$t, tolerance = 1e-12)
})

test_that("local seeds restore state and NULL advances the normal stream", {
  fit <- small_fit()
  set.seed(15); before <- .Random.seed
  a <- bboot.gir(fit, R = 2, seed = 6, inference = FALSE, keep_weights = TRUE)
  expect_identical(.Random.seed, before)
  b <- bboot.gir(fit, R = 2, seed = 6, inference = FALSE)
  expect_equal(a$t, b$t)
  c <- bboot.gir(fit, R = 2, inference = FALSE, keep_weights = TRUE)
  after <- .Random.seed
  assign(".Random.seed", before, envir = .GlobalEnv)
  e1 <- rexp(8); e2 <- rexp(8)
  expect_identical(.Random.seed, after)
  expect_equal(unname(c$bootstrap_weights), cbind(e1 / mean(e1), e2 / mean(e2)))
})

test_that("supplied nuisances are never silently treated as estimated refits", {
  d <- small_gir_data(); d$m <- naive_mid(d$y, d$d, rep(1, 8))
  fit <- gir(y ~ d, data = d, midcdf = m)
  expect_error(bboot.gir(fit, R = 2), "supplied mid-CDF")
  w <- small_weights()
  fixed <- bboot.gir(fit, bootstrap_weights = w, normalize_weights = FALSE,
                      first_stage_boot = "fixed", inference = FALSE)
  cb <- function(object, weights, replicate) object$midcdf
  external <- bboot.gir(fit, bootstrap_weights = w, normalize_weights = FALSE,
                         midcdf_boot = cb, inference = FALSE)
  expect_equal(fixed$t, external$t)
  from_formula <- bboot.gir(y ~ d, data = d, midcdf = m,
    bootstrap_weights = w, normalize_weights = FALSE,
    first_stage_boot = "fixed", inference = FALSE)
  expect_equal(fixed$t, from_formula$t, tolerance = 1e-12)
  fail <- function(object, weights, replicate) {
    if (replicate == 2L) stop("intentional failure")
    object$midcdf
  }
  expect_error(bboot.gir(fit, bootstrap_weights = w, midcdf_boot = fail,
                         inference = FALSE), "replicate 2 failed")
  omitted <- expect_gir_warnings(bboot.gir(fit, bootstrap_weights = w,
                         midcdf_boot = fail, failure = "omit"),
    c("1 failed replicate(s) were omitted, not redrawn.",
      "Fewer than ten expected draws per confidence-interval tail."))
  expect_identical(omitted$success, c(TRUE, FALSE, TRUE))
  expect_true(all(is.na(omitted$t[2, ])))
  expect_equal(omitted$R.success, 2)
})

test_that("formula convenience preserves subset, NA and weight evaluation", {
  set.seed(16)
  d <- data.frame(y = runif(100, 0.1, 2), x = runif(100), a = runif(100, 1, 2),
                   id = 1:100, group = factor(rep(1:2, 50)))
  d$x[3] <- NA_real_
  fit <- expect_gir_warnings(gir(y ~ x + group, data = d, subset = id <= 80,
    weights = a, na.action = na.omit, first_stage = ~x + group,
    continuous = "x", discrete = "group",
    sieve_size = 2L, model = FALSE, x = FALSE, y = FALSE, singular.ok = FALSE),
    "1 observation(s) with missing values were removed by na.omit.")
  a <- expect_gir_warnings(bboot.gir(fit, R = 3L, seed = 9L, inference = FALSE),
    "Bootstrap diagnostics: 1 draw(s) with mid-CDF values outside [0,1]")
  b <- expect_gir_warnings(bboot.gir(y ~ x + group, data = d, subset = id <= 80,
    weights = a, na.action = na.omit, first_stage = ~x + group,
    continuous = "x", discrete = "group",
    sieve_size = 2L, model = FALSE, x = FALSE, y = FALSE, singular.ok = FALSE,
    R = 3L, seed = 9L, inference = FALSE),
    c("1 observation(s) with missing values were removed by na.omit.",
      "Bootstrap diagnostics: 1 draw(s) with mid-CDF values outside [0,1]"))
  expect_equal(a$t0, b$t0, tolerance = 1e-12)
  expect_equal(a$t, b$t, tolerance = 1e-12)
  expect_equal(fit$nobs, 79L)
  expect_null(fit$model); expect_null(fit$x); expect_null(fit$y)
  expect_equal(predict(fit), fitted(fit))
  expect_error(bboot.gir(fit, R = 2L, basis = "spline"), "already specifies")
})

test_that("prediction confidence intervals use joint coefficient draws", {
  fit <- small_fit()
  bb <- expect_gir_warnings(bboot.gir(fit, bootstrap_weights = small_weights(),
                                   normalize_weights = FALSE),
    "Fewer than ten expected draws per confidence-interval tail.")
  nd <- data.frame(d = factor(c("a", "b"), levels = c("a", "b")))
  set.seed(17); before <- .Random.seed
  pred <- expect_gir_warnings(predict(fit, nd, interval = "confidence",
                                     se.fit = TRUE, bootstrap = bb),
    "Few bootstrap draws per prediction-interval tail")
  expect_identical(.Random.seed, before)
  X <- model.matrix(~d, nd)
  pb <- bb$t %*% t(X); p0 <- drop(X %*% coef(fit))
  expect_equal(unname(pred$fit[, "fit"]), unname(p0), tolerance = 1e-12)
  expect_equal(unname(pred$se.fit), unname(apply(pb, 2, sd)), tolerance = 1e-12)
  roots <- sweep(pb, 2, p0, "-")
  lo <- apply(roots, 2, quantile, probs = .025, names = FALSE)
  hi <- apply(roots, 2, quantile, probs = .975, names = FALSE)
  expect_equal(unname(pred$fit[, c("lwr", "upr")]),
                unname(cbind(p0 - hi, p0 - lo)), tolerance = 1e-12)
  expect_error(predict(fit, nd, interval = "prediction"))
  expect_error(predict(fit, data.frame(d = "c")), "new level")
  expect_error(predict(fit, nd, bootstrap = bb, se.fit = TRUE, R = 3), "Do not supply")
  fit2 <- gir(y ~ d, data = transform(small_gir_data(), y = y + 1),
                discrete = "d", minimum_cell_size = 1)
  expect_error(predict(fit2, nd, bootstrap = bb, se.fit = TRUE), "different")
})

test_that("newdata NA policy and scalar designs keep dimensions", {
  fit <- small_fit()
  nd <- data.frame(d = c("a", NA, "b"))
  p <- predict(fit, nd)
  expect_length(p, 3); expect_true(is.na(p[2]))
  expect_length(predict(fit, nd, na.action = na.omit), 2)
  px <- predict(fit, nd, na.action = na.exclude)
  expect_length(px, 3); expect_true(is.na(px[2]))
  expect_error(predict(fit, nd, na.action = na.fail))
  d <- data.frame(y = 1:20)
  one <- gir(y ~ 1, data = d)
  b <- expect_gir_warnings(bboot.gir(one, R = 3, seed = 2),
    "Fewer than ten expected draws per confidence-interval tail.")
  expect_equal(dim(b$t), c(3L, 1L))
  expect_equal(dim(vcov(b)), c(1L, 1L))
  pr <- expect_gir_warnings(predict(one, data.frame(dummy = 1:2),
                          bootstrap = b, interval = "confidence", se.fit = TRUE),
    "Few bootstrap draws per prediction-interval tail")
  expect_equal(dim(pr$fit), c(2L, 3L)); expect_length(pr$se.fit, 2)
})

test_that("safe formula predictions reuse polynomial parameters", {
  set.seed(24)
  d <- data.frame(x = runif(100), y = runif(100, .1, 2))
  pp <- poly(d$x, 2)
  fit <- expect_gir_warnings(gir(y ~ poly(x, 2), data = d,
              first_stage = ~x, continuous = "x", sieve_size = 3))
  nd <- data.frame(x = c(.2, .5, .8))
  expected <- drop(cbind(1, predict(pp, nd$x)) %*% coef(fit))
  expect_equal(unname(expect_gir_warnings(predict(fit, nd))), unname(expected),
                tolerance = 1e-12)
})

test_that("aliased fits are protected against invalid inference and prediction", {
  set.seed(25)
  d <- data.frame(x = runif(80), y = runif(80, .1, 2)); d$xx <- 2 * d$x
  fit <- expect_gir_warnings(gir(y ~ x + xx, data = d,
                               first_stage = ~x, continuous = "x", sieve_size = 2),
    "The second-stage design is rank deficient; 1 coefficient(s) are aliased")
  expect_true(any(fit$aliased))
  expect_error(bboot.gir(fit, R = 2), "full-rank")
  good <- expect_gir_warnings(predict(fit, data.frame(x = .3, xx = .6)))
  expect_true(is.finite(good))
  bad <- expect_gir_warnings(predict(fit, data.frame(x = .3, xx = .9)),
    "1 newdata row(s) are non-estimable for this aliased fit")
  expect_true(is.na(bad))
})

# ---- clustered multipliers --------------------------------------------------
cluster_data <- function() {
  set.seed(8)
  d <- data.frame(ind = factor(rep(c("a", "b", "c", "d"), each = 30)),
                  t = factor(rep(rep(c("t1", "t2"), each = 15), 4)))
  d$x <- rnorm(nrow(d))
  d$y <- rgamma(nrow(d), 2, 1)
  d
}

test_that("cluster multipliers are shared within clusters and match a manual draw", {
  d <- cluster_data()
  fit <- gir(y ~ ind + t, data = d, first_stage = ~ ind + t,
             discrete = c("ind", "t"))
  bb <- bboot.gir(fit, R = 3L, seed = 5L, cluster = ~ ind, keep_weights = TRUE,
                  keep_midcdf = TRUE, inference = FALSE)
  expect_equal(bb$cluster$n_clusters, 4L)
  expect_true(bb$cluster$first_stage_invariant)
  for (b in 1:3) {
    per_cluster <- tapply(bb$bootstrap_weights[, b], d$ind,
                          function(v) length(unique(v)))
    expect_equal(as.vector(per_cluster), rep(1L, 4))
  }
  # cells nested in clusters: the first stage does not move
  expect_equal(unname(bb$midcdf[, 1]), unname(midcdf(fit)), tolerance = 1e-12)
  set.seed(5L)
  e <- rexp(4); w <- e[as.integer(d$ind)]; w <- w / mean(w)
  manual <- bboot.gir(fit, bootstrap_weights = w, normalize_weights = FALSE,
                      inference = FALSE)
  expect_equal(unname(bb$t[1, ]), unname(manual$t[1, ]), tolerance = 1e-12)
  by_vector <- bboot.gir(fit, R = 3L, seed = 5L, cluster = d$ind,
                         inference = FALSE)
  expect_equal(bb$t, by_vector$t)
  by_character <- bboot.gir(fit, R = 3L, seed = 5L, cluster = as.character(d$ind),
                            inference = FALSE)
  expect_equal(bb$t, by_character$t)
  two_way <- bboot.gir(fit, R = 2L, seed = 2L, cluster = ~ ind + t,
                       inference = FALSE)
  expect_equal(two_way$cluster$n_clusters, 8L)
})

test_that("cluster vectors align to the original data when rows were dropped", {
  d <- cluster_data()
  d$y[c(3, 50)] <- NA
  fit <- expect_gir_warnings(gir(y ~ ind + t, data = d, first_stage = ~ ind + t,
                              discrete = c("ind", "t"), na.action = na.omit),
    "2 observation(s) with missing values were removed by na.omit.")
  expect_equal(fit$nobs, nrow(d) - 2L)
  a <- bboot.gir(fit, R = 2L, seed = 1L, cluster = d$ind, inference = FALSE)
  b <- bboot.gir(fit, R = 2L, seed = 1L, cluster = d$ind[fit$row.index],
                 inference = FALSE)
  expect_equal(a$t, b$t)
  expect_error(bboot.gir(fit, R = 2L, cluster = d$ind[-1], inference = FALSE),
               "one value per")
})

test_that("cluster arguments are validated", {
  d <- cluster_data()
  fit <- gir(y ~ ind + t, data = d, first_stage = ~ ind + t,
             discrete = c("ind", "t"))
  expect_error(bboot.gir(fit, R = 2L, cluster = ~ ind, inference = FALSE,
                         bootstrap_weights = matrix(1, fit$nobs, 2)),
               "either cluster or bootstrap_weights")
  expect_error(bboot.gir(fit, R = 2L, cluster = ~ nothing, inference = FALSE),
               "not a retained")
  expect_error(bboot.gir(fit, R = 2L, cluster = rep(1, fit$nobs),
                         inference = FALSE), "at least two")
  expect_error(bboot.gir(fit, R = 2L, cluster = y ~ ind, inference = FALSE),
               "one-sided")
  d$ind_na <- d$ind; d$ind_na[1] <- NA
  expect_error(bboot.gir(fit, R = 2L, cluster = d$ind_na, inference = FALSE),
               "missing")
})

test_that("a first stage pooled across clusters varies with cluster multipliers", {
  d <- cluster_data()
  fit <- gir(y ~ x + ind, data = d, first_stage = ~ x, continuous = "x",
             sieve_size = 3)
  expect_warning(
    bb <- bboot.gir(fit, R = 2L, seed = 3L, cluster = ~ ind, keep_midcdf = TRUE,
                    inference = FALSE),
    "Bootstrap diagnostics: 1 draw(s) with mid-CDF values outside [0,1]",
    fixed = TRUE)
  expect_false(bb$cluster$first_stage_invariant)
  expect_true(any(abs(bb$midcdf[, 1] - midcdf(fit)) > 1e-8))
  expect_equal(sum(bb$diagnostics$midcdf_outside > 0), 1L)
})

test_that("clustered inference prints its cluster structure and feeds predict", {
  d <- cluster_data()
  fit <- gir(y ~ ind + t, data = d, first_stage = ~ ind + t,
             discrete = c("ind", "t"))
  bb <- expect_gir_warnings(bboot.gir(fit, R = 40L, seed = 7L, cluster = ~ ind),
    "Fewer than ten expected draws per confidence-interval tail.")
  expect_output(print(bb), "4 clusters")
  expect_output(print(expect_gir_warnings(summary(bb),
    "Fewer than ten expected draws per confidence-interval tail.")), "nested in clusters")
  pr <- expect_gir_warnings(predict(fit, newdata = d[c(1, 31), ], se.fit = TRUE,
                                 R = 40L, seed = 7L, cluster = ~ ind),
    "Fewer than ten expected draws per confidence-interval tail.")
  direct <- expect_gir_warnings(predict(fit, newdata = d[c(1, 31), ], se.fit = TRUE,
                                     bootstrap = bb))
  expect_equal(pr$se.fit, direct$se.fit)
})
