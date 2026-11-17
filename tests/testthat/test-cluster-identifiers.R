test_that("cluster formulas distinguish exact value combinations", {
  cases <- list(
    dots = list(a = c("a.b", "a", "z"), b = c("c", "b.c", "q")),
    separators = list(a = c("a|b", "a", "z"), b = c("c", "b|c", "q")),
    factors = list(
      a = factor(c("a.b", "a", "z"), levels = c("unused", "z", "a", "a.b")),
      b = ordered(c("c", "b.c", "q"), levels = c("q", "b.c", "c", "unused"))),
    numeric = list(a = c(1, 1 + 4 * .Machine$double.eps, 2),
                   b = c("same", "same", "other")),
    three_columns = list(a = c("a.b", "a", "z"), b = c("c", "b.c", "q"),
                         c = c("u", "u", "v")),
    renamed = list(a = c("z", "a", "m"), b = c("third", "first", "second"))
  )
  reference <- NULL
  for (columns in cases) {
    d <- as.data.frame(columns)[rep(1:3, 4), , drop = FALSE]
    d$y <- c(1, 2, 4, 3, 7, 6, 2, 5, 3, 4, 1, 8)
    grouping <- reformulate(names(columns))
    fit <- gir(y ~ 1, data = d, first_stage = grouping,
               discrete = names(columns), minimum_cell_size = 1L)
    bb <- bboot.gir(fit, R = 4L, seed = 23L, cluster = grouping,
                    inference = FALSE, keep_weights = TRUE)
    expect_identical(bb$cluster$n_clusters, 3L)
    expect_identical(bb$cluster$index, rep(1:3, 4))
    if (is.null(reference)) reference <- bb
    # Relabeling the same partition must preserve its seeded multipliers/draws.
    expect_identical(bb$bootstrap_weights, reference$bootstrap_weights)
    expect_equal(bb$t, reference$t, tolerance = 1e-12)
  }
})

test_that("colliding formula labels match explicit clusters and independent refits", {
  pairs <- expand.grid(a = c("a.b", "a", "z"), b = c("c", "b.c", "q"),
                       stringsAsFactors = FALSE)
  index <- rep(1:9, 3)
  d <- pairs[index, , drop = FALSE]
  d$y <- rep(c(1, 4, 2, 3, 5, 2, 7, 6, 4), 3) + rep(c(0, 1, 3), each = 9)
  # Pool the first stage across clusters so both stages change in each draw.
  fit <- gir(y ~ a + b, data = d, first_stage = ~ 1, singular.ok = FALSE)
  bb <- bboot.gir(fit, R = 4L, seed = 37L, cluster = ~ a + b,
                  inference = FALSE, keep_weights = TRUE, keep_midcdf = TRUE)
  explicit <- bboot.gir(fit, R = 4L, seed = 37L, cluster = index,
                        inference = FALSE, keep_weights = TRUE)
  expect_identical(bb$cluster$n_clusters, 9L)
  expect_identical(bb$cluster$index, index)
  expect_false(bb$cluster$first_stage_invariant)
  expect_identical(bb$bootstrap_weights, explicit$bootstrap_weights)
  expect_identical(bb$t, explicit$t)
  expect_true(any(abs(bb$midcdf[, 1] - midcdf(fit)) > 1e-8))

  X <- model.matrix(y ~ a + b, d)
  set.seed(37L)
  for (b in 1:4) {
    # Expected weights come from the known nine groups, not the resolver.
    w <- rexp(9)[index]
    w <- w / mean(w)
    m <- vapply(seq_along(d$y), function(i) {
      sum(w * ((d$y < d$y[i]) + 0.5 * (d$y == d$y[i]))) / sum(w)
    }, numeric(1))
    theta <- solve(crossprod(X, X * (w * d$y)),
                   colSums(X * (w * d$y * (2 * m - 1))))
    expect_equal(unname(bb$bootstrap_weights[, b]), w, tolerance = 1e-12)
    expect_equal(unname(bb$midcdf[, b]), m, tolerance = 1e-12)
    expect_equal(unname(bb$t[b, ]), unname(theta), tolerance = 1e-12)
  }
})

test_that("coding formula columns does not turn missing values into valid clusters", {
  object <- list(nobs = 4L, original.nobs = 4L, row.index = 1:4,
                 model = data.frame(a = c("a", "b", "a", "b"),
                                    b = factor(c("c", "d", "c", "d"))),
                 estimation = list(first_stage_data = NULL))
  for (column in c("a", "b")) {
    missing <- object
    missing$model[[column]][1] <- NA
    expect_error(GIR:::.gir_cluster_index(~ a + b, missing),
                 "cluster contains missing values", fixed = TRUE)
  }
  object$model$a[] <- "a"
  object$model$b[] <- "c"
  expect_error(GIR:::.gir_cluster_index(~ a + b, object),
               "cluster must define at least two clusters", fixed = TRUE)
})
