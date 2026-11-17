cluster_order_data <- function() {
  data.frame(y = c(1, 4, 2, 7, 3, 9), id = c("A", "A", "B", "B", "C", "C"))
}

cluster_order_fit <- function(d, rows = seq_len(nrow(d)), model = TRUE,
                              na.action = na.fail) {
  gir(y ~ 1, data = d, subset = rows, first_stage = ~ id, discrete = "id",
      minimum_cell_size = 1L, model = model, na.action = na.action)
}

test_that("equal-length reordered cluster vectors require an explicit order", {
  d <- cluster_order_data()
  rows <- c(1L, 3L, 5L, 2L, 4L, 6L)
  fit <- cluster_order_fit(d, rows)
  expect_identical(fit$row.index, rows)
  set.seed(10L)
  before <- .Random.seed
  expect_error(bboot.gir(fit, R = 4L, cluster = d$id, inference = FALSE),
               "cluster vector order is ambiguous", fixed = TRUE)
  expect_error(bboot.gir(fit, R = 4L, cluster = d$id, cluster_order = "auto",
                         inference = FALSE), "Set cluster_order", fixed = TRUE)
  expect_identical(.Random.seed, before)

  original <- bboot.gir(fit, R = 4L, seed = 17L, cluster = d$id,
                        cluster_order = "original", inference = FALSE,
                        keep_weights = TRUE, keep_midcdf = TRUE)
  estimation <- bboot.gir(fit, R = 4L, seed = 17L, cluster = d$id[rows],
                          cluster_order = "estimation", inference = FALSE,
                          keep_weights = TRUE, keep_midcdf = TRUE)
  formula <- bboot.gir(fit, R = 4L, seed = 17L, cluster = ~ id,
                       inference = FALSE, keep_weights = TRUE, keep_midcdf = TRUE)
  expect_identical(original$cluster$index, rep(1:3, 2))
  expect_identical(original$cluster$n_clusters, 3L)
  expect_identical(original$t0, coef(fit))
  expect_identical(original$bootstrap_weights, formula$bootstrap_weights)
  expect_identical(estimation$bootstrap_weights, formula$bootstrap_weights)
  expect_identical(original$t, formula$t)
  expect_identical(estimation$t, formula$t)
  expect_identical(original$midcdf, formula$midcdf)

  # Calculate weights and both stages directly in the known estimation order.
  set.seed(17L)
  index <- rep(1:3, 2)
  y <- d$y[rows]
  for (b in 1:4) {
    w <- rexp(3)[index]
    w <- w / mean(w)
    m <- vapply(seq_along(y), function(i) {
      g <- index == index[i]
      sum(w[g] * ((y[g] < y[i]) + 0.5 * (y[g] == y[i]))) / sum(w[g])
    }, numeric(1))
    theta <- sum(w * y * (2 * m - 1)) / sum(w * y)
    expect_equal(unname(original$bootstrap_weights[, b]), w, tolerance = 1e-12)
    expect_equal(unname(original$midcdf[, b]), m, tolerance = 1e-12)
    expect_equal(unname(original$t[b, 1]), theta, tolerance = 1e-12)
  }
})

test_that("automatic ordering handles unchanged, shorter and expanded samples", {
  d <- cluster_order_data()
  selections <- list(1:6, c(6L, 2L, 1L, 4L), c(6L, 2L, 1L, 4L, 3L, 5L, 1L, 3L))
  for (rows in selections) {
    fit <- cluster_order_fit(d, rows, model = FALSE)
    formula <- bboot.gir(fit, R = 3L, seed = 19L, cluster = ~ id,
                         inference = FALSE, keep_weights = TRUE)
    for (order in c("auto", "original", "estimation")) {
      values <- if (order == "estimation") d$id[rows] else d$id
      bb <- bboot.gir(fit, R = 3L, seed = 19L, cluster = values,
                      cluster_order = order, inference = FALSE, keep_weights = TRUE)
      expect_identical(bb$cluster$index, match(d$id[rows], unique(d$id[rows])))
      expect_identical(bb$bootstrap_weights, formula$bootstrap_weights)
      expect_identical(bb$t, formula$t)
    }
    # Also infer an already-aligned vector when only its length matches.
    aligned <- bboot.gir(fit, R = 3L, seed = 19L, cluster = d$id[rows],
                         inference = FALSE, keep_weights = TRUE)
    expect_identical(aligned$bootstrap_weights, formula$bootstrap_weights)
    expect_identical(aligned$t, formula$t)
  }
})

test_that("equal-length selections with repeated rows are also ambiguous", {
  d <- cluster_order_data()
  rows <- c(1L, 3L, 5L, 2L, 4L, 1L)
  fit <- cluster_order_fit(d, rows)
  expect_error(bboot.gir(fit, R = 3L, cluster = d$id, inference = FALSE),
               "cluster vector order is ambiguous", fixed = TRUE)
  original <- bboot.gir(fit, R = 3L, seed = 21L, cluster = d$id,
                        cluster_order = "original", inference = FALSE)
  estimation <- bboot.gir(fit, R = 3L, seed = 21L, cluster = d$id[rows],
                          cluster_order = "estimation", inference = FALSE)
  formula <- bboot.gir(fit, R = 3L, seed = 21L, cluster = ~ id, inference = FALSE)
  expect_identical(original$cluster$index, c(1L, 2L, 3L, 1L, 2L, 1L))
  expect_identical(original$t, formula$t)
  expect_identical(estimation$t, formula$t)
})

test_that("alignment precedes missing-cluster checks after subset and NA removal", {
  d <- cluster_order_data()
  d$y[c(2, 5)] <- NA_real_
  expect_warning(fit <- cluster_order_fit(d, c(6L, 2L, 1L, 5L, 4L, 3L),
                                          na.action = na.omit),
                  "2 observation(s) with missing values", fixed = TRUE)
  expect_identical(fit$row.index, c(6L, 1L, 4L, 3L))
  values <- d$id
  values[c(2, 5)] <- NA_character_
  formula <- bboot.gir(fit, R = 3L, seed = 23L, cluster = ~ id,
                       inference = FALSE, keep_weights = TRUE)
  for (order in c("auto", "original", "estimation")) {
    v <- if (order == "estimation") values[fit$row.index] else values
    bb <- bboot.gir(fit, R = 3L, seed = 23L, cluster = v, cluster_order = order,
                    inference = FALSE, keep_weights = TRUE)
    expect_identical(bb$bootstrap_weights, formula$bootstrap_weights)
    expect_identical(bb$t, formula$t)
  }
  values[1] <- NA_character_
  for (order in c("auto", "original", "estimation")) {
    v <- if (order == "estimation") values[fit$row.index] else values
    expect_error(bboot.gir(fit, R = 3L, cluster = v, cluster_order = order,
                           inference = FALSE), "cluster contains missing values",
                 fixed = TRUE)
  }
})

test_that("cluster order and its required vector lengths are validated", {
  d <- cluster_order_data()
  fit <- cluster_order_fit(d, c(6L, 2L, 1L, 4L))
  expect_error(bboot.gir(fit, R = 3L, cluster = d$id[fit$row.index],
                         cluster_order = "original", inference = FALSE),
               "cluster_order = 'original' requires a cluster vector of length 6",
               fixed = TRUE)
  expect_error(bboot.gir(fit, R = 3L, cluster = d$id,
                         cluster_order = "estimation", inference = FALSE),
               "cluster_order = 'estimation' requires a cluster vector of length 4",
               fixed = TRUE)
  expect_error(bboot.gir(fit, R = 3L, cluster = d$id[-1], inference = FALSE),
               "one value per", fixed = TRUE)
  expect_error(bboot.gir(fit, R = 3L, cluster = d$id,
                         cluster_order = "unknown", inference = FALSE), "arg.*should be")
  for (order in c("original", "estimation")) {
    expect_error(bboot.gir(fit, R = 3L, cluster = ~ id, cluster_order = order,
                           inference = FALSE), "applies only to cluster vectors",
                 fixed = TRUE)
    expect_error(bboot.gir(fit, R = 3L, cluster_order = order, inference = FALSE),
                 "cluster is NULL", fixed = TRUE)
  }
})

test_that("formula fitting and prediction forward cluster_order correctly", {
  d <- cluster_order_data()
  rows <- c(1L, 3L, 5L, 2L, 4L, 6L)
  fit <- cluster_order_fit(d, rows)
  bb <- bboot.gir(fit, R = 40L, level = 0.5, seed = 17L, cluster = d$id,
                  cluster_order = "original")
  from_formula <- bboot.gir(y ~ 1, data = d, subset = rows, first_stage = ~ id,
                            discrete = "id", minimum_cell_size = 1L,
                            R = 40L, level = 0.5, seed = 17L, cluster = d$id,
                            cluster_order = "original")
  expect_identical(from_formula$t, bb$t)
  automatic <- predict(fit, se.fit = TRUE, R = 40L, level = 0.5, seed = 17L,
                        cluster = d$id, cluster_order = "original")
  reused <- predict(fit, se.fit = TRUE, level = 0.5, bootstrap = bb)
  expect_identical(automatic$fit, reused$fit)
  expect_identical(automatic$se.fit, reused$se.fit)
  expect_error(predict(fit, se.fit = TRUE, cluster = d$id),
               "cluster vector order is ambiguous", fixed = TRUE)
})
