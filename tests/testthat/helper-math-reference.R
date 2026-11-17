# Small deterministic fixtures for the paper's finite-sample identities.
# Reference calculations never call GIR's basis, solver or cumulative-score
# helpers. Keeping the response support small makes direct threshold fits cheap.
math_data <- function() {
  make_cell <- function(nx, nz, group) {
    d <- expand.grid(i = seq_len(nx), j = seq_len(nz), repeat_id = 1:5)
    d$x <- -2 + 5 * (d$i - 1) / (nx - 1)
    d$z <- 10 + 8 * (d$j - 1) / (nz - 1)
    category <- ifelse(d$repeat_id <= 4, d$repeat_id,
                       1 + (d$i + 2 * d$j + group) %% 4)
    data.frame(y = c(0, 1, 3, 6)[category], x = d$x, z = d$z,
               g = if (group == 1) "a" else "b")
  }
  d <- rbind(make_cell(7, 6, 1), make_cell(6, 5, 2))
  d$g <- factor(d$g, levels = c("a", "b"))
  d
}

math_multipliers <- function(n, R = 40L) {
  # Positive, heterogeneous columns; neither random draws nor a coverage study.
  i <- seq_len(n)
  b <- seq_len(R)
  W <- 1 + 0.35 * sin(outer(i, b, function(i, b) i * (b + 1) * 0.37)) +
    0.2 * cos(outer(i, b, function(i, b) (i + 3) * b * 0.19))
  sweep(W, 2, colMeans(W), "/")
}

math_basis <- function(d, basis) {
  u <- (d$x + 2) / 5
  v <- (d$z - 10) / 8
  margin <- function(s) {
    if (basis == "spline") {
      return(splines::bs(s, knots = 0.5, degree = 3, intercept = TRUE,
                         Boundary.knots = c(0, 1)))
    }
    t <- 2 * s - 1
    cbind(1, sqrt(3) * t, sqrt(5) * (3 * t^2 - 1) / 2)
  }
  A <- margin(u)
  B <- margin(v)
  pairs <- expand.grid(a = seq_len(ncol(A)), b = seq_len(ncol(B)))
  vapply(seq_len(nrow(pairs)), function(k) {
    A[, pairs$a[k]] * B[, pairs$b[k]]
  }, numeric(nrow(d)))
}

math_fit <- function(d, basis, a = rep(1, nrow(d)), GIR = TRUE) {
  gir(y ~ x + z + g, data = d, weights = a, first_stage = ~ x + z + g,
      continuous = c("x", "z"), discrete = "g", basis = basis,
      sieve_size = if (basis == "legendre") 9L else c(1L, 1L),
      sieve_size_type = if (basis == "legendre") "dimension" else "knots",
      basis_bounds = list(x = c(-2, 3), z = c(10, 18)),
      first_stage_singular = "error", clip = FALSE, singular.ok = FALSE,
      GIR = GIR)
}

math_reference <- function(d, P, a = rep(1, nrow(d))) {
  m <- numeric(nrow(d))
  for (rows in split(seq_len(nrow(d)), d$g)) {
    for (threshold in unique(d$y[rows])) {
      response <- (d$y[rows] < threshold) + 0.5 * (d$y[rows] == threshold)
      regression <- stats::lm.wfit(P[rows, , drop = FALSE], response, w = a[rows])
      at_threshold <- d$y[rows] == threshold
      m[rows[at_threshold]] <- regression$fitted.values[at_threshold]
    }
  }
  X <- model.matrix(~ x + z + g, d)
  theta <- stats::lm.wfit(X, 2 * m - 1, w = a * d$y)$coefficients
  list(midcdf = m, coefficients = theta, fitted = drop(X %*% theta))
}

math_boot <- function(fit, W) {
  # Forty draws and level 0.5 keep >= 10 expected draws per tail. This small
  # deterministic calculation tests interval arithmetic, not statistical coverage.
  bboot.gir(fit, bootstrap_weights = W, normalize_weights = FALSE,
            level = 0.5, keep_midcdf = TRUE)
}
