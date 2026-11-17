# GIR: linear Gini Index Regression
# ---------------------------------
# Point estimation (gir), two-stage Bayesian bootstrap (bboot.gir) and
# linear-index prediction (predict.gir). Runtime dependencies are the
# recommended packages splines, stats and utils. Help files are generated
# from the roxygen comments in this file.


# =============================================================================
# Internal utilities
# =============================================================================

.gir_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

.gir_positive <- function(x, name, zero_ok = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  ok <- ok && if (zero_ok) x >= 0 else x > 0
  if (!ok) {
    stop(
      name, " must be a finite ",
      if (zero_ok) "nonnegative" else "strictly positive",
      " scalar.", call. = FALSE
    )
  }
  invisible(x)
}

.gir_names <- function(x, name) {
  if (is.null(x)) return(character(0))
  if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
    stop(name, " must be NULL or a character vector of column names.",
         call. = FALSE)
  }
  unique(x)
}

.gir_subset_index <- function(value, n, row_names = NULL) {
  if (is.null(value)) return(seq_len(n))
  if (is.character(value)) {
    if (is.null(row_names)) {
      stop("Character subset indices require row names.", call. = FALSE)
    }
    index <- match(value, row_names)
  } else if (is.logical(value) || is.numeric(value)) {
    index <- seq_len(n)[value]
  } else {
    stop("subset must evaluate to a logical, numeric, or character index.",
         call. = FALSE)
  }
  index <- index[!is.na(index)]
  if (!length(index)) stop("subset selects no observations.", call. = FALSE)
  as.integer(index)
}

.gir_align <- function(value, selected, n_original, n_selected, name) {
  if (length(value) == 1L) return(rep(value, n_selected))
  if (length(value) == n_original) return(value[selected])
  if (length(value) == n_selected) return(value)
  stop(
    name, " must have length 1, nrow(data) = ", n_original,
    ", or the number of selected rows = ", n_selected, ".",
    call. = FALSE
  )
}

.gir_na_name <- function(fun) {
  if (identical(fun, stats::na.fail)) return("na.fail")
  if (identical(fun, stats::na.omit)) return("na.omit")
  if (identical(fun, stats::na.exclude)) return("na.exclude")
  if (identical(fun, stats::na.pass)) return("na.pass")
  "custom na.action"
}

.gir_apply_na <- function(model_frame, fs_frame, weights, supplied_midcdf,
                          na.action) {
  n <- nrow(model_frame)
  probe <- data.frame(.GIR_row = seq_len(n))
  for (j in seq_along(model_frame)) {
    probe[[paste0(".GIR_model_", j)]] <- model_frame[[j]]
  }
  for (j in seq_along(fs_frame)) {
    probe[[paste0(".GIR_fs_", j)]] <- fs_frame[[j]]
  }
  probe$.GIR_weights <- weights
  if (!is.null(supplied_midcdf)) probe$.GIR_midcdf <- supplied_midcdf

  if (is.null(na.action)) na.action <- getOption("na.action", "na.fail")
  fun <- match.fun(na.action)
  action_name <- .gir_na_name(fun)

  acted <- tryCatch(
    fun(probe),
    error = function(e) {
      prefix <- if (identical(fun, stats::na.fail)) {
        "Missing values are present: "
      } else paste0(action_name, " failed: ")
      stop(prefix, conditionMessage(e), call. = FALSE)
    }
  )
  if (!is.data.frame(acted) || !(".GIR_row" %in% names(acted))) {
    stop("na.action must return a data frame retaining its columns.",
         call. = FALSE)
  }
  keep <- acted$.GIR_row
  if (!is.numeric(keep) || anyNA(keep) || any(!is.finite(keep)) ||
      any(keep != floor(keep)) || anyDuplicated(keep) ||
      any(keep < 1L | keep > n)) {
    stop("na.action returned invalid row indices.", call. = FALSE)
  }
  keep <- as.integer(keep)
  if (!length(keep)) {
    stop("No observations remain after ", action_name, ": all ", n,
         " selected observations were removed.", call. = FALSE)
  }
  incomplete <- which(!stats::complete.cases(probe))
  if (length(intersect(keep, incomplete))) {
    stop(
      "Missing values remain after na.action. Use na.fail, na.omit, or ",
      "na.exclude, or a custom action that removes incomplete rows.",
      call. = FALSE
    )
  }
  removed <- setdiff(seq_len(n), keep)
  omitted <- NULL
  if (length(removed)) {
    omitted <- removed
    names(omitted) <- rownames(model_frame)[removed]
    class(omitted) <- if (identical(fun, stats::na.exclude)) "exclude" else "omit"
  }
  list(keep = keep, removed = removed, object = omitted, name = action_name,
       removed_incomplete = sum(removed %in% incomplete))
}

.gir_fs_variables <- function(terms_object, first_stage, continuous, discrete,
                              data) {
  variables <- if (is.null(first_stage)) {
    all.vars(stats::delete.response(terms_object))
  } else if (inherits(first_stage, "formula")) {
    fs_terms <- stats::terms(first_stage, data = data)
    if (attr(fs_terms, "response") > 0L) {
      stop("first_stage must be one-sided, for example ~ x1 + x2.",
           call. = FALSE)
    }
    all.vars(fs_terms)
  } else if (is.character(first_stage)) {
    first_stage
  } else {
    stop("first_stage must be NULL, a one-sided formula, or column names.",
         call. = FALSE)
  }
  variables <- unique(variables)
  absent <- setdiff(variables, names(data))
  if (length(absent)) {
    stop("First-stage columns not found in data: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }

  continuous <- .gir_names(continuous, "continuous")
  discrete <- .gir_names(discrete, "discrete")
  if (length(setdiff(continuous, variables))) {
    stop("continuous contains columns not included in first_stage.",
         call. = FALSE)
  }
  if (length(setdiff(discrete, variables))) {
    stop("discrete contains columns not included in first_stage.",
         call. = FALSE)
  }
  overlap <- intersect(continuous, discrete)
  if (length(overlap)) {
    stop("Columns cannot be both continuous and discrete: ",
         paste(overlap, collapse = ", "), call. = FALSE)
  }

  remaining <- setdiff(variables, union(continuous, discrete))
  defaulted_continuous <- character(0)
  if (length(remaining)) {
    is_discrete <- vapply(
      data[remaining],
      function(z) is.factor(z) || is.ordered(z) || is.character(z) ||
        is.logical(z),
      logical(1)
    )
    discrete <- unique(c(discrete, remaining[is_discrete]))
    defaulted_continuous <- remaining[!is_discrete]
    continuous <- unique(c(continuous, defaulted_continuous))
  }

  for (v in continuous) {
    if (!is.numeric(data[[v]])) {
      stop("Continuous first-stage column '", v, "' must be numeric.",
           call. = FALSE)
    }
  }

  list(all = variables, continuous = continuous, discrete = discrete,
       defaulted_continuous = defaulted_continuous)
}

.gir_check_fs_support <- function(spec, data) {
  # Check support only after missing-value handling and sample validation.
  # The classification is the user's decision. Numeric columns left
  # undeclared default to continuous, and that default is reported. The
  # one default that can never be right, a two-point support carrying a
  # continuous basis, is refused; an explicit declaration is still accepted.
  for (v in spec$defaulted_continuous) {
    n_distinct <- length(unique(data[[v]][!is.na(data[[v]])]))
    if (n_distinct <= 2L) {
      stop("Numeric first-stage column '", v, "' has only ", n_distinct,
           " distinct value(s) and was not declared. A two-point support ",
           "cannot carry a continuous basis: list it in 'discrete', or in ",
           "'continuous' if that is really intended.", call. = FALSE)
    }
  }
  if (length(spec$defaulted_continuous)) {
    message("First stage: numeric column(s) ",
            paste(spec$defaulted_continuous, collapse = ", "),
            " were not declared and are treated as continuous. List them ",
            "in 'continuous' or 'discrete' to make the choice explicit.")
  }

  invisible(NULL)
}

.gir_bounds <- function(data, continuous, basis_bounds) {
  if (!length(continuous)) return(list())
  out <- vector("list", length(continuous)); names(out) <- continuous

  if (is.null(basis_bounds)) {
    for (v in continuous) out[[v]] <- range(data[[v]], finite = TRUE)
  } else if (is.numeric(basis_bounds) && length(continuous) == 1L &&
             length(basis_bounds) == 2L) {
    out[[continuous]] <- as.numeric(basis_bounds)
  } else if (is.list(basis_bounds)) {
    if (is.null(names(basis_bounds))) {
      if (length(basis_bounds) != length(continuous)) {
        stop("Unnamed basis_bounds must have one entry per continuous column.",
             call. = FALSE)
      }
      names(basis_bounds) <- continuous
    }
    missing_names <- setdiff(continuous, names(basis_bounds))
    if (length(missing_names)) {
      stop("basis_bounds is missing: ", paste(missing_names, collapse = ", "),
           call. = FALSE)
    }
    out <- basis_bounds[continuous]
  } else if (is.matrix(basis_bounds) || is.data.frame(basis_bounds)) {
    bb <- as.matrix(basis_bounds)
    if (ncol(bb) != 2L) stop("basis_bounds must have two columns.",
                             call. = FALSE)
    if (!is.null(rownames(bb)) && all(continuous %in% rownames(bb))) {
      bb <- bb[continuous, , drop = FALSE]
    } else if (nrow(bb) != length(continuous)) {
      stop("basis_bounds must have one row per continuous column.",
           call. = FALSE)
    }
    out <- lapply(seq_along(continuous), function(j) as.numeric(bb[j, ]))
    names(out) <- continuous
  } else {
    stop("Invalid basis_bounds specification.", call. = FALSE)
  }

  tol <- sqrt(.Machine$double.eps)
  for (v in continuous) {
    b <- as.numeric(out[[v]])
    if (length(b) != 2L || any(!is.finite(b)) || b[1L] >= b[2L]) {
      stop("Bounds for '", v, "' must be finite with lower < upper.",
           call. = FALSE)
    }
    z <- data[[v]]
    if (any(z < b[1L] - tol | z > b[2L] + tol)) {
      stop("Observed values of '", v, "' lie outside basis_bounds.",
           call. = FALSE)
    }
    out[[v]] <- b
  }
  out
}

.gir_scale <- function(data, continuous, bounds) {
  out <- data.frame(row.names = seq_len(nrow(data)))
  for (v in continuous) {
    b <- bounds[[v]]
    out[[v]] <- (data[[v]] - b[1L]) / (b[2L] - b[1L])
  }
  out
}

.gir_cells <- function(data, discrete) {
  n <- nrow(data)
  if (!length(discrete)) {
    return(factor(rep("all observations", n), levels = "all observations"))
  }
  # Codes are built with match(), so grouping does not rely on display rounding.
  codes <- lapply(discrete, function(v) match(data[[v]], unique(data[[v]])))
  key <- do.call(paste, c(codes, sep = "|"))
  level_keys <- unique(key)

  display <- lapply(discrete, function(v) {
    z <- data[[v]]
    if (is.factor(z) || is.ordered(z)) as.character(z) else format(z, trim = TRUE)
  })
  display_key <- do.call(paste, c(lapply(seq_along(discrete), function(j) {
    paste0(discrete[j], "=", display[[j]])
  }), sep = ", "))
  level_labels <- display_key[match(level_keys, key)]
  if (anyDuplicated(level_labels)) {
    level_labels <- paste0(level_labels, " [cell ", level_keys, "]")
  }
  factor(key, levels = level_keys, labels = level_labels)
}


# =============================================================================
# Basis construction
# =============================================================================

.gir_legendre <- function(z, number_terms, name) {
  z <- as.numeric(z); number_terms <- as.integer(number_terms)
  if (length(number_terms) != 1L || is.na(number_terms) || number_terms < 1L) {
    stop("Legendre dimension must be a positive integer.", call. = FALSE)
  }
  if (any(!is.finite(z)) || any(z < 0 | z > 1)) {
    stop("Scaled continuous regressors must lie in [0,1].", call. = FALSE)
  }
  t <- 2 * z - 1
  P <- matrix(0, nrow = length(z), ncol = number_terms)
  P[, 1L] <- 1
  if (number_terms >= 2L) P[, 2L] <- t
  if (number_terms >= 3L) {
    for (degree in 2:(number_terms - 1L)) {
      P[, degree + 1L] <- ((2 * degree - 1) * t * P[, degree] -
        (degree - 1) * P[, degree - 1L]) / degree
    }
  }
  P <- sweep(P, 2L, sqrt(2 * (0:(number_terms - 1L)) + 1), "*")
  colnames(P) <- paste0(name, ":leg", 0:(number_terms - 1L))
  P
}

.gir_spline <- function(z, order, knots, name) {
  z <- as.numeric(z); order <- as.integer(order); knots <- as.integer(knots)
  if (order < 2L) stop("spline_order must be at least 2.", call. = FALSE)
  if (knots < 0L) stop("Spline knot count cannot be negative.",
                        call. = FALSE)
  if (any(!is.finite(z)) || any(z < 0 | z > 1)) {
    stop("Scaled continuous regressors must lie in [0,1].", call. = FALSE)
  }
  interior <- if (knots == 0L) numeric(0) else seq_len(knots) / (knots + 1)
  P <- as.matrix(splines::bs(
    z, knots = interior, degree = order - 1L, intercept = TRUE,
    Boundary.knots = c(0, 1)
  ))
  colnames(P) <- paste0(name, ":spl", seq_len(ncol(P)))
  P
}

.gir_tensor <- function(mats, max_columns) {
  if (!length(mats)) stop("Internal error: empty tensor basis.", call. = FALSE)
  out <- as.matrix(mats[[1L]])
  if (length(mats) == 1L) return(out)
  for (j in 2:length(mats)) {
    B <- as.matrix(mats[[j]])
    new_p <- ncol(out) * ncol(B)
    if (!is.finite(new_p) || new_p > max_columns) {
      stop("Tensor-product basis exceeds max_basis_columns = ", max_columns,
           ".", call. = FALSE)
    }
    li <- rep(seq_len(ncol(out)), each = ncol(B))
    ri <- rep(seq_len(ncol(B)), times = ncol(out))
    left_names <- colnames(out); right_names <- colnames(B)
    out <- out[, li, drop = FALSE] * B[, ri, drop = FALSE]
    colnames(out) <- paste(left_names[li], right_names[ri], sep = " * ")
  }
  out
}

.gir_balanced_dims <- function(target, d, minimum, maximum_product) {
  if (length(minimum) == 1L) minimum <- rep(as.integer(minimum), d)
  dims <- minimum
  target <- min(max(as.integer(floor(target)), prod(minimum)), maximum_product)
  repeat {
    candidate_list <- lapply(seq_len(d), function(j) {
      z <- dims; z[j] <- z[j] + 1L; z
    })
    products <- vapply(candidate_list, prod, numeric(1))
    feasible <- which(products <= target & products <= maximum_product)
    if (!length(feasible)) break
    current <- dims[feasible]
    candidates <- feasible[current == min(current)]
    chosen <- candidates[which.max(products[candidates])]
    dims <- candidate_list[[chosen]]
  }
  dims
}

.gir_reduce_knots <- function(knots, order, maximum_product) {
  knots <- as.integer(knots)
  while (prod(knots + order) > maximum_product && any(knots > 0L)) {
    j <- which.max(knots + order)
    knots[j] <- knots[j] - 1L
  }
  knots
}

.gir_size_spec <- function(n_cell, n_positive, n_total, d, basis, kappa,
                           basis_constant, sieve_size, sieve_size_type,
                           sieve_size_reference, spline_order,
                           max_basis_columns) {
  if (d == 0L) {
    return(list(dims = integer(0), knots = integer(0), requested = 1L,
                actual = 1L, capped = FALSE, floored = FALSE,
                type = "constant"))
  }
  n_ref <- if (sieve_size_reference == "cell") n_cell else n_total
  raw <- if (is.null(sieve_size)) floor(basis_constant * n_ref^kappa) else sieve_size
  # "auto" targets the total tensor-product dimension, the K of the paper,
  # for both bases. A knots-per-margin alternative is available as
  # sieve_size_type = "knots".
  type <- if (sieve_size_type == "auto") "dimension" else sieve_size_type
  max_product <- max(1L, min(as.integer(max_basis_columns), n_positive - 2L))

  if (basis == "legendre") {
    if (type == "knots") stop("'knots' is unavailable for Legendre bases.",
                               call. = FALSE)
    if (length(raw) != 1L || !is.finite(raw) || raw < 1) {
      stop("Legendre sieve_size must be a positive scalar.", call. = FALSE)
    }
    # Every continuous variable enters at least linearly: the minimum
    # Legendre tensor dimension is 2^d, as the spline minimum is s0^d.
    minimum_product <- 2L^d
    if (max_product < minimum_product) {
      stop("Cell is too small for a tensor-product Legendre basis with a ",
           "linear term in every continuous variable.", call. = FALSE)
    }
    requested <- max(1L, as.integer(floor(raw)))
    dims <- .gir_balanced_dims(requested, d, 2L, max_product)
    actual <- prod(dims)
    return(list(dims = dims, knots = integer(0), requested = requested,
                actual = actual,
                capped = actual < max(requested, minimum_product),
                floored = requested < minimum_product,
                type = "dimension"))
  }

  if (type == "dimension") {
    if (length(raw) != 1L || !is.finite(raw) || raw < 1) {
      stop("Dimension sieve_size must be a positive scalar.", call. = FALSE)
    }
    minimum_product <- spline_order^d
    if (max_product < minimum_product) {
      stop("Cell is too small for a tensor-product spline basis of this order.",
           call. = FALSE)
    }
    requested <- max(1L, as.integer(floor(raw)))
    dims <- .gir_balanced_dims(requested, d, spline_order, max_product)
    actual <- prod(dims)
    return(list(dims = dims, knots = dims - spline_order,
                requested = requested, actual = actual,
                capped = actual < max(requested, minimum_product),
                floored = requested < minimum_product,
                type = "dimension"))
  }

  if (length(raw) == 1L) knots <- rep(as.integer(floor(raw)), d)
  else if (length(raw) == d) knots <- as.integer(floor(raw))
  else stop("For knot sizing, sieve_size must be scalar or length(continuous).",
            call. = FALSE)
  if (any(!is.finite(knots)) || any(knots < 0L)) {
    stop("Spline knot counts must be finite nonnegative integers.",
         call. = FALSE)
  }
  requested <- knots
  if (max_product < spline_order^d) {
    stop("Cell is too small for a tensor-product spline basis of this order.",
         call. = FALSE)
  }
  knots <- .gir_reduce_knots(knots, spline_order, max_product)
  dims <- knots + spline_order
  list(dims = dims, knots = knots,
       requested = if (length(requested) == 1L) requested else
         paste(requested, collapse = ","),
       actual = prod(dims), capped = any(knots < requested),
       floored = FALSE, type = "knots")
}

.gir_make_basis <- function(scaled, continuous, rows, spec, basis,
                            spline_order, max_basis_columns) {
  if (!length(continuous)) {
    P <- matrix(1, nrow = length(rows), ncol = 1L)
    colnames(P) <- "constant"
    return(P)
  }
  mats <- lapply(seq_along(continuous), function(j) {
    v <- continuous[j]; z <- scaled[[v]][rows]
    if (basis == "legendre") .gir_legendre(z, spec$dims[j], v)
    else .gir_spline(z, spline_order, spec$knots[j], v)
  })
  .gir_tensor(mats, max_basis_columns)
}


# =============================================================================
# Fast conditional mid-CDF first stage
# =============================================================================

.gir_col_cumsum <- function(A) {
  A <- as.matrix(A)
  out <- matrix(0, nrow(A), ncol(A))
  for (j in seq_len(ncol(A))) out[, j] <- cumsum(A[, j])
  out
}

.gir_diagnostics <- function(A, design, tolerance) {
  rc <- tryCatch(base::rcond(A), error = function(e) 0)
  rank <- qr(design, tol = tolerance)$rank
  list(rcond = rc, rank = rank,
       singular = !is.finite(rc) || rc < tolerance || rank < ncol(design))
}

.gir_ridge_adjustment <- function(gram, original, ridge, label, multiplier = 1) {
  # A common rescaling of the weights must also rescale the penalty.
  scale <- mean(diag(original))
  if (!is.finite(scale) || scale <= 0) {
    stop("Cannot apply ridge in cell ", label, ": the average Gram diagonal ",
         "must be finite and positive.", call. = FALSE)
  }
  increment <- (ridge * scale) * multiplier
  if (!is.finite(increment) || increment <= 0) {
    stop("Ridge penalty underflowed or overflowed in cell ", label,
         "; rescale the weights or adjust ridge.", call. = FALSE)
  }
  adjusted <- gram + diag(increment, nrow(gram))
  if (any(!is.finite(adjusted))) {
    stop("Ridge-adjusted Gram matrix contains non-finite values in cell ",
         label, "; rescale the weights or adjust ridge.", call. = FALSE)
  }
  list(gram = adjusted, increment = increment)
}

.gir_fit_cell <- function(y, P, weights, label, singular_method, tolerance,
                          ridge, order_info = NULL) {
  gram_original <- crossprod(P, P * weights)
  weighted_design <- P * sqrt(weights)
  diag <- .gir_diagnostics(gram_original, weighted_design, tolerance)

  gram <- gram_original
  P_used <- P
  dropped <- character(0)
  action <- "none"
  lambda <- 0

  if (diag$singular) {
    if (singular_method == "error") {
      return(list(ok = FALSE, label = label, reason = "singular or ill-conditioned",
                  rank = diag$rank, dimension = ncol(P), rcond = diag$rcond))
    }
    if (singular_method == "drop") {
      q <- qr(weighted_design, tol = tolerance)
      keep <- q$pivot[seq_len(q$rank)]
      omit <- setdiff(seq_len(ncol(P)), keep)
      if (!length(omit)) {
        return(list(
          ok = FALSE, label = label,
          reason = paste0("ill-conditioned but QR found no dependent column; ",
                          "use first_stage_singular='ridge' or 'none'"),
          rank = q$rank, dimension = ncol(P), rcond = diag$rcond
        ))
      }
      dropped <- colnames(P)[omit]
      P_used <- P[, keep, drop = FALSE]
      gram <- crossprod(P_used, P_used * weights)
      action <- "drop"
    }
    if (singular_method == "ridge") {
      adjusted <- .gir_ridge_adjustment(gram, gram_original, ridge, label)
      lambda <- adjusted$increment
      gram <- adjusted$gram
      action <- "ridge"
    }
  }

  if (is.null(order_info)) {
    ord <- order(y)
    y_ord <- y[ord]
    tie <- rle(y_ord)
    tie_ends <- cumsum(tie$lengths)
    tie_starts <- c(1L, utils::head(tie_ends, -1L) + 1L)
    first_position <- rep(tie_starts, tie$lengths)
    last_position <- rep(tie_ends, tie$lengths)
  } else {
    # ADDED: Bayesian draws reuse the same ordering and exact-tie groups.
    ord <- order_info$order
    first_position <- order_info$first
    last_position <- order_info$last
  }
  P_ord <- P_used[ord, , drop = FALSE]

  cumulative <- .gir_col_cumsum(P_ord * weights[ord])
  with_zero <- rbind(rep(0, ncol(P_ord)), cumulative)
  strictly_below <- with_zero[first_position, , drop = FALSE]
  weakly_below <- cumulative[last_position, , drop = FALSE]
  threshold_scores <- 0.5 * (strictly_below + weakly_below)

  coefs <- tryCatch(
    solve(gram, t(threshold_scores)),
    error = function(e) {
      if (singular_method != "ridge") {
        return(structure(list(message = conditionMessage(e)),
                         class = "gir_solve_error"))
      }
      adjusted <- .gir_ridge_adjustment(gram, gram_original, ridge, label,
                                        multiplier = 100)
      result <- solve(adjusted$gram, t(threshold_scores))
      lambda <<- lambda + adjusted$increment
      action <<- "ridge"
      result
    }
  )
  if (inherits(coefs, "gir_solve_error")) {
    return(list(ok = FALSE, label = label, reason = coefs$message,
                rank = diag$rank, dimension = ncol(P), rcond = diag$rcond))
  }

  fitted_ord <- rowSums(P_ord * t(coefs))
  fitted <- numeric(length(y)); fitted[ord] <- fitted_ord
  list(
    ok = TRUE, fitted = fitted, rank = if (action == "drop") ncol(P_used) else diag$rank,
    dimension = ncol(P), effective_dimension = ncol(P_used),
    rcond = diag$rcond, singular = diag$singular, action = action,
    dropped = dropped, ridge_lambda = lambda
  )
}

.gir_estimate_midcdf <- function(y, fs_data, weights, continuous, discrete,
                                 basis, kappa, basis_constant, sieve_size,
                                 sieve_size_type, sieve_size_reference,
                                 spline_order, bounds, first_stage_singular,
                                 first_stage_tol, ridge, minimum_cell_size,
                                 max_basis_columns) {
  n <- length(y)
  cells <- .gir_cells(fs_data, discrete)
  scaled <- .gir_scale(fs_data, continuous, bounds)
  fitted <- rep(NA_real_, n)
  diagnostics <- vector("list", nlevels(cells))
  failures <- character(0)
  modified <- character(0)
  capped <- character(0)
  degenerate <- character(0)

  for (g in seq_len(nlevels(cells))) {
    label <- levels(cells)[g]
    rows <- which(cells == label)
    n_positive <- sum(weights[rows] > 0)
    if (n_positive < minimum_cell_size) {
      failures <- c(failures, paste0(label, " [", length(rows), " rows; ",
        n_positive, " positive-weight rows; minimum ", minimum_cell_size, "]"))
      next
    }

    spec <- tryCatch(.gir_size_spec(
      n_cell = length(rows), n_positive = n_positive, n_total = n,
      d = length(continuous), basis = basis, kappa = kappa,
      basis_constant = basis_constant, sieve_size = sieve_size,
      sieve_size_type = sieve_size_type,
      sieve_size_reference = sieve_size_reference,
      spline_order = spline_order, max_basis_columns = max_basis_columns
    ), error = function(e) e)
    if (inherits(spec, "error")) {
      failures <- c(failures, paste0(label, " [", conditionMessage(spec), "]"))
      next
    }
    if (isTRUE(spec$capped)) capped <- c(capped, label)

    if (length(continuous)) {
      # A marginal basis with more columns than the variable has distinct
      # values in the cell is rank deficient by construction. Name the
      # cause instead of reporting a generic singular Gram matrix.
      positive <- rows[weights[rows] > 0]
      n_distinct <- vapply(continuous, function(v) {
        length(unique(scaled[[v]][positive]))
      }, integer(1))
      short <- which(spec$dims > n_distinct)
      if (length(short)) {
        detail <- paste0(continuous[short], " has ", n_distinct[short],
                         " distinct value(s) but a marginal basis of ",
                         "dimension ", spec$dims[short], collapse = "; ")
        if (first_stage_singular == "error") {
          failures <- c(failures, paste0(label, " [", detail,
            "; declare the variable discrete or reduce the basis]"))
          next
        }
        degenerate <- c(degenerate, paste0(detail, " in ", label))
      }
    }

    P <- tryCatch(.gir_make_basis(
      scaled, continuous, rows, spec, basis, spline_order, max_basis_columns
    ), error = function(e) e)
    if (inherits(P, "error")) {
      failures <- c(failures, paste0(label, " [", conditionMessage(P), "]"))
      next
    }

    cell_fit <- .gir_fit_cell(
      y[rows], P, weights[rows], label, first_stage_singular,
      first_stage_tol, ridge
    )
    if (!isTRUE(cell_fit$ok)) {
      failures <- c(failures, paste0(
        label, " [rank ", cell_fit$rank, "/", cell_fit$dimension,
        "; rcond=", format(cell_fit$rcond, digits = 4),
        "; ", cell_fit$reason, "]"
      ))
      next
    }

    fitted[rows] <- cell_fit$fitted
    if (cell_fit$action %in% c("ridge", "drop")) modified <- c(modified, label)
    marginal_dimensions <- stats::setNames(
      as.integer(spec$dims), continuous
    )
    marginal_knot_counts <- if (basis == "spline") {
      stats::setNames(as.integer(spec$knots), continuous)
    } else {
      stats::setNames(rep(NA_integer_, length(continuous)), continuous)
    }

    diagnostics[[g]] <- data.frame(
      cell = label,
      observations = length(rows),
      positive_weight_observations = n_positive,
      requested_size = as.character(spec$requested),
      size_type = spec$type,
      marginal_dimensions = I(list(marginal_dimensions)),
      marginal_knot_counts = I(list(marginal_knot_counts)),
      basis_dimension = ncol(P),
      effective_dimension = cell_fit$effective_dimension,
      rank = cell_fit$rank,
      rcond = cell_fit$rcond,
      singular = cell_fit$singular,
      action = cell_fit$action,
      ridge_lambda = cell_fit$ridge_lambda,
      dropped_columns = if (length(cell_fit$dropped))
        paste(cell_fit$dropped, collapse = "; ") else "",
      capped = spec$capped,
      floored = isTRUE(spec$floored),
      stringsAsFactors = FALSE
    )
  }

  if (length(failures)) {
    stop(
      "First-stage conditional-mid-CDF estimation failed in:\n  - ",
      paste(failures, collapse = "\n  - "),
      "\nReduce the sieve size, combine sparse categories, change basis, or ",
      "choose first_stage_singular='ridge' or 'drop'.",
      call. = FALSE
    )
  }
  diagnostics <- do.call(rbind, diagnostics); rownames(diagnostics) <- NULL
  if (length(capped)) {
    warning("The first-stage basis was reduced in: ",
            paste(unique(capped), collapse = "; "),
            ". Inspect fit$first_stage$diagnostics.", call. = FALSE)
  }
  if (length(degenerate)) {
    warning("A marginal basis exceeds the number of distinct values of its ",
            "variable, so the cell Gram matrix is singular by construction ",
            "and was handled with '", first_stage_singular, "': ",
            paste(degenerate, collapse = "; "),
            ". Declare the variable discrete or reduce the basis.",
            call. = FALSE)
  }
  if (length(modified)) {
    warning("First-stage singularity/ill-conditioning was handled with '",
            first_stage_singular, "' in: ",
            paste(unique(modified), collapse = "; "),
            ". Inspect fit$first_stage$diagnostics.", call. = FALSE)
  }
  list(values = fitted, cells = cells, diagnostics = diagnostics)
}

.gir_maybe_clip <- function(x, clip, source) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop("The first-stage mid-CDF vector contains non-finite values.",
         call. = FALSE)
  }
  outside <- which(x < 0 | x > 1)
  original_range <- range(x)
  if (length(outside)) {
    if (clip) {
      warning(length(outside), " mid-CDF value(s) from ", source,
              " were outside [0,1] and were clipped; original range [",
              format(original_range[1], digits = 5), ", ",
              format(original_range[2], digits = 5), "].", call. = FALSE)
      x <- pmin(1, pmax(0, x))
    } else {
      warning(length(outside), " mid-CDF value(s) from ", source,
              " are outside [0,1]; range [",
              format(original_range[1], digits = 5), ", ",
              format(original_range[2], digits = 5),
              "]. Set clip=TRUE to clip them.", call. = FALSE)
    }
  }
  list(values = x, outside = length(outside),
       original_range = original_range, clipped = clip && length(outside) > 0L)
}


# =============================================================================
# Second-stage point estimator and GIR-R^2
# =============================================================================

.gir_row_groups <- function(frame) {
  # Index the distinct rows of a regressor frame in order of first appearance,
  # together with the position of the first row of each group. Observations
  # with the same frame row receive the same design row, so the second stage
  # can accumulate its sums over the groups. Matrix columns (poly(), spline
  # bases) are split into their columns, factors use their codes, and other
  # columns are matched on their values.
  n <- nrow(frame)
  key <- rep.int(1L, n)
  for (v in seq_along(frame)) {
    column <- frame[[v]]
    parts <- if (is.factor(column)) list(as.integer(column)) else {
      column <- unclass(column)
      if (is.matrix(column)) {
        lapply(seq_len(ncol(column)), function(j) column[, j])
      } else list(column)
    }
    for (part in parts) {
      code <- match(part, unique(part))
      # key and code are at most n, so the combined key is exact in doubles.
      combined <- as.numeric(key - 1L) * max(code) + (code - 1L)
      key <- match(combined, unique(combined))
    }
  }
  list(index = key, first = match(seq_len(max(key)), key))
}

.gir_expand_design <- function(X, group, row_names) {
  full <- X[group, , drop = FALSE]
  dimnames(full) <- list(row_names, colnames(X))
  attr(full, "assign") <- attr(X, "assign")
  attr(full, "contrasts") <- attr(X, "contrasts")
  full
}

.gir_second_stage <- function(X, y, z, weights, singular.ok, tolerance,
                              group = NULL) {
  # X holds one row per distinct regressor row and `group` maps each of the n
  # observations to its row; without `group` the rows are the observations.
  n <- length(y); p <- ncol(X)
  if (is.null(group)) group <- seq_len(n)
  wy <- weights * y
  # The moments depend on X_i only through its value, so the products
  # a_i Y_i and a_i Y_i (2 M_i - 1) are first summed within each group and
  # the moments are then formed from the distinct rows. This is exact; only
  # the order of summation differs from the observation-level formulas.
  sums <- rowsum(cbind(wy, wy * z), group)
  S <- sums[, 1L]
  Tz <- sums[, 2L]
  A <- crossprod(X, X * S) / n
  b <- colSums(X * Tz) / n
  rc <- tryCatch(base::rcond(A), error = function(e) 0)

  nonnegative <- all(wy >= 0)
  design_w <- if (nonnegative) X * sqrt(S) else NULL
  q <- if (nonnegative) qr(design_w, tol = tolerance) else qr(A, tol = tolerance)
  rank <- q$rank
  rank_deficient <- rank < p
  ill_conditioned <- is.finite(rc) && rc < tolerance
  singular <- rank_deficient || !is.finite(rc)

  cn <- colnames(X)
  coefficients <- setNames(rep(NA_real_, p), cn)
  complete <- setNames(rep(0, p), cn)
  aliased <- setNames(rep(FALSE, p), cn)

  if (!singular) {
    solved <- tryCatch(solve(A, b), error = function(e) e)
    if (inherits(solved, "error")) {
      singular <- TRUE; rank_deficient <- TRUE
    } else {
      coefficients[] <- as.numeric(solved)
      complete[] <- as.numeric(solved)
    }
  }

  if (singular) {
    if (!singular.ok) {
      stop(
        "The second-stage matrix sum_i w_i Y_i X_i X_i' is singular. ",
        "Set singular.ok=TRUE for an lm-like aliased fit.", call. = FALSE
      )
    }
    if (nonnegative) {
      q <- qr(design_w, tol = tolerance); rank <- q$rank
      if (rank < 1L) stop("The second stage has zero numerical rank.",
                          call. = FALSE)
      keep <- q$pivot[seq_len(rank)]
      omit <- setdiff(seq_len(p), keep)
      reduced_qr <- qr(design_w[, keep, drop = FALSE], tol = tolerance)
      # Weighted least squares of the group means of 2M-1 with the summed
      # weights solves the same normal equations as the observation-level fit.
      group_mean <- ifelse(S > 0, Tz / S, 0)
      beta <- qr.coef(reduced_qr, sqrt(S) * group_mean)
      coefficients[keep] <- as.numeric(beta)
      complete[keep] <- as.numeric(beta)
      aliased[omit] <- TRUE
    } else {
      q <- qr(A, tol = tolerance); rank <- q$rank
      if (rank < 1L) stop("The second-stage moment matrix has zero rank.",
                          call. = FALSE)
      keep <- q$pivot[seq_len(rank)]
      omit <- setdiff(seq_len(p), keep)
      beta <- tryCatch(
        solve(A[keep, keep, drop = FALSE], b[keep]),
        error = function(e) qr.solve(A[keep, keep, drop = FALSE], b[keep],
                                     tol = tolerance)
      )
      coefficients[keep] <- as.numeric(beta)
      complete[keep] <- as.numeric(beta)
      aliased[omit] <- TRUE
    }
    warning("The second-stage design is rank deficient; ", sum(aliased),
            " coefficient(s) are aliased and reported as NA.", call. = FALSE)
  } else if (ill_conditioned) {
    warning("The second-stage matrix is full rank but ill-conditioned (rcond=",
            format(rc, digits = 4), "). Estimates were not regularized.",
            call. = FALSE)
  }

  fitted <- unname(drop(X %*% complete))[group]
  list(coefficients = coefficients, complete = complete,
       fitted = fitted, residuals = z - fitted, rank = rank,
       df.residual = n - rank, aliased = aliased, singular = singular,
       ill.conditioned = ill_conditioned, rcond = rc, A = A, b = b)
}

.gir_r2 <- function(y, z, fitted, weights, has_intercept) {
  wy <- weights * y
  total_weight <- sum(wy)
  if (!is.finite(total_weight) || abs(total_weight) <= .Machine$double.eps) {
    warning("GIR-R^2 is undefined because sum(weights*Y) is zero or non-finite.",
            call. = FALSE)
    return(list(value = NA_real_, restricted_coefficient = NA_real_,
                unrestricted_loss = NA_real_, restricted_loss = NA_real_))
  }
  restricted_coefficient <- sum(wy * z) / total_weight
  unrestricted_loss <- sum(wy * (z - fitted)^2)
  restricted_loss <- sum(wy * (z - restricted_coefficient)^2)
  value <- if (!is.finite(restricted_loss) ||
               abs(restricted_loss) <= .Machine$double.eps) {
    warning("GIR-R^2 is undefined because the restricted loss is zero.",
            call. = FALSE)
    NA_real_
  } else 1 - unrestricted_loss / restricted_loss
  if (!has_intercept) {
    warning("GIR-R^2 was requested without an intercept. The constant-only ",
            "benchmark is not nested, so GIR-R^2 need not lie in [0,1].",
            call. = FALSE)
  }
  list(value = value, restricted_coefficient = restricted_coefficient,
       unrestricted_loss = unrestricted_loss, restricted_loss = restricted_loss)
}


# =============================================================================
# Public estimator
# =============================================================================

#' Linear Gini Index Regression
#'
#' Fit the linear GIR point estimator with a tie-aware least-squares-series
#' conditional mid-distribution first stage and optional GIR goodness of fit.
#'
#' `gir()` computes point estimates for the linear GIR model. Its formula/data
#' interface is patterned after `lm()` and `quantreg::rq()`. The conditional
#' mid-CDF is either estimated with least-squares series or supplied by the
#' user.
#'
#' @aliases gir.object
#' @param formula Two-sided second-stage formula. An intercept is included
#'   unless removed with `0 +` or `-1`. Offsets are not implemented.
#' @param data Data frame containing the needed variables. Required.
#' @param subset Optional row selector evaluated in `data`, before
#'   missing-value handling.
#' @param weights Nonnegative finite case weights entering both stages. At
#'   least one must be positive in the retained estimation sample, after
#'   subsetting and missing-value handling. A scalar is repeated; vectors align
#'   to original or selected rows. Missing weights are handled by `na.action`
#'   together with the other model variables.
#' @param na.action Missing-value action for variables in either stage and
#'   supplied mid-CDFs. Defaults to `stats::na.fail`; `na.omit` or
#'   `na.exclude` removes incomplete rows with a warning. Remaining missing
#'   values error. If no rows or no positive weights remain, fitting stops
#'   with a sample-specific error before the routine removal warning. Custom
#'   actions must return a data frame selecting original rows with valid,
#'   unique row identifiers. Their removal warning distinguishes missing and
#'   complete rows, and their errors are reported as custom-action failures.
#' @param first_stage One-sided formula or character vector selecting raw
#'   conditioning columns. `NULL` uses underlying variables from the
#'   second-stage right-hand side.
#' @param continuous Character vector of numeric continuous first-stage
#'   columns. The full tensor product of their marginal bases is used.
#'   Numeric first-stage columns listed in neither `continuous` nor
#'   `discrete` are treated as continuous, and a message names them. An
#'   undeclared numeric column with at most two distinct values in the retained
#'   estimation sample is an error, because a two-point support cannot carry
#'   a continuous basis; declare it in `discrete`, or in `continuous` if that
#'   is intended.
#' @param discrete Character vector of categorical first-stage columns.
#'   Observed joint values form separate cells. Unclassified factors, ordered
#'   factors, character and logical columns are automatic discrete variables;
#'   numeric binary columns must be declared explicitly.
#' @param midcdf Optional supplied numeric vector or data-column name for
#'   conditional mid-distribution values. Bypasses internal first-stage
#'   estimation. A scalar, original-row vector, or selected-row vector is
#'   aligned before missing-value removal.
#' @param basis `"legendre"` or `"spline"`. The spline sieve uses
#'   `splines::bs()` internally; there is no separate B-spline or wavelet
#'   option.
#' @param kappa Exponent in the automatic sieve-size rule. `NULL` selects
#'   0.245 for Legendre and 0.300 for splines.
#' @param basis_constant Positive multiplier in
#'   `floor(basis_constant * n_reference^kappa)`.
#' @param sieve_size Optional fixed total-dimension target, or scalar/vector
#'   interior-knot count according to `sieve_size_type`. A knot vector follows
#'   `continuous` ordering. `NULL` uses the automatic rule. A
#'   marginal basis with more columns than the variable has distinct values
#'   in a cell is singular by construction: an error under
#'   `first_stage_singular = "error"`, otherwise a warning naming the
#'   variable while the chosen policy is applied.
#' @param sieve_size_type `"auto"` is `"dimension"` for both bases: the size
#'   rule gives the total tensor-product dimension \eqn{K} of the paper,
#'   allocated greedily across balanced integer marginal dimensions, never
#'   below two columns per Legendre margin or `spline_order` columns per
#'   spline margin. `"knots"` (splines only) applies the rule, or
#'   `sieve_size`, to the interior-knot count of every continuous margin
#'   instead.
#' @param sieve_size_reference Use each cell size (`"cell"`) or full
#'   estimation-sample size (`"full"`) in the automatic rule. Regression
#'   observations remain cell-specific in either case.
#' @param spline_order Order of marginal polynomial splines. Default 4 means
#'   cubic, with degree 3. Each marginal basis has \eqn{J+s_0} columns for
#'   \eqn{J} interior knots and order \eqn{s_0}.
#' @param basis_bounds Bounds for rescaling continuous columns to
#'   \eqn{[0,1]}, preferably a named list of two-number vectors. Default uses
#'   final estimation-sample minima and maxima. Bounds are common across
#'   cells.
#' @param clip Clip first-stage values to \eqn{[0,1]}? Default `FALSE`.
#'   Out-of-range values always warn. Clipping changes the estimator; it does
#'   not constrain the second-stage linear index.
#' @param first_stage_singular `"error"`, `"ridge"`, `"drop"`, or `"none"`.
#'   Error is the conservative default; ridge adds a scaled diagonal
#'   correction, drop removes QR-dependent columns, and none attempts the
#'   unmodified solve.
#' @param first_stage_tol Positive first-stage numerical
#'   rank/reciprocal-condition tolerance.
#' @param ridge Relative diagonal regularization constant, default `1e-10`,
#'   used only with `first_stage_singular = "ridge"`. The initial penalty is
#'   `ridge * mean(diag(G))` for the weighted cell Gram matrix `G`, without an
#'   absolute floor. If the solve fails, a further 100 times this base penalty
#'   is added. Diagnostics report the total applied penalty. The default is a
#'   numerical heuristic, not a data-tuned or theoretically optimal choice.
#'   Unusable scales, penalties, or adjusted matrices cause an informative
#'   error. Ridge changes the first-stage estimator; assess its sensitivity.
#' @param minimum_cell_size Minimum number of positive-weight observations per
#'   observed cell. Default 8.
#' @param max_basis_columns Upper limit on constructed tensor-product columns.
#'   Cell-size constraints can reduce the requested allocation.
#' @param singular.ok Allow second-stage aliasing? `TRUE` retains an
#'   identified subset and marks aliased coefficients `NA`. `FALSE` errors.
#'   The second stage is not ridge-regularized.
#' @param second_stage_tol Positive second-stage numerical rank tolerance and
#'   conditioning warning threshold.
#' @param GIR Compute the GIR-specific goodness-of-fit measure using the same
#'   first-stage ranks in unrestricted and constant-only losses? Default
#'   `FALSE`.
#' @param model,x,y Retain public model-frame, design-matrix and response
#'   slots. Defaults `TRUE`, `FALSE`, `FALSE`. Independently, an internal
#'   resolved estimation snapshot is always retained for bootstrap and
#'   prediction; these flags do not anonymize the object.
#' @param contrasts Optional contrasts passed to `stats::model.matrix()`.
#'
#' @details
#' The first-stage response at threshold \eqn{t} is
#' \eqn{1(Y<t)+0.5\,1(Y=t)}. Sorting and cumulative basis sums calculate all
#' observed-threshold scores exactly, including ties. With case weights
#' \eqn{a_i}, the second stage solves
#' \deqn{\left\{\sum_i a_iY_iX_iX_i'\right\}\widehat\theta
#'        =\sum_i a_iY_iX_i(2\widehat M_i-1).}
#' For splines, equally spaced interior knots are
#' specified explicitly and a stable B-spline representation is used. Under
#' total-dimension sizing, integer products need not equal the target; the
#' minimum spline tensor dimension \eqn{s_0^d} and the minimum Legendre
#' tensor dimension \eqn{2^d} can exceed a small nominal target, in which
#' case the minimum is used and `floored` is recorded in the diagnostics.
#' Cell-size capping and QR dropping are reported.
#'
#' Observations with the same regressor values share a design row. The
#' second-stage sums are accumulated over those groups and the moments are
#' formed from the distinct rows, which is exact (only the order of
#' summation differs). With cell-level regressors, such as industry and year
#' effects in firm data, the retained design and every bootstrap draw are
#' therefore the size of the distinct rows rather than of the sample.
#'
#' Use `summary(fit)` for the point fit and `summary(bboot.gir(fit))` for
#' bootstrap inference. Fitted values target the linear conditional Gini
#' index. See `vignette("GIR-guide", package = "GIR")` for examples,
#' the underlying calculations and guidance on numerical settings.
#'
#' @return An object with classes `gir` and `gir.object`, containing named
#'   coefficients, complete coefficients for fitted values, fitted indices and
#'   generated-rank residuals, weights, first-stage mid-CDF values, formula
#'   metadata, row alignment, and diagnostics. `first_stage$diagnostics`
#'   contains named list-columns `marginal_dimensions` and
#'   `marginal_knot_counts`, total `basis_dimension`, and
#'   `effective_dimension` after any column removal, and `floored` records
#'   whether the minimum tensor dimension exceeded the requested size. Knots
#'   are `NA` for
#'   Legendre; an intercept-only first stage has empty marginal vectors. The
#'   `estimation` component retains the distinct design rows with the index
#'   mapping observations to them, the outcome, weights, first-stage columns
#'   and controls for later refits; `x = TRUE` stores the full design.
#'
#' @references
#' Bhardwaj, V. and Diaz-Villarejo, D. (2026). Gini Index Regression. Working
#' paper, Universidad Carlos III de Madrid.
#' \url{https://ddiazvil.com/assets/papers/JMP_GIR_DiazVillarejo.pdf}.
#'
#' @seealso [bboot.gir()], [predict.gir()], [midcdf()],
#'   \code{\link{GIR-validation}} for mathematical checks and their scope.
#'
#' @examples
#' d <- data.frame(y = rep(c(1, 2, 2, 4, 5), 8),
#'                 group = factor(rep(c("a", "b"), each = 20)))
#' f <- gir(y ~ group, data = d, discrete = "group", GIR = TRUE)
#' summary(f)
#' summary(midcdf(f))
#' f$first_stage$diagnostics
#'
#' @importFrom stats coef confint fitted formula model.frame model.matrix nobs
#'   predict residuals setNames terms vcov weights
#' @export

gir <- function(
    formula,
    data,
    subset,
    weights,
    na.action = stats::na.fail,
    first_stage = NULL,
    continuous = NULL,
    discrete = NULL,
    midcdf = NULL,
    basis = c("legendre", "spline"),
    kappa = NULL,
    basis_constant = 1,
    sieve_size = NULL,
    sieve_size_type = c("auto", "dimension", "knots"),
    sieve_size_reference = c("cell", "full"),
    spline_order = 4L,
    basis_bounds = NULL,
    clip = FALSE,
    first_stage_singular = c("error", "ridge", "drop", "none"),
    first_stage_tol = 1e-12,
    ridge = 1e-10,
    minimum_cell_size = 8L,
    max_basis_columns = 1000L,
    singular.ok = TRUE,
    second_stage_tol = 1e-7,
    GIR = FALSE,
    model = TRUE,
    x = FALSE,
    y = FALSE,
    contrasts = NULL) {

  call <- match.call()
  caller <- parent.frame()
  subset_expr <- substitute(subset)
  weights_expr <- substitute(weights)
  midcdf_expr <- substitute(midcdf)

  for (arg in c("clip", "singular.ok", "GIR", "model", "x", "y")) {
    .gir_flag(get(arg), arg)
  }
  basis <- match.arg(basis)
  sieve_size_type <- match.arg(sieve_size_type)
  sieve_size_reference <- match.arg(sieve_size_reference)
  first_stage_singular <- match.arg(first_stage_singular)

  if (!is.null(sieve_size) && (!is.numeric(sieve_size) || anyNA(sieve_size) ||
                               any(!is.finite(sieve_size)))) {
    stop("sieve_size must be NULL or a finite numeric vector.", call. = FALSE)
  }
  if (is.null(kappa)) kappa <- if (basis == "legendre") 0.245 else 0.300
  .gir_positive(kappa, "kappa", zero_ok = TRUE)
  .gir_positive(basis_constant, "basis_constant")
  .gir_positive(first_stage_tol, "first_stage_tol")
  .gir_positive(ridge, "ridge", zero_ok = TRUE)
  .gir_positive(second_stage_tol, "second_stage_tol")

  spline_order <- as.integer(spline_order)
  minimum_cell_size <- as.integer(minimum_cell_size)
  max_basis_columns <- as.integer(max_basis_columns)
  if (length(spline_order) != 1L || is.na(spline_order) || spline_order < 2L) {
    stop("spline_order must be an integer of at least 2.", call. = FALSE)
  }
  if (length(minimum_cell_size) != 1L || is.na(minimum_cell_size) ||
      minimum_cell_size < 1L) {
    stop("minimum_cell_size must be a positive integer.", call. = FALSE)
  }
  if (length(max_basis_columns) != 1L || is.na(max_basis_columns) ||
      max_basis_columns < 1L) {
    stop("max_basis_columns must be a positive integer.", call. = FALSE)
  }
  if (first_stage_singular == "ridge" && ridge <= 0) {
    stop("ridge must be positive when first_stage_singular='ridge'.",
         call. = FALSE)
  }
  if (!inherits(formula, "formula")) stop("formula must be a formula.",
                                           call. = FALSE)
  if (missing(data)) stop("data must be supplied.", call. = FALSE)

  data_full <- as.data.frame(data)
  n_original <- nrow(data_full)
  if (!n_original) stop("data has no rows.", call. = FALSE)

  selected <- if (missing(subset)) seq_len(n_original) else {
    value <- eval(subset_expr, data_full, caller)
    .gir_subset_index(value, n_original, rownames(data_full))
  }
  data_selected <- data_full[selected, , drop = FALSE]
  n_selected <- nrow(data_selected)

  obs_weights <- if (missing(weights)) rep(1, n_selected) else {
    raw <- eval(weights_expr, data_full, caller)
    if (!is.numeric(raw) || !is.null(dim(raw))) {
      stop("weights must evaluate to a numeric vector.", call. = FALSE)
    }
    as.numeric(.gir_align(raw, selected, n_original, n_selected, "weights"))
  }
  # Missing weights are left to na.action, like any other model variable;
  # only infinite weights are rejected here.
  if (any(is.infinite(obs_weights))) {
    stop("weights contains infinite values.", call. = FALSE)
  }
  if (any(obs_weights < 0, na.rm = TRUE)) {
    stop("weights must be nonnegative.", call. = FALSE)
  }
  supplied <- NULL
  supplied_flag <- !missing(midcdf) && !identical(midcdf_expr, quote(NULL))
  if (supplied_flag) {
    raw <- eval(midcdf_expr, data_full, caller)
    if (is.character(raw) && length(raw) == 1L && raw %in% names(data_full)) {
      raw <- data_full[[raw]]
    }
    if (!is.numeric(raw) || !is.null(dim(raw))) {
      stop("midcdf must evaluate to a numeric vector or a numeric column name.",
           call. = FALSE)
    }
    supplied <- as.numeric(.gir_align(raw, selected, n_original, n_selected,
                                      "midcdf"))
  }

  mf0 <- stats::model.frame(formula, data = data_selected,
                            na.action = stats::na.pass,
                            drop.unused.levels = TRUE)
  tt <- attr(mf0, "terms")
  if (is.null(tt) || attr(tt, "response") != 1L) {
    stop("formula must contain exactly one response.", call. = FALSE)
  }
  if (length(attr(tt, "offset"))) {
    stop("Offsets are not implemented in gir().", call. = FALSE)
  }

  fs_spec <- if (supplied_flag) {
    list(all = character(0), continuous = character(0), discrete = character(0))
  } else .gir_fs_variables(tt, first_stage, continuous, discrete, data_selected)
  fs0 <- if (length(fs_spec$all)) data_selected[fs_spec$all] else
    data.frame(row.names = seq_len(n_selected))

  na_result <- .gir_apply_na(mf0, fs0, obs_weights, supplied, na.action)
  keep <- na_result$keep
  data_used <- data_selected[keep, , drop = FALSE]
  obs_weights <- obs_weights[keep]
  if (!is.null(supplied)) supplied <- supplied[keep]
  row_index <- selected[keep]

  if (!any(obs_weights > 0)) {
    stop("At least one weight must be positive in the retained estimation ",
         "sample after ", na_result$name, ".", call. = FALSE)
  }
  if (!supplied_flag) .gir_check_fs_support(fs_spec, data_used)
  if (length(na_result$removed)) {
    if (na_result$name == "custom na.action") {
      warning(length(na_result$removed), " observation(s) were removed by ",
              "custom na.action (", na_result$removed_incomplete,
              " with missing values, ",
              length(na_result$removed) - na_result$removed_incomplete,
              " complete).", call. = FALSE)
    } else {
      warning(length(na_result$removed), " observation(s) with missing values ",
              "were removed by ", na_result$name, ".", call. = FALSE)
    }
  }

  mf <- stats::model.frame(formula, data = data_used,
                           na.action = stats::na.pass,
                           drop.unused.levels = TRUE)
  tt <- attr(mf, "terms")
  outcome <- stats::model.response(mf)
  if (!is.numeric(outcome) || is.matrix(outcome) || is.data.frame(outcome)) {
    stop("The GIR response must be a numeric vector.", call. = FALSE)
  }
  outcome <- as.numeric(outcome)
  if (any(!is.finite(outcome))) stop("The response contains non-finite values.",
                                      call. = FALSE)
  if (any(outcome < 0)) {
    warning(sum(outcome < 0), " response value(s) are negative. The Gini ",
            "interpretation and usual GIR-R^2 properties require Y >= 0; ",
            "the moment estimator is nevertheless computed.", call. = FALSE)
  }

  # Observations with the same regressor values share a design row. The
  # design is built from the distinct rows only and the second stage sums
  # over the groups, which is exact and, with cell-level regressors, far
  # smaller than one row per observation.
  groups <- .gir_row_groups(mf[-1L])
  mf_unique <- mf[groups$first, , drop = FALSE]
  attr(mf_unique, "terms") <- tt
  X <- as.matrix(stats::model.matrix(tt, mf_unique, contrasts.arg = contrasts))
  if (!ncol(X)) stop("formula creates no second-stage regressors.",
                      call. = FALSE)
  if (any(!is.finite(X))) stop("The second-stage matrix is non-finite.",
                                call. = FALSE)

  if (supplied_flag) {
    processed <- .gir_maybe_clip(supplied, clip, "the supplied vector")
    mhat <- processed$values
    fs_output <- list(
      method = "supplied", variables = character(0), continuous = character(0),
      discrete = character(0), basis = NA_character_, kappa = NA_real_,
      basis_constant = NA_real_, sieve_size = NULL,
      sieve_size_type = NA_character_, sieve_size_reference = NA_character_,
      spline_order = NA_integer_, bounds = list(),
      diagnostics = data.frame(cell = "supplied mid-CDF",
                               observations = length(mhat)),
      outside_count = processed$outside, original_range = processed$original_range,
      clipped = processed$clipped
    )
  } else {
    fs_data <- data_used[fs_spec$all]
    for (v in fs_spec$continuous) {
      if (any(!is.finite(fs_data[[v]]))) {
        stop("Continuous first-stage column '", v, "' is non-finite.",
             call. = FALSE)
      }
    }
    for (v in fs_spec$discrete) {
      if (is.numeric(fs_data[[v]]) && any(!is.finite(fs_data[[v]]))) {
        stop("Discrete first-stage column '", v, "' is non-finite.",
             call. = FALSE)
      }
    }
    resolved_bounds <- .gir_bounds(fs_data, fs_spec$continuous, basis_bounds)
    estimated <- .gir_estimate_midcdf(
      y = outcome, fs_data = fs_data, weights = obs_weights,
      continuous = fs_spec$continuous, discrete = fs_spec$discrete,
      basis = basis, kappa = kappa, basis_constant = basis_constant,
      sieve_size = sieve_size, sieve_size_type = sieve_size_type,
      sieve_size_reference = sieve_size_reference,
      spline_order = spline_order, bounds = resolved_bounds,
      first_stage_singular = first_stage_singular,
      first_stage_tol = first_stage_tol, ridge = ridge,
      minimum_cell_size = minimum_cell_size,
      max_basis_columns = max_basis_columns
    )
    processed <- .gir_maybe_clip(estimated$values, clip, "series estimation")
    mhat <- processed$values
    fs_output <- list(
      method = "least-squares series conditional mid-CDF",
      variables = fs_spec$all, continuous = fs_spec$continuous,
      discrete = fs_spec$discrete, basis = basis, kappa = kappa,
      basis_constant = basis_constant, sieve_size = sieve_size,
      sieve_size_type = sieve_size_type,
      sieve_size_reference = sieve_size_reference,
      spline_order = if (basis == "spline") spline_order else NA_integer_,
      bounds = resolved_bounds, cells = estimated$cells,
      diagnostics = estimated$diagnostics,
      outside_count = processed$outside,
      original_range = processed$original_range,
      clipped = processed$clipped,
      singular_method = first_stage_singular,
      tolerance = first_stage_tol, ridge = ridge
    )
  }

  centered <- 2 * mhat - 1
  stage2 <- .gir_second_stage(X, outcome, centered, obs_weights,
                              singular.ok, second_stage_tol,
                              group = groups$index)
  names(stage2$fitted) <- rownames(mf)
  names(stage2$residuals) <- rownames(mf)
  has_intercept <- attr(tt, "intercept") == 1L
  r2 <- if (GIR) .gir_r2(outcome, centered, stage2$fitted, obs_weights,
                          has_intercept) else NULL

  factor_cols <- vapply(mf, is.factor, logical(1))
  xlevels <- lapply(mf[factor_cols], levels)

  fit <- list(
    coefficients = stage2$coefficients,
    coefficients.complete = stage2$complete,
    residuals = stage2$residuals,
    fitted.values = stage2$fitted,
    rank = stage2$rank,
    df.residual = stage2$df.residual,
    aliased = stage2$aliased,
    call = call, formula = formula, terms = tt,
    contrasts = attr(X, "contrasts"), xlevels = xlevels,
    assign = attr(X, "assign"), weights = obs_weights,
    midcdf = setNames(mhat, rownames(mf)),
    centered.midrank = setNames(centered, rownames(mf)),
    GIR.R2 = r2,
    first_stage = fs_output,
    second_stage = list(
      singular = stage2$singular,
      ill.conditioned = stage2$ill.conditioned,
      rcond = stage2$rcond, tolerance = second_stage_tol,
      moment_matrix = stage2$A, moment_vector = stage2$b
    ),
    na.action = na_result$object,
    nobs = length(outcome), original.nobs = n_original,
    subset.nobs = n_selected, row.index = row_index,
    model = if (model) mf else NULL,
    x = if (x) .gir_expand_design(X, groups$index, rownames(mf)) else NULL,
    y = if (y) outcome else NULL,
    outcome.name = names(mf)[1L], has.intercept = has_intercept,
    # ADDED: a resolved, immutable-by-convention snapshot for bboot.gir/predict.
    # This avoids re-evaluating a call against data that may since have changed.
    # X holds the distinct design rows; group maps observations to them.
    estimation = list(
      version = 2L, X = X, group = groups$index, y = outcome,
      weights = obs_weights,
      first_stage_data = if (supplied_flag) NULL else fs_data,
      predictor_variables = intersect(
        all.vars(stats::delete.response(tt)), names(data_used)
      ),
      predictor_levels = lapply(mf[vapply(mf, function(z) {
        is.factor(z) || is.character(z)
      }, logical(1))], function(z) {
        if (is.factor(z)) levels(z) else levels(factor(z))
      }),
      control = list(
        clip = clip, spline_order = spline_order,
        first_stage_singular = first_stage_singular,
        first_stage_tol = first_stage_tol, ridge = ridge,
        minimum_cell_size = minimum_cell_size,
        max_basis_columns = max_basis_columns,
        singular.ok = singular.ok, second_stage_tol = second_stage_tol
      )
    )
  )
  class(fit) <- c("gir", "gir.object")
  fit
}


# =============================================================================
# S3 methods
# =============================================================================

#' @rdname summary.gir
#' @export
print.gir <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nPoint coefficients:\n")
  print.default(format(x$coefficients, digits = digits), quote = FALSE,
                print.gap = 2L)
  if (any(x$aliased)) {
    cat("\n", sum(x$aliased), " coefficient(s) are aliased.\n", sep = "")
  }
  if (!is.null(x$GIR.R2)) {
    cat("\nGIR-R^2: ", format(x$GIR.R2$value, digits = digits), "\n", sep = "")
  }
  cat("\nFirst stage: ", x$first_stage$method,
      "\nObservations used: ", x$nobs, "\n", sep = "")
  invisible(x)
}

#' Summaries and Extraction for GIR Point Fits
#'
#' Print and summarize a GIR point fit and extract its model components.
#' Bootstrap inference remains a separate operation.
#'
#' @param object,x,formula A gir fit, or the summary object for its print
#'   method. In the `model.frame` method the argument named `formula` is the
#'   fit.
#' @param digits Print precision.
#' @param complete Return all coefficients including aliased `NA` values?
#'   Default `TRUE`.
#' @param ... Further method arguments; currently unused except supported
#'   print formatting.
#'
#' @return `summary` returns a point-estimator summary with a one-column
#'   Estimate table, rank and alias information, optional GIR-R2, mid-CDF
#'   summary and diagnostics. Standard extraction methods return the indicated
#'   components. `fitted` returns the linear conditional Gini index;
#'   `residuals` returns centered estimated mid-rank minus that index.
#'   `model.frame` requires the public `model` slot, and `model.matrix`
#'   requires the public `x` or `model` slot. The internal estimation snapshot
#'   is reserved for bootstrap and prediction.
#'
#' @seealso [gir()], [midcdf()], [bboot.gir()], [predict.gir()].
#' @export
summary.gir <- function(object, ...) {
  tab <- matrix(object$coefficients, ncol = 1L,
                dimnames = list(names(object$coefficients), "Estimate"))
  out <- list(
    call = object$call, coefficients = tab, aliased = object$aliased,
    rank = object$rank, df.residual = object$df.residual, nobs = object$nobs,
    GIR.R2 = if (is.null(object$GIR.R2)) NULL else object$GIR.R2$value,
    GIR.R2.details = object$GIR.R2,
    midcdf.summary = summary(object$midcdf),
    first_stage = object$first_stage, second_stage = object$second_stage
  )
  class(out) <- "summary.gir"; out
}

#' @rdname summary.gir
#' @export
print.summary.gir <- function(x, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nPoint estimates (use bboot.gir() for two-stage inference):\n")
  stats::printCoefmat(x$coefficients, digits = digits, na.print = "NA")
  if (!is.null(x$GIR.R2)) {
    cat("\nGIR-R^2: ", format(x$GIR.R2, digits = digits), "\n", sep = "")
  }
  cat("\nConditional mid-CDF fitted values:\n")
  print(x$midcdf.summary, digits = digits)
  cat("\nFirst stage:\n  Method: ", x$first_stage$method, "\n", sep = "")
  if (!identical(x$first_stage$method, "supplied")) {
    cat("  Basis: ", x$first_stage$basis, "\n", sep = "")
    cat("  Continuous variables: ",
        if (length(x$first_stage$continuous))
          paste(x$first_stage$continuous, collapse = ", ") else "none",
        "\n", sep = "")
    cat("  Discrete variables: ",
        if (length(x$first_stage$discrete))
          paste(x$first_stage$discrete, collapse = ", ") else "none",
        "\n", sep = "")
    cat("  Joint discrete cells: ", nrow(x$first_stage$diagnostics), "\n",
        sep = "")
    cat("  Basis dimension range: ",
        min(x$first_stage$diagnostics$basis_dimension), " to ",
        max(x$first_stage$diagnostics$basis_dimension), "\n", sep = "")
  }
  if (any(x$aliased)) {
    cat("\nAliased coefficients: ",
        paste(names(x$aliased)[x$aliased], collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname summary.gir
#' @export
coef.gir <- function(object, complete = TRUE, ...) {
  if (isTRUE(complete)) object$coefficients else object$coefficients[!object$aliased]
}
#' @rdname summary.gir
#' @export
fitted.gir <- function(object, ...) object$fitted.values
#' @rdname summary.gir
#' @export
residuals.gir <- function(object, ...) object$residuals
#' @rdname summary.gir
#' @export
weights.gir <- function(object, ...) object$weights
#' @rdname summary.gir
#' @export
nobs.gir <- function(object, ...) object$nobs
#' @rdname summary.gir
#' @export
formula.gir <- function(x, ...) x$formula
#' @rdname summary.gir
#' @export
terms.gir <- function(x, ...) x$terms
#' @rdname summary.gir
#' @export
model.frame.gir <- function(formula, ...) {
  if (is.null(formula$model)) stop("Refit with model=TRUE.", call. = FALSE)
  formula$model
}
#' @rdname summary.gir
#' @export
model.matrix.gir <- function(object, ...) {
  if (!is.null(object$x)) return(object$x)
  if (is.null(object$model)) stop("Neither x nor model was retained.",
                                  call. = FALSE)
  stats::model.matrix(object$terms, object$model,
                      contrasts.arg = object$contrasts)
}

#' Extract Conditional Mid-Distribution Estimates
#'
#' Extract the named first-stage conditional mid-distribution vector from a
#' GIR fit.
#'
#' @param object A fitted GIR model.
#' @param ... Further arguments; currently unused.
#'
#' @return A numeric vector in retained estimation-sample order. Standard
#'   numeric summary methods give minimum, quartiles, median, mean and
#'   maximum.
#'
#' @seealso [gir()], [bboot.gir()].
#' @export
midcdf <- function(object, ...) UseMethod("midcdf")
#' @rdname midcdf
#' @export
midcdf.gir <- function(object, ...) object$midcdf


# =============================================================================
# Bayesian-bootstrap and prediction utilities
# =============================================================================

.gir_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x != floor(x) || x > .Machine$integer.max) {
    stop(name, " must be an integer of at least ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.gir_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be strictly between 0 and 1.", call. = FALSE)
  }
  invisible(level)
}

.gir_no_dots <- function(dots, where) {
  if (length(dots)) {
    stop("Unused argument(s) in ", where, ": ",
         paste(names(dots), collapse = ", "), ".", call. = FALSE)
  }
}

.gir_state <- function(object) {
  if (!inherits(object, "gir")) stop("object must be a gir fit.", call. = FALSE)
  s <- object$estimation
  if (is.null(s) || !isTRUE(s$version %in% c(1L, 2L))) {
    stop("This fit has no current estimation snapshot. Refit it with the ",
         "updated gir() before using bootstrap inference.", call. = FALSE)
  }
  if (identical(s$version, 1L)) {
    # An earlier snapshot stored one design row per observation.
    s$group <- seq_len(object$nobs)
    s$version <- 2L
  }
  if (!is.matrix(s$X) || !is.integer(s$group) ||
      length(s$group) != object$nobs || nrow(s$X) != max(s$group) ||
      length(s$y) != object$nobs || length(s$weights) != object$nobs ||
      !identical(colnames(s$X), names(object$coefficients))) {
    stop("The retained estimation snapshot is inconsistent with the fit.",
         call. = FALSE)
  }
  s
}

.gir_boot_cache <- function(object) {
  s <- .gir_state(object)
  fs <- object$first_stage
  if (identical(fs$method, "supplied")) return(NULL)
  scaled <- .gir_scale(s$first_stage_data, fs$continuous, fs$bounds)
  dg <- fs$diagnostics
  cells <- fs$cells
  if (is.null(cells) || length(cells) != object$nobs) {
    stop("The fit does not retain its first-stage cell assignment.", call. = FALSE)
  }
  lapply(seq_len(nrow(dg)), function(g) {
    label <- as.character(dg$cell[g])
    rows <- which(as.character(cells) == label)
    spec <- list(dims = dg$marginal_dimensions[[g]],
                 knots = dg$marginal_knot_counts[[g]])
    P <- .gir_make_basis(
      scaled, fs$continuous, rows, spec, fs$basis,
      s$control$spline_order, s$control$max_basis_columns
    )
    if (ncol(P) != dg$basis_dimension[g]) {
      stop("Reconstructed basis does not match the fit in cell ", label, ".",
           call. = FALSE)
    }
    ord <- order(s$y[rows])
    tie <- rle(s$y[rows][ord])
    ends <- cumsum(tie$lengths)
    starts <- c(1L, utils::head(ends, -1L) + 1L)
    list(label = label, rows = rows, P = P,
         order_info = list(order = ord,
                           first = rep(starts, tie$lengths),
                           last = rep(ends, tie$lengths)))
  })
}

.gir_boot_rep <- function(object, cache, multiplier, first_stage_boot,
                           midcdf_boot, replicate) {
  s <- object$estimation
  a <- s$weights * multiplier
  if (any(!is.finite(a)) || any(a[s$weights > 0] <= 0)) {
    stop("Combined observation/bootstrap weights overflowed or underflowed.",
         call. = FALSE)
  }
  ctl <- s$control
  modified <- character(0)
  min_rc <- NA_real_
  if (first_stage_boot == "fixed") {
    m <- unname(object$midcdf)
  } else if (!is.null(midcdf_boot)) {
    m <- midcdf_boot(object = object, weights = a, replicate = replicate)
    if (!is.numeric(m) || !is.null(dim(m)) || length(m) != object$nobs) {
      stop("midcdf_boot must return a numeric vector in estimation-row order.",
           call. = FALSE)
    }
  } else {
    m <- rep(NA_real_, object$nobs)
    rc <- numeric(length(cache))
    for (g in seq_along(cache)) {
      cg <- cache[[g]]
      fit <- .gir_fit_cell(
        s$y[cg$rows], cg$P, a[cg$rows], cg$label,
        ctl$first_stage_singular, ctl$first_stage_tol, ctl$ridge,
        order_info = cg$order_info
      )
      if (!isTRUE(fit$ok)) {
        stop("First stage, cell ", cg$label, ": ", fit$reason,
             " [rank ", fit$rank, "/", fit$dimension,
             "; rcond=", format(fit$rcond, digits = 5), "].", call. = FALSE)
      }
      m[cg$rows] <- fit$fitted
      rc[g] <- fit$rcond
      if (fit$action != "none" || fit$ridge_lambda > 0) {
        modified <- c(modified, cg$label)
      }
    }
    min_rc <- min(rc)
  }
  clipped <- .gir_maybe_clip(m, ctl$clip, "a bootstrap first stage")
  stage2 <- .gir_second_stage(
    s$X, s$y, 2 * clipped$values - 1, a,
    singular.ok = FALSE, tolerance = ctl$second_stage_tol, group = s$group
  )
  if (any(!is.finite(stage2$coefficients))) {
    stop("The bootstrap coefficient vector is non-finite.", call. = FALSE)
  }
  list(coefficients = stage2$coefficients, midcdf = clipped$values,
       first_stage_modified = paste(modified, collapse = "; "),
       first_stage_rcond = min_rc, midcdf_outside = clipped$outside,
       second_stage_rcond = stage2$rcond)
}

.gir_transformer <- function(transform, theta) {
  p <- length(theta); cn <- names(theta)
  if (is.null(transform)) {
    return(list(fun = function(b) b, names = cn, t0 = theta, identity = TRUE))
  }
  if (is.function(transform)) {
    f <- transform
  } else if (is.numeric(transform)) {
    L <- transform
    if (is.null(dim(L))) {
      nm <- names(L)
      if (!is.null(nm)) {
        if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) ||
            any(!nm %in% cn)) {
          stop("Names in a contrast must be unique coefficient names.", call. = FALSE)
        }
        v <- stats::setNames(numeric(p), cn); v[nm] <- L
        L <- matrix(v, nrow = 1L, dimnames = list("contrast", cn))
      } else {
        if (length(L) != p) {
          stop("An unnamed contrast must have one entry per coefficient.", call. = FALSE)
        }
        L <- matrix(L, nrow = 1L, dimnames = list("contrast", cn))
      }
    } else if (is.matrix(L)) {
      if (!nrow(L) || !ncol(L)) stop("The contrast matrix is empty.", call. = FALSE)
      nm <- colnames(L)
      if (is.null(nm)) {
        if (ncol(L) != p) {
          stop("An unnamed contrast matrix must have p columns.", call. = FALSE)
        }
        colnames(L) <- cn
      } else {
        if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) ||
            any(!nm %in% cn)) {
          stop("Contrast matrix columns must name distinct coefficients.", call. = FALSE)
        }
        full <- matrix(0, nrow(L), p, dimnames = list(rownames(L), cn))
        full[, nm] <- L; L <- full
      }
      if (is.null(rownames(L))) rownames(L) <- paste0("contrast", seq_len(nrow(L)))
    } else stop("transform must be a function, a vector, or a matrix.", call. = FALSE)
    if (any(!is.finite(L))) stop("Contrasts must be finite.", call. = FALSE)
    f <- function(b) stats::setNames(as.numeric(L %*% b), rownames(L))
  } else {
    stop("transform must be NULL, a numeric contrast, or a function.", call. = FALSE)
  }
  t0 <- f(theta)
  if (!is.numeric(t0) || !is.null(dim(t0)) || !length(t0) ||
      any(!is.finite(t0))) {
    stop("transform must return a nonempty finite numeric vector.", call. = FALSE)
  }
  nm <- names(t0)
  if (is.null(nm)) nm <- paste0("statistic", seq_along(t0))
  if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("Transformed statistics must have unique nonempty names.", call. = FALSE)
  }
  original_names <- names(t0)
  checked <- function(b) {
    value <- f(b)
    if (!is.numeric(value) || !is.null(dim(value)) ||
        length(value) != length(nm) || any(!is.finite(value)) ||
        !identical(names(value), original_names)) {
      stop("transform returned non-finite values or changed its length/names ",
           "on a bootstrap draw. Undefined transformations cannot be omitted.",
           call. = FALSE)
    }
    stats::setNames(as.numeric(value), nm)
  }
  list(fun = checked, names = nm,
       t0 = stats::setNames(as.numeric(t0), nm), identity = FALSE)
}

.gir_basic_ci <- function(t0, draws, level) {
  .gir_level(level)
  draws <- as.matrix(draws)
  if (nrow(draws) < 2L) {
    stop("At least two successful draws are required for confidence intervals.",
         call. = FALSE)
  }
  alpha <- 1 - level
  ans <- matrix(NA_real_, length(t0), 2L,
                dimnames = list(names(t0), c("lwr", "upr")))
  for (j in seq_along(t0)) {
    delta <- draws[, j] - t0[j]
    q <- stats::quantile(delta, c(alpha / 2, 1 - alpha / 2),
                         names = FALSE, type = 7L)
    ans[j, ] <- c(t0[j] - q[2L], t0[j] - q[1L])
  }
  ans
}

.gir_boot_stats <- function(t0, draws, level) {
  draws <- as.matrix(draws)
  if (nrow(draws) < 2L || any(!is.finite(draws))) {
    stop("Inference requires at least two finite successful bootstrap draws.",
         call. = FALSE)
  }
  V <- stats::cov(draws)
  dimnames(V) <- list(names(t0), names(t0))
  se <- stats::setNames(sqrt(pmax(0, diag(V))), names(t0))
  bias <- stats::setNames(colMeans(draws) - t0, names(t0))
  ci <- .gir_basic_ci(t0, draws, level)
  tab <- cbind(Estimate = t0, `Boot. SE` = se, `Boot. bias` = bias,
               lwr = ci[, 1L], upr = ci[, 2L])
  list(table = tab, covariance = V, std.error = se, bias = bias, conf.int = ci)
}

.gir_null <- function(null, nm) {
  if (!is.numeric(null) || !is.null(dim(null)) || any(!is.finite(null))) {
    stop("null must be a finite numeric vector.", call. = FALSE)
  }
  if (length(null) == 1L && is.null(names(null))) {
    return(stats::setNames(rep(null, length(nm)), nm))
  }
  if (!is.null(names(null))) {
    if (anyDuplicated(names(null)) || !setequal(names(null), nm)) {
      stop("Named null values must match the transformed statistic names.", call. = FALSE)
    }
    null <- null[nm]
  }
  if (length(null) != length(nm)) {
    stop("null must be a scalar or have one value per transformed statistic.",
         call. = FALSE)
  }
  stats::setNames(as.numeric(null), nm)
}

.gir_boot_test <- function(t0, draws, level, null, alternative, joint) {
  sm <- .gir_boot_stats(t0, draws, level)
  null <- .gir_null(null, names(t0))
  delta <- sweep(draws, 2L, t0, "-")
  observed <- t0 - null
  B <- nrow(delta)
  pvalue <- vapply(seq_along(t0), function(j) {
    low <- (1 + sum(delta[, j] <= observed[j])) / (B + 1)
    high <- (1 + sum(delta[, j] >= observed[j])) / (B + 1)
    switch(alternative, less = low, greater = high,
           two.sided = min(1, 2 * min(low, high)))
  }, numeric(1))
  names(pvalue) <- names(t0)
  degenerate <- sm$std.error == 0
  if (any(degenerate)) {
    pvalue[degenerate] <- NA_real_
    warning("No bootstrap variation for: ",
            paste(names(t0)[degenerate], collapse = ", "),
            ". Their p-values are NA; a degenerate interval is not evidence ",
            "of regular, perfectly precise inference.", call. = FALSE)
  }
  joint_result <- NULL
  if (joint) {
    # Common-covariance bootstrap Wald statistic; not a nested bootstrap-t.
    # Scaling to a correlation matrix avoids a units-dependent rank decision.
    V <- sm$covariance
    se <- sm$std.error
    good <- all(se > 0)
    C <- if (good) V / outer(se, se) else NULL
    rc <- if (good) tryCatch(rcond(C), error = function(e) 0) else 0
    if (!is.finite(rc) || rc < 1e-12) {
      joint_result <- list(statistic = NA_real_, p.value = NA_real_,
                           df = length(t0),
                           reason = "Singular bootstrap covariance of restrictions")
      warning("The joint test is undefined: bootstrap covariance is singular. ",
              "Use nonredundant restrictions and more replicates.", call. = FALSE)
    } else {
      obs <- observed / se
      roots <- sweep(delta, 2L, se, "/")
      statistic <- as.numeric(crossprod(obs, solve(C, obs)))
      boot_stat <- rowSums(roots * t(solve(C, t(roots))))
      joint_result <- list(
        statistic = statistic,
        p.value = (1 + sum(boot_stat >= statistic)) / (B + 1),
        df = length(t0), bootstrap.statistics = boot_stat,
        method = "Centered, common-covariance bootstrap Wald test"
      )
    }
  }
  list(
    table = cbind(sm$table, Null = null, Difference = observed,
                   `Pr(bootstrap)` = pvalue),
    estimate = t0, null = null, difference = observed,
    std.error = sm$std.error, covariance = sm$covariance,
    conf.int = sm$conf.int, p.value = pvalue,
    alternative = alternative, joint = joint_result
  )
}

.gir_boot_compatible <- function(object, bootstrap) {
  if (!inherits(bootstrap, "gir_bootstrap")) {
    stop("bootstrap must be a result from bboot.gir().", call. = FALSE)
  }
  .gir_state(object)
  if (!identical(object$coefficients, bootstrap$t0) ||
      !identical(object$estimation, bootstrap$fit$estimation) ||
      !identical(object$first_stage, bootstrap$fit$first_stage) ||
      !identical(object$midcdf, bootstrap$fit$midcdf)) {
    stop("The bootstrap result belongs to a different or modified gir fit.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.gir_cluster_index <- function(cluster, object, cluster_order = "auto") {
  # Resolve the cluster argument of bboot.gir() to an integer index over the
  # estimation rows: a one-sided formula is evaluated in the retained model
  # frame or first-stage columns and is already in estimation order. Vector
  # order is explicit, or inferred only when the row mapping is unambiguous.
  n <- object$nobs
  values <- cluster
  is_formula <- inherits(cluster, "formula")
  if (is_formula) {
    if (cluster_order != "auto") {
      stop("cluster_order applies only to cluster vectors; formulas already ",
           "use estimation-row order. Leave cluster_order = 'auto'.", call. = FALSE)
    }
    if (length(cluster) != 2L || !length(all.vars(cluster))) {
      stop("cluster must be a one-sided formula naming variables, such as ~ id.",
           call. = FALSE)
    }
    sources <- list(object$model, object$estimation$first_stage_data)
    columns <- lapply(all.vars(cluster), function(v) {
      for (src in sources) {
        if (!is.null(src) && v %in% names(src)) return(src[[v]])
      }
      stop("cluster variable '", v, "' is not a retained model or first-stage ",
           "column. Refit with model = TRUE, or pass the cluster values as a ",
           "vector.", call. = FALSE)
    })
    if (any(vapply(columns, anyNA, logical(1L)))) {
      stop("cluster contains missing values.", call. = FALSE)
    }
    values <- if (length(columns) == 1L) columns[[1L]] else {
      # Group exact value combinations, not potentially ambiguous display labels.
      # Integer codes keep punctuation and numeric formatting out of the keys.
      codes <- lapply(columns, function(x) match(x, unique(x)))
      do.call(paste, c(codes, sep = "|"))
    }
  }
  if (!is.null(dim(values)) || is.list(values)) {
    stop("cluster must be a vector or a one-sided formula.", call. = FALSE)
  }
  if (!is_formula) {
    if (cluster_order == "auto") {
      original_length <- length(values) == object$original.nobs
      estimation_length <- length(values) == n
      if (!original_length && !estimation_length) {
        stop("cluster must have one value per estimation row (", n, ") or per ",
             "row of the original data (", object$original.nobs, ").", call. = FALSE)
      }
      if (original_length && estimation_length &&
          !identical(object$row.index, seq_len(n))) {
        stop("cluster vector order is ambiguous: the original data and ",
             "estimation sample have the same length but a different row mapping. ",
             "Set cluster_order = 'original' or 'estimation'.", call. = FALSE)
      }
      cluster_order <- if (original_length && !estimation_length) "original" else
        "estimation"
    }
    expected <- if (cluster_order == "original") object$original.nobs else n
    if (length(values) != expected) {
      stop("cluster_order = '", cluster_order, "' requires a cluster vector of ",
           "length ", expected, ".", call. = FALSE)
    }
    if (cluster_order == "original") values <- values[object$row.index]
  } else if (length(values) != n) {
    stop("cluster formula columns must have one value per estimation row (",
         n, ").", call. = FALSE)
  }
  if (anyNA(values)) stop("cluster contains missing values.", call. = FALSE)
  index <- match(values, unique(values))
  if (max(index) < 2L) {
    stop("cluster must define at least two clusters.", call. = FALSE)
  }
  list(index = index, n_clusters = max(index))
}

# =============================================================================
# Public Bayesian-bootstrap function
# =============================================================================

#' Bayesian Bootstrap for Gini Index Regression
#'
#' Refit both GIR stages using common normalized exponential multipliers, and
#' obtain centered/basic intervals, bootstrap standard errors, covariance, and
#' tests of smooth coefficient transformations.
#'
#' @aliases gir_bootstrap
#' @param object A current `gir` fit, preferably. Alternatively a two-sided
#'   formula, with fully named `gir()` arguments in `...`.
#' @param R Number of replicates, default 999. Omitted with a supplied
#'   multiplier matrix, use its number of columns. `R=1` is allowed only with
#'   `inference=FALSE`.
#' @param level Confidence level strictly between zero and one. Basic
#'   intervals remain two-sided regardless of the test alternative.
#' @param transform `NULL` for raw coefficient tests; a numeric contrast
#'   vector or matrix; or a deterministic function of named coefficients
#'   returning a fixed-length finite numeric vector. Named contrast columns
#'   can select a subset, with others zero. Smooth, nondegenerate
#'   transformations are intended.
#' @param null Scalar or vector null values for transformed statistics. Named
#'   vectors must match all transformed statistic names.
#' @param alternative `"two.sided"` for equal-tail marginal tests; `"less"`
#'   or `"greater"` for directional alternatives. Marginal p-values are
#'   unadjusted for multiplicity.
#' @param joint Also compute a bootstrap-calibrated common-covariance joint
#'   Wald test? Requires two-sided, nonredundant regular restrictions.
#'   Singular restriction covariance gives an undefined joint test with
#'   warning.
#' @param seed Optional nonnegative integer local seed. Restore the caller RNG
#'   state on exit. `NULL` advances the existing stream. Cannot be combined
#'   with supplied multipliers.
#' @param cluster Optional clustering of the multipliers: a one-sided formula
#'   naming retained model or first-stage columns, whose combination defines
#'   the clusters, or a vector with one value per estimation row or per row
#'   of the original data (see `cluster_order`). One Exp(1) multiplier is drawn
#'   per cluster and shared by all its observations in both stages, then
#'   normalized as usual.
#'   Use it when observations are dependent within groups, for example firms
#'   within industries under an industry-level treatment. Cannot be combined
#'   with `bootstrap_weights`.
#' @param cluster_order Row order of a `cluster` vector. `"original"` requires
#'   one value per original-data row, in the order at fitting time, and aligns
#'   it using `object$row.index`. `"estimation"` requires one value per final
#'   estimation row and uses the vector unchanged. The default `"auto"` infers
#'   the order from length. If both lengths match, automatic ordering is
#'   allowed only when `object$row.index` is the unchanged row sequence;
#'   otherwise specify `"original"` or `"estimation"` explicitly. Missing
#'   cluster values are checked after alignment. Leave `"auto"` for formulas,
#'   whose retained columns are already aligned, or when `cluster = NULL`.
#' @param bootstrap_weights Positive finite nobs-by-R multipliers, or an
#'   nobs-vector for a single draw, in final estimation-row order. They
#'   multiply original case weights once.
#' @param normalize_weights Divide each supplied multiplier column by its
#'   mean? For the internal estimator, including relative ridge handling,
#'   common rescaling preserves the draws up to floating-point error within
#'   usable numerical ranges. `FALSE` uses and stores the supplied values
#'   unchanged. External first-stage callbacks must respect weight scaling
#'   themselves. Generated weights always use the paper normalization.
#' @param first_stage_boot `"refit"` re-estimates the first stage for each
#'   draw. `"fixed"` explicitly freezes it for known/oracle nuisances; it is
#'   not full two-stage inference for an estimated nuisance.
#' @param midcdf_boot Optional callback `function(object, weights, replicate)`
#'   for a point fit with a supplied first stage. Return a finite nobs-vector
#'   of refitted mid-CDFs. The callback weights already combine original
#'   weights and multipliers. Required to refit an external first stage;
#'   disallowed for internal series fits or fixed mode.
#' @param inference Compute summaries and tests? `FALSE` returns raw draws and
#'   diagnostics, enabling a one-draw warp-speed experiment. At least two
#'   successes are needed for inference.
#' @param keep_weights,keep_midcdf Store nobs-by-R multipliers as used (raw
#'   when supplied with `normalize_weights = FALSE`) and first-stage fitted
#'   draw values, respectively. Defaults `FALSE` to limit memory.
#' @param failure `"error"` stops on any failed draw. `"omit"` records a
#'   whole `NA` coefficient row and warns; no replacement draw is generated.
#'   Conditioning on successful draws may distort inference.
#' @param ... Fully named `gir()` arguments only when `object` is a formula.
#'   A fitted object already fixes its point-estimation specification, so
#'   extra arguments are an error. Bootstrap controls after the dots require
#'   exact names, preventing `first_stage` and `midcdf` from partially
#'   matching bootstrap controls.
#'
#' @details
#' For each replicate generate iid Exp(1) values and divide by their
#' full-sample mean. Multiply the same vector into the first-stage
#' Gram/scores and second-stage GIR moments. Cells, outcome ties, original
#' basis allocation and support bounds are held fixed; first-stage numerical
#' handling and clipping rules are inherited. Bases and tie orders are cached
#' within the call, and each draw accumulates the second-stage sums over the
#' distinct regressor rows, as `gir()` does. Original estimation data are
#' retained in the fit and are not re-evaluated from a possibly changed data
#' frame.
#'
#' With `cluster`, the Exp(1) values are drawn per cluster and repeated
#' within it, a Bayesian cluster bootstrap. When every first-stage cell lies
#' inside a single cluster, the cell fits are invariant to a common
#' multiplier, so the draws reflect between-cluster variation only and the
#' estimated first stage enters through the residual variation across
#' clusters, as in a clustered regression of cell Ginis; the returned
#' `cluster$first_stage_invariant` reports whether this is the case. The
#' paper's theorem covers independent observations; clustered inference
#' rests on the usual many-clusters argument.
#'
#' A formula such as `cluster = ~ industry + region` defines one cluster per
#' distinct industry-region combination: this is one-way clustering at their
#' intersections, not two-way clustering by industry and region separately.
#' Group membership uses the values in each column, so punctuation in labels
#' cannot merge distinct combinations. Numeric identifiers are distinguished
#' by their stored values, not their printed formatting. Clusters are numbered
#' in order of first occurrence in the estimation rows; changing labels
#' without changing group membership preserves seeded multipliers.
#'
#' For errors \eqn{D_{bj}=\widehat\theta_j^{*(b)}-\widehat\theta_j},
#' Algorithm 2 gives
#' \deqn{[\widehat\theta_j-q_{1-\alpha/2}(D_j),\quad
#'         \widehat\theta_j-q_{\alpha/2}(D_j)].}
#' Quantiles use R's default type 7. These are basic
#' confidence intervals, not percentile credible intervals. SEs are column
#' sample standard deviations of successful draws with B-1 denominator, not
#' divided by sqrt(B) or sqrt(n). Bootstrap bias is the draw mean minus the
#' point estimate and does not replace that point estimate.
#'
#' For a transformed null, compare observed \eqn{h(t_0)-null} with centered
#' draws \eqn{h(t^*)-h(t_0)}. Tail estimates use (1+exceedances)/(B+1);
#' two-sided p-values are twice the smaller tail, capped at one. This
#' finite-simulation safeguard does not confer exact finite-sample validity.
#' The joint statistic uses a common sample covariance of transformed draws
#' for observed and centered-bootstrap Wald quadratic forms, with an
#' upper-tail bootstrap p-value. There is no OLS F calibration, nested
#' bootstrap or replicate-specific studentization. Zero bootstrap variation
#' produces an undefined marginal p-value with warning.
#'
#' Smooth transformation validity requires nondegenerate derivatives. Joint
#' restriction derivatives must have full row rank at the null; numerical
#' rank checks are not a proof of that condition. Weak bootstrap convergence
#' alone does not imply bootstrap variance moment consistency; additional
#' moment/uniform-integrability conditions are needed for interpreting
#' covariance estimates consistently. The paper's inference theory covers
#' continuous or finite-discrete response regimes, not all mixed laws.
#' Nonuniform case weights, clustered multipliers, survey designs,
#' regularization and external first stages require their own justification.
#'
#' A supplied mid-CDF vector requires a refitting callback or an explicit
#' choice of fixed mode. Bootstrap inference requires a full-rank second stage.
#' Omitting failed draws produces a warning. Use `R=1`, `inference=FALSE` as
#' a draw generator, for
#' example in a warp-speed Monte Carlo study that pools centered draws across
#' simulated datasets.
#'
#' @return An object of class `gir_bootstrap`, not class `boot`. Components
#'   `t0` and `t` contain raw coefficients and the R-by-p coefficient draw
#'   matrix, even when a transformation is tested. It also contains the point
#'   `fit`, requested/successful counts, success indicators, per-draw
#'   diagnostics, RNG metadata, the cluster structure (`cluster`: number of
#'   clusters, index and first-stage invariance) when clustering was
#'   requested, and optional multiplier/mid-CDF matrices. With
#'   inference, `std.error`, `covariance`, `conf.int`, `bias` and
#'   `coefficient.table` concern raw coefficients; `test` contains transformed
#'   inference and an optional joint test. These summary fields are `NULL`
#'   for draw-only output.
#'
#' @references
#' Bhardwaj, V. and Diaz-Villarejo, D. (2026). Gini Index Regression. Working
#' paper, Universidad Carlos III de Madrid.
#' \url{https://ddiazvil.com/assets/papers/JMP_GIR_DiazVillarejo.pdf}.
#'
#' Cheng, G. (2014, revised preprint). Moment Consistency of the Exchangeably
#' Weighted Bootstrap for Semiparametric M-Estimation.
#' \url{https://arxiv.org/abs/1109.4204}.
#'
#' @seealso [gir()], [predict.gir()], [summary.gir_bootstrap()].
#'
#' @examples
#' d <- data.frame(y = rep(c(1, 2, 2, 4, 5), 8),
#'                 g = factor(rep(c("a", "b"), each = 20)))
#' f <- gir(y ~ g, data = d, discrete = "g", singular.ok = FALSE)
#' bb <- bboot.gir(f, R = 499L, seed = 42L)
#' summary(bb)
#' confint(bb)
#' vcov(bb)
#' summary(bb, transform = c(gb = 1), null = 0)
#'
#' # Inspect cluster alignment with a reordered estimation sample.
#' cluster_dat <- data.frame(y = c(1, 4, 2, 7, 3, 9),
#'                           id = c("A", "A", "B", "B", "C", "C"))
#' cluster_fit <- gir(y ~ 1, data = cluster_dat, subset = c(1, 3, 5, 2, 4, 6),
#'                    first_stage = ~ id, discrete = "id", minimum_cell_size = 1L)
#' cluster_boot <- bboot.gir(cluster_fit, cluster = cluster_dat$id,
#'                           cluster_order = "original", R = 3L, seed = 17L,
#'                           inference = FALSE)
#' cluster_boot$cluster$index  # 1 2 3 1 2 3
#'
#' @export bboot.gir
bboot.gir <- function(
    object,
    R = 999L,
    level = 0.95,
    ...,
    transform = NULL,
    null = 0,
    alternative = c("two.sided", "less", "greater"),
    joint = FALSE,
    seed = NULL,
    cluster = NULL,
    cluster_order = c("auto", "original", "estimation"),
    bootstrap_weights = NULL,
    normalize_weights = TRUE,
    first_stage_boot = c("refit", "fixed"),
    midcdf_boot = NULL,
    inference = TRUE,
    keep_weights = FALSE,
    keep_midcdf = FALSE,
    failure = c("error", "omit")) {
  call <- match.call()
  unevaluated <- match.call(expand.dots = FALSE)$...
  caller <- parent.frame()
  if (inherits(object, "formula")) {
    args <- as.list(unevaluated)
    nm <- names(args)
    if (length(args) && (is.null(nm) || any(!nzchar(nm)) ||
                        anyDuplicated(nm) || any(!nm %in% names(formals(gir))))) {
      stop("All extra arguments must be fully named gir() arguments.", call. = FALSE)
    }
    if ("formula" %in% nm) stop("The formula is already supplied as object.", call. = FALSE)
    supplied_formula <- object
    env <- new.env(parent = caller)
    env$.GIR_formula <- supplied_formula
    env$gir <- gir
    fc <- as.call(c(list(as.name("gir"), formula = as.name(".GIR_formula")), args))
    object <- eval(fc, envir = env)
    object$call$formula <- supplied_formula
  } else if (length(unevaluated)) {
    stop("A fitted gir object already specifies all point-estimation options. ",
         "Refit gir() to change them; extra arguments are not used.", call. = FALSE)
  }
  s <- .gir_state(object)
  n <- object$nobs
  if (anyNA(object$coefficients) || any(object$aliased) || object$rank < ncol(s$X)) {
    stop("Bootstrap coefficient inference requires an identified, full-rank ",
         "second stage. Refit an identified specification.", call. = FALSE)
  }
  .gir_level(level)
  alternative <- match.arg(alternative)
  cluster_order <- match.arg(cluster_order)
  if (is.null(cluster) && cluster_order != "auto") {
    stop("cluster_order requires a cluster vector; cluster is NULL.", call. = FALSE)
  }
  first_stage_boot <- match.arg(first_stage_boot)
  failure <- match.arg(failure)
  for (a in c("joint", "normalize_weights", "inference", "keep_weights", "keep_midcdf")) {
    .gir_flag(get(a), a)
  }
  if (joint && alternative != "two.sided") {
    stop("joint=TRUE requires alternative='two.sided'.", call. = FALSE)
  }
  if (!is.null(midcdf_boot) && !is.function(midcdf_boot)) {
    stop("midcdf_boot must be NULL or a function.", call. = FALSE)
  }
  supplied_fs <- identical(object$first_stage$method, "supplied")
  if (first_stage_boot == "fixed" && !is.null(midcdf_boot)) {
    stop("midcdf_boot cannot be combined with first_stage_boot='fixed'.", call. = FALSE)
  }
  if (!supplied_fs && !is.null(midcdf_boot)) {
    stop("midcdf_boot is only for a point fit with a supplied first stage.", call. = FALSE)
  }
  if (first_stage_boot == "refit" && supplied_fs && is.null(midcdf_boot)) {
    stop("A supplied mid-CDF cannot be re-estimated from its vector alone. ",
         "Supply midcdf_boot, or explicitly use first_stage_boot='fixed' ",
         "for a known/oracle first stage.", call. = FALSE)
  }
  if (inference && first_stage_boot == "fixed") {
    warning("The first stage is treated as fixed. These are not full two-stage ",
            "GIR intervals when that first stage was estimated from the data.", call. = FALSE)
  }
  if (inference && (any(s$weights == 0) ||
                    length(unique(s$weights[s$weights > 0])) > 1L)) {
    warning("Observation weights are multiplied by iid bootstrap multipliers. ",
            "This is not design-based survey or cluster inference; the supplied ",
            "paper's iid, unweighted theorem does not by itself cover it.", call. = FALSE)
  }
  if (inference && any(s$y < 0)) {
    warning("Negative outcomes lie outside the paper's Gini inference assumptions.",
            call. = FALSE)
  }

  cl <- NULL
  if (!is.null(cluster)) {
    if (!is.null(bootstrap_weights)) {
      stop("Supply either cluster or bootstrap_weights, not both.", call. = FALSE)
    }
    cl <- .gir_cluster_index(cluster, object, cluster_order)
  }
  W <- NULL
  if (!is.null(bootstrap_weights)) {
    if (!is.null(seed)) {
      stop("Supply either seed or bootstrap_weights, not both.", call. = FALSE)
    }
    if (!is.numeric(bootstrap_weights) || length(dim(bootstrap_weights)) > 2L) {
      stop("bootstrap_weights must be a numeric vector or matrix.", call. = FALSE)
    }
    W <- if (is.null(dim(bootstrap_weights))) {
      matrix(bootstrap_weights, ncol = 1L)
    } else as.matrix(bootstrap_weights)
    if (nrow(W) != n || ncol(W) < 1L || any(!is.finite(W)) || any(W <= 0)) {
      stop("bootstrap_weights must have nobs(object) rows, positive finite ",
           "entries, and one column per replicate.", call. = FALSE)
    }
    if (missing(R)) R <- ncol(W)
    R <- .gir_integer(R, "R")
    if (R != ncol(W)) stop("R does not match the multiplier matrix.", call. = FALSE)
    # Internal refits, including relative ridge, preserve common weight scaling
    # within numerical precision. External callbacks receive the weights as used.
    if (normalize_weights) {
      for (b in seq_len(ncol(W))) {
        average <- mean(W[, b])
        if (!is.finite(average) || average <= 0) stop("Invalid multiplier mean.", call. = FALSE)
        W[, b] <- W[, b] / average
      }
    }
    if (any(!is.finite(W)) || any(W <= 0)) {
      stop("Multiplier normalization underflowed or overflowed.", call. = FALSE)
    }
  }
  R <- .gir_integer(R, "R")
  if (inference && R < 2L) {
    stop("R=1 yields a draw, not an interval or standard error. Use ",
         "inference=FALSE for a warp-speed Monte Carlo replication.", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .gir_integer(seed, "seed", minimum = 0L)
    existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (existed) get(".Random.seed", envir = .GlobalEnv) else NULL
    old_kind <- RNGkind()
    on.exit({
      do.call(RNGkind, as.list(old_kind))
      if (existed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  rng_kind <- RNGkind()
  rng_start <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else NULL
  # Validate transformation and null before undertaking the expensive refits.
  tr <- .gir_transformer(transform, object$coefficients)
  if (joint && length(tr$t0) > length(object$coefficients)) {
    stop("A regular joint test cannot have more independent restrictions than ",
         "coefficients. Reduce the transformation.", call. = FALSE)
  }
  .gir_null(null, tr$names)
  cache <- if (first_stage_boot == "refit" && !supplied_fs) .gir_boot_cache(object) else NULL
  cluster_info <- NULL
  if (!is.null(cl)) {
    invariant <- NA
    if (!is.null(cache)) {
      # When every first-stage cell lies inside one cluster, the cell fits
      # are invariant to the common multiplier of that cluster and only the
      # second stage varies across draws.
      cells <- as.integer(object$first_stage$cells)
      key <- as.numeric(cells) * (cl$n_clusters + 1) + cl$index
      invariant <- length(unique(key)) == length(unique(cells))
    }
    cluster_info <- list(n_clusters = cl$n_clusters, index = cl$index,
                         first_stage_invariant = invariant)
  }
  draws <- matrix(NA_real_, R, length(object$coefficients),
                  dimnames = list(NULL, names(object$coefficients)))
  saved_W <- if (keep_weights) matrix(NA_real_, n, R) else NULL
  saved_M <- if (keep_midcdf) matrix(NA_real_, n, R) else NULL
  dg <- data.frame(
    replicate = seq_len(R), success = rep(FALSE, R),
    first_stage_rcond = rep(NA_real_, R),
    first_stage_modified_cells = rep("", R),
    midcdf_outside = rep(NA_integer_, R), second_stage_rcond = rep(NA_real_, R),
    message = rep("", R), stringsAsFactors = FALSE
  )
  warning_log <- vector("list", R)
  for (b in seq_len(R)) {
    w <- if (is.null(W)) {
      # One Exp(1) value per observation, or per cluster shared by all its
      # observations, normalized by the observation-level mean.
      e <- stats::rexp(if (is.null(cl)) n else cl$n_clusters, rate = 1)
      if (!is.null(cl)) e <- e[cl$index]
      e / mean(e)
    } else W[, b]
    if (any(!is.finite(w)) || any(w <= 0)) {
      stop("A generated multiplier is nonpositive or non-finite; no redraw ",
           "has been made.", call. = FALSE)
    }
    if (keep_weights) saved_W[, b] <- w
    messages <- character(0)
    value <- tryCatch(withCallingHandlers(
      .gir_boot_rep(object, cache, w, first_stage_boot, midcdf_boot, b),
      warning = function(warn) {
        messages <<- c(messages, conditionMessage(warn))
        invokeRestart("muffleWarning")
      }
    ), error = function(e) e)
    warning_log[[b]] <- unique(messages)
    if (inherits(value, "error")) {
      dg$message[b] <- conditionMessage(value)
      if (failure == "error") {
        stop("Bayesian bootstrap replicate ", b, " failed: ", dg$message[b],
             " No failed draw was redrawn or silently discarded.", call. = FALSE)
      }
      next
    }
    draws[b, ] <- value$coefficients
    if (keep_midcdf) saved_M[, b] <- value$midcdf
    dg$success[b] <- TRUE
    dg$first_stage_rcond[b] <- value$first_stage_rcond
    dg$first_stage_modified_cells[b] <- value$first_stage_modified
    dg$midcdf_outside[b] <- value$midcdf_outside
    dg$second_stage_rcond[b] <- value$second_stage_rcond
  }
  dg$warnings <- I(warning_log)
  B <- sum(dg$success)
  if (!B || (inference && B < 2L)) {
    stop("Too few successful replicates for the requested output.", call. = FALSE)
  }
  if (B < R) {
    warning(R - B, " failed replicate(s) were omitted, not redrawn. Inference ",
            "is conditional on successful fits and can be distorted; inspect ",
            "diagnostics before using it.", call. = FALSE)
  }
  nout <- sum(dg$midcdf_outside > 0L, na.rm = TRUE)
  nmod <- sum(nzchar(dg$first_stage_modified_cells))
  nwarn <- sum(lengths(warning_log) > 0L)
  if (nout || nmod || nwarn) {
    warning("Bootstrap diagnostics: ", nout, " draw(s) with mid-CDF values ",
            "outside [0,1], ", nmod, " draw(s) with first-stage numerical ",
            "adjustments, and ", nwarn, " draw(s) with captured warnings. ",
            "Inspect the returned diagnostics; clipping/singularity policies ",
            "are inherited from gir().", call. = FALSE)
  }
  if (keep_weights) dimnames(saved_W) <- list(names(object$midcdf), paste0("b", seq_len(R)))
  if (keep_midcdf) dimnames(saved_M) <- list(names(object$midcdf), paste0("b", seq_len(R)))
  out <- list(
    call = call, fit = object, t0 = object$coefficients, t = draws,
    R = R, R.success = B, success = dg$success,
    level = level, quantile.type = 7L, transform = transform, null = null,
    alternative = alternative, joint = joint, inference = inference,
    first_stage_boot = if (!is.null(midcdf_boot)) "external refit" else first_stage_boot,
    cluster = cluster_info,
    seed = seed, rng.kind = rng_kind, rng.start = rng_start,
    bootstrap_weights = saved_W, midcdf = saved_M, diagnostics = dg,
    std.error = NULL, covariance = NULL, conf.int = NULL, bias = NULL,
    coefficient.table = NULL, test = NULL
  )
  class(out) <- "gir_bootstrap"
  if (inference) {
    sm <- summary(out)
    out$std.error <- sm$std.error
    out$covariance <- sm$covariance
    out$conf.int <- sm$conf.int
    out$bias <- sm$bias
    out$coefficient.table <- sm$coefficients
    out$test <- sm$test
  }
  out
}

#' Summaries and Extraction for GIR Bootstrap Results
#'
#' Reanalyze stored raw bootstrap draws at a new level or for smooth
#' coefficient transformations, and extract raw coefficients, covariance or
#' confidence intervals.
#'
#' @aliases summary.gir_bootstrap.object
#' @param object,x A `gir_bootstrap` object, or its summary for the summary
#'   print method.
#' @param level,transform,null,alternative,joint As in `bboot.gir()`. Summary
#'   defaults inherit the original bootstrap call. A new summary reuses the
#'   stored successful raw draws and does not modify the bootstrap object.
#' @param digits Print precision.
#' @param parm Optional raw coefficient names or positive integer indices for
#'   `confint`.
#' @param ... Print methods forward formatting arguments. Summary, covariance
#'   and interval methods reject unused arguments; coefficient extraction
#'   ignores extras.
#'
#' @details
#' At least two successful draws are needed for inferential summaries. A low
#' expected number of observations in an interval tail warns. Marginal tests
#' are unadjusted and confidence intervals remain two-sided. Standard errors,
#' moment-consistency qualifications and the common-covariance Wald test are
#' described under [bboot.gir()].
#'
#' @return `summary` returns class `summary.gir_bootstrap` with raw
#'   coefficient table, SEs, covariance, intervals, bias, transformed test
#'   results and diagnostics. `print` invisibly returns its argument. `coef`
#'   returns the original raw coefficient estimate, not bootstrap means.
#'   `vcov` returns raw coefficient covariance; `confint` returns basic
#'   intervals for raw coefficients. Obtain transformed intervals from
#'   `summary(...)$test$conf.int`.
#'
#' @seealso [bboot.gir()], [predict.gir()].
#' @export
summary.gir_bootstrap <- function(
    object, level = object$level, transform = object$transform,
    null = object$null, alternative = object$alternative,
    joint = object$joint, ...) {
  .gir_no_dots(list(...), "summary.gir_bootstrap()")
  .gir_level(level); .gir_flag(joint, "joint")
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  if (joint && alternative != "two.sided") stop("A joint test is two-sided.", call. = FALSE)
  draws <- object$t[object$success, , drop = FALSE]
  B <- nrow(draws)
  if (B < 2L) stop("A single draw has no bootstrap standard error or interval.", call. = FALSE)
  if (B * (1 - level) / 2 < 10) {
    warning("Fewer than ten expected draws per confidence-interval tail. ",
            "Use more replicates for stable tail inference.", call. = FALSE)
  }
  raw <- .gir_boot_stats(object$t0, draws, level)
  tr <- .gir_transformer(transform, object$t0)
  if (joint && length(tr$t0) > length(object$t0)) {
    stop("A regular joint test cannot have more independent restrictions than ",
         "coefficients. Reduce the transformation.", call. = FALSE)
  }
  transformed <- if (tr$identity) draws else t(vapply(
    seq_len(B), function(b) tr$fun(stats::setNames(draws[b, ], names(object$t0))),
    numeric(length(tr$t0))
  ))
  # vapply() returns a vector for a scalar statistic; force B-by-q orientation.
  if (!tr$identity && length(tr$t0) == 1L) {
    transformed <- matrix(transformed, nrow = B, ncol = 1L)
  }
  colnames(transformed) <- names(tr$t0)
  test <- .gir_boot_test(tr$t0, transformed, level, null, alternative, joint)
  out <- list(
    call = object$call, fit.call = object$fit$call, nobs = object$fit$nobs,
    R = object$R, R.success = B, level = level,
    coefficients = raw$table, std.error = raw$std.error,
    covariance = raw$covariance, conf.int = raw$conf.int, bias = raw$bias,
    test = test, transformed = !tr$identity,
    first_stage_boot = object$first_stage_boot, cluster = object$cluster,
    diagnostics = object$diagnostics
  )
  class(out) <- "summary.gir_bootstrap"
  out
}

#' @rdname summary.gir_bootstrap
#' @export
print.gir_bootstrap <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!isTRUE(x$inference)) {
    cat("\nBayesian bootstrap draws for GIR\n")
    cat("Successful replicates:", x$R.success, "of", x$R, "\n")
    cat("Point estimate:\n"); print(x$t0, digits = digits)
    cat("\nDraws are in $t. No inference was requested.\n")
  } else {
    # Stored inference is printed without another warning-producing analysis.
    view <- list(fit.call = x$fit$call, R = x$R, R.success = x$R.success,
                 level = x$level, coefficients = x$coefficient.table,
                 test = x$test, transformed = !is.null(x$transform),
                 first_stage_boot = x$first_stage_boot, cluster = x$cluster)
    class(view) <- "summary.gir_bootstrap"
    print(view, digits = digits, ...)
  }
  invisible(x)
}

#' @rdname summary.gir_bootstrap
#' @export
print.summary.gir_bootstrap <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nBayesian bootstrap for Gini Index Regression\n\nCall:\n")
  print(x$fit.call)
  cat("\nSuccessful replicates:", x$R.success, "of", x$R,
      "\nFirst stage:", x$first_stage_boot)
  if (!is.null(x$cluster)) {
    cat("\nMultipliers:", x$cluster$n_clusters, "clusters",
        if (isTRUE(x$cluster$first_stage_invariant)) {
          "(first-stage cells nested in clusters, so cell fits do not vary)"
        } else "")
  }
  cat("\nConfidence level:", x$level, "(basic, type-7 quantiles)\n")
  if (!x$transformed) {
    cat("\nCoefficient inference and marginal null tests:\n")
    print(x$test$table, digits = digits, ...)
  } else {
    cat("\nCoefficient estimates and intervals:\n")
    print(x$coefficients, digits = digits, ...)
    cat("\nTransformed statistics and marginal null tests:\n")
    print(x$test$table, digits = digits, ...)
  }
  cat("\nMarginal alternative:", x$test$alternative,
      "\nMarginal p-values are not multiplicity adjusted.\n")
  if (!is.null(x$test$joint)) {
    cat("\nJoint bootstrap Wald statistic:",
        format(x$test$joint$statistic, digits = digits),
        "\nJoint bootstrap p-value:",
        format(x$test$joint$p.value, digits = digits), "\n")
    if (!is.null(x$test$joint$reason)) cat(x$test$joint$reason, "\n")
  }
  invisible(x)
}

#' @rdname summary.gir_bootstrap
#' @export
coef.gir_bootstrap <- function(object, ...) object$t0

#' @rdname summary.gir_bootstrap
#' @export
vcov.gir_bootstrap <- function(object, ...) {
  .gir_no_dots(list(...), "vcov.gir_bootstrap()")
  if (!is.null(object$covariance)) return(object$covariance)
  draws <- object$t[object$success, , drop = FALSE]
  if (nrow(draws) < 2L) stop("At least two draws are required for covariance.", call. = FALSE)
  stats::cov(draws)
}

#' @rdname summary.gir_bootstrap
#' @export
confint.gir_bootstrap <- function(object, parm, level = object$level, ...) {
  .gir_no_dots(list(...), "confint.gir_bootstrap()")
  ci <- .gir_basic_ci(object$t0, object$t[object$success, , drop = FALSE], level)
  if (!missing(parm)) {
    if (is.character(parm)) {
      if (anyNA(parm) || any(!parm %in% rownames(ci))) stop("Unknown coefficient in parm.", call. = FALSE)
    } else if (!is.numeric(parm) || anyNA(parm) || any(parm != floor(parm)) ||
               any(parm < 1L | parm > nrow(ci))) {
      stop("parm must select valid coefficient names or positive indices.", call. = FALSE)
    }
    ci <- ci[parm, , drop = FALSE]
  }
  colnames(ci) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), trim = TRUE), "%")
  ci
}

# =============================================================================
# Prediction of the conditional Gini index (not prediction of Y)
# =============================================================================

.gir_new_matrix <- function(object, newdata, na.action) {
  s <- .gir_state(object)
  # For the estimation sample, X holds the distinct regressor rows and
  # `group` maps the observations to them.
  if (is.null(newdata)) return(list(X = s$X, action = NULL, group = s$group))
  data <- as.data.frame(newdata)
  absent <- setdiff(s$predictor_variables, names(data))
  if (length(absent)) {
    stop("newdata is missing predictor column(s): ", paste(absent, collapse = ", "),
         ".", call. = FALSE)
  }
  tt <- stats::delete.response(object$terms)
  lev <- s$predictor_levels
  # A factor response is not permitted by gir(), but remove its name defensively.
  lev[[object$outcome.name]] <- NULL
  mf <- stats::model.frame(tt, data = data, na.action = stats::na.pass,
                           xlev = lev, drop.unused.levels = FALSE)
  fun <- if (is.null(na.action)) stats::na.pass else match.fun(na.action)
  acted <- fun(mf)
  if (!is.data.frame(acted)) stop("na.action must return a data frame.", call. = FALSE)
  action <- attr(acted, "na.action")
  X <- stats::model.matrix(tt, acted, contrasts.arg = object$contrasts)
  cn <- names(object$coefficients)
  if (!setequal(colnames(X), cn) || anyDuplicated(colnames(X))) {
    stop("newdata produces a design incompatible with the original fit. ",
         "Check variable classes, contrasts and factor levels.", call. = FALSE)
  }
  X <- X[, cn, drop = FALSE]
  if (any(is.infinite(X))) stop("newdata produces infinite design values.", call. = FALSE)
  list(X = X, action = action)
}

.gir_estimable_rows <- function(object, X) {
  if (!any(object$aliased)) return(rep(TRUE, nrow(X)))
  s <- .gir_state(object)
  wy <- s$weights * s$y
  design <- if (all(wy >= 0)) {
    # The distinct rows scaled by the square root of their summed weights
    # span the same space as the observation-level weighted design.
    s$X * sqrt(rowsum(wy, s$group)[, 1L])
  } else object$second_stage$moment_matrix
  V <- svd(design, nu = 0L, nv = ncol(design))$v
  V <- V[, seq_len(object$rank), drop = FALSE]
  difference <- X - (X %*% V) %*% t(V)
  sqrt(rowSums(difference^2)) <=
    object$second_stage$tolerance * pmax(1, sqrt(rowSums(X^2)))
}

#' Predict the Linear Conditional Gini Index
#'
#' Evaluate the GIR linear index on the estimation sample or new covariates,
#' with optional bootstrap standard errors and pointwise basic confidence
#' intervals.
#'
#' @param object A current `gir` fit.
#' @param newdata Optional data frame of raw second-stage predictors. `NULL`
#'   uses retained estimation rows in the order of `fitted(object)`. No
#'   response or first-stage-only variables are required.
#' @param se.fit Return a list with prediction(s) and bootstrap SEs? Default
#'   `FALSE`.
#' @param interval `"none"` or `"confidence"`. The latter gives pointwise
#'   basic bootstrap confidence intervals for the linear Gini index, not
#'   individual-outcome prediction intervals.
#' @param level Confidence level strictly between zero and one.
#' @param bootstrap Optional matching `gir_bootstrap` object. Its raw
#'   coefficient draws are reused without refitting or new randomness. A
#'   different estimator or sample is rejected.
#' @param R,seed Replicate count and local seed for an automatic bootstrap
#'   only. Do not supply them together with an existing bootstrap object.
#' @param na.action For `newdata` only: default `na.pass` preserves missing
#'   rows as `NA`; omit removes them; exclude reinserts them; fail stops.
#'   Empty `newdata`, or omission of every prediction row, returns a
#'   length-zero point vector or a zero-row, three-column confidence-interval
#'   matrix. Exclude restores `NA` positions, including when every row is missing.
#' @param ... Additional `bboot.gir()` controls for an automatic bootstrap
#'   only, including `cluster` and `cluster_order`. Cluster vectors refer to
#'   the original fitting data or estimation sample, never `newdata`.
#'   These controls do not change `gir()` settings.
#'
#' @details
#' The target is \eqn{x'\theta}, not the response Y, its conditional mean, or
#' the first-stage conditional CDF. Under misspecification it is the fitted
#' linear-approximation target. Predictor terms, safe transformation
#' attributes, original contrasts and factor levels are retained;
#' unrecognized factor levels or missing predictor columns error.
#'
#' Given raw coefficient draws, predict each \eqn{x'\widehat\theta^{*(b)}}
#' and take the Algorithm 2 basic quantiles of its centered errors.
#' Prediction SE is \eqn{\sqrt{x'\widehat V_{boot}x}}, equivalently their
#' sample standard deviation. Do not combine marginal coefficient interval
#' endpoints or add an OLS outcome-noise variance. Uncertainty is pointwise,
#' not simultaneous. Point predictions outside \eqn{[0,1]} warn and are not
#' clipped; intervals are not clipped either.
#'
#' A pure point-prediction call uses no bootstrap. Supplied bootstrap draws
#' are preferred for repeated requests; otherwise uncertainty triggers
#' `bboot.gir()`. At least two successful draws are required. Original-sample
#' output contains estimation rows only, not rows omitted from fitting;
#' `row.index` maps them back. `newdata` `na.exclude` is reinserted in the
#' returned objects.
#'
#' Point predictions from aliased fits are returned only for estimable linear
#' functions, with nonestimable new rows `NA`. Bootstrap prediction
#' uncertainty requires a full-rank coefficient specification in this
#' version. The supplied inference theory and covariance moment
#' qualifications apply.
#'
#' @return A named prediction vector when `interval` is `"none"`, or a matrix
#'   with columns `fit`, `lwr` and `upr` for confidence intervals. With
#'   `se.fit=TRUE` return a list with `fit`, `se.fit`, `level`, `method` and
#'   `R`. The result has no `residual.scale` or Student degrees of freedom,
#'   since uncertainty is calculated from the bootstrap draws.
#'
#' @seealso [gir()], [bboot.gir()], [stats::predict.lm()].
#'
#' @examples
#' d <- data.frame(y = rep(c(1, 2, 2, 4, 5), 8),
#'                 g = factor(rep(c("a", "b"), each = 20)))
#' f <- gir(y ~ g, data = d, discrete = "g")
#' predict(f)
#' nd <- data.frame(g = c("a", "b"))
#' bb <- bboot.gir(f, R = 499L, seed = 9L)
#' predict(f, newdata = nd, bootstrap = bb,
#'         interval = "confidence", se.fit = TRUE)
#'
#' @method predict gir
#' @export
#' @rawNamespace export(predict.gir)
predict.gir <- function(
    object,
    newdata = NULL,
    se.fit = FALSE,
    interval = c("none", "confidence"),
    level = 0.95,
    bootstrap = NULL,
    R = 999L,
    seed = NULL,
    na.action = stats::na.pass,
    ...) {
  s <- .gir_state(object)
  .gir_flag(se.fit, "se.fit"); .gir_level(level)
  interval <- match.arg(interval)
  extra <- list(...)
  need_boot <- se.fit || interval == "confidence"
  if (!need_boot && (length(extra) || !missing(R) || !is.null(seed))) {
    stop("Bootstrap controls require se.fit=TRUE or interval='confidence'.", call. = FALSE)
  }
  mm <- .gir_new_matrix(object, newdata, na.action)
  X <- mm$X
  # For the estimation sample, X holds the distinct regressor rows and
  # `expand` maps the observations to them; everything is computed once per
  # distinct row and expanded. New data are used row by row.
  expand <- mm$group
  point_rows <- drop(X %*% object$coefficients.complete)
  valid_rows <- stats::complete.cases(X)
  if (any(object$aliased) && !is.null(newdata)) {
    est <- rep(FALSE, nrow(X))
    if (any(valid_rows)) {
      est[valid_rows] <- .gir_estimable_rows(object, X[valid_rows, , drop = FALSE])
    }
    bad <- valid_rows & !est
    if (any(bad)) {
      point_rows[bad] <- NA_real_
      warning(sum(bad), " newdata row(s) are non-estimable for this aliased fit; ",
              "predictions are NA.", call. = FALSE)
    }
    valid_rows <- valid_rows & est
  }
  point_rows[!valid_rows] <- NA_real_
  point <- if (is.null(expand)) point_rows else point_rows[expand]
  names(point) <- if (is.null(expand)) rownames(X) else names(object$midcdf)
  valid <- if (is.null(expand)) valid_rows else valid_rows[expand]
  if (any(!is.finite(point[valid]))) stop("Linear predictions overflowed.", call. = FALSE)
  outside <- sum(point[valid] < 0 | point[valid] > 1)
  if (outside) {
    warning(outside, " fitted Gini index value(s) are outside [0,1]. ",
            "The linear predictions have not been clipped.", call. = FALSE)
  }
  se <- rep(NA_real_, length(point)); names(se) <- names(point)
  result <- point
  B <- NULL
  if (need_boot) {
    if (any(object$aliased)) {
      stop("Bootstrap prediction inference currently requires a full-rank ",
           "coefficient specification.", call. = FALSE)
    }
    if (is.null(bootstrap)) {
      forbidden <- intersect(names(extra), c("object", "R", "level", "seed"))
      if (length(forbidden)) stop("Duplicated bootstrap controls.", call. = FALSE)
      bootstrap <- do.call(bboot.gir, c(
        list(object = object, R = R, level = level, seed = seed), extra
      ))
    } else {
      if (length(extra) || !missing(R) || !is.null(seed)) {
        stop("Do not supply R, seed, or extra bootstrap controls with an ",
             "existing bootstrap result.", call. = FALSE)
      }
    }
    .gir_boot_compatible(object, bootstrap)
    draws <- bootstrap$t[bootstrap$success, , drop = FALSE]
    B <- nrow(draws)
    if (B < 2L) stop("Prediction uncertainty requires at least two draws.", call. = FALSE)
    if (interval == "confidence" && B * (1 - level) / 2 < 10) {
      warning("Few bootstrap draws per prediction-interval tail; use more replicates.",
              call. = FALSE)
    }
    V <- stats::cov(draws)
    se_rows <- rep(NA_real_, nrow(X))
    if (any(valid_rows)) {
      XV <- X[valid_rows, , drop = FALSE]
      variance <- rowSums((XV %*% V) * XV)
      se_rows[valid_rows] <- sqrt(pmax(0, variance))
    }
    se[] <- if (is.null(expand)) se_rows else se_rows[expand]
    if (interval == "confidence") {
      ci_rows <- matrix(NA_real_, nrow(X), 2L)
      rows <- which(valid_rows)
      # Bound temporary memory: never allocate n_predictions-by-R in full.
      block <- max(1L, min(1000L, floor(2e6 / B)))
      if (length(rows)) {
        for (start in seq.int(1L, length(rows), by = block)) {
          index <- rows[start:min(length(rows), start + block - 1L)]
          predictions_star <- draws %*% t(X[index, , drop = FALSE])
          ci_rows[index, ] <- .gir_basic_ci(point_rows[index], predictions_star,
                                            level)
        }
      }
      result <- matrix(NA_real_, length(point), 3L,
                       dimnames = list(names(point), c("fit", "lwr", "upr")))
      result[, "fit"] <- point
      result[, c("lwr", "upr")] <- if (is.null(expand)) ci_rows else
        ci_rows[expand, , drop = FALSE]
    }
  }
  # Match na.exclude reinsertion for newdata; original-sample output always
  # follows the retained estimation rows, consistently with fitted.gir().
  if (!is.null(mm$action) && inherits(mm$action, "exclude")) {
    result <- stats::napredict(mm$action, result)
    se <- stats::napredict(mm$action, se)
  }
  if (se.fit) {
    return(list(fit = result, se.fit = se, level = level,
                method = "Bayesian bootstrap; pointwise linear-index inference", R = B))
  }
  result
}
