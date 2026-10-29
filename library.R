# HF-GNAR empirical functions
# Start with fit_gnar(): the full estimation flow is written first.
# Helpers then define least squares, initial networks, the two MIP models,
# and reporting. Sourcing this file defines functions only.
# Empirical implementation: all9_v4. The runner supplies the paper's settings.
# Reporting uses combined specification labels below the membership figures.

# ============================================================================
# 1. One empirical specification: initial network -> LS -> MIP -> MIP-it
# ============================================================================

fit_gnar <- function(y, method, budget, seed, threads = 10L, time_limit = 200,
                     max_iter = 2L, min_improvement = .001) {
  y <- unname(as.matrix(y))
  method <- match.arg(method, c("SIS", "BS", "AL"))
  d <- nrow(y)
  if (!is.numeric(y) || ncol(y) < 3L || d < 2L || any(!is.finite(y))) {
    stop("y must be a finite numeric matrix with markets in rows.")
  }
  if (length(budget) != 1L || !is.finite(budget) || budget != floor(budget) ||
      budget < 1L || budget > d * (d - 1L)) stop("Invalid off-diagonal edge budget.")
  if (max_iter < 1L || max_iter != floor(max_iter)) stop("max_iter must be positive.")
  centered_y <- sweep(y, 1L, rowMeans(y), "-")
  total_edges <- as.integer(budget + d)

  # These are the two original RNG streams. The function intentionally sets
  # the RNG state directly; it does not install or replace global functions.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed + 1000L))
  message("  Initial network: ", method)
  sis_B <- sis_initial(y, total_edges)
  initial_B <- switch(method,
    SIS = sis_B,
    BS = bs_initial(y, total_edges, sis_B, threads, time_limit)$B,
    AL = al_initial(y)
  )
  initial_A <- adjacency_matrix(initial_B)
  initial_group <- initial_groups(initial_B, y)
  LS <- grouped_ls(centered_y, initial_A, initial_group)
  post_ls_group <- initial_groups(LS$B, y)
  message("  Initial fixed-network group MIP")
  MIP <- fixed_network_fit(y, initial_A, threads, time_limit)

  # Preserve the historical mixed initial state. It is not replaced by MIP
  # unless a network/group proposal passes the raw-data RSS acceptance rule.
  B <- as.matrix(initial_B)
  group <- as.integer(post_ls_group)
  beta <- as.numeric(MIP$beta)
  rss_raw <- gnar_rss(y, B)
  accepted_iterations <- 0L
  stop_reason <- "max_iter"
  history <- vector("list", as.integer(max_iter))

  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed + 30000L))
  for (iteration in seq_len(as.integer(max_iter))) {
    message("  Refinement ", iteration, ": full network MIP")
    network <- network_mip(y, group, beta, budget, threads, time_limit)
    message("  Refinement ", iteration, ": fixed-network group MIP")
    proposal <- fixed_network_fit(y, network$A, threads, time_limit)
    proposal_rss <- gnar_rss(y, proposal$B)
    change <- (proposal_rss - rss_raw) / rss_raw
    accepted <- is.finite(change) && change <= -min_improvement
    history[[iteration]] <- data.frame(
      iteration = iteration, accepted = accepted,
      RSS_before_raw = rss_raw, RSS_proposal_raw = proposal_rss,
      relative_change_raw = change,
      network_status = network$status, network_runtime = network$runtime,
      network_gap = network$gap,
      group_status = proposal$status, group_runtime = proposal$solver$runtime,
      group_gap = proposal$solver$relative_gap,
      edges_before = sum(adjacency_matrix(B)),
      edges_candidate = sum(adjacency_matrix(proposal$B))
    )
    message(sprintf("  Raw RSS change: %.4f%%; %s", 100 * change,
                    if (accepted) "accepted" else "rejected"))
    if (!accepted) {
      stop_reason <- "RSS_rejected"
      break
    }
    B <- proposal$B
    group <- as.integer(proposal$group)
    beta <- as.numeric(proposal$beta)
    rss_raw <- proposal_rss
    accepted_iterations <- accepted_iterations + 1L
  }
  MIP_it <- list(
    B = B, A = adjacency_matrix(B), group = group, beta = beta,
    RSS = gnar_rss(centered_y, B), n_iter = accepted_iterations,
    n_attempted = iteration, stop_reason = stop_reason,
    status = if (accepted_iterations == 0L) "HISTORICAL_INITIAL_STATE_RETAINED"
             else "LAST_ACCEPTED_MIP",
    coefficient_group = group,
    coefficient_source = if (accepted_iterations == 0L) "mixed_initial_Sep30"
                         else "last_accepted_Step4"
  )
  # LS reporting uses post-LS labels without changing the fitted coefficients.
  LS$coefficient_group <- LS$group
  LS$group <- post_ls_group
  LS$coefficient_source <- "LS_before_regroup"
  MIP$coefficient_group <- MIP$group
  MIP$coefficient_source <- "Step2_MIP"
  list(LS = LS, MIP = MIP, MIP_it = MIP_it,
       history = do.call(rbind, history[seq_len(iteration)]),
       initial_B = initial_B, initial_A = initial_A,
       initial_group = post_ls_group, initial_beta = MIP$beta)
}

# ============================================================================
# 2. GNAR coefficients, residuals, and unrestricted grouped least squares
#    Data matrices have markets in rows and observations in columns.
#    Internal group labels are 1 and 0, in that order in beta.
# ============================================================================

adjacency_matrix <- function(B) {
  A <- 1L * (abs(as.matrix(B)) > 0)
  diag(A) <- 0L
  A
}

gnar_matrix <- function(A, group, beta) {
  A <- adjacency_matrix(A)
  block <- 2L - as.integer(group)
  phi <- matrix(beta[1:4], 2L, 2L, byrow = TRUE)
  B <- (A / pmax(1, rowSums(A))) * phi[block, block, drop = FALSE]
  diag(B) <- beta[4L + block]
  dimnames(B) <- dimnames(A)
  B
}

gnar_rss <- function(y, B) {
  residual <- y[, -1L, drop = FALSE] - B %*% y[, -ncol(y), drop = FALSE]
  sum(residual^2)
}

grouped_ls <- function(y, A, group) {
  A <- adjacency_matrix(A)
  group <- as.integer(group)
  d <- nrow(y)
  nt <- ncol(y) - 1L
  X <- y[, seq_len(nt), drop = FALSE]
  Y <- y[, seq.int(2L, nt + 1L), drop = FALSE]
  block <- 2L - group
  degree <- rowSums(A)
  design <- matrix(0, d * nt, 6L)
  for (i in seq_len(d)) {
    rows <- seq.int((i - 1L) * nt + 1L, i * nt)
    if (degree[i] > 0) {
      for (j in 1:2) {
        nodes <- which(A[i, ] != 0L & block == j)
        if (length(nodes)) {
          design[rows, (block[i] - 1L) * 2L + j] <-
            colSums(X[nodes, , drop = FALSE]) / degree[i]
        }
      }
    }
    design[rows, 4L + block[i]] <- X[i, ]
  }
  response <- as.vector(t(Y))
  decomposition <- svd(design, nu = min(nrow(design), 6L), nv = 6L)
  beta <- numeric(6L)
  rank <- 0L
  if (length(decomposition$d) && max(decomposition$d) > 0) {
    tolerance <- max(dim(design)) * max(decomposition$d) * .Machine$double.eps
    keep <- which(decomposition$d > tolerance)
    rank <- length(keep)
    if (rank > 0L) {
      beta <- as.numeric(decomposition$v[, keep, drop = FALSE] %*%
        (crossprod(decomposition$u[, keep, drop = FALSE], response) /
           decomposition$d[keep]))
    }
  }
  names(beta) <- c("phi11", "phi12", "phi21", "phi22", "nu1", "nu2")
  B <- gnar_matrix(A, group, beta)
  list(B = B, A = A, group = group, beta = beta, RSS = gnar_rss(y, B),
       status = "LS", ls_rank = rank)
}

# ============================================================================
# 3. Initial networks and groups
#    SIS and BS retain the original initializer and its numerical operations.
# ============================================================================

marginal_var <- function(y) {
  y <- sweep(y, 1L, rowMeans(y), "-")
  X <- y[, -ncol(y), drop = FALSE]
  Y <- y[, -1L, drop = FALSE]
  denominator <- rowSums(X^2)
  B <- tcrossprod(Y, X)
  valid <- denominator > 0
  B[, valid] <- sweep(B[, valid, drop = FALSE], 2L, denominator[valid], "/")
  B[, !valid] <- 0
  B
}

sis_initial <- function(y, total_edges) {
  y <- y - rowMeans(y)
  d <- nrow(y)
  y <- diag(d) %*% y
  X <- y[, seq_len(ncol(y) - 1L), drop = FALSE]
  Y <- y[, 2:ncol(y), drop = FALSE]
  beta <- as.vector(tcrossprod(Y, X)) / rep(rowSums(X^2), each = d)
  ranking <- matrix(beta, d)
  diag(ranking) <- max(abs(beta))
  selected <- order(abs(as.vector(ranking)), decreasing = TRUE)[seq_len(total_edges)]
  beta[-selected] <- 0
  matrix(beta, d)
}

initial_groups <- function(B, y) {
  d <- nrow(B)
  nstart <- floor(d / 2)
  fixed <- rowMeans(y[, 2:ncol(y)]) - B %*% rowMeans(y[, -ncol(y)])
  values <- as.vector(B)
  nonzero <- which(values != 0)
  network <- NULL
  if (length(unique(values[nonzero])) > 4L) {
    clusters <- stats::kmeans(values[nonzero], centers = 4L, nstart = nstart)
    mapped <- 0 * values
    mapped[nonzero] <- clusters$cluster
    mapped <- matrix(mapped, d, d)
    network <- matrix(0, d, 4L)
    for (i in seq_len(d)) {
      for (j in 1:4) {
        nodes <- which(mapped[i, ] == j)
        if (length(nodes)) network[i, j] <- sum(B[i, nodes]) / length(nodes)
      }
    }
  }
  decomposition <- svd(B + t(B), nu = 2L, nv = 2L)
  spectral <- matrix(0, d, 2L)
  for (i in seq_len(d)) {
    spectral[i, ] <- decomposition$u[i, ] /
      (base::norm(decomposition$u[i, ], "2") + 1e-9)
  }
  features <- cbind(diag(B), fixed, network)
  # All four PAM calls are retained in their original order. The unused calls
  # consume random numbers and therefore affect subsequent initialization.
  community <- cluster::pam(spectral, k = 2L)
  grouping <- cluster::pam(features, k = 2L, nstart = nstart)
  combined <- cluster::pam(cbind(features, spectral), k = 2L, nstart = nstart)
  reduced <- cluster::pam(cbind(fixed, network), k = 2L, nstart = nstart)
  2 - grouping$clustering
}

detect_communities <- function(B) {
  support <- 1L * (abs(B) > 1e-8)
  d <- nrow(support)
  decomposition <- svd(support + t(support), nu = min(2L, d), nv = 0)
  spectral <- matrix(0, d, 2L)
  spectral[, seq_len(ncol(decomposition$u))] <- decomposition$u
  spectral <- spectral / (sqrt(rowSums(spectral^2)) + 1e-9)
  if (nrow(unique(as.data.frame(spectral))) < 2L) {
    return(2L - c(rep(1L, ceiling(d / 2)), rep(2L, floor(d / 2))))
  }
  2L - cluster::pam(spectral, k = 2L, nstart = 1L)$clustering
}

# BS includes its original bounded coefficient fit. This is distinct from the
# unrestricted grouped_ls() used for the LS benchmark and group-MIP starts.
bs_coefficient_objective <- function(beta, y, group, A, return_matrix = FALSE) {
  if (any(beta < 0) || any(beta > 1)) return(1e10)
  d <- nrow(y)
  n_group1 <- sum(group)
  order <- c(which(group == 1), which(group == 0))
  ordered_y <- rbind(y[group == 1, ], y[group == 0, ])
  ordered_A <- (1 * (A != 0))[order, order]
  diag(ordered_A) <- 0
  weights <- diag(1 / rowSums(ordered_A))
  weights[is.infinite(weights)] <- 0
  normalized_A <- weights %*% ordered_A
  block_beta <- rbind(
    cbind(matrix(beta[1], n_group1, n_group1),
          matrix(beta[2], n_group1, d - n_group1)),
    cbind(matrix(beta[3], d - n_group1, n_group1),
          matrix(beta[4], d - n_group1, d - n_group1))
  )
  B <- block_beta * normalized_A +
    diag(c(rep(beta[5], n_group1), rep(beta[6], d - n_group1)))
  if (return_matrix) return(B)
  gnar_rss(ordered_y, B)
}

bs_coefficient_fit <- function(y, A, group) {
  A <- 1 * (A != 0)
  diag(A) <- 0
  start <- c(0.2, 0, 0, 0.2, 0.1, 0.1)
  bounded <- tryCatch(
    stats::optim(start, bs_coefficient_objective, y = y, group = group, A = A,
                 method = "L-BFGS-B", lower = 0, upper = 1,
                 control = list(maxit = 200, trace = FALSE,
                                factr = 1e12, pgtol = 1e-4)),
    error = function(e) list(par = start, value = Inf, convergence = 999)
  )
  simplex <- stats::optim(
    pmax(0, pmin(1, bounded$par)), bs_coefficient_objective,
    y = y, group = group, A = A, method = "Nelder-Mead",
    control = list(maxit = 200, trace = FALSE, reltol = 1e-8)
  )
  beta <- pmax(0, pmin(1, simplex$par))
  ordered_B <- bs_coefficient_objective(beta, y, group, A, return_matrix = TRUE)
  order <- c(which(group == 1), which(group == 0))
  B <- 0 * ordered_B
  B[order, order] <- ordered_B
  list(B = B, beta = beta, group = group)
}

bs_initial <- function(y, total_edges, sis_B, threads, time_limit) {
  y <- y - rowMeans(y)
  d <- nrow(y)
  adjusted_y <- diag(d) %*% y
  X <- adjusted_y[, seq_len(ncol(y) - 1L), drop = FALSE]
  Y <- adjusted_y[, 2:ncol(y), drop = FALSE]
  nB <- d * d
  off <- which(row(matrix(0, d, d)) != col(matrix(0, d, d)))
  selector <- Matrix::sparseMatrix(i = seq_along(off), j = off, x = 1,
                                   dims = c(length(off), nB))
  big_m <- 2 * max(abs(sis_B))
  start_A <- 1 * (sis_B != 0)
  diag(start_A) <- 1
  lower_A <- matrix(0, d, d)
  diag(lower_A) <- 1
  gram <- Matrix::Matrix(kronecker(tcrossprod(X), diag(d)), sparse = TRUE)
  model <- list(
    A = rbind(cbind(selector, -big_m * selector),
              cbind(-selector, -big_m * selector),
              c(rep(0, nB), rep(1, nB))),
    sense = c(rep("<", 2L * length(off)), "="),
    rhs = c(rep(0, 2L * length(off)), total_edges),
    ub = c(rep(big_m, nB), rep(1, nB)),
    lb = c(rep(-big_m, nB), as.vector(lower_A)),
    obj = c(-2 * as.vector(tcrossprod(Y, X)), rep(0, nB)),
    Q = Matrix::bdiag(gram, Matrix::Matrix(0, nB, nB)),
    start = c(as.vector(sis_B), as.vector(start_A))
  )
  # VAR coefficients are continuous; support indicators are explicitly binary.
  model[["vtype"]] <- c(rep("C", nB), rep("B", nB))
  # Keep the BS defaults: no explicit MIPGap, Seed, Cuts, or Presolve setting.
  solution <- gurobi::gurobi(model, list(TimeLimit = time_limit,
                                        Threads = threads, OutputFlag = 0L))
  B <- if (is.null(solution$x)) sis_B else matrix(solution$x[seq_len(nB)], d)
  group <- initial_groups(B, y)
  fitted <- bs_coefficient_fit(y, B, group)
  list(B = fitted$B, status = solution$status)
}

al_initial <- function(y) {
  if (as.character(utils::packageVersion("sparseVAR")) != "0.4.4") {
    stop("AL initialization requires sparseVAR 0.4.4.")
  }
  fit <- sparseVAR::sVAR_adalasso_fista(
    Yt = sweep(y, 1L, rowMeans(y), "-"), p = 1, fold = 5L,
    penalize_diag = FALSE, standardize = TRUE, whiten = "none",
    lambda_rule = "cv", nlambda = 100L, max_iter = 1500L, tol = 1e-6,
    adaptive_zero = "finite", adaptive_cv = "fixed_weights",
    eps_w = 1e-8, cap_w = 1e6, debiasTF = FALSE, pilot = "lasso"
  )
  if (!identical(fit$pilot, "lasso") || !identical(fit$penalize_diag, FALSE)) {
    stop("Unexpected adaptive-lasso pilot or diagonal penalty.")
  }
  B <- as.matrix(fit$Phi_hat_Alasso)
  if (!identical(dim(B), c(nrow(y), nrow(y))) || any(!is.finite(B))) {
    stop("Invalid adaptive-lasso coefficient matrix.")
  }
  1L * (abs(B) > 1e-8)
}

# ============================================================================
# 4. Fixed-network group MIP
#    The variable order and constraint order match the empirical computation.
# ============================================================================

mip_feasible <- function(model, x, tolerance = 1e-6) {
  if (is.null(x) || length(x) != length(model$obj) || any(!is.finite(x))) {
    return(FALSE)
  }
  if (any(x < model$lb - tolerance | x > model$ub + tolerance)) return(FALSE)
  integer <- which(model$vtype %in% c("B", "I"))
  if (any(abs(x[integer] - round(x[integer])) > tolerance)) return(FALSE)
  difference <- as.numeric(model$A %*% x) - model$rhs
  if (any(difference[model$sense == "<"] > tolerance) ||
      any(difference[model$sense == ">"] < -tolerance) ||
      any(abs(difference[model$sense == "="]) > tolerance)) return(FALSE)
  for (indicator in model$genconind) {
    if (round(x[indicator$binvar]) == indicator$binval) {
      difference <- sum(as.numeric(indicator$a) * x) - indicator$rhs
      if (abs(difference) > tolerance) return(FALSE)
    }
  }
  TRUE
}

group_model <- function(y, A, start) {
  A <- adjacency_matrix(A)
  d <- nrow(y)
  warm <- list(A = A, group = start$group, beta = pmax(.05, pmin(1, start$beta)))
  names(warm$beta) <- c("phi11", "phi12", "phi21", "phi22", "nu1", "nu2")
  warm$B <- gnar_matrix(A, warm$group, warm$beta)
  warm$RSS <- gnar_rss(y, warm$B)
  degree <- rowSums(A)
  edges <- which(A != 0L, arr.ind = TRUE)
  nedge <- nrow(edges)
  nB <- d * d
  base <- nB + 6L + 2L * d
  nvars <- base + 4L * nedge
  coefficient_index <- function(i, j) i + (j - 1L) * d
  phi_index <- function(k, l) nB + k + (l - 1L) * 2L
  group_index <- function(i, k) nB + 6L + (i - 1L) * 2L + k

  # Assemble sparse linear constraints, in the same order as their equations.
  columns <- values <- list()
  rhs <- numeric()
  senses <- character()
  add_row <- function(index, value, bound, sense = "<") {
    row <- length(rhs) + 1L
    keep <- value != 0
    columns[[row]] <<- as.integer(index[keep])
    values[[row]] <<- value[keep]
    rhs[row] <<- bound
    senses[row] <<- sense
  }
  indicators <- list()
  for (i in seq_len(d)) {
    add_row(group_index(i, 1:2), c(1, 1), 1, "=")
    for (k in 1:2) {
      indicators[[length(indicators) + 1L]] <- list(
        binvar = as.integer(group_index(i, k)), binval = 1L,
        a = Matrix::sparseVector(i = c(coefficient_index(i, i), nB + 4L + k),
                                  x = c(1, -1), length = nvars),
        sense = "=", rhs = 0
      )
    }
  }
  for (e in seq_len(nedge)) {
    i <- edges[e, 1L]
    j <- edges[e, 2L]
    for (k in 1:2) {
      for (l in 1:2) {
        w <- base + (e - 1L) * 4L + (k - 1L) * 2L + l
        add_row(c(w, group_index(i, k)), c(1, -1), 0)
        add_row(c(w, group_index(j, l)), c(1, -1), 0)
        add_row(c(w, group_index(i, k), group_index(j, l)), c(-1, 1, 1), 1)
        indicators[[length(indicators) + 1L]] <- list(
          binvar = as.integer(w), binval = 1L,
          a = Matrix::sparseVector(i = c(coefficient_index(i, j), phi_index(k, l)),
                                    x = c(degree[i], -1), length = nvars),
          sense = "=", rhs = 0
        )
      }
    }
  }
  X <- y[, -ncol(y), drop = FALSE]
  Y <- y[, -1L, drop = FALSE]
  model <- list(
    A = Matrix::sparseMatrix(i = rep(seq_along(columns), lengths(columns)),
                             j = as.integer(unlist(columns, use.names = FALSE)),
                             x = as.numeric(unlist(values, use.names = FALSE)),
                             dims = c(length(rhs), nvars)),
    rhs = rhs, sense = senses,
    Q = Matrix::bdiag(
      Matrix::Matrix(kronecker(tcrossprod(X), diag(d)), sparse = TRUE),
      Matrix::Matrix(0, nvars - nB, nvars - nB, sparse = TRUE)),
    obj = c(-2 * as.vector(tcrossprod(Y, X)), rep(0, nvars - nB)),
    objcon = sum(Y^2), modelsense = "min",
    lb = c(rep(-Inf, nB + 6L), rep(0, nvars - nB - 6L)),
    ub = c(rep(Inf, nB + 6L), rep(1, nvars - nB - 6L)),
    vtype = c(rep("C", nB + 6L), rep("B", nvars - nB - 6L)),
    genconind = indicators, modelname = "GNAR88_group"
  )
  model$lb[nB + seq_len(6L)] <- .05
  model$ub[nB + seq_len(6L)] <- 1
  absent <- which(A == 0L & row(A) != col(A))
  model$lb[absent] <- model$ub[absent] <- 0
  group <- as.integer(start$group)
  group_id <- 2L - group
  conjunction <- numeric(4L * nedge)
  for (e in seq_len(nedge)) {
    conjunction[(e - 1L) * 4L + (group_id[edges[e, 1L]] - 1L) * 2L +
                  group_id[edges[e, 2L]]] <- 1
  }
  model$start <- c(as.vector(warm$B), as.numeric(warm$beta)[c(1, 3, 2, 4, 5, 6)],
                   as.vector(rbind(group, 1L - group)), conjunction)
  list(model = model, warm = warm)
}

group_mip <- function(y, A, start, threads, time_limit) {
  specification <- group_model(y, A, start)
  model <- specification$model
  warm <- specification$warm
  if (!mip_feasible(model, model$start)) stop("Infeasible group-MIP starting value.")
  parameters <- list(MIPGap = .01, TimeLimit = time_limit, Threads = threads,
                     Cuts = 2L, Presolve = 2L, Heuristics = .5, OutputFlag = 0L,
                     Seed = 20261009L, FeasibilityTol = 1e-8, IntFeasTol = 1e-8)
  solution <- gurobi::gurobi(model, parameters)
  if (length(solution$status) != 1L ||
      !solution$status %in% c("OPTIMAL", "TIME_LIMIT") ||
      !mip_feasible(model, solution$x)) {
    stop("Group MIP has no valid incumbent: ", solution$status)
  }
  d <- nrow(y)
  nB <- d * d
  group <- as.integer(matrix(solution$x[nB + 6L + seq_len(2L * d)], 2L, d)[1L, ] >= .5)
  beta <- as.numeric(solution$x[nB + c(1, 3, 2, 4, 5, 6)])
  names(beta) <- c("phi11", "phi12", "phi21", "phi22", "nu1", "nu2")
  A <- adjacency_matrix(A)
  B <- gnar_matrix(A, group, beta)
  if (max(abs(matrix(solution$x[seq_len(nB)], d, d) - B)) > 1e-5) {
    stop("Group-MIP coefficients disagree with the fitted group structure.")
  }
  rss <- gnar_rss(y, B)
  objective <- sum(solution$x * as.numeric(model$Q %*% solution$x)) +
    sum(model$obj * solution$x) + model$objcon
  if (abs(objective - rss) > 1e-5 * max(1, abs(rss))) {
    stop("Group-MIP objective and residual sum of squares disagree.")
  }
  if (rss > warm$RSS + 1e-9 * max(1, abs(warm$RSS))) {
    stop("Group-MIP incumbent is worse than its bounded LS starting value.")
  }
  bound <- if (length(solution$objbound) == 1L && is.finite(solution$objbound)) {
    max(0, as.numeric(solution$objbound))
  } else NA_real_
  if (is.finite(bound) && bound > rss + 1e-5 * max(1, abs(rss))) bound <- NA_real_
  gap <- if (is.finite(bound)) max(0, rss - bound) / max(1e-12, abs(rss)) else NA_real_
  storage.mode(A) <- "integer"
  B <- gnar_matrix(A, group, beta)
  list(B = B, A = A, group = group, beta = beta, RSS = gnar_rss(y, B),
       status = solution$status, is_MIP = TRUE, coefficient_source = "group_MIP",
       solver = list(runtime = solution$runtime, relative_gap = gap))
}

fixed_network_fit <- function(y, A, threads, time_limit) {
  centered_y <- sweep(y, 1L, rowMeans(y), "-")
  A <- adjacency_matrix(A)
  sis_B <- sis_initial(y, sum(A) + nrow(y))
  group <- initial_groups(sis_B, y)
  start <- grouped_ls(centered_y, A, group)
  group_mip(centered_y, A, start, threads, time_limit)
}

# ============================================================================
# 5. Full quadratic network update, with the original cold solver start
# ============================================================================

network_model <- function(y, group, beta, budget) {
  centered_y <- sweep(y, 1L, rowMeans(y), "-")
  d <- nrow(y)
  X <- centered_y[, -ncol(centered_y), drop = FALSE]
  block <- 2L - as.integer(group)
  residual <- centered_y[, -1L, drop = FALSE] - beta[4L + block] * X
  nC <- d * d
  nvars <- 2L * nC
  off <- which(row(matrix(0L, d, d)) != col(matrix(0L, d, d)))
  # Alternating rows impose C_ij <= A_ij and -C_ij <= 0.1 A_ij.
  row_ids <- rep(seq_len(2L * length(off)), each = 2L)
  column_ids <- as.vector(rbind(rep(off, each = 2L), nC + rep(off, each = 2L)))
  row_values <- rep(c(1, -1, -1, -.1), length(off))
  row_ids <- c(row_ids, rep(2L * length(off) + 1L, length(off)))
  column_ids <- c(column_ids, nC + off)
  row_values <- c(row_values, rep(1, length(off)))
  model <- list(
    A = Matrix::sparseMatrix(i = row_ids, j = column_ids, x = row_values,
                             dims = c(2L * length(off) + 1L, nvars)),
    rhs = c(rep(0, 2L * length(off)), as.integer(budget)),
    sense = rep("<", 2L * length(off) + 1L)
  )
  q_i <- q_j <- q_value <- vector("list", d)
  objective <- numeric(nvars)
  for (i in seq_len(d)) {
    nodes <- setdiff(seq_len(d), i)
    indices <- i + (nodes - 1L) * d
    Z <- t(X[nodes, , drop = FALSE])
    q_i[[i]] <- rep(indices, times = length(indices))
    q_j[[i]] <- rep(indices, each = length(indices))
    q_value[[i]] <- as.vector(crossprod(Z))
    objective[indices] <- -2 * as.numeric(crossprod(Z, residual[i, ]))
  }
  model$Q <- Matrix::sparseMatrix(
    i = as.integer(unlist(q_i, use.names = FALSE)),
    j = as.integer(unlist(q_j, use.names = FALSE)),
    x = as.numeric(unlist(q_value, use.names = FALSE)), dims = c(nvars, nvars))
  model$obj <- objective
  model$modelsense <- "min"
  model$lb <- c(rep(-.1, nC), rep(0, nC))
  model$ub <- rep(1, nvars)
  diagonal <- seq.int(1L, nC, by = d + 1L)
  model$lb[diagonal] <- model$ub[diagonal] <- 0
  model$lb[nC + diagonal] <- model$ub[nC + diagonal] <- 1
  model$vtype <- c(rep("C", nC), rep("B", nC))
  model$modelname <- "GNAR_Sep30_rules_full"
  # The objective omits the constant ||residual||^2, as in the original solver.
  # No model$start is supplied: this is a cold network optimization.
  list(model = model, X = X, residual = residual, budget = as.integer(budget))
}

network_mip <- function(y, group, beta, budget, threads, time_limit) {
  specification <- network_model(y, group, beta, budget)
  model <- specification$model
  parameters <- list(MIPGap = .01, Cuts = 2L, Presolve = 2L, Heuristics = .5,
                     TimeLimit = as.numeric(time_limit), Threads = as.integer(threads),
                     OutputFlag = 0L)
  solution <- gurobi::gurobi(model, parameters)
  if (!mip_feasible(model, solution$x)) {
    stop("Network MIP has no feasible incumbent: ", solution$status)
  }
  d <- nrow(y)
  nC <- d * d
  A <- matrix(as.integer(round(solution$x[nC + seq_len(nC)])), d, d)
  diag(A) <- 0L
  C <- matrix(solution$x[seq_len(nC)], d, d)
  C[A == 0L] <- 0
  diag(C) <- 0
  if (sum(A) > budget) stop("Network proposal exceeds the edge budget.")
  free_rss <- sum((specification$residual - C %*% specification$X)^2)
  objective <- as.numeric(solution$objval) + sum(specification$residual^2)
  if (abs(objective - free_rss) > 1e-5 * max(1, abs(free_rss))) {
    stop("Full network objective and residual sum of squares disagree.")
  }
  list(A = A, C = C, free_RSS = free_rss, status = as.character(solution$status),
       runtime = solution$runtime, gap = solution$mipgap)
}


# ============================================================================
# 6. Tables and membership figures
# ============================================================================

# ---- Membership figures --------------------------------------------------------------

plot_memberships <- function(membership, title, fill, filename) {
    nr <- nrow(membership)
    nc <- ncol(membership)
    stopifnot(!is.null(rownames(membership)), !is.null(colnames(membership)),
              all(membership %in% c(1L, 2L)))
    # Landscape layout with embedded fonts when Cairo is available.
    if (capabilities("cairo")) {
        grDevices::cairo_pdf(filename, width = 11.4359375, height = 7.1523958,
                            pointsize = 12, family = "sans", bg = "white")
    } else {
        grDevices::pdf(filename, width = 11.4359375, height = 7.1523958,
                       pointsize = 12, family = "Helvetica", paper = "special",
                       useDingbats = FALSE, title = title)
    }
    on.exit(grDevices::dev.off())
    graphics::par(mai = c(28.54, 99.07, 22.04, 2.88) / 72,
                  xaxs = "i", yaxs = "i", xpd = FALSE, family = "sans")
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, nc), ylim = c(nr, 0))
    for (i in seq_len(nr)) for (j in seq_len(nc)) {
        graphics::rect(j - 1, i - 1, j, i, border = NA,
                       col = if (membership[i, j] == 1L) fill else "white")
    }
    graphics::abline(h = 0:nr, v = 0:nc, col = "grey80", lwd = 0.73)
    graphics::abline(h = 6, col = "black", lwd = 1.53)
    graphics::box(lwd = 0.93)
    # Use all labels explicitly; axis() can omit labels that it regards as close.
    graphics::text(-0.08, seq_len(nr) - 0.5, labels = rownames(membership),
                   adj = c(1, 0.5), cex = 12.5 / 12, xpd = NA)
    graphics::text(seq_len(nc) - 0.5, nr + 0.75,
                   labels = colnames(membership), adj = c(0.5, 1),
                   cex = 12.5 / 12, family = "mono", xpd = NA)
    graphics::mtext(title, side = 3, line = 0.25, font = 2, cex = 16 / 12,
                    xpd = NA)
}

# ---- Final estimates, group assignments, and iteration records ------------

write_results <- function(results, labels, outdir,
                          community_seed = 20261011L, anchor = "us.SPX") {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    coefficients <- market_groups <- block_effects <- flow <- list()
    groups <- communities <- matrix(NA_integer_, length(labels), 9)
    column_names <- character(9)
    i <- 0L

    for (rv in c("RV10", "RV15", "RV20")) {
        for (method in c("SIS", "BS", "AL")) {
            i <- i + 1L
            result <- results[[paste(rv, method, sep = "_")]]
            f <- result$fit$MIP_it
            history <- result$fit$history
            group <- 2L - as.integer(f$group)
            beta <- as.numeric(f$beta)

            # Display group 1 contains the anchor. Switch its beta blocks too.
            if (group[match(anchor, labels)] == 2L) {
                group <- 3L - group
                beta <- beta[c(4, 3, 2, 1, 6, 5)]
            }
            names(beta) <- c("phi11", "phi12", "phi21", "phi22", "nu1", "nu2")

            # The original community calculation resets this stream for each
            # support matrix and includes the diagonal in abs(B) > 1e-8.
            RNGkind("L'Ecuyer-CMRG")
            set.seed(as.integer(community_seed))
            community <- 2L - detect_communities(f$B)
            if (sum(group != 3L - community) < sum(group != community))
                community <- 3L - community

            # Record an inconsistent retained initial state without refitting.
            W <- sweep(f$A, 1L, pmax(rowSums(f$A), 1), "/")
            Bcheck <- W * matrix(beta[1:4], 2, byrow = TRUE)[group, group]
            diag(Bcheck) <- beta[4L + group]
            state_error <- max(abs(f$B - Bcheck))
            coefficients[[i]] <- data.frame(
                RV = rv, method = method,
                n_G1 = sum(group == 1L), n_G2 = sum(group == 2L),
                as.list(beta), RSS = f$RSS, A_edges = sum(f$A),
                budget = result$budget,
                rho = max(Mod(eigen(f$B, only.values = TRUE)$values)),
                status = as.character(f$status),
                accepted_iterations = sum(history$accepted),
                state_error = state_error, state_consistent = state_error < 1e-5,
                stringsAsFactors = FALSE, check.names = FALSE)
            market_groups[[i]] <- data.frame(
                RV = rv, method = method, market = labels,
                group = group, community = community, stringsAsFactors = FALSE)

            # B[receiver, source]: exclude own lags and average by receivers.
            Boff <- f$B
            diag(Boff) <- 0
            E <- matrix(0, 2, 2)
            for (r in 1:2) for (s in 1:2) {
                E[r, s] <- sum(Boff[group == r, group == s, drop = FALSE]) /
                    sum(group == r)
            }
            block_effects[[i]] <- data.frame(
                RV = rv, method = method,
                E11 = E[1, 1], E12 = E[1, 2], E21 = E[2, 1], E22 = E[2, 2],
                E21_minus_E12 = E[2, 1] - E[1, 2],
                E11_minus_E22 = E[1, 1] - E[2, 2],
                E21_over_E12 = if (E[1, 2] == 0) NA_real_ else E[2, 1] / E[1, 2],
                stringsAsFactors = FALSE)
            flow[[i]] <- data.frame(RV = rv, method = method, history,
                                    stringsAsFactors = FALSE, check.names = FALSE)
            groups[, i] <- group
            communities[, i] <- community
            column_names[i] <- paste(method, tolower(rv), sep = "-")
        }
    }
    coefficients <- do.call(rbind, coefficients)
    market_groups <- do.call(rbind, market_groups)
    block_effects <- do.call(rbind, block_effects)
    flow <- do.call(rbind, flow)

    # Cross-classification displayed in the paper: AL at RV20.
    selected <- market_groups[market_groups$RV == "RV20" &
                                  market_groups$method == "AL", ]
    cross <- list()
    i <- 0L
    for (community in 1:2) for (group in 1:2) {
        i <- i + 1L
        members <- sort(selected$market[selected$group == group &
                                            selected$community == community],
                        method = "radix")
        cross[[i]] <- data.frame(community = community, group = group,
                                 n_indices = length(members),
                                 markets = paste(members, collapse = ", "),
                                 stringsAsFactors = FALSE)
    }
    cross <- do.call(rbind, cross)

    utils::write.csv(coefficients, file.path(outdir, "coefficients.csv"), row.names = FALSE)
    utils::write.csv(market_groups, file.path(outdir, "market_groups.csv"), row.names = FALSE)
    utils::write.csv(block_effects, file.path(outdir, "block_effects.csv"), row.names = FALSE)
    utils::write.csv(cross, file.path(outdir, "AL_RV20_crossclassification.csv"),
                     row.names = FALSE)
    utils::write.csv(flow, file.path(outdir, "iteration_flow.csv"), row.names = FALSE)

    # Keep the paper's display order; this does not change estimation order.
    leading <- c("us.DJI", "us.IXIC", "us.RUT", "us.SPX", "pk.KSE", "it.FTMIB",
                 "ca.GSPTSE", "ch.SSMI", "sg.STI", "br.BVSP", "jp.N225")
    market_order <- c(leading, labels[!labels %in% leading])
    groups <- groups[match(market_order, labels), , drop = FALSE]
    communities <- communities[match(market_order, labels), , drop = FALSE]
    dimnames(groups) <- dimnames(communities) <- list(market_order, column_names)
    plot_memberships(groups, "Group structure", "#00008B",
                     file.path(outdir, "groupstructure.pdf"))
    plot_memberships(communities, "Community detection", "#FF8C00",
                     file.path(outdir, "communitydetection.pdf"))
    invisible(list(coefficients = coefficients, market_groups = market_groups,
                   block_effects = block_effects, crossclassification = cross,
                   iteration_flow = flow, groups = groups, communities = communities))
}
