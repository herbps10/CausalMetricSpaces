sqrtm2 <- function(X) {
  A <- X[1, 1]
  B <- X[1, 2]
  C <- X[2, 1]
  D <- X[2, 2]

  tau <- A + D
  delta <- A*D - B*C
  s <- sqrt(delta)
  t <- sqrt(tau + 2 * s)
  matrix(1 / t * c(
    A + s, B,
    C, D + s
  ), 2, 2, byrow = TRUE)
}

tr <- function(x) sum(diag(x))

frechet_mean_gaussian <- function(mu, Sigma, weights = rep(1 / length(Sigma), length(Sigma))) {
  list(
    mu = Reduce(`+`, map2(mu, weights, function(mu, w) w * mu)),
    Sigma = frechet_mean_gaussian_cpp(Sigma, weights, 100, 1e-5)
  )
}

distance_gaussian <- function(Q, S) {
  sqrtS <- sqrtm2(S$Sigma)
  sum((Q$mu - S$mu)^2) + tr(Q$Sigma) + tr(S$Sigma) - 2 * tr(sqrtm2(sqrtS %*% Q$Sigma %*% sqrtS))
}

population_effect_gaussian <- function(data, mu, Sigma, A, W, fit_glm = TRUE, permutation_test = TRUE, permutation_iter = 1e3, bootstrap_se = TRUE, bootstrap_iter = 1e3, progress = TRUE) {
  # Naive estimator
  Qstar0_naive <- frechet_mean_gaussian(data[[mu]][data[[A]] == 0], data[[Sigma]][data[[A]] == 0])
  Qstar1_naive <- frechet_mean_gaussian(data[[mu]][data[[A]] == 1], data[[Sigma]][data[[A]] == 1])
  naive_effect <- distance_gaussian(Qstar0_naive, Qstar1_naive)

  if(fit_glm == TRUE) {
    g_formula <- as.formula(paste0(A, " ~ ", paste0(W, collapse = " + ")))
    fit <- glm(g_formula, family = binomial(link = "logit"), data = data)
    gn <- predict(fit, type = "response")
  }
  else {
    gn <- rep(0.5, nrow(data))
  }
  w0 <- unname(1 / (1 - gn[data[[A]] == 0]))
  w0 <- w0 / sum(w0)
  w1 <- unname(1 / gn[data$A == 1])
  w1 <- w1 / sum(w1)

  Qstar0_iptw <- frechet_mean_gaussian(data$mu[data[[A]] == 0], data[[Sigma]][data[[A]] == 0], weights = w0)
  Qstar1_iptw <- frechet_mean_gaussian(data$mu[data[[A]] == 1], data[[Sigma]][data[[A]] == 1], weights = w1)

  iptw_effect <- distance_gaussian(Qstar0_iptw, Qstar1_iptw)

  ps <- NULL
  p_value <- NA
  if(permutation_test) {
    if(progress) pb <- progress_bar$new(total = permutation_iter, format = "Permutation test [:bar] :percent")
    ps <- numeric(permutation_iter)
    for(index in 1:permutation_iter) {
      As <- data[[A]][sample(1:nrow(data), nrow(data), replace = FALSE)]
      data2 <- data
      data2[[A]] <- As
      ps[index] <- population_effect_gaussian(data2, mu, Sigma, A, W, bootstrap_se = FALSE, permutation_test = FALSE, progress = FALSE)$iptw$psi
      if(progress) pb$tick()
    }
    p_value <- mean(ps > iptw_effect)
  }

  bs <- NULL
  if(bootstrap_se) {
    if(progress) pb <- progress_bar$new(total = bootstrap_iter, format = "Bootstrap [:bar] :percent")
    bs <- numeric(bootstrap_iter)
    for(index in 1:bootstrap_iter) {
      indices <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
      bs[index] <- population_effect_gaussian(data[indices,], mu, Sigma, A, W, bootstrap_se = FALSE, permutation_test = FALSE, progress = FALSE)$iptw$psi
      if(progress) pb$tick()
    }
  }

  res <- list(
    naive = list(psi = naive_effect, Qstar0 = Qstar0_naive, Qstar1 = Qstar1_naive),
    iptw = list(psi = iptw_effect, Qstar0 = Qstar0_iptw, Qstar1 = Qstar1_iptw, bs = bs, ps = ps, gn = gn, p_value = p_value)
  )

  class(res) <- "GaussianMetricEffect"

  res
}

