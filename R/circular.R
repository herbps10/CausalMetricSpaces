objective_circular <- Vectorize(function(X, mu, weights = rep(1 / length(X), length(X))) {
  sum(weights * distance_circular(rep(mu, length(X)), X)^2)
}, vectorize.args = "mu")

frechet_mean_circular <- function(X, weights = rep(1 / length(X), length(X))) {
  mu <- (sum(weights * X) + 2*pi*(1:length(X)) / length(X)) %% (2 * pi) - pi
  mu[which.min(objective_circular(X, mu, weights = weights))]
}

distance_circular <- function(x, y) {
  pi - abs(pi - abs(x - y))
}

circular_recenter <- \(x) (x / (2 * pi) - floor(x / (2 * pi))) * 2 * pi

#' @importFrom stats glm predict
#' @importFrom progress progress_bar
#' @export
population_effect_circular <- function(data, Y, A, W, fit_glm = TRUE, permutation_test = TRUE, permutation_iter = 1e3, bootstrap_se = TRUE, bootstrap_iter = 1e3, progress = TRUE) {
  # Naive estimator
  Qstar0_naive <- frechet_mean_circular(data[[Y]][data[[A]] == 0])
  Qstar1_naive <- frechet_mean_circular(data[[Y]][data[[A]] == 1])
  naive_bw <- distance_circular(Qstar0_naive, Qstar1_naive)

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
  w1 <- unname(1 / gn[data[[A]] == 1])
  w1 <- w1 / sum(w1)

  Qstar0_iptw <- frechet_mean_circular(data[[Y]][data[[A]] == 0], weights = w0)
  Qstar1_iptw <- frechet_mean_circular(data[[Y]][data[[A]] == 1], weights = w1)

  iptw_bw <- distance_circular(Qstar0_iptw, Qstar1_iptw)

  ps <- NULL
  p_value <- NA
  if(permutation_test) {
    if(progress) pb <- progress_bar$new(total = permutation_iter, format = "Permutation test [:bar] :percent")
    ps <- numeric(permutation_iter)
    for(index in 1:permutation_iter) {
      As <- data[[A]][sample(1:nrow(data), nrow(data), replace = FALSE)]
      data2 <- data
      data2[[A]] <- As
      ps[index] <- population_effect_circular(data2, A = A, Y = Y, W = W, bootstrap_se = FALSE, permutation_test = FALSE, progress = FALSE)$iptw$psi
      if(progress) pb$tick()
    }
    p_value <- mean(ps > iptw_bw)
  }

  bs <- NULL
  if(bootstrap_se) {
    if(progress) pb <- progress_bar$new(total = bootstrap_iter, format = "Bootstrap [:bar] :percent")
    bs <- numeric(bootstrap_iter)
    for(index in 1:bootstrap_iter) {
      indices <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
      bs[index] <- population_effect_circular(data[indices,], A = A, Y = Y, W = W, bootstrap_se = FALSE, permutation_test = FALSE, progress = FALSE)$iptw$psi
      if(progress) pb$tick()
    }
  }

  res <- list(
    data = data,
    Y = Y,
    A = A,
    W = W,
    naive = list(psi = naive_bw, Qstar0 = Qstar0_naive, Qstar1 = Qstar1_naive),
    iptw = list(psi = iptw_bw, Qstar0 = Qstar0_iptw, Qstar1 = Qstar1_iptw, bs = bs, ps = ps, gn = gn, p_value = p_value)
  )
  class(res) <- "CircularMetricEffect"
  res
}
