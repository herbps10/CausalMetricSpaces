perturb_gaussian <- function(frechet_mean, s = 0.1) {
  x <- diag(runif(2, 1 - s, 1 + s))
  x %*% frechet_mean %*% x
}

#'
#' @import tibble dplyr purrr
#' @export
simulate_gaussian <- function(
    N = 1e2,
    p = 2,
    e = function(W) 0.5,

    m_mu = function(W) 0,
    m_sigma = function(W) 0.5,
    m_rho = function(W) 0,

    tau_mu = function(W) 0,
    tau_sigma = function(W) 0,
    tau_rho = function(W) 0,

    sigma_perturb = 0.1,
    mu_perturb = 0.1,
    seed = NA
) {
  if(!is.na(seed)) {
    set.seed(seed)
  }

  W <- matrix(runif(N * p), ncol = p)
  colnames(W) <- paste0("W", 1:p)

  A <- rbinom(N, 1, e(W))

  rho0 <- map_dbl(1:N, function(index) m_rho(W[index, ]) - 0 * tau_rho(W[index, ]))
  rho1 <- map_dbl(1:N, function(index) m_rho(W[index, ]) + 1 * tau_rho(W[index, ]))

  # Generate conditional barycenters
  mu0_star <- map(1:N, function(index) m_mu(W[index, ]) - 0 * tau_mu(W[index, ]))
  mu1_star <- map(1:N, function(index) m_mu(W[index, ]) + 1 * tau_mu(W[index, ]))

  sigma0_star <- map(1:N, function(index) m_sigma(W[index, ]) - 0 * tau_sigma(W[index, ]))
  sigma1_star <- map(1:N, function(index) m_sigma(W[index, ]) + 1 * tau_sigma(W[index, ]))

  Sigma0_star <- map2(rho0, sigma0_star, function(rho, sigma) {
    diag(sigma^2, 2) %*% matrix(c(1, rho, rho, 1), ncol = 2) %*% diag(sigma^2, 2)
  })

  Sigma1_star <- map2(rho1, sigma1_star, function(rho, sigma) {
    diag(sigma^2, 2) %*% matrix(c(1, rho, rho, 1), ncol = 2) %*% diag(sigma^2, 2)
  })

  Sigma0 <- map(Sigma0_star, perturb_gaussian, s = sigma_perturb)
  Sigma1 <- map(Sigma1_star, perturb_gaussian, s = sigma_perturb)
  mu0 <- map(mu0_star, function(mu) mu + rnorm(2, 0, mu_perturb))
  mu1 <- map(mu1_star, function(mu) mu + rnorm(2, 0, mu_perturb))

  Sigma <- pmap(list(Sigma0, Sigma1, A), function(Sigma0, Sigma1, A) {
    if(A == 0) return(Sigma0)
    return(Sigma1)
  })

  mu <- pmap(list(mu0, mu1, A), function(mu0, mu1, A) {
    if(A == 0) return(mu0)
    return(mu1)
  })

  as_tibble(W) %>%
    mutate(A = A,
           e = e(W),
           mu0_star = mu0_star,
           mu1_star = mu1_star,
           mu0 = mu0,
           mu1 = mu1,
           sigma0_star = sigma0_star,
           sigma1_star = sigma1_star,
           Sigma0_star = Sigma0_star,
           Sigma1_star = Sigma1_star,
           Sigma0 = Sigma0,
           Sigma1 = Sigma1,
           Sigma = Sigma,
           mu = mu)
}

plot_gaussian_simulated_data <- function(data) {
  pmap(list(data$A, data$mu, data$Sigma), function(A, mu, Sigma) {
    ellipse::ellipse(Sigma) %>%
      as.data.frame() %>%
      mutate(x = x + mu[1],
             y = y + mu[2],
             A = A)
  }) %>%
    bind_rows(.id = "index") %>%
    ggplot(aes(x = x, y = y, group = index, color = factor(A))) +
    geom_path() +
    labs(x = expression(Y[1]), y = expression(Y[2]))
}
