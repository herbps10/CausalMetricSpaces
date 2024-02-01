#'
#' @import tibble dplyr purrr
#' @importFrom circular rvonmises circular
#' @export
simulate_circular <- function(
    seed,
    N = 1e2,
    p = 2,
    e = function(W) 0.5,

    m = function(W) 0,
    tau = function(W) 0,
    sigma = 0.1
) {
  set.seed(seed)

  W <- matrix(runif(N * p), ncol = p)
  colnames(W) <- paste0("W", 1:p)

  A <- rbinom(N, 1, e(W))

  # Generate conditional barycenters
  mu0 <- circular_recenter(m(W))
  mu1 <- circular_recenter(m(W) + tau(W))

  sigma <- rep(sigma, N)

  Y0 <- map2_dbl(mu0, sigma, \(mu, sigma) circular::rvonmises(1, mu = circular(mu), kappa = 1 / sigma))
  Y1 <- map2_dbl(mu1, sigma, \(mu, sigma) circular::rvonmises(1, mu = circular(mu), kappa = 1 / sigma))

  Y <- ifelse(A == 1, Y1, Y0)

  as_tibble(W) %>%
    mutate(A = A,
           sigma = sigma,
           e = e(W),
           mu0 = mu0,
           mu1 = mu1,
           Y0 = Y0,
           Y1 = Y1,
           Y = Y)
}
