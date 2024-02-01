#'
#' @import tibble dplyr purrr
#' @importFrom circular rvonmises circular
#' @export
simulate_circular <- function(
    N = 1e2,
    p = 2,
    e = function(W) 0.5,

    m = function(W) 0,
    tau = function(W) 0,
    sigma = 0.1,
    seed = NA
) {
  if(!is.na(seed)) {
    set.seed(seed)
  }

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

#'
#' @import ggplot2 dplyr
#' @export
plot_circular_simulated_data <- function(data) {
  breaks <- c(0, pi / 2, pi, 3 * pi / 2)
  labels <- c(expression(0), expression(frac(pi, 2)), expression(pi), expression(frac(3*pi, 2)))

  data %>%
    mutate(radius = ifelse(A == 0, 1, 1.2),
           A = factor(A)) %>%
    ggplot(aes(x = Y, y = radius)) +
    geom_point(aes(color = A)) +
    coord_polar(theta = "x") +
    scale_y_continuous(limits = c(0, 1.2)) +
    scale_x_continuous(limits = c(0, 2 * pi), breaks = breaks, labels = labels) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
}
