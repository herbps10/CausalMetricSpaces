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
