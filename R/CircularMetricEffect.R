#' @importFrom glue glue
#' @export
print.CircularMetricEffect <- function(x) {
  cat("Circular Metric Effect\n")
  cat(glue::glue("Naive estimator: \t{round(x$naive$psi, 3)}"))
  cat("\n")
  cat(glue::glue("IPTW estimator: \t{round(x$iptw$psi, 3)}"))
}

#' @importFrom glue glue
#' @export
summary.CircularMetricEffect <- function(x) {
  res <- list(
    naive = list(
      psi = x$naive$psi,
      Qstar0 = x$naive$Qstar0,
      Qstar1 = x$naive$Qstar1
    ),
    iptw = list(
      psi = x$iptw$psi,
      Qstar0 = x$iptw$Qstar0,
      Qstar1 = x$iptw$Qstar1
    )
  )

  if(!is.na(x$iptw$p_value)) {

  }

  if(!is.null(x$iptw$bs)) {
    ci <- quantile(x$iptw$bs, c(0.025, 0.975))
  }
}

#' @importFrom glue glue
#' @export
print.summaryCircularMetricEffect <- function(x, digits = 3) {
  cat(glue::glue("Naive estimator\n\n"))
  cat(glue::glue("Population effect: \t{format(x$naive$psi, digits = digits)}\n\n"))
  cat(glue::glue("Control: \t\t{format(x$naive$Qstar0, digits = digits)}\n\n"))
  cat(glue::glue("Treatment: \t{format(x$naive$Qstar1, digits = digits)}\n\n"))

  cat("\n")
  cat(glue::glue("IPTW estimator\n\n"))
  cat(glue::glue("Population effect: \t{format(x$iptw$psi, digits = digits)}\n\n"))
  cat(glue::glue("Control: \t\t{format(x$iptw$Qstar0, digits = digits)}\n\n"))
  cat(glue::glue("Treatment: \t{format(x$iptw$Qstar1, digits = digits)}\n\n"))
  if(!is.na(x$iptw$p_value)) {
    cat(glue::glue("Permutation P-value: \t{format(x$iptw$p_value, digits = digits)}\n\n"))
  }
  if(!is.null(x$iptw$bs)) {
    ci <- quantile(x$iptw$bs, c(0.025, 0.975))
    cat(glue::glue("Bootstrap 95% CI: \t({format(ci[1], digits = digits)}, {format(ci[2], digits = digits)})\n\n"))
  }
}

#' @import ggplot2 dplyr tibble
#' @export
plot.CircularMetricEffect <- function(x) {
  breaks <- c(0, pi / 2, pi, 3 * pi / 2)
  labels <- c(expression(0), expression(frac(pi, 2)), expression(pi), expression(frac(3*pi, 2)))

  frechet_means <- tibble(
    A = factor(c(0, 1)),
    radius = c(1.1, 1.3),
    Y = c(x$iptw$Qstar0, x$iptw$Qstar1)
  )

  x$data %>%
    mutate(radius = ifelse(x$data[[x$A]] == 0, 1, 1.2),
           A = factor(x$data[[x$A]])) %>%
    ggplot(aes(x = Y, y = radius)) +
    geom_point(aes(color = A)) +
    geom_text(label = "\u25B2", data = frechet_means, aes(color = A, angle = -Y * 180 / pi)) +
    coord_polar(theta = "x") +
    scale_y_continuous(limits = c(0, 1.3)) +
    scale_x_continuous(limits = c(0, 2 * pi), breaks = breaks, labels = labels) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
}
