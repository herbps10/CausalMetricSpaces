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
  cat(glue::glue("Naive estimator\n\n"))
  cat(glue::glue("Population effect: \t{round(x$naive$psi, 3)}\n\n"))
  cat(glue::glue("Control: \t\t{round(x$naive$Qstar0, 3)}\n\n"))
  cat(glue::glue("Treatment: \t{round(x$naive$Qstar1, 3)}\n\n"))

  cat("\n")
  cat(glue::glue("IPTW estimator\n\n"))
  cat(glue::glue("Population effect: \t{round(x$iptw$psi, 3)}\n\n"))
  cat(glue::glue("Control: \t\t{round(x$iptw$Qstar0, 3)}\n\n"))
  cat(glue::glue("Treatment: \t{round(x$iptw$Qstar1, 3)}\n\n"))
  if(!is.na(x$iptw$p_value)) {
    cat(glue::glue("Permutation P-value: \t{signif(x$iptw$p_value, 3)}\n\n"))
  }
  if(!is.null(fit$iptw$bs)) {
    ci <- quantile(fit$iptw$bs, c(0.025, 0.975))
    cat(glue::glue("Bootstrap 95% CI: \t({signif(ci[1], 3)}, {signif(ci[2], 3)})\n\n"))
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
