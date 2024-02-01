#' @importFrom glue glue
#' @export
print.GaussianMetricEffect <- function(x) {
  cat("Gaussian Metric Effect\n")
  cat(glue::glue("Naive estimator: \t{round(x$naive$psi, 3)}"))
  cat("\n")
  cat(glue::glue("IPTW estimator: \t{round(x$iptw$psi, 3)}"))
}

format_mu <- function(x, indent = 0, digits = 3) {
  tabs <- paste0(rep("\t", indent), collapse = "")
  glue::glue("{tabs}({paste0(format(x$mu, digits = digits), collapse = ', ')})")
}

format_Sigma <- function(x, indent = 0, digits = 3) {
  tabs <- paste0(rep("\t", indent), collapse = "")
  glue::glue("({paste0(format(x$Sigma, digits = digits), collapse = ', ')})")
}

#' @importFrom glue glue
#' @export
summary.GaussianMetricEffect <- function(x) {
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
    res$iptw$p_value <- x$iptw$p_value
  }

  if(!is.null(x$iptw$bs)) {
    ci <- quantile(x$iptw$bs, c(0.025, 0.975))
    res$iptw$ci <- ci
  }
  class(res) <- "summaryGaussianMetricEffect"
  res
}

#' @importFrom glue glue
#' @export
print.summaryGaussianMetricEffect <- function(x, digits = 3) {
  cat(glue::glue("Naive estimator\n\n"))
  cat(glue::glue("Population effect: \t{format(x$naive$psi, digits = digits)}\n\n"))
  cat(glue::glue("Control: \t\tmu = {format_mu(x$naive$Qstar0, digits = digits)}\n\n"))
  cat("\t\t\tSigma = ")
  cat(glue::glue("{format_Sigma(x$naive$Qstar0, 4, digits = digits)}\n\n"))
  cat(glue::glue("Treatment: \t\tmu = {format_mu(x$naive$Qstar1, digits = digits)}\n\n"))
  cat("\t\t\tSigma = ")
  cat(glue::glue("{format_Sigma(x$naive$Qstar1, 4, digits = digits)}\n\n"))
  cat("\n")

  cat(glue::glue("IPTW estimator\n\n"))
  cat(glue::glue("Population effect: \t{format(x$iptw$psi, digits = digits)}\n\n"))
  cat(glue::glue("Control: \t\tmu = {format_mu(x$iptw$Qstar0, digits = digits)}\n\n"))
  cat("\t\t\tSigma = ")
  cat(glue::glue("{format_Sigma(x$iptw$Qstar0, 4, digits = digits)}\n\n"))
  cat(glue::glue("Treatment: \t\tmu = {format_mu(x$iptw$Qstar1, digits = digits)}\n\n"))
  cat("\t\t\tSigma = ")
  cat(glue::glue("{format_Sigma(x$iptw$Qstar1, 4, digits = digits)}\n\n"))
  cat("\n")

  if(!is.na(x$iptw$p_value)) {
    cat(glue::glue("Permutation P-value: \t{format(x$iptw$p_value, digits = digits)}\n\n"))
  }
  if(!is.null(x$iptw$ci)) {
    cat(glue::glue("Bootstrap 95% CI: \t({format(x$iptw$ci[1], digits = digits)}, {format(x$iptw$ci[2], digits = digits)})\n\n"))
  }
}

#' @import ggplot2
#' @export
plot.GaussianMetricEffect <- function(x) {
  data <- x$data

  get_ellipse <- function(A, mu, Sigma) {
    ellipse::ellipse(Sigma) %>%
      as.data.frame() %>%
      mutate(x = x + mu[1],
             y = y + mu[2],
             A = A)
  }

  effects <- pmap(list(
      list(0, 1),
      list(x$iptw$Qstar0$mu, x$iptw$Qstar1$mu),
      list(x$iptw$Qstar0$Sigma, x$iptw$Qstar1$Sigma)
    ),
    get_ellipse
  ) %>%
    bind_rows(.id = "index")

  pmap(list(data[[x$A]], data[[x$mu]], data[[x$Sigma]]), get_ellipse) %>%
    bind_rows(.id = "index") %>%
    ggplot(aes(x = x, y = y, group = index, color = factor(A))) +
    geom_path(alpha = 0.25) +
    geom_path(data = effects, linewidth = 1) +
    labs(x = expression(Y[1]), y = expression(Y[2])) +
    guides(color = guide_legend(title = "A"))
}
