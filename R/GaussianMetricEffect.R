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
  cat(glue::glue("Naive estimator\n\n"))
  cat(glue::glue("Population effect: \t{round(x$naive$psi, 3)}\n\n"))
  cat(glue::glue("Control: \t\tmu = {format_mu(x$naive$Qstar0)}\n\n"))
  cat("\t\t\tSigma = ")
  cat(glue::glue("{format_Sigma(x$naive$Qstar0, 4)}\n\n"))
  cat(glue::glue("Treatment: \t\tmu = {format_mu(x$naive$Qstar1)}\n\n"))
  cat("\t\t\tSigma = ")
  cat(glue::glue("{format_Sigma(x$naive$Qstar1, 4)}\n\n"))
  cat("\n")

  cat(glue::glue("IPTW estimator\n\n"))
  cat(glue::glue("Population effect: \t{round(x$iptw$psi, 3)}\n\n"))
  cat(glue::glue("Control: \t\tmu = {format_mu(x$iptw$Qstar0)}\n\n"))
  cat("\t\t\tSigma = ")
  cat(glue::glue("{format_Sigma(x$iptw$Qstar0, 4)}\n\n"))
  cat(glue::glue("Treatment: \t\tmu = {format_mu(x$iptw$Qstar1)}\n\n"))
  cat("\t\t\tSigma = ")
  cat(glue::glue("{format_Sigma(x$iptw$Qstar1, 4)}\n\n"))
  cat("\n")

  if(!is.na(x$iptw$p_value)) {
    cat(glue::glue("Permutation P-value: \t{signif(x$iptw$p_value, 3)}\n\n"))
  }
  if(!is.null(fit$iptw$bs)) {
    ci <- quantile(fit$iptw$bs, c(0.025, 0.975))
    cat(glue::glue("Bootstrap 95% CI: \t({signif(ci[1], 3)}, {signif(ci[2], 3)})\n\n"))
  }
}
