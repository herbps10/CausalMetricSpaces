#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat sqrtm(arma::mat X) {
  double A = X(0, 0);
  double B = X(0, 1);
  double C = X(1, 0);
  double D = X(1, 1);

  double tau = A + D;
  double delta = A * D - B * C;
  double s = sqrt(delta);
  double t = sqrt(tau + 2 * s);

  arma::mat::fixed<2, 2> M;

  M(0, 0) = (A + s) / t;
  M(0, 1) = B / t;
  M(1, 0) = C / t;
  M(1, 1) = (D + s) / t;
  return std::move(M);
}

// [[Rcpp::export]]
double distance_gaussian_cpp(arma::mat Sigma0, arma::mat Sigma1) {
  arma::mat sqrtSigma0 = sqrtm(Sigma0);
  return trace(Sigma0) + trace(Sigma1) - 2 * trace(sqrtm(sqrtSigma0 * Sigma1* sqrtSigma0));
}

// [[Rcpp::export]]
arma::mat frechet_mean_gaussian_cpp(std::vector<arma::mat> Sigma, arma::vec weights, int max_iter, double tol) {
  arma::mat b = arma::mat(2, 2, arma::fill::eye);
  arma::mat sqrtb = arma::mat(2, 2);
  arma::mat bnew = arma::mat(2, 2);
  int N = Sigma.size();

  for(int i = 0; i < max_iter; i++) {
    sqrtb = sqrtm(b);
    bnew.fill(0);
    for(int n = 0; n < N; n++) {
      bnew += weights[n] * sqrtm(sqrtb * Sigma[n] * sqrtb);
    }
    if(arma::abs(b - bnew).max() < tol) {
      return bnew;
    }
    else {
      b = bnew;
    }
  }

  return b;
}

arma::vec frechet_mean_gaussian_mu_cpp(std::vector<arma::vec> mu, arma::vec weights) {
  int N = mu.size();
  arma::vec res(mu[0].size());
  for(int i = 0; i < N; i++) {
    res += mu[i] * weights[i];
  }
  return res;
}

// [[Rcpp::export]]
double bw_local(std::vector<arma::vec> mu, std::vector<arma::mat> Sigma, arma::mat K0, arma::mat K1, arma::vec A, int max_iter, double tol) {
  int N = mu.size();
  double dist = 0;
  for(int i = 0; i < N; i++) {
    arma::vec w0 = K0.row(i).t() / arma::sum(K0.row(i));
    arma::vec w1 = K1.row(i).t() / arma::sum(K1.row(i));

    arma::vec loo_mu0 = frechet_mean_gaussian_mu_cpp(mu, w0);
    arma::mat loo_Sigma0 = frechet_mean_gaussian_cpp(Sigma, w0, max_iter, tol);

    arma::vec loo_mu1 = frechet_mean_gaussian_mu_cpp(mu, w1);
    arma::mat loo_Sigma1 = frechet_mean_gaussian_cpp(Sigma, w1, max_iter, tol);

    dist += 1/(double)N * (arma::sum(arma::pow(loo_mu0 - loo_mu1, 2)) + distance_gaussian_cpp(loo_Sigma0, loo_Sigma1));
  }
  return dist;
}

// [[Rcpp::export]]
double bw_local_loo(std::vector<arma::vec> mu, std::vector<arma::mat> Sigma, arma::mat K0, arma::mat K1, arma::vec A, int max_iter, double tol) {
  int N = mu.size();
  double error = 0;
  for(int i = 0; i < N; i++) {
    if(A[i] == 1) {
      double tmp = K1(i, i);
      K1(i, i) = 0;

      arma::vec w = K1.row(i).t() / arma::sum(K1.row(i));
      arma::vec loo_mu = frechet_mean_gaussian_mu_cpp(mu, w);
      arma::mat loo_Sigma = frechet_mean_gaussian_cpp(Sigma, w, max_iter, tol);

      error += 1/(double)N * (arma::sum(arma::pow(loo_mu - mu[i], 2)) + distance_gaussian_cpp(loo_Sigma, Sigma[i]));

      K1(i, i) = tmp;
    }
    else {
      double tmp = K0(i, i);
      K0(i, i) = 0;

      arma::vec w = K0.row(i).t() / arma::sum(K0.row(i));
      arma::vec loo_mu = frechet_mean_gaussian_mu_cpp(mu, w);
      arma::mat loo_Sigma = frechet_mean_gaussian_cpp(Sigma, w, max_iter, tol);

      error += 1/(double)N * (arma::sum(arma::pow(loo_mu - mu[i], 2)) + distance_gaussian_cpp(loo_Sigma, Sigma[i]));

      K0(i, i) = tmp;
    }
  }
  return error;
}
