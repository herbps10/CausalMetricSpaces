#include <RcppArmadillo.h>
#include <numbers>
using namespace Rcpp;

const double pi = 3.14159265358979323846;

// [[Rcpp::depends(RcppArmadillo)]]

double circular_distance(double x, double y) {
  return pi - abs(pi - abs(x - y));
}

double circular_recenter(double x) {
  return (x / (2.0 * pi) - floor(x / (2.0 * pi))) * 2.0 * pi;
}

// [[Rcpp::export]]
double circular_objective_cpp(arma::vec X, double mu, arma::vec weights) {
  double x = 0;
  int N = X.size();
  for(int i = 0; i < N; i++) {
    x += weights(i) * pow(circular_distance(mu, X(i)), 2);
  }
  return x;
}

// [[Rcpp::export]]
double circular_frechet_mean_cpp(arma::vec X, arma::vec weights, bool sorted) {
  double s = sum(weights % X);
  int N = X.size();
  double lowest_obj = std::numeric_limits<double>::max();
  double lowest_mu;

  int start_index = std::floor(s * N / (2 * pi));

  for(int i = 0; i < N; i++) {
    double mu = std::fmod((s + 2 * pi * (i + 1) / N), (2 * pi)) - pi;
    if(sorted && mu <= s - 2 * pi * start_index / N) {
      continue;
    }
    double obj = circular_objective_cpp(X, mu, weights);

    if(obj < lowest_obj) {
      lowest_obj = obj;
      lowest_mu = mu;
    }
  }
  return lowest_mu;
}


// [[Rcpp::export]]
arma::mat circular_local(arma::vec y, arma::mat K0, arma::mat K1, arma::vec A) {
  int N = y.size();
  double dist = 0;
  arma::mat res = arma::mat(N, 2, arma::fill::eye);
  for(int i = 0; i < N; i++) {
    arma::vec w0 = K0.row(i).t() / arma::sum(K0.row(i));
    arma::vec w1 = K1.row(i).t() / arma::sum(K1.row(i));

    double mu0 = circular_frechet_mean_cpp(y, w0, false);
    double mu1 = circular_frechet_mean_cpp(y, w1, false);

    res(i, 0) = mu0;
    res(i, 1) = mu1;
  }
  return res;
}

// [[Rcpp::export]]
arma::mat dist_cpp(arma::mat w1, arma::mat w2) {
  int N = w1.n_rows;

  arma::mat res = arma::mat(N, N, arma::fill::zeros);
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < i; j++) {
      res(i, j) = sum(w1.row(j) - w2.row(i));
      res(j, i) = -res(i, j);
    }
  }
  return res;
}

// [[Rcpp::export]]
arma::mat circular_local_linear_cpp(arma::vec y0, arma::vec y1, arma::mat K0, arma::mat K1, arma::vec K0_means, arma::vec K1_means, arma::mat dist1_0, arma::mat dist2_0,  arma::mat dist1_1, arma::mat dist2_1, arma::vec A) {
  int N = y0.size() + y1.size();
  double dist = 0;
  arma::mat res = arma::mat(N, 2, arma::fill::zeros);
  arma::vec w0 = arma::vec(N);
  arma::vec w1 = arma::vec(N);

  for(int i = 0; i < N; i++) {
    double mu2_0 = as_scalar(K0.row(i) * dist2_0.col(i)) / (double) N;
    double mu1_0 = as_scalar(K0.row(i) * dist1_0.col(i)) / (double) N;
    //double mu0_0 = arma::mean(K0.row(i));
    double mu0_0 = K0_means(i);
    double sigmasq0 = mu0_0 * mu2_0 - pow(mu1_0, 2);

    double mu2_1 = as_scalar(K1.row(i) * dist2_1.col(i)) / (double) N;
    double mu1_1 = as_scalar(K1.row(i) * dist1_1.col(i)) / (double) N;
    //double mu0_1 = arma::mean(K1.row(i));
    double mu0_1 = K1_means(i);
    double sigmasq1 = mu0_1 * mu2_1 - pow(mu1_1, 2);

    w0 = 1 / sigmasq0 * (K0.row(i).t() % (mu2_0 - mu1_0 * dist1_0.col(i)));
    w0 = w0 / arma::sum(w0);

    w1 = 1 / sigmasq1 * (K1.row(i).t() % (mu2_1 - mu1_1 * dist1_1.col(i)));
    w1 = w1 / arma::sum(w1);

    double mu0 = circular_frechet_mean_cpp(y0, w0, true);
    double mu1 = circular_frechet_mean_cpp(y1, w1, true);

    res(i, 0) = circular_recenter(mu0);
    res(i, 1) = circular_recenter(mu1);
  }
  return res;
}


// [[Rcpp::export]]
double circular_local_linear_loo_cpp(arma::vec y, arma::mat K0, arma::mat K1, arma::mat dist1, arma::mat dist2, arma::vec A) {
  int N = y.size();
  double dist = 0;
  arma::mat res = arma::mat(N, 2, arma::fill::eye);
  arma::vec w0 = arma::vec(N);
  arma::vec w1 = arma::vec(N);

  double error = 0;
  for(int i = 0; i < N; i++) {
    if(A(i) == 0) {
      double tmp = K0(i, i);
      K0(i, i) = 0;
      double mu2_0 = arma::mean(K0.row(i) % dist2.row(i));
      double mu1_0 = arma::mean(K0.row(i) % dist1.row(i));
      double mu0_0 = arma::mean(K0.row(i));
      double sigmasq0 = mu0_0 * mu2_0 - pow(mu1_0, 2);

      w0 = 1 / sigmasq0 * (K0.row(i).t() % (mu2_0 - mu1_0 * dist1.row(i).t()));
      w0 = w0 / arma::sum(w0);

      double mu0 = circular_recenter(circular_frechet_mean_cpp(y, w0, false));
      error += circular_distance(y(i), mu0);
      K0(i, i) = tmp;
    }
    else {
      double tmp = K1(i, i);
      K1(i, i) = 0;
      double mu2_1 = arma::mean(K1.row(i) % dist2.row(i));
      double mu1_1 = arma::mean(K1.row(i) % dist1.row(i));
      double mu0_1 = arma::mean(K1.row(i));
      double sigmasq1 = mu0_1 * mu2_1 - pow(mu1_1, 2);

      w1 = 1 / sigmasq1 * (K1.row(i).t() % (mu2_1 - mu1_1 * dist1.row(i).t()));
      w1 = w1 / arma::sum(w1);
      double mu1 = circular_recenter(circular_frechet_mean_cpp(y, w1, false));
      error += circular_distance(y(i), mu1);
      K1(i, i) = tmp;
    }
  }
  return error / (double) N;
}

// [[Rcpp::export]]
double circular_local_loo(arma::vec y, arma::mat K0, arma::mat K1, arma::vec A) {
  int N = y.size();
  double error = 0;
  for(int i = 0; i < N; i++) {
    if(A(i) == 1) {
      double tmp = K1(i, i);
      K1(i, i) = 0;

      arma::vec w = K1.row(i).t() / arma::sum(K1.row(i));
      double loo_mu = circular_frechet_mean_cpp(y, w, false);

      error += circular_distance(loo_mu, y(i));

      K1(i, i) = tmp;
    }
    else {
      double tmp = K0(i, i);
      K0(i, i) = 0;

      arma::vec w = K0.row(i).t() / arma::sum(K0.row(i));
      double loo_mu = circular_frechet_mean_cpp(y, w, false);

      error += circular_distance(loo_mu, y(i));

      K0(i, i) = tmp;
    }
  }
  return error / (double)N;
}
