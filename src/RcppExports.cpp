// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// circular_objective_cpp
double circular_objective_cpp(arma::vec X, double mu, arma::vec weights);
RcppExport SEXP _CausalMetricSpaces_circular_objective_cpp(SEXP XSEXP, SEXP muSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(circular_objective_cpp(X, mu, weights));
    return rcpp_result_gen;
END_RCPP
}
// circular_frechet_mean_cpp
double circular_frechet_mean_cpp(arma::vec X, arma::vec weights, bool sorted);
RcppExport SEXP _CausalMetricSpaces_circular_frechet_mean_cpp(SEXP XSEXP, SEXP weightsSEXP, SEXP sortedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type sorted(sortedSEXP);
    rcpp_result_gen = Rcpp::wrap(circular_frechet_mean_cpp(X, weights, sorted));
    return rcpp_result_gen;
END_RCPP
}
// circular_local
arma::mat circular_local(arma::vec y, arma::mat K0, arma::mat K1, arma::vec A);
RcppExport SEXP _CausalMetricSpaces_circular_local(SEXP ySEXP, SEXP K0SEXP, SEXP K1SEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(circular_local(y, K0, K1, A));
    return rcpp_result_gen;
END_RCPP
}
// dist_cpp
arma::mat dist_cpp(arma::mat w1, arma::mat w2);
RcppExport SEXP _CausalMetricSpaces_dist_cpp(SEXP w1SEXP, SEXP w2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w2(w2SEXP);
    rcpp_result_gen = Rcpp::wrap(dist_cpp(w1, w2));
    return rcpp_result_gen;
END_RCPP
}
// circular_local_linear_cpp
arma::mat circular_local_linear_cpp(arma::vec y0, arma::vec y1, arma::mat K0, arma::mat K1, arma::vec K0_means, arma::vec K1_means, arma::mat dist1_0, arma::mat dist2_0, arma::mat dist1_1, arma::mat dist2_1, arma::vec A);
RcppExport SEXP _CausalMetricSpaces_circular_local_linear_cpp(SEXP y0SEXP, SEXP y1SEXP, SEXP K0SEXP, SEXP K1SEXP, SEXP K0_meansSEXP, SEXP K1_meansSEXP, SEXP dist1_0SEXP, SEXP dist2_0SEXP, SEXP dist1_1SEXP, SEXP dist2_1SEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K0_means(K0_meansSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K1_means(K1_meansSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist1_0(dist1_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist2_0(dist2_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist1_1(dist1_1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist2_1(dist2_1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(circular_local_linear_cpp(y0, y1, K0, K1, K0_means, K1_means, dist1_0, dist2_0, dist1_1, dist2_1, A));
    return rcpp_result_gen;
END_RCPP
}
// circular_local_linear_loo_cpp
double circular_local_linear_loo_cpp(arma::vec y, arma::mat K0, arma::mat K1, arma::mat dist1, arma::mat dist2, arma::vec A);
RcppExport SEXP _CausalMetricSpaces_circular_local_linear_loo_cpp(SEXP ySEXP, SEXP K0SEXP, SEXP K1SEXP, SEXP dist1SEXP, SEXP dist2SEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist1(dist1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist2(dist2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(circular_local_linear_loo_cpp(y, K0, K1, dist1, dist2, A));
    return rcpp_result_gen;
END_RCPP
}
// circular_local_loo
double circular_local_loo(arma::vec y, arma::mat K0, arma::mat K1, arma::vec A);
RcppExport SEXP _CausalMetricSpaces_circular_local_loo(SEXP ySEXP, SEXP K0SEXP, SEXP K1SEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(circular_local_loo(y, K0, K1, A));
    return rcpp_result_gen;
END_RCPP
}
// distance_gaussian_cpp
double distance_gaussian_cpp(arma::mat Sigma0, arma::mat Sigma1);
RcppExport SEXP _CausalMetricSpaces_distance_gaussian_cpp(SEXP Sigma0SEXP, SEXP Sigma1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Sigma0(Sigma0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma1(Sigma1SEXP);
    rcpp_result_gen = Rcpp::wrap(distance_gaussian_cpp(Sigma0, Sigma1));
    return rcpp_result_gen;
END_RCPP
}
// frechet_mean_gaussian_cpp
arma::mat frechet_mean_gaussian_cpp(std::vector<arma::mat> Sigma, arma::vec weights, int max_iter, double tol);
RcppExport SEXP _CausalMetricSpaces_frechet_mean_gaussian_cpp(SEXP SigmaSEXP, SEXP weightsSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(frechet_mean_gaussian_cpp(Sigma, weights, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}
// bw_local
double bw_local(std::vector<arma::vec> mu, std::vector<arma::mat> Sigma, arma::mat K0, arma::mat K1, arma::vec A, int max_iter, double tol);
RcppExport SEXP _CausalMetricSpaces_bw_local(SEXP muSEXP, SEXP SigmaSEXP, SEXP K0SEXP, SEXP K1SEXP, SEXP ASEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::vec> >::type mu(muSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(bw_local(mu, Sigma, K0, K1, A, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}
// bw_local_loo
double bw_local_loo(std::vector<arma::vec> mu, std::vector<arma::mat> Sigma, arma::mat K0, arma::mat K1, arma::vec A, int max_iter, double tol);
RcppExport SEXP _CausalMetricSpaces_bw_local_loo(SEXP muSEXP, SEXP SigmaSEXP, SEXP K0SEXP, SEXP K1SEXP, SEXP ASEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::vec> >::type mu(muSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(bw_local_loo(mu, Sigma, K0, K1, A, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CausalMetricSpaces_circular_objective_cpp", (DL_FUNC) &_CausalMetricSpaces_circular_objective_cpp, 3},
    {"_CausalMetricSpaces_circular_frechet_mean_cpp", (DL_FUNC) &_CausalMetricSpaces_circular_frechet_mean_cpp, 3},
    {"_CausalMetricSpaces_circular_local", (DL_FUNC) &_CausalMetricSpaces_circular_local, 4},
    {"_CausalMetricSpaces_dist_cpp", (DL_FUNC) &_CausalMetricSpaces_dist_cpp, 2},
    {"_CausalMetricSpaces_circular_local_linear_cpp", (DL_FUNC) &_CausalMetricSpaces_circular_local_linear_cpp, 11},
    {"_CausalMetricSpaces_circular_local_linear_loo_cpp", (DL_FUNC) &_CausalMetricSpaces_circular_local_linear_loo_cpp, 6},
    {"_CausalMetricSpaces_circular_local_loo", (DL_FUNC) &_CausalMetricSpaces_circular_local_loo, 4},
    {"_CausalMetricSpaces_distance_gaussian_cpp", (DL_FUNC) &_CausalMetricSpaces_distance_gaussian_cpp, 2},
    {"_CausalMetricSpaces_frechet_mean_gaussian_cpp", (DL_FUNC) &_CausalMetricSpaces_frechet_mean_gaussian_cpp, 4},
    {"_CausalMetricSpaces_bw_local", (DL_FUNC) &_CausalMetricSpaces_bw_local, 7},
    {"_CausalMetricSpaces_bw_local_loo", (DL_FUNC) &_CausalMetricSpaces_bw_local_loo, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_CausalMetricSpaces(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
