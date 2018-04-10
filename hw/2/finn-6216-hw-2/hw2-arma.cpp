#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
arma::colvec put_option_pricer_arma(arma::colvec s, double k, double r, double y, double t, double sigma) {
  
  arma::colvec d1 = (log(s / k) + (r - y + sigma * sigma / 2.0) * t) / (sigma * sqrt(t));
  arma::colvec d2 = d1 - sigma * sqrt(t);
  
  arma::colvec V = arma::normcdf(-d2) * k * exp(-r * t) - s * exp(-y * t) % arma::normcdf(-d1);
  return as<NumericVector>(wrap(V));
}

// [[Rcpp::export]]
NumericVector put_option_pricer_rcpp(NumericVector s, double k, double r, double y, double t, double sigma) {
  NumericVector d1 = (log(s / k) + (r - y + sigma * sigma / 2.0) * t) / (sigma * sqrt(t));
  NumericVector d2 = d1 - sigma * sqrt(t);
  
  NumericVector V = Rcpp::pnorm(-d2) * k * exp(-r * t) - s * exp(-y * t) * Rcpp::pnorm(-d1);
  return V;
}




/*** R
s <- 50:55
put_option_pricer_arma(s, 60, .01, .02, 1, .05)
put_option_pricer_rcpp(s, 60, .01, .02, 1, .05)
*/
