#include <Rcpp.h>
#include <mvtnorm.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sort_cpp(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
NumericVector cdf_empirical_cpp(NumericVector x) {
  NumericVector x_sort = sort_cpp(x);
  double n = x.length();
  NumericVector rank = as<NumericVector>(match(x, x_sort));
  return rank / n;
}

// [[Rcpp::export]]
double e_prod_cpp(NumericVector x, NumericVector y, double rho, Function gauss_box) {

  NumericVector x_sort = sort_cpp(x);
  NumericVector y_sort = sort_cpp(y);

  NumericVector cdf_sorted_x = cdf_empirical_cpp(x_sort);
  NumericVector cdf_sorted_y = cdf_empirical_cpp(y_sort);

  double e_product = 0;
  double sum_p     = 0;

  // 0 based indexing
  for (int i = 1; i < 503; ++i) {
    
    double x_i   = x_sort[i];
    double px_i  = cdf_sorted_x[i];
    double px_i1 = cdf_sorted_x[i-1];
    
    for (int j = 1; j < 503; ++j) {
      
      double y_j   = y_sort[j];
      double py_j  = cdf_sorted_y[j];
      double py_j1 = cdf_sorted_y[j-1];
      
      NumericVector p_temp = gauss_box(px_i, py_j, px_i1, py_j1, rho);
      double p = as<double>(p_temp);
      
      e_product = e_product + x_i * y_j * p;
      sum_p     = sum_p + p;
      
    }
  }
  
  return e_product;
}

/*** R
#3.74271e-05
e_prod_cpp(portfolio_shifts$AAPL, portfolio_shifts$SPY, rho, gaussian_copula_mvtnorm_box)
*/
