#include "misc.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
// Add function declarations here.

// FUNCTION DEFINITIONS
// --------------------
// This implements update_q_theta with version = "Rcpp".
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List update_q_theta_rcpp (const arma::vec& x, const arma::vec& s, 
			  const arma::vec& mu, const arma::vec& bias, 
			  const arma::vec& c2, const arma::vec& psi2, 
			  double w, const arma::mat& U, 
			  const arma::vec& m, const arma::mat& V, 
			  unsigned int maxiter, double tol, double lwr, 
			  double upr) {
  
}

void update_q_theta (const vec& x, const vec& s, const vec& mu, 
		     const vec& bias, const vec& c2, const vec& psi2, 
		     double w, const mat& U, vec& m, mat& V, vec& a, 
		     unsigned int maxiter, double tol, double lwr, 
		     double upr) {
  mat mnew(m);
  mat Vnew(V);
  compute_poisson_rates(s,mu,bias,m,V.diag(),a);
  for (int iter = 0; iter < maxiter; iter++) {
    
  }
}
