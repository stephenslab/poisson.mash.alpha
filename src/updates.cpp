#include <RcppArmadillo.h>

using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void update_rho (const vec& x, const mat& Fuv, double s, const  vec& mu,
		 const vec& l, vec& rho, unsigned int maxiter, double tol, 
		 double maxrho);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat & update_rho_rcpp (const arma::mat& X, const arma::mat& Fuv,
			     const arma::vec& s, const arma::mat& mu,
			     const arma::mat& L, const arma::mat& rho,
			     unsigned int maxiter, double tol, double maxrho) {
  unsigned int D = rho.n_rows;
  unsigned int R = rho.n_cols;
  mat rhonew(D,R);
  vec t(D);
  for (unsigned int r = 0; r < R; r++) {
    t = rho.col(r);
    update_rho(X.col(r),Fuv,s(r),mu.col(r),L.col(r),t,maxiter,tol,maxrho);
    rhonew.col(r) = t;
  }
  return rhonew;
}

// TO DO: Explain here what this function does, and how to use it.
void update_rho (const vec& x, const mat& Fuv, double s, const  vec& mu,
		 const vec& l, vec& rho, unsigned int maxiter, double tol, 
		 double maxrho) {

}
