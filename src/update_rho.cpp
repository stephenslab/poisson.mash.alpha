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
arma::mat update_rho_all_rcpp (const arma::mat& X, const arma::mat& Fuv,
			       const arma::vec& s, const arma::mat& mu,
			       const arma::mat& L, const arma::mat& rho,
			       unsigned int maxiter, double tol, 
			       double maxrho) {
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

// Update the d-vector rho for a given condition r.
void update_rho (const vec& x, const mat& Fuv, double s, const  vec& mu,
		 const vec& l, vec& rho, unsigned int maxiter, double tol, 
		 double maxrho) {
  unsigned int J = Fuv.n_rows;
  unsigned int D = Fuv.n_cols;
  vec rhonew(D);
  vec bias(J);
  vec d1F(D);
  mat d2F(D,D);

  // Repeat until we reach the maximum number of iterations, or until
  // the convergence criterion is satisfied.
  for (unsigned int iter = 0; iter < maxiter; iter++) {
    bias = Fuv * rho;
    for (unsigned int d = 0; d < D; d++)
      d1F(d) = dot(x,Fuv.col(d)) - s * dot(Fuv.col(d),l % exp(mu + bias));
    for (unsigned int d = 0; d < D; d++)
      for (unsigned int t = 0; t < D; t++)
        d2F(d,t) = -s * dot(Fuv.col(d) % Fuv.col(t),l % exp(mu + bias));
    rhonew = rho - solve(d2F,d1F);
    rhonew.clamp(-maxrho,maxrho);
    if (max(abs(rhonew - rho)) < tol)
      break;
    rho = rhonew;
  }
}
