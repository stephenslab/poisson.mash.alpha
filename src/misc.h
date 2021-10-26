#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
void compute_poisson_rates (const arma::vec& s, const arma::vec& mu, 
			    const arma::vec& bias, const arma::vec& gamma, 
			    const arma::vec& V, arma::vec& a);

#endif
