#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Compute the expected Poisson rates under the variational
// approximation of the Poisson mash parameters, and store the result
// in "a".
void compute_poisson_rates (const vec& s, const vec& mu, const vec& bias,
			    const vec& gamma, const vec& V, vec& a) {
  a = s % exp(mu + bias + gamma + V/2);
}
