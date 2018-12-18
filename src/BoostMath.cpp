/*******************************************************************************
 * BoostMath.cpp
 *******************************************************************************
 * 2018-12-06
 * Curtis Miller
 *******************************************************************************
 * C++ functions providing Boost mathematical functionality to R.
 ******************************************************************************/

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>
#include <vector>
#include <iterator>

using namespace Rcpp;
// [[Rcpp::export]]
NumericVector besselJ_zeros_cpp(const double& nu, const unsigned& a,
                                const unsigned& b) {
    std::vector<double> roots;
    const unsigned int m = b - a + 1;

    boost::math::cyl_bessel_j_zero(nu, a, m, std::back_inserter(roots));

    NumericVector result = wrap(roots);
    return(result);
}
