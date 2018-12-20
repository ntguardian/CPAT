/*******************************************************************************
* XOUTExperimental.cpp
********************************************************************************
* 2018-08-30
* Curtis Miller
********************************************************************************
* C++ functions for experimental functionality and accompanying function in
* XOUTExprimental.R.
*******************************************************************************/

// #include <Rcpp.h>
#include <RcppArmadillo.h>

// Constants used for marking derivative locations
#define PO 0
#define PA 1
#define PB 2

using namespace Rcpp;

/* Function used for computing partial derivatives of conditional variances;
   see R function get_gradient_hessian() */

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List cond_var_gradient_hessian_cpp(const NumericVector& var,
                                   const NumericVector& eps,
                                   const double& omega, const double& alpha,
                                   const double& beta,
                                   const NumericVector& init_vals) {
  int n = var.size();
  NumericVector gradient = NumericVector(Dimension(3,n));  /* 3xn for three
                                                              parameters and n
                                                              gradients */
  arma::cube hessian(3,3,n);  /* 3x3xn for three parameters and n matrices */
  
  /* Starting partial derivatives; these are temporary and will update through
     the loop, while their results are saved in the loop */
  // Read: p equates to partial, o equates to omega, a to alpha, b to beta
  double po = init_vals[0];
  double pa = init_vals[1];
  double pb = init_vals[2];
  
  double popo = 0;
  double papo = init_vals[3];
  double pbpo = init_vals[4];
  double papa = init_vals[5];
  double pbpa = init_vals[6];
  double pbpb = init_vals[7];
  
  for (int i = 0; i < n; ++i) {
    // First entries of gradient and hessian should be init values
    // Otherwise, iterate as needed
    if (i > 0) {
      /* Do second derivatives first since first derivatives' old values show
         up in calculation of second derivative values */
      papo = beta * papo;
      pbpo = po + beta * pbpo;
      papa = beta * papa;
      pbpa = pa + beta * pbpa;
      pbpb = 2 * pb + beta * pbpb;
      
      // Now first derivatives
      po = 1 + beta * po;
      pa = eps[i - 1] + beta * pa;
      pb = var[i - 1] + beta * pb;
    }
    
    // Assign to matrices
    // First, gradient
    gradient(PO, i) = po;
    gradient(PA, i) = pa;
    gradient(PB, i) = pb;
    
    // Now, Hessian
    hessian(PO, PO, i) = popo;
    hessian(PA, PO, i) = papo;
    hessian(PO, PA, i) = papo;
    hessian(PB, PO, i) = pbpo;
    hessian(PO, PB, i) = pbpo;
    hessian(PA, PA, i) = papa;
    hessian(PB, PA, i) = pbpa;
    hessian(PA, PB, i) = pbpa;
    hessian(PB, PB, i) = pbpb;
  }
  
  return List::create(Named("gradient") = gradient,
                      Named("hessian") = hessian);
}
