/*******************************************************************************
* ChangePointTests.cpp
********************************************************************************
* 2018-08-30
* Curtis Miller
********************************************************************************
* C++ functions accompanying functions in R/ChangePointTests.R.
*******************************************************************************/

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

/* Function to be wrapped; to be used only in stat_Vn(); see a description
 * there. */
// [[Rcpp::export]]
List stat_Vn_cpp(const NumericVector& dat, const double& kn, const double& tau,
                 const bool& use_kernel_var, const NumericVector& lrv_est,
                 const bool& get_all_vals) {
    double n = dat.size();
    // Get total sum and sum of squares; this will be the "upper sum"
    // (i.e. the sum above k)
    double s_upper = 0, s_square_upper = 0;
    // The "lower sums" (i.e. those below k)
    double s_lower = 0, s_square_lower = 0;
    // Get lower sums
    // Go to kn - 1 to prevent double-counting in main
    // loop
    for (int i = 0; i < kn - 1; ++i) {
        s_lower += dat[i];
        s_square_lower += dat[i] * dat[i];
    }
    // Get upper sum
    for (int i = kn - 1; i < n; ++i) {
        s_upper += dat[i];
        s_square_upper += dat[i] * dat[i];
    }
    // The maximum, which will be returned
    double M = 0;
    // A candidate for the new maximum, used in a loop+
    double M_candidate;
    // Get the complete upper sum
    double s_total = s_upper + s_lower;
    // Estimate for change location
    int est = n;
    /* A vector that will contain the value of the statistic at each n checked;
     *  relevant only if get_all_vals == TRUE */
    NumericVector all_vals = NumericVector::create();
    
    // Compute the test statistic
    for (int k = kn; k <= (n - kn); ++k) {
        // Update s and s_square for both lower and upper
        s_lower += dat[k - 1];
        s_square_lower += dat[k - 1] * dat[k - 1];
        s_upper -= dat[k - 1];
        s_square_upper -= dat[k - 1] * dat[k - 1];
        
        // Get estimate of sd for this k
        double sdk = 0;
        if (use_kernel_var) {
            sdk = std::sqrt(lrv_est[k - 1]);
        } else {
            sdk = std::sqrt((s_square_lower - s_lower * s_lower / k +
                s_square_upper -
                s_upper * s_upper / (n - k))/n);
        }
        M_candidate = std::abs(s_lower - (k / n) * s_total) / sdk /
            std::pow(k / n * (n - k) / n, tau);
        
            // If we are getting all values, add another value to all_vals
            if (get_all_vals) {
                all_vals.push_back(M_candidate / std::sqrt(n));
            }
        // Choose new maximum
        if (M_candidate > M) {
            M = M_candidate;
            est = k;
        }
    }
    // Final step; the maximum M should be normalized
    M /= std::sqrt(n);
    
    return List::create(Named("statistic") = M,
                        Named("estimate") = est,
                        Named("stat_vals") = all_vals);
}

/* Function to be wrapped; to be used only in stat_Zn(); see a description
 * there. */
// [[Rcpp::export]]
List stat_Zn_cpp(const NumericVector& dat, const double& kn,
                 const bool& use_kernel_var, const NumericVector& lrv_est,
                 const bool& get_all_vals) {
    double n = dat.size();
    // Get total sum and sum of squares; this will be the "upper sum"
    // (i.e. the sum above k)
    double s_upper = 0, s_square_upper = 0;
    // The "lower sums" (i.e. those below k)
    double s_lower = 0, s_square_lower = 0;
    // Get lower sums
    // Go to kn - 1 to prevent double-counting in main
    // loop
    for (int i = 0; i < kn - 1; ++i) {
        s_lower += dat[i];
        s_square_lower += dat[i] * dat[i];
    }
    // Get upper sum
    for (int i = kn - 1; i < n; ++i) {
        s_upper += dat[i];
        s_square_upper += dat[i] * dat[i];
    }
    // The maximum, which will be returned
    double M = 0;
    // A candidate for the new maximum, used in a loop
    double M_candidate;
    // Estimate for change location
    int est = n;
    /* A vector that will contain the value of the statistic at each n checked;
     * relevant only if get_all_vals == TRUE */
    NumericVector all_vals = NumericVector::create();
    
    // Compute the test statistic
    for (int k = kn; k <= (n - kn); ++k) {
        // Update s and s_square for both lower and upper
        s_lower += dat[k - 1];
        s_square_lower += dat[k - 1] * dat[k - 1];
        s_upper -= dat[k - 1];
        s_square_upper -= dat[k - 1] * dat[k - 1];
        
        // Get estimate of sd for this k
        double sdk = 0;
        if (use_kernel_var) {
            sdk = std::sqrt(lrv_est[k - 1]);
        } else {
            sdk = std::sqrt((s_square_lower - s_lower * s_lower / k +
                               s_square_upper -
                               s_upper * s_upper / (n - k))/n);
        }
        M_candidate = std::abs(s_lower / k - s_upper / (n - k)) / sdk;
        
        // If we are getting all values, add another value to all_vals
        if (get_all_vals) {
            all_vals.push_back(M_candidate * std::sqrt(kn));
        }
        // Choose new maximum
        if (M_candidate > M) {
            M = M_candidate;
            est = k;
        }
    }
    
    // One final step to get the test statistic;
    // Multiply the maximum by kn
    M *= std::sqrt(kn);
    
    return List::create(Named("statistic") = M,
                        Named("estimate") = est,
                        Named("stat_vals") = all_vals);
}

/* Short and unsafe function; won't check that x and A have proper dimensions.
 * But we get better speed when we can make that assumption. */
inline double norm_A_square(arma::vec x, arma::mat A) {
    arma::mat prod_mat = x.t() * A * x;
    return prod_mat(0,0);
}

inline double norm_A(arma::vec x, arma::mat A) {
    return sqrt(norm_A_square(x, A));
}

/* XXX: This may be the faster function but it will also be numerically
 * unstable since it uses the normal equations for finding least-squares
 * solutions. Perhaps should consider creating a function that is slower but
 * values numerical accuracy more. */
// [[Rcpp::export]]
List stat_Zn_lm_cpp(const NumericMatrix& X_input, const NumericMatrix& y_input,
                    const double& kn, const List& lrv_est,
                    const bool& get_all_vals) {
    unsigned int n = X_input.rows();
    unsigned int d = X_input.cols();
    if (y_input.rows() != n) {
        throw std::range_error("Bad y passed; must have one column and same "
                               "number of rows as data matrix X")
    }

    // Create Armadillo matrices
    const arma::mat X(X_input.begin(), n, d, false);
    const arma::vec y(y_input.begin(), n, false);

    /* Call X the data matrix and X' its transpose (I usually don't do this);
     * then X'X and X'y are sums. I want sums; these will be "upper sums" (i.e.
     * the sum above k). */
    arma::mat sxx_upper(d, d, arma::fill::zeros);
    arma::mat sxy_upper(d, 1, arma::fill::zeros);
    // The "lower sums" (i.e. those below k)
    arma::mat sxx_lower(d, d, arma::fill::zeros);
    arma::mat sxy_lower(d, 1, arma::fill::zeros);
    /* Get lower sums; go to kn - 1 to prevent double-counting in main loop */
    for (int i = 0; i < kn - 1; ++i) {
        sxx_lower += X.row(i).t() * X.row(i);
        sxy_lower += X.row(i).t() * y(i);
    }
    // Get upper sum
    for (int i = kn - 1; i < n; ++i) {
        sxx_upper += X.row(i).t() * X.row(i);
        sxy_upper += X.row(i).t() * y(i);
    }

    // The maximum, which will be returned
    double M = 0;
    // A candidate for the new maximum, used in a loop
    double M_candidate;
    // Estimate for change location
    int est = n;
    /* A vector that will contain the value of the statistic at each n checked;
     * relevant only if get_all_vals == TRUE */
    NumericVector all_vals = NumericVector::create();

    // Compute the test statistic
    for (int k = kn; k <= (n - kn); ++k) {
        sxx_lower += X.row(k).t() * X.row(k);
        sxy_lower += X.row(k).t() * y(i);
        sxx_upper -= X.row(k).t() * X.row(k);
        sxy_upper -= X.row(k).t() * y(i);
    }
}

// Function used for computing long-run variance; see R function get_lrv_vec()
// [[Rcpp::export]]
NumericVector get_lrv_vec_cpp(const NumericMatrix& Y, const NumericVector& kern,
                              const int& max_l) {
    // Number of data points, inferred from Y
    double n = Y.nrow();
    /* Vector that will contain estimated variances at points t; 2 <= t <= n - 2
     * (initialize with -1, an impossible value that indicates an error) */
    NumericVector sigma = NumericVector(n - 1, -1);
    
    // Start computing variances
    for (int t = 1; t <= n - 1; ++t) {
        /* Will be added over through the loop, and eventually added to sigma
         * vector */
        double sum = 0;     
        // Iterate through lags
        for (int l = 0; l <= std::min(double(max_l), n - 1); ++l) {
            int m = 2;      // A multiplier used in the sum
            if (l == 0) {
                m = 1;
            }
            
            // Iterate through sample
            for (int i = 1; i <= n - l; ++i) {
                // Next summand
                sum += m * kern[l] * Y(i - 1, t - 1) * Y(i + l - 1, t - 1) /
                       (n - l);
            }
        }
        
        sigma[t - 1] = sum;
    }
    
    return(sigma);
}
