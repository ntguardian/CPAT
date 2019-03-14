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
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

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

inline double norm_inv_A_square(arma::vec x, arma::mat A) {
    arma::mat prod_mat = x.t() * arma::solve(A, x);
    return prod_mat(0, 0);
}

inline double norm_inv_A(arma::vec x, arma::mat A) {
    return sqrt(norm_inv_A_square(x, A));
}

// TODO: curtis: THIS HAS BAD STATISTICAL PROPERTIES; WHAT IS WRONG? -- Fri 08 Mar 2019 12:00:51 AM MST
// [[Rcpp::export]]
List stat_Zn_reg_cpp(const NumericMatrix& X_input, const NumericVector& y_input,
                     const double& kn, const bool& use_kernel_var,
                     NumericVector lrv_est, const bool& get_all_vals,
                     const bool& fast = false) {
    const unsigned int n = X_input.rows();
    const unsigned int d = X_input.cols();
    if (y_input.size() != n) {
        throw std::range_error("Bad y passed; must have one column and same "
                               "number of rows as data matrix X");
    }

    /* NOTE: lrv_est needs to be a 3D array */
    const IntegerVector lrv_est_dims = lrv_est.attr("dim");
    if ((lrv_est_dims[0] != d) || (lrv_est_dims[1] != d) ||
        (lrv_est_dims[2] != n)) {
        throw std::range_error("Bad lrv_est passed");
    }
    // Create Armadillo objects
    const arma::mat X = as<arma::mat>(X_input);
    const arma::vec y = as<arma::vec>(y_input);
    arma::cube lrv_est_cube(lrv_est.begin(), d, d, n, false);
    const arma::vec eps = y - X * arma::solve(X, y);
    const arma::mat eX = X.each_col() % eps;

    /* Call X the data matrix and X' its transpose (I usually don't do this);
     * then X'X and X'y are sums. I want sums; these will be "upper sums" (i.e.
     * the sum above k). */
    arma::mat sxx_upper(d, d, arma::fill::zeros);
    arma::mat sxy_upper(d, 1, arma::fill::zeros);
    // The "lower sums" (i.e. those below k)
    arma::mat sxx_lower(d, d, arma::fill::zeros);
    arma::mat sxy_lower(d, 1, arma::fill::zeros);

    arma::mat C(d, d, arma::fill::zeros);   // A scaling matrix used later
    // Normal equations will be used for fast computation
    /* Get lower sums; go to kn - 1 to prevent double-counting in main
     * loop */
    for (int i = 0; i < kn - 1; ++i) {
        sxx_lower += X.row(i).t() * X.row(i);
        if (fast) {
            sxy_lower += X.row(i).t() * y(i);
        }
    }
    // Get upper sum
    for (int i = kn - 1; i < n; ++i) {
        sxx_upper += X.row(i).t() * X.row(i);
        if (fast) {
            sxy_upper += X.row(i).t() * y(i);
        }
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

    // Lower/upper coefficient estimates
    arma::vec beta_lower(d, arma::fill::zeros);
    arma::vec beta_upper(d, arma::fill::zeros);

    // If not using pre-computed LRV, prepare for computing it
    arma::mat Q_lower(d, d, arma::fill::zeros);
    arma::mat Q_upper(d, d, arma::fill::zeros);
    if (!(use_kernel_var)) {
        for (int i = 0; i < n / 2; ++i) {
            Q_lower += eX.row(i).t() * eX.row(i);
            lrv_est_cube.slice(i) = Q_lower / (double(i) + 1);
        }
        for (int i = (n - 1); i >= n / 2; --i) {
            Q_upper += eX.row(i).t() * eX.row(i);
            lrv_est_cube.slice(i) = Q_upper / (double(n) - i);
        }
    }

    // Compute the test statistic
    for (int k = kn; k <= (n - kn); ++k) {
        sxx_lower += X.row(k - 1).t() * X.row(k - 1);
        sxx_upper -= X.row(k - 1).t() * X.row(k - 1);
        if (fast) {
            sxy_lower += X.row(k - 1).t() * y(k - 1);
            sxy_upper -= X.row(k - 1).t() * y(k - 1);

            beta_lower = arma::solve(sxx_lower, sxy_lower);
            beta_upper = arma::solve(sxx_upper, sxy_upper);
        } else {
            beta_lower = arma::solve(X.head_rows(k), y.head(k));
            beta_upper = arma::solve(X.tail_rows(n - k), y.tail(n - k));
        }

        if (k <= n/2) {
            C = sxx_lower / k;
        } else {
            C = sxx_upper / (n - k + 1);
        }
        
        M_candidate = norm_inv_A_square(C * (beta_lower - beta_upper),
                                        lrv_est_cube.slice(k - 1));
        // If we are getting all values, add another value to all_vals
        if (get_all_vals) {
            all_vals.push_back(sqrt(M_candidate * kn));
        }
        // Choose new maximum
        if (M_candidate > M) {
            M = M_candidate;
            est = k;
        }
    }

    /* One final step to get the test statistic; multiply the maximum by kn */
    M = sqrt(M);
    M *= sqrt(kn);

    return List::create(Named("statistic") = M,
                        Named("estimate") = est,
                        Named("stat_vals") = all_vals);
}

// XXX: curtis: THIS CODE WORKS CORRECTLY BUT THE RESULTS ARE GARBAGE FROM A
// STATISTICAL PERSPECTIVE -- Fri 11 Jan 2019 06:12:52 PM MST
List stat_de_reg_cpp(const NumericMatrix& X_input, const NumericVector& y_input,
                     const double& kn, const double& a_n, const double& b_n,
                     const bool& get_all_vals, const bool& fast = false) {
    const unsigned int n = X_input.rows();
    const unsigned int d = X_input.cols();
    if (y_input.size() != n) {
        throw std::range_error("Bad y passed; must have one column and same "
                               "number of rows as data matrix X");
    }
    // Create Armadillo objects
    const arma::mat X = as<arma::mat>(X_input);
    const arma::vec y = as<arma::vec>(y_input);
    const arma::mat Sigma = X.t() * X;
    
    /* Call X the data matrix and X' its transpose; then X'X and X'y are sums. I
     * want sums; these will be "upper sums" (i.e. the sum above k). */
    arma::mat sxx_upper(d, d, arma::fill::zeros);
    arma::mat sxy_upper(d, 1, arma::fill::zeros);
    // The "lower sums" (i.e. those below k)
    arma::mat sxx_lower(d, d, arma::fill::zeros);
    arma::mat sxy_lower(d, 1, arma::fill::zeros);
    if (fast) {
        // Normal equations will be used for fast computation
        /* Get lower sums; go to kn - 1 to prevent double-counting in main 
         * loop */
        for (int i = 0; i < kn - 1; ++i) {
            sxx_lower += X.row(i).t() * X.row(i);
            sxy_lower += X.row(i).t() * y(i);
        }
        // Get upper sum
        for (int i = kn - 1; i < n; ++i) {
            sxx_upper += X.row(i).t() * X.row(i);
            sxy_upper += X.row(i).t() * y(i);
        }
    }

    // The maximum, which will be returned
    double M;
    // A candidate for the new maximum, used in a loop
    double M_candidate;
    // Estimate for change location
    int est = n;
    // Omega matrix involved in covariance estimator
    arma::mat Omega(d, d, arma::fill::zeros);
    /* A vector that will contain the value of the statistic at each n checked;
     * relevant only if get_all_vals == TRUE */
    NumericVector all_vals = NumericVector::create();

    // Lower/upper coefficient estimates
    arma::vec beta_lower(d, arma::fill::zeros);
    arma::vec beta_upper(d, arma::fill::zeros);

    double s2 = 0;  // Used in the loop (part of estimating Omega)

    // Compute the test statistic
    for (int k = kn + 1; k < (n - kn); ++k) {
        if (fast) {
            sxx_lower += X.row(k - 1).t() * X.row(k - 1);
            sxy_lower += X.row(k - 1).t() * y(k - 1);
            sxx_upper -= X.row(k - 1).t() * X.row(k - 1);
            sxy_upper -= X.row(k - 1).t() * y(k - 1);

            beta_lower = arma::solve(sxx_lower, sxy_lower);
            beta_upper = arma::solve(sxx_upper, sxy_upper);
        } else {
            beta_lower = arma::solve(X.head_rows(k), y.head(k));
            beta_upper = arma::solve(X.tail_rows(n - k), y.tail(n - k));
        }

        Omega *= 0;
        for (int j = 0; j < n; ++j) {
            if (j < k) {
                s2 = y(j) - arma::dot(beta_lower, X.row(j));
            } else {
                s2 = y(j) - arma::dot(beta_upper, X.row(j));
            }
            s2 *= s2;   // Square it
            Omega += s2 * (X.row(j).t() * X.row(j));
        }

        // Compute candidate
        M_candidate = norm_inv_A_square(Sigma * (beta_lower - beta_upper),
                                        Omega);
        M_candidate *= static_cast<double>(k * (n - k)) /
            static_cast<double>(n * n);
        M_candidate -= b_n;
        M_candidate /= a_n;

        // If we are getting all values, add another value to all_vals
        if (get_all_vals) {
            all_vals.push_back(M_candidate);
        }
        // Choose new maximum
        if ((M_candidate > M) || (k == kn)) {
            M = M_candidate;
            est = k;
        }
    }

    return List::create(Named("statistic") = M,
                        Named("estimate") = est,
                        Named("stat_vals") = all_vals);
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

// XXX: curtis: INCOMPLETE -- Mon 21 Jan 2019 11:55:38 PM MST
// [[Rcpp::export]]
NumericVector get_reg_lrv_arr_cpp(const NumericMatrix& X_input,
                                  const NumericVector& y_input,
                                  const NumericVector& kern,
                                  const int& max_l, const bool& fast = false) {
    unsigned int n = X_input.rows();
    unsigned int d = X_input.cols();
    if (y_input.size() != n) {
        throw std::range_error("Bad y passed; must have one column and same "
                               "number of rows as data matrix X");
    }

    const arma::mat X = as<arma::mat>(X_input);
    const arma::vec y = as<arma::vec>(y_input);
    // n - 2 * d is because you need at least d observations to not be
    // underdetermined
    arma::cube lrv_est_cube(d, d, n - 2 * d, arma::fill::zeros);

    /* Call X the data matrix and X' its transpose; then X'X and X'y are sums. I
     * want sums; these will be "upper sums" (i.e. the sum above k). */
    arma::mat sxx_upper(d, d, arma::fill::zeros);
    arma::mat sxy_upper(d, 1, arma::fill::zeros);
    // The "lower sums" (i.e. those below k)
    arma::mat sxx_lower(d, d, arma::fill::zeros);
    arma::mat sxy_lower(d, d, arma::fill::zeros);
    if (fast) {
        // Normal equations will be used for fast computation
        /* Get lower sums; go to d - 1 to prevent double-counting in main
         * loop */
        for (int i = 0; i < d - 1; ++i) {
            sxx_lower += X.row(i).t() * X.row(i);
            sxy_lower += X.row(i).t() * y(i);
        }
        // Get upper sum
        for (int i = d - 1; i < n; ++i) {
            sxx_upper += X.row(i).t() * X.row(i);
            sxy_upper += X.row(i).t() * y(i);
        }
    }

    // Lower/upper coefficient estimates
    arma::vec beta_lower(d, arma::fill::zeros);
    arma::vec beta_upper(d, arma::fill::zeros);

    // Estimated error vectors
    arma::vec eps_lower(n, arma::fill::zeros);
    arma::vec eps_upper(n, arma::fill::zeros);

    for (int k = d; k <= (n - d); ++k) {
        if (fast) {
            sxx_lower += X.row(k - 1).t() * X.row(k - 1);
            sxy_lower += X.row(k - 1).t() * y(k - 1);
            sxx_upper -= X.row(k - 1).t() * X.row(k - 1);
            sxy_upper -= X.row(k - 1).t() * y(k - 1);

            beta_lower = arma::solve(sxx_lower, sxy_lower);
            beta_upper = arma::solve(sxx_upper, sxy_upper);
        } else {
            beta_lower = arma::solve(X.head_rows(k), y.head(k));
            beta_upper = arma::solve(X.tail_rows(n - k), y.tail(n - k));
        }

        eps_lower = y - arma::vec(X * beta_lower);
        eps_upper = y - arma::vec(X * beta_upper);

        // TODO: curtis: REST OF ALGORITHM -- Mon 07 Jan 2019 05:38:36 PM MST
    }

    // XXX: curtis: INCOMPLETE -- Mon 21 Jan 2019 11:54:44 PM MST
    // Following line added to suppress compiler warning
    return(y_input);
}
