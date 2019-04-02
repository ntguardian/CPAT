/*******************************************************************************
 * ChangePointTests.cpp
 *******************************************************************************
 * 2018-08-30
 * Curtis Miller
 *******************************************************************************
 * C++ functions accompanying functions in R/ChangePointTests.R.
 ******************************************************************************/

// [[Rcpp::depends(BH)]]
#include <boost/math/constants/constants.hpp>
#include <list>
#include <iterator>
// #include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Constants
// For kernels
const unsigned char KERN_CUSTOM   = 0;
const unsigned char KERN_TRUNC    = 1;
const unsigned char KERN_BARTLETT = 2;
const unsigned char KERN_PARZEN   = 3;
const unsigned char KERN_TUKHAN   = 4;
const unsigned char KERN_QS       = 5;
// Math constants
const double L_PI = boost::math::constants::pi<double>();

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
    return std::sqrt(norm_A_square(x, A));
}

inline double norm_inv_A_square(arma::vec x, arma::mat A) {
    arma::mat prod_mat = x.t() * arma::solve(A, x);
    return prod_mat(0, 0);
}

inline double norm_inv_A(arma::vec x, arma::mat A) {
    return std::sqrt(norm_inv_A_square(x, A));
}

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
            all_vals.push_back(std::sqrt(M_candidate * kn));
        }
        // Choose new maximum
        if (M_candidate > M) {
            M = M_candidate;
            est = k;
        }
    }

    /* One final step to get the test statistic; multiply the maximum by kn */
    M = std::sqrt(M);
    M *= std::sqrt(kn);

    return List::create(Named("statistic") = M,
                        Named("estimate") = est,
                        Named("stat_vals") = all_vals);
}

// Kernel functions used in LRV estimation
// Truncated kernel
inline double tr_kernel(const double& x) {
    if ((x < -1) || (x > 1)) {
        return 0;
    } else {
        return 1;
    }
}

// Bartlett kernel
inline double ba_kernel(const double& x) {
    const double absx = std::abs(x);
    if (absx > 1) {
        return absx;
    } else {
        return 1 - absx;
    }
}

// Parzen kernel
inline double pa_kernel(const double& x) {
    const double absx = std::abs(x);
    if (absx <= 0.5) {
        return 1 - 6 * std::pow(absx, 2) + 6 * std::pow(absx, 3);
    } else if (absx <= 1) {
        return 2 * std::pow(1 - absx, 3);
    } else {
        return 0;
    }
}

// Tukey-Hanning kernel
inline double th_kernel(const double& x) {
    if ((x < -1) || (x > 1)) {
        return 0;
    } else {
        return (1 + cos(L_PI * x))/2;
    }
}

// Quadratic spectral kernel
inline double qs_kernel(const double& x) {
    const double spxd5 = 6 * L_PI * x / 5;
    if (x == 0) {
        return 1;
    }
    return (sin(spxd5)/spxd5 - cos(spxd5)) / (std::pow(spxd5, 2));
}

// Bandwidth functions, for getting the bandwidth
// Bartlett kernel
inline double ba_bandwidth(const double& param, const unsigned int& n) {
    // In Andrews (1991), the parameter was 1.1447; round up, to be conservative
    return 1.1448 * param * std::pow(n, 1.0/3.0);
}

// Parzen kernel
inline double pa_bandwidth(const double& param, const unsigned int& n) {
    // In Andrews (1991), the parameter was 2.6614; round up, to be conservative
    return 2.6615 * param * std::pow(n, 1.0/5.0);
}

// Tukey-Hanning kernel
inline double th_bandwidth(const double& param, const unsigned int& n) {
    // In Andrews (1991), the parameter was 1.7462; round up, to be conservative
    return 1.7463 * param * std::pow(n, 1.0/5.0);
}

// Quadratic spectral kernel
inline double qs_bandwidth(const double& param, const unsigned int& n) {
    // In Andrews (1991), the parameter was 1.3221; round up, to be conservative
    return 1.3222 * param * std::pow(n, 1.0/5.0);
}

// Truncated kernel
inline double tr_bandwidth(const double& param, const unsigned int& n) {
    // As recommended by Lin and Sataka (2013)
    // (DOI: https://doi.org/10.1007/978-1-4614-1653-1_15)
    return qs_bandwidth(param, n) / 3.0;
}

// Short function for outer product of vectors
inline arma::mat outer_product(arma::vec x, arma::vec y) {
    return x * y.t();
}

inline arma::mat lag_product(arma::mat X, int s, int u) {
    arma::mat slice = outer_product(X.row(s).t(), X.row(s + u).t());
    return slice + slice.t();
}

double get_bandwidth(const double& param, const unsigned int& n,
                     const unsigned char& kernel) {
    switch (kernel) {
        case KERN_TRUNC:
            return tr_bandwidth(param, n);
            break;  // Safety first!
        case KERN_BARTLETT:
            return ba_bandwidth(param, n);
            break;
        case KERN_PARZEN:
            return pa_bandwidth(param, n);
            break;
        case KERN_TUKHAN:
            return th_bandwidth(param, n);
            break;
        case KERN_QS:
            return qs_bandwidth(param, n);
            break;
        default:
            throw std::domain_error("get_bandwidth() given bad kernel"
                                    " identifier; should not be in default"
                                    " case!");
            return 1;
            break;
    }

    throw std::runtime_error("get_bandwidth() should not have left its switch!"
                             " How did I get here?");
    return 1;
}

// Generic kernel function
double kernel_function(const double& x, const unsigned char& kernel) {
    switch (kernel) {
        case KERN_TRUNC:
            return tr_kernel(x);
            break;  // Safety first!
        case KERN_BARTLETT:
            return ba_kernel(x);
            break;
        case KERN_PARZEN:
            return pa_kernel(x);
            break;
        case KERN_TUKHAN:
            return th_kernel(x);
            break;
        case KERN_QS:
            return qs_kernel(x);
            break;
        default:
            throw std::domain_error("kernel_function() given bad kernel"
                                    " identifier; should not be in default"
                                    " case!");
            return 1;
            break;
    }

    throw std::runtime_error("kernel_function() should not have left its"
                             " switch! How did I get here?");
    return 1;
}

/* 
 * Function that computes a sequence of long-run variance matrices and returns
 * them in an armadillo cube.
 */
arma::cube lrv_matrix_cube_computer(const arma::mat X,
                                    const unsigned char& kernel,
                                    const double& param,
                                    const NumericVector& custom_bw,
                                    const Function& custom_kernel,
                                    const bool& use_custom_bw = false) {
    typedef std::list<arma::mat> mat_list;
    typedef mat_list::size_type mat_list_size;

    const unsigned int n = X.n_rows;
    const unsigned int d = X.n_cols;
    arma::cube covs(d, d, n, arma::fill::zeros);
    arma::mat slice_cov(d, d, arma::fill::zeros);
    mat_list mat_cum_sums;
    double bandwidth;
    mat_list_size maxlag;
    // Loop variables
    unsigned int u = 0;
    mat_list::iterator li;

    // Determine number of elements in list, and initialize
    if ((kernel == KERN_QS) || (kernel == KERN_CUSTOM)) {
        maxlag = n - 1;
    } else {
        if (use_custom_bw) {
            bandwidth = custom_bw[n - 1];
        } else {
            bandwidth = get_bandwidth(param, n - 1, kernel);
        }
        maxlag = std::min(int(bandwidth), int(n - 1));
    }
    mat_cum_sums.resize(maxlag + 1, arma::mat(d, d, arma::fill::zeros));

    for (int k = 0; k < n; ++k) {
        slice_cov *= 0;
        if (use_custom_bw) {
            bandwidth = custom_bw[k];
        } else {
            bandwidth = get_bandwidth(param, k, kernel);
        }
        if ((kernel == KERN_QS) || (kernel == KERN_CUSTOM)) {
            maxlag = k;
        } else {
            maxlag = std::min(int(bandwidth), k);
        }

        for (u = 0, li = mat_cum_sums.begin();
                (li != mat_cum_sums.end()) && (u <= k);
                ++u, ++li) {
            *li += lag_product(X, k - u, u);
            if (u == 0) {
                slice_cov += *li / double(2 * (k + 1));
            } else if (u <= maxlag) {
                if (kernel == KERN_CUSTOM) {
                    NumericVector kval = custom_kernel(double(u)/bandwidth);
                    slice_cov += kval[0] * (*li) / 
                        double(k + 1 - u);
                } else {
                    slice_cov += kernel_function(double(u)/bandwidth, kernel) *
                        (*li) / double(k + 1 - u);
                }
            }
        }

        covs.slice(k) = slice_cov;
    }

    return covs;
}

arma::cube flipfb(const arma::cube& X) {
    const unsigned int m = X.n_rows;
    const unsigned int n = X.n_cols;
    const unsigned int p = X.n_slices;
    arma::cube res(m, n, p);

    for (int i = 0; i < p; ++i) {
        res.slice(p - i - 1) = X.slice(i);
    }

    return res;
}

/* Function used for computing long-run covariance matrices; see R function
 * get_lrv_arr() */
// [[Rcpp::export]]
NumericVector get_lrv_arr_cpp(const NumericMatrix& X_input,
                              const unsigned char& kernel,
                              const double& bandwidth_param,
                              const NumericVector& custom_bw,
                              const Function& custom_kernel,
                              const bool& use_custom_bw = false) {
    const unsigned int n = X_input.rows();
    const unsigned int d = X_input.cols();
    const arma::mat X = as<arma::mat>(X_input);
    arma::cube lrv_est(d, d, n, arma::fill::zeros);

    lrv_est.slices(0, n/2 - 1) = lrv_matrix_cube_computer(X.rows(0, n/2 - 1),
                                                          kernel,
                                                          bandwidth_param,
                                                          custom_bw,
                                                          custom_kernel,
                                                          use_custom_bw);
    arma::cube last_cube_half = lrv_matrix_cube_computer(
            arma::flipud(X.rows(n/2, n - 1)), kernel, bandwidth_param,
            custom_bw, custom_kernel, use_custom_bw);
    lrv_est.slices(n/2, n - 1) = flipfb(last_cube_half);

    return wrap(lrv_est);
}

// Type definitions used in get_lrv_vec_cpp()

typedef struct tempXsum {
    long double lower_nolag;
    long double lower_lag;
    long double upper_nolag;
    double upper_lag;
} Xsum;
typedef struct tempXmean {
    long double lower_mean;
    long double upper_mean;
} Xmean;

/* Function used for computing long-run variance of univariate data; see
 * get_lrv_vec() for documentation. */
// [[Rcpp::export]]
NumericVector get_lrv_vec_cpp(const NumericVector& X, const NumericVector& kern,
                              const unsigned int& max_l) {
    typedef std::list<Xsum> Xsum_list;
    typedef Xsum_list::size_type Xsum_list_size;
    typedef std::list<Xmean> Xmean_list;
    typedef Xmean_list::size_type Xmean_list_size;

    const unsigned int n = X.length();
    /* Vector that will contain estimated variances at points t; 2 <= t <= n - 2
     * (initialize with -1, an impossible value that indicates an error) */
    NumericVector sigma = NumericVector(n - 1, 0);
    // Storage for t sums and means
    Xsum_list Xsums_data;
    Xmean_list Xmeans_data;

    // Intermediate variables
    long double sum = 0;
    long double lower_sum = 0;
    long double upper_sum = 0;
    long double lag_sum = 0;
    long double gamma = 0;
    long double m;
    Xsum_list_size vecsize = n - 1;
    // Loop variables
    unsigned int u;
    Xsum_list::iterator li;
    Xmean_list::iterator mi;

    // Initialize list
    const Xsum Xsum_zero = {0, 0, 0, 0};
    Xsums_data.resize(vecsize, Xsum_zero);

    // Prepare upper sum
    for (int i = 0; i < n; ++i) {
        upper_sum += X[i];
    }

    for (int l = 0; l <= std::min(int(max_l), int(n - 1)); ++l) {
        if (l == 0) {
            m = 1;
        } else {
            m = 2 * kern[l];
        }
        
        // Compute lag_sum Σ X_s X_{s + l}
        lag_sum = 0;
        for (int i = 0; i < n - l; ++i) {
            lag_sum += X[i] * X[i + l];
        }

        /*
        std::cout << "\nl = " << l << "; m = " << m << "; lag_su = " <<
            lag_sum << '\n';
        */

        // Other initialization
        mi = Xmeans_data.begin();
        // TODO: curtis: FIX THIS LOOP -- Sun 31 Mar 2019 01:32:53 AM MDT
        for (u = 0, li = Xsums_data.begin();
                (li != Xsums_data.end()) && (u < n); ++u, ++li) {
            if (l == 0) {
                // First time through, this basically needs initialization
                lower_sum += X[u];
                upper_sum -= X[u];
                li->lower_nolag = lower_sum;
                li->lower_lag = lower_sum;
                li->upper_nolag = upper_sum;
                li->upper_lag = upper_sum;
                Xmeans_data.push_back((Xmean){
                        .lower_mean = (lower_sum / (u + 1)),
                        .upper_mean = (upper_sum / (n - u - 1))
                        });
                gamma = (lag_sum - (u + 1) * std::pow(lower_sum / (u + 1), 2) -
                        (n - u - 1) * std::pow(upper_sum / (n - u - 1), 2)) /
                        (n);

            } else {
                // TODO: curtis: FIX ME -- Sun 31 Mar 2019 01:32:59 AM MDT
                li->lower_lag -= X[l - 1];
                li->upper_nolag -= X[n - l];
                if (u + 1 >= l) {
                    li->lower_nolag -= X[u + 1 - l];
                    li->upper_nolag += X[u + 1 - l];
                }
                if (u <= n - l - 1) {
                    li->lower_lag += X[u + l];
                    li->upper_lag -= X[u + l];
                } else {
                    li->upper_lag = 0;
                }
                if ((l <= u) && (u <= n - l - 1)) {
                    gamma = (lag_sum - mi->lower_mean * li->lower_lag -
                            mi->upper_mean * li->upper_lag -
                            mi->lower_mean * li->lower_nolag -
                            mi->upper_mean * li->upper_nolag + (u + 1 - l) *
                            std::pow(mi->lower_mean, 2) + l * mi->lower_mean *
                            mi->upper_mean + (n - l - u - 1) *
                            std::pow(mi->upper_mean, 2)) / (n - l);
                } else if ((u < l) && (u <= n - l - 1)) {
                    gamma = (lag_sum - mi->lower_mean * li->lower_lag -
                            mi->upper_mean * li->upper_nolag -
                            mi->upper_mean * li->upper_lag + (u + 1) *
                            mi->upper_mean * mi->lower_mean + (n - l - u - 1) *
                            std::pow(mi->upper_mean, 2)) / (n - l);
                } else if ((u >= l) && (u > n - l - 1)) {
                    gamma = (lag_sum - mi->lower_mean * li->lower_lag -
                            mi->lower_mean * li->lower_nolag - mi->upper_mean *
                            li->upper_nolag + (u + 1 - l) *
                            std::pow(mi->lower_mean, 2) + (n - u - 1) *
                            mi->lower_mean * mi->upper_mean) / (n - l);
                } else if ((n - l - 1 < u) && (u < l)) {
                    gamma = (lag_sum - mi->lower_mean * li->lower_lag -
                            mi->upper_mean * li->upper_nolag + (n - l) *
                            mi->lower_mean * mi->upper_mean) / (n - l);
                } else {
                    // THROW EXCEPTION; SHOULD NOT BE HERE!
                    throw std::runtime_error("get_lrv_vec_cpp() should not have"
                                             " reached the 'else' in the long "
                                             "if/else statement; how did I get "
                                             "here?");
                    break;
                }
            }

            /*
            std::cout << "u = " << u << "; gamma = " << gamma << '\n';
            std::cout << "lower_mean = " << mi->lower_mean << 
                "; upper_mean = " << mi->upper_mean << '\n';
            std::cout << "lower_lag = " << li->lower_lag <<
                "; upper_lag = " << li->upper_lag << "; lower_nolag = " <<
                li->lower_nolag << "; upper_nolag = " << li->upper_nolag <<
                '\n';
            */

            sigma[u] += m * gamma;
            ++mi;
        }
    }

    return sigma;
}

/* Function used for computing long-run variance; see R function
 * get_lrv_vec_old1(); now deprecated in favor of get_lrv_vec_cpp() */
// [[Rcpp::export]]
NumericVector get_lrv_vec_old1_cpp(const NumericMatrix& Y,
                                   const NumericVector& kern,
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
    
    return sigma;
}

