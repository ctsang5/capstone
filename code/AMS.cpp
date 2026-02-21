#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <limits>
#include <cmath>
using namespace Rcpp;
using namespace arma;
/**
 * simulate_AMS
 * Monte-Carlo simulation of price paths under:
 *   1 = Black–Scholes (exact solution)
 *   2 = Heston (Euler discretisation)
 *   3 = Heston (Milstein discretisation)
 *   4 = Heston (Quadratic–Exponential scheme, Andersen 2008)
 *
 * @param model  Integer {1,2,3,4} selecting the model.
 * @param n      Number of simulated paths (>0).
 * @param t      Maturity in years (>0).
 * @param p      Total time steps (>0).
 * @param r      Risk–free rate.
 * @param sigma  Black–Scholes volatility (>=0, used only when model==1).
 * @param S0     Initial spot price (>0).
 * @param rho    Correlation between asset and variance Brownian motions
 *               (required for Heston models, must be finite and in [-1,1]).
 * @param rim    Left–trim: discard the first `rim` time steps
 *               (0 <= rim < p). Returned matrices keep p - rim + 1 columns
 *               including the initial time.
 * @param v0     Initial variance for Heston models (>=0).
 *
 * @return Rcpp::List
 *   - model 1:  S  (n x (p - rim + 1))
 *   - models 2–3: S, V
 *   - model 4:    S, V, scheme = "QE"
 *
 * Notes
 * -----
 *  * Column-major loops follow Armadillo’s memory layout.
 *  * Rcpp::RNGScope is not needed when using Armadillo RNGs,
 *    but can be added in calling code to keep R and C++ RNGs in sync.
 *  * Complexity: O(n * (p - rim)) time, O(n * (p - rim)) memory.
 */
// [[Rcpp::export]]
Rcpp::List simulate_AMS(int model,
                        int n,
                        double t,
                        int p,
                        double r,
                        double sigma,
                        double S0,
                        Nullable<double> rho = R_NilValue,
                        int rim = 0,
                        double v0 = 0.04)
{ Rcpp::RNGScope scope; // keep RNG coherent with R
  // -------- Input checks --------
  if (n <= 0)                  stop("n must be > 0");
  if (p <= 0)                  stop("p must be > 0");
  if (t <= 0.0)                stop("t must be > 0");
  if (rim < 0 || rim >= p)     stop("rim must be in [0, p)");
  if (model < 1 || model > 4)  stop("model must be in {1,2,3,4}");
  if (S0 <= 0.0)               stop("S0 must be > 0");
  if (sigma < 0.0)             stop("sigma must be >= 0");
  if (model >= 2) {
    if (rho.isNull()) stop("rho is required for Heston models");
    double rh = as<double>(rho);
    if (!std::isfinite(rh) || rh < -1.0 || rh > 1.0)
      stop("rho must be finite and in [-1,1]");
  }

  const int cols = p - rim;                 // number of simulation steps to generate
  const double dt = t / p;                  // step size

  // -------- Black–Scholes (exact) --------
  if (model == 1) {
    arma::mat Z   = randn<arma::mat>(n, cols);
    const double drift = (r - 0.5 * sigma * sigma) * dt;
    const double sd    = sigma * std::sqrt(dt);
    arma::mat L = arma::cumsum(drift + sd * Z, 1);
    arma::mat S = S0 * arma::exp(L);
    arma::mat Sfull(n, cols + 1);
    Sfull.col(0).fill(S0);
    Sfull.cols(1, cols) = std::move(S);
    return List::create(Named("S") = Sfull);
  }

  // Common settings for Heston models
  const double rho_val = as<double>(rho);
  const double sqrt_dt = std::sqrt(dt);
  const double sqrt_1m_rho2 = std::sqrt(1.0 - rho_val * rho_val);
  const double kappa = 2.0;
  const double theta = 0.04;
  const double volvol = 0.3;  // vol of vol
  mat S(n, cols + 1, fill::zeros);
  mat V(n, cols + 1, fill::zeros);
  S.col(0).fill(S0);
  V.col(0).fill(v0);

  // -------- Heston Euler --------
  if (model == 2) {
    mat Z = randn(n * cols, 2);
    for (int c = 1; c <= cols; ++c) {
      for (int i = 0; i < n; ++i) {
        int idx = (c - 1) * n + i;
        double z1 = Z(idx, 0);
        double z2 = Z(idx, 1);
        double dW_s = sqrt_dt * z1;
        double dW_v = sqrt_dt * (rho_val * z1 + sqrt_1m_rho2 * z2);

        double v_prev = V(i, c - 1);
        double s_prev = S(i, c - 1);
        double sqrt_v = std::sqrt(std::max(0.0, v_prev));

        double v = v_prev + kappa * (theta - v_prev) * dt
                   + volvol * sqrt_v * dW_v;
        v = std::max(0.0, v);

        double inc = (r - 0.5 * v_prev) * dt + sqrt_v * dW_s;
        S(i, c) = s_prev * std::exp(inc);
        V(i, c) = v;
      }
    }
    return List::create(Named("S") = S, Named("V") = V);
  }

  // -------- Heston Milstein --------
  if (model == 3) {
    mat Z = randn(n * cols, 2);
    for (int c = 1; c <= cols; ++c) {
      for (int i = 0; i < n; ++i) {
        int idx = (c - 1) * n + i;
        double z1 = Z(idx, 0);
        double z2 = Z(idx, 1);
        double dW_s = sqrt_dt * z1;
        double dW_v = sqrt_dt * (rho_val * z1 + sqrt_1m_rho2 * z2);

        double v_prev = V(i, c - 1);
        double s_prev = S(i, c - 1);
        double sqrt_v = std::sqrt(std::max(0.0, v_prev));

        double drift_v = kappa * (theta - v_prev);
        double diff_v  = volvol * sqrt_v;

        double v = v_prev + drift_v * dt + diff_v * dW_v
                   + 0.25 * volvol * volvol * (dW_v * dW_v - dt);
        v = std::max(0.0, v);

        double inc = (r - 0.5 * v_prev) * dt + sqrt_v * dW_s;
        S(i, c) = s_prev * std::exp(inc);
        V(i, c) = v;
      }
    }
    return List::create(Named("S") = S, Named("V") = V);
  }

  // -------- Heston Quadratic–Exponential --------
  if (model == 4) {
    if (2.0 * kappa * theta <= volvol * volvol)
      stop("Feller condition not satisfied: 2*kappa*theta must be > volvol^2");

    const double psi_c  = 1.5;
    const double gamma1 = 0.5;
    const double gamma2 = 0.5;
    const double ekdt   = std::exp(-kappa * dt);
    const double K0 = -rho_val * kappa * theta * dt / volvol;
    const double K1 = gamma1 * dt * (kappa * rho_val / volvol - 0.5) - rho_val / volvol;
    const double K2 = gamma2 * dt * (kappa * rho_val / volvol - 0.5) + rho_val / volvol;
    const double K3 = gamma1 * dt * (1.0 - rho_val * rho_val);
    const double K4 = gamma2 * dt * (1.0 - rho_val * rho_val);

    mat Z_norm = randn(n * cols, 2);
    vec U      = randu(n * cols);

    for (int c = 1; c <= cols; ++c) {
      for (int i = 0; i < n; ++i) {
        int idx = (c - 1) * n + i;
        double v_prev = V(i, c - 1);
        double s_prev = S(i, c - 1);
         
        const double one_m_ekdt = 1.0 - ekdt;
        const double one_m_ekdt_sq = one_m_ekdt * one_m_ekdt;

        // conditional mean and variance of V_{t+dt} | V_t
        double m  = theta + (v_prev - theta) * ekdt;
        m = std::max(m, 1e-14);               // guard contro underflow
        double s2 = (v_prev * volvol * volvol * ekdt / kappa) * one_m_ekdt
          + (theta * volvol * volvol / (2.0 * kappa)) * one_m_ekdt_sq;
        double psi = s2 / (m * m);

        double v_new;
        if (psi <= psi_c) {
          double invpsi = 1.0 / std::max(psi, 1e-14);
          double tmp    = 2.0 * invpsi;
          double b2     = tmp - 1.0 + 2.0 * std::sqrt(tmp * (tmp - 1.0));
          b2 = std::max(b2, 1e-14);
          double a      = m / (1.0 + b2);
          double Zv     = Z_norm(idx, 0);
          v_new = a * std::pow(std::sqrt(b2) + Zv, 2.0);
        } else {
          double p_star = (psi - 1.0) / (psi + 1.0);
          p_star = std::min(std::max(p_star, 0.0), 1.0 - 1e-14);
          double beta = 2.0 / (m * (1.0 + psi));
          double Uv   = U[idx];
          if (Uv <= p_star) {
            v_new = 0.0;
          } else {
            v_new = std::log((1.0 - p_star) / std::max(1.0 - Uv, 1e-14)) / beta;
          }
        }
        v_new = std::max(0.0, v_new);

        double variance_term = std::max(0.0, K3 * v_prev + K4 * v_new);
        double eps = Z_norm(idx, 1);
        double log_inc = r * dt + K0 + K1 * v_prev + K2 * v_new
                       + std::sqrt(variance_term) * eps;
        S(i, c) = s_prev * std::exp(log_inc);
        V(i, c) = v_new;
      }
    }
    return List::create(Named("S") = S,
                        Named("V") = V,
                        Named("scheme") = "QE");
  }

  // Defensive fallback (should be unreachable)
  stop("Invalid model code");
}

/**
 * Binary payoffs from simulated paths.
 *
 * Returns 0/1 indicators for path events. Computes only what is needed per option.
 *
 * @param paths  n x m matrix of simulated prices (rows = paths, cols = time steps)
 * @param K      strike or threshold (finite)
 * @param option integer in 1..6:
 *               1 = I{ S_T > K }
 *               2 = I{ S_T < K }
 *               3 = I{ mean_t S_t > K }
 *               4 = I{ mean_t S_t < K }
 *               5 = I{ max_t S_t > K }
 *               6 = I{ min_t S_t < K }
 * @return Numeric vector of length n with values in {0,1}.
 *
 * Notes:
 *  - Comparisons with NaN yield 0 (false). No NA handling is performed.
 *  - Complexity: O(n) for options 1–2; O(n*m) for 3–6.
 */

arma::vec payoff(const arma::mat& paths, double K, int option) {
  const arma::uword n = paths.n_rows;
  const arma::uword m = paths.n_cols;
  if (n == 0 || m == 0) Rcpp::stop("paths must have at least one row and one column");
  if (!std::isfinite(K)) Rcpp::stop("K must be finite");
  if (option < 1 || option > 6) Rcpp::stop("option must be an integer in 1..6");

  const arma::uword last = m - 1;

  switch (option) {
    case 1: // I{ S_T > K }
      return arma::conv_to<arma::vec>::from(paths.col(last) > K);
    case 2: // I{ S_T < K }
      return arma::conv_to<arma::vec>::from(paths.col(last) < K);
    case 3: { // I{ mean_t S_t > K }
      arma::vec A = arma::mean(paths, 1);
      return arma::conv_to<arma::vec>::from(A > K);
    }
    case 4: { // I{ mean_t S_t < K }
      arma::vec A = arma::mean(paths, 1);
      return arma::conv_to<arma::vec>::from(A < K);
    }
    case 5: { // I{ max_t S_t > K }
      arma::vec M = arma::max(paths, 1);
      return arma::conv_to<arma::vec>::from(M > K);
    }
    case 6: { // I{ min_t S_t < K }
      arma::vec mvec = arma::min(paths, 1);
      return arma::conv_to<arma::vec>::from(mvec < K);
    }
    default:
      Rcpp::stop("unreachable");
  }
}

/**
 * function_AMS_Cpp
 * Score matrix for AMS algorithms.
 *
 * Feature x_{i,j} per path i and step j:
 *   1: spot S_{i,j}            ->  I{S_T > K} (funz=1) or +S_{i,j} (funz=2)
 *   2: spot S_{i,j}            ->  I{S_T < K} (funz=1) or -S_{i,j} (funz=2)
 *   3: running average         ->  I{avg > K} (funz=1) or +avg     (funz=2)
 *   4: running average         ->  I{avg < K} (funz=1) or -avg     (funz=2)
 *   5: running maximum         ->  I{max > K} (funz=1) or +max     (funz=2)
 *   6: running minimum         ->  I{min < K} (funz=1) or -min     (funz=2)
 *
 * If funz==1:
 *   a(i,j) = exp(-r * T_j) * Phi( sgn * d2 ),
 *   where d2 = [log(x/K) + (r - 0.5*sigma^2)*T_j] / [sigma * sqrt(T_j)],
 *   sgn = +1 for options {1,3,5}, -1 for {2,4,6}.
 *
 * If funz==2:
 *   a(i,j) = ± x_{i,j}  (sign as above).
 *
 * @param S_paths  n x p matrix of simulated prices (row=path, col=time step).
 * @param option   integer in {1..6}.
 * @param funz     1 for BS digital proxy, 2 for signed feature.
 * @param strike   strike K (>0 when funz==1).
 * @param r        risk-free rate.
 * @param sigma    BS volatility (>0 when funz==1).
 * @param time     maturity in years (>0); T_j = (p - j) / (time * 252).
 *
 * @return n x p NumericMatrix of scores.
 *
 * Notes:
 *  - Column-major outer loop for cache locality.
 *  - Clamp x to avoid log(<=0).
 *  - Complexity O(n*p).
 */
NumericMatrix function_AMS_Cpp(const NumericMatrix& S_paths,
                               int option,
                               int funz,
                               double strike,
                               double r,
                               double sigma,
                               double time)
{
  const int n = S_paths.nrow();
  const int p = S_paths.ncol();

  if (n <= 0 || p <= 0) stop("S_paths must be non-empty");
  if (option < 1 || option > 6) stop("option must be in 1..6");
  if (funz != 1 && funz != 2) stop("funz must be 1 or 2");
  if (time <= 0.0) stop("time must be > 0");
  if (funz == 1) {
    if (sigma <= 0.0) stop("sigma must be > 0 when funz == 1");
    if (strike <= 0.0) stop("strike must be > 0 when funz == 1");
  }

  NumericMatrix a(n, p);

  const double inv_time252 = 1.0 / (time * 252.0);
  const double mu = r - 0.5 * sigma * sigma;

  std::vector<double> Tj(p), disc(p), denom(p);
  for (int j = 0; j < p; ++j) {
    const double T = (p - j) * inv_time252;
    Tj[j]   = T;
    disc[j] = std::exp(-r * T);
    denom[j]= sigma * std::sqrt(T > 0.0 ? T : 0.0);
  }

  const double logK = (funz == 1) ? std::log(strike) : 0.0;
  const double EPS  = 1e-300;

  std::vector<double> run_sum, run_max, run_min;
  if (option == 3 || option == 4) run_sum.assign(n, 0.0);
  if (option == 5) {
    run_max.assign(n, -std::numeric_limits<double>::infinity());
  }
  if (option == 6) {
    run_min.assign(n,  std::numeric_limits<double>::infinity());
  }

  const int sgn_comp = (option == 1 || option == 3 || option == 5) ? +1 : -1;
  const double sgn_feat = (option == 1 || option == 3 || option == 5) ? +1.0 : -1.0;

  for (int j = 0; j < p; ++j) {
    const double T   = Tj[j];
    const double dsc = disc[j];
    const double den = denom[j];
    if (funz == 1 && den == 0.0)
      stop("denominator is zero; adjust sigma or time/grid");

    for (int i = 0; i < n; ++i) {
      const double s = S_paths(i, j);
      double x;

      switch (option) {
        case 1: case 2:
          x = s; break;
        case 3: case 4:
          run_sum[i] += s;
          x = run_sum[i] / static_cast<double>(j + 1);
          break;
        case 5:
          if (s > run_max[i]) run_max[i] = s;
          x = run_max[i];
          break;
        case 6:
          if (s < run_min[i]) run_min[i] = s;
          x = run_min[i];
          break;
        default:
          x = s;
      }

      if (funz == 1) {
        const double lx  = std::log(std::max(x, EPS)) - logK;
        const double d2  = (lx + mu * T) / den;
        const double arg = static_cast<double>(sgn_comp) * d2;
        const double cdf = R::pnorm(arg, 0.0, 1.0, 1, 0);
        a(i, j) = dsc * cdf;
      } else {
        a(i, j) = sgn_feat * x;
      }
    }
  }

  return a;
}

/**
 * AMS
 * Adaptive Multilevel Splitting estimator for rare-event option payoffs.
 *
 * Pipeline per iteration:
 *   1) Simulate n paths under the chosen model (BS/Heston-family).
 *   2) Compute continuation scores a_{i,j} via function_AMS_Cpp.
 *   3) Set level L = K-th order statistic of max_j a_{i,j}.
 *   4) Identify survivors (top n-K) and parents (K indices that cleared the level).
 *   5) For each parent, cut at first index that exceeds L and resimulate suffix.
 *   6) Repeat until L >= Lmax. Then compute discounted payoff on final population.
 *
 * @param model   1=Black–Scholes; 2,3,4=Heston variants (as in simulate_AMS).
 * @param type    payoff type passed to payoff() and function_AMS_Cpp (1..6).
 * @param funz    1=BS digital proxy in continuation; 2=raw feature (signed).
 * @param n       population size (>K).
 * @param t       maturity in years (>0).
 * @param p       total time steps (>0).
 * @param r       risk-free rate.
 * @param sigma   BS volatility (used by continuation; >0 if funz==1).
 * @param S0      initial spot.
 * @param rho     correlation for Heston models (required for model>=2, in [-1,1]).
 * @param rim     left-trim for simulation (keep last p-rim steps; 0 <= rim < p).
 * @param Lmax    stopping level: iterate while L < Lmax.
 * @param strike  strike K used by continuation and final payoff.
 * @param K       number of resampled offspring per iteration (1..n-1).
 *
 * @return List with:
 *   - price : AMS estimator of discounted payoff E[e^{-rt} * payoff]
 *   - std   : Monte Carlo standard error of the mean estimator
 *
 * @details
 *  - Uses column count cols = p - rim + 1 from simulate_AMS (includes column 0).
 *  - Parents are chosen uniformly at random among paths that exceed L at some time.
 *    If fewer than K exceeders exist, the K parents are taken from the top-K maxima.
 *  - Copies between Armadillo and R are minimized but unavoidable when calling
 *    function_AMS_Cpp (operates on R matrices).
 *  - Complexity per iteration dominated by sorting O(n log n) and K suffix resimulations.
 *
 * @throws std::runtime_error on invalid inputs or inconsistent dimensions.
 */
// [[Rcpp::export]]
List AMS(int model,
         int type,
         int funz,
         int n,
         double t,
         int p,
         double r,
         double sigma,
         double S0,
         Nullable<double> rho = R_NilValue,  // Heston only
         int rim      = 0,
         double v0 = 0.04,
         double Lmax  = 0.0,
         double strike = 1.0,
         int K        = 1) {

  Rcpp::RNGScope scope; // keep RNG coherent with R

  /* ---- input validation ----------------------------------------- */
  if (n <= 1)               stop("n must be > 1");
  if (K <= 0 || K >= n)     stop("K must be in [1, n-1]");
  if (p <= 0)               stop("p must be > 0");
  if (rim < 0 || rim >= p)  stop("rim must satisfy 0 <= rim < p");
  if (t <= 0.0)             stop("t must be > 0");
  if (strike <= 0.0 || !R_finite(strike)) stop("strike must be finite and > 0");
  if (funz == 1 && (sigma <= 0.0 || !R_finite(sigma)))
    stop("sigma must be finite and > 0 when funz==1");
  if (model < 1 || model > 4) stop("model must be in {1,2,3,4}");
  if (type < 1 || type > 6)   stop("type must be in {1..6}");
  if (model >= 2) {
    if (rho.isNull()) stop("rho is required for Heston models (2..4)");
    const double rh = Rcpp::as<double>(rho);
    if (!std::isfinite(rh) || rh < -1.0 || rh > 1.0)
      stop("rho must be finite and in [-1,1]");
  }

  const int cols = p - rim + 1; // simulate_AMS returns [0..cols-1], column 0 is initial

  /* ---- 1) initial simulation ------------------------------------ */
  arma::mat S, V; // prices and variance (V only for Heston)
  Rcpp::List sim0 = simulate_AMS(model, n, t, p, r, sigma, S0, rho, rim, v0);
  S = as<arma::mat>(sim0["S"]);            // n x cols
  if (S.n_cols != static_cast<uword>(cols))
    stop("simulate_AMS returned unexpected number of columns");
  if (model >= 2) V = as<arma::mat>(sim0["V"]);

  /* ---- 2) initial continuation ---------------------------------- */
  // function_AMS_Cpp operates on R matrices; wrap S once per evaluation.
  Rcpp::NumericMatrix cont0 = function_AMS_Cpp(Rcpp::wrap(S), type, funz, strike, r, sigma, t);
  arma::mat A = as<arma::mat>(cont0);      // n x cols continuation values

  arma::vec amax = arma::max(A, 1);        // max over time per path
  arma::uvec order = arma::stable_sort_index(amax); // ascending indices
  double L = amax(order(K - 1));           // K-th order statistic
  double product = 1.0;                    // survival product ∏ (n-K)/n across iterations

  // Working buffers reused across iterations
  arma::ivec first_hit(n, fill::zeros);

  /* ---- 3) AMS loop ---------------------------------------------- */
  while (L < Lmax) {

    // 3.1 first exceedance index per path, -1 means never exceeds
    for (int i = 0; i < n; ++i) {
      int idx = -1;
      // scan left-to-right; columns are contiguous in memory
      for (int j = 0; j < cols; ++j) {
        if (A(i, j) > L) { idx = j; break; }
      }
      first_hit(i) = idx;
    }

    // 3.2 survivors (top n-K by amax) and parents
    //      parents: choose K uniformly among exceeders; if fewer than K,
    //      fall back to the K best overall (rightmost K in 'order').
    arma::uvec survivors = order.tail(n - K);

    arma::uvec exceeders = arma::find(first_hit >= 0);
    arma::uvec parents;
    if (exceeders.n_elem >= static_cast<uword>(K)) {
      // sample without replacement among exceeders
      parents = exceeders( arma::randperm(exceeders.n_elem).head(K) );
    } else {
      // fallback: use top-K overall maxima
      parents = order.tail(K);
    }

    // 3.3 generate offspring: copy prefix up to 'cut' and resimulate suffix
    arma::mat S_child(K, cols, arma::fill::zeros);
    arma::mat V_child; if (model >= 2) V_child.set_size(K, cols);

    for (int ii = 0; ii < K; ++ii) {
      const int id  = static_cast<int>(parents[ii]);
      const int cut = first_hit(id); // -1 possible if fallback used

      const int cut_clamped = std::max(cut, -1); // ensure -1..cols-1
      const int pref_last   = std::max(cut_clamped, -1); // last prefix column
      // Copy prefix [0..cut] if cut>=0; else empty prefix
      if (pref_last >= 0) {
        S_child.row(ii).cols(0, pref_last) = S.row(id).cols(0, pref_last);
        if (model >= 2) V_child.row(ii).cols(0, pref_last) = V.row(id).cols(0, pref_last);
      }

      // If cut < cols-1, resimulate suffix starting from state at 'cut'
      if (cut_clamped < cols - 1) {
        const double newS0  = (cut_clamped >= 0) ? S(id, cut_clamped) : S(id, 0);
        if (model == 1) {
          // Resimulate a single BS path with same grid, trimming left at 'cut'
          Rcpp::List child = simulate_AMS(model, 1, t, p, r, sigma, newS0, rho, /*rim=*/ (rim + std::max(cut_clamped,0)));
          arma::mat S1 = child["S"];             // 1 x (cols - cut_clamped)
          // Child S includes the initial column at its col 0; skip it.
          S_child.row(ii).cols(cut_clamped + 1, cols - 1) =
            S1.row(0).cols(1, S1.n_cols - 1);
        } else {
          const double newV0 = (cut_clamped >= 0) ? V(id, cut_clamped) : V(id, 0);
          Rcpp::List child = simulate_AMS(model, 1, t, p, r, sigma, newS0, rho,
                                          /*rim=*/ (rim + std::max(cut_clamped,0)), newV0);
          arma::mat S1 = child["S"];
          arma::mat V1 = child["V"];
          S_child.row(ii).cols(cut_clamped + 1, cols - 1) =
            S1.row(0).cols(1, S1.n_cols - 1);
          V_child.row(ii).cols(cut_clamped + 1, cols - 1) =
            V1.row(0).cols(1, V1.n_cols - 1);
        }
      }
    }

    // 3.4 rebuild population: K children + survivors (n-K)
    S = arma::join_cols(S_child, S.rows(survivors));
    if (model >= 2) V = arma::join_cols(V_child, V.rows(survivors));

    // 3.5 update continuation, level L, and survival product
    Rcpp::NumericMatrix cont = function_AMS_Cpp(Rcpp::wrap(S), type, funz, strike, r, sigma, t);
    A = as<arma::mat>(cont);

    amax  = arma::max(A, 1);
    order = arma::stable_sort_index(amax);
    L     = amax(order(K - 1));
    product *= static_cast<double>(n - K) / static_cast<double>(n);
  }

  /* ---- 4) final discounted payoff ------------------------------- */
  arma::vec pay = product * payoff(S, strike, type) * std::exp(-r * t);
  // MC standard error of the sample mean: sqrt(Var / n), with unbiased sample variance.
  const double se = std::sqrt( arma::var(pay) / static_cast<double>(n) );

  return List::create(
    _["price"] = arma::mean(pay),
    _["std"]   = se
  );
}