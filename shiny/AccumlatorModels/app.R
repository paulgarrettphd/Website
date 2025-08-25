################################################################################
#  ONE-STOP ACCUMULATOR SIMULATOR (snapshot-based reactivity + kill switch)
#  Models:
#   - RW (random walk)
#   - LBA (fixed bound)
#   - LBA-CB (linear collapsing upper bound)
#   - LCA (leaky competing accumulator)
#   - DDM (fixed bounds)
#   - DDM-HCB (hyperbolic collapsing bound; upper-only or ±a(t))
#   - WIENER (Ratcliff DDM with across-trial variability)
#   - CDM  (vector 2D diffusion to circle)
#   - SCDM (sinusoidal field over angle; crest hits bound)
################################################################################

library(shiny)
library(rtdists)     # DDM & Wiener
library(circular)    # polar density helpers
library(Rcpp)        # fast CDM core
library(shinybusy)   # spinner
library(later)       # chunked, cancel-able runs

cache_dir <- file.path(tempdir(), "accum_rcpp_cache")
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

.onLoad <- local({
  function(...) {
    # no-op call that forces JIT compile
    try({ invisible(sim_cdm_cpp(1, 0.1, 0.1, 1.0, 0.01, 1.0, 0.0)) }, silent = TRUE)
  }
})

Rcpp::sourceCpp(code='
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// -------------------------------------------------------------------
// CDM
// -------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericMatrix sim_cdm_cpp(int    n,
                                double v,
                                double s,
                                double R,
                                double dt,
                                double cosT,
                                double sinT) {
  Rcpp::NumericMatrix out(n, 2);              // rt, angle
  double vx = v * cosT,  vy = v * sinT;
  for (int i = 0; i < n; ++i) {
    double x = 0.0, y = 0.0, t = 0.0;
    while (true) {
      t += dt;
      x += vx * dt + R::rnorm(0.0, s * std::sqrt(dt));
      y += vy * dt + R::rnorm(0.0, s * std::sqrt(dt));
      if (x*x + y*y >= R*R) {
        out(i,0) = t;
        out(i,1) = std::atan2(y, x);
        break;
      }
      if (t > 20.0) { out(i,0)=NA_REAL; out(i,1)=NA_REAL; break; }
    }
  }
  return out;
}

// -------------------------------------------------------------------
// LCA (path + batch simulate)
// -------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List lca_one_trial_path_cpp(double v, double s, double b, double a,
                                  double kappa, double beta, double dt,
                                  double max_t, int downsample_every) {
  int max_steps = (int)std::ceil(max_t / dt);
  std::vector<double> t;     t.reserve(max_steps + 1);
  std::vector<double> x1_v;  x1_v.reserve(max_steps + 1);
  std::vector<double> x2_v;  x2_v.reserve(max_steps + 1);

  double x1 = R::runif(0.0, a);
  double x2 = R::runif(0.0, a);
  double time = 0.0;

  t.push_back(0.0); x1_v.push_back(x1); x2_v.push_back(x2);

  for (int k = 0; k < max_steps; ++k) {
    time += dt;
    // Euler–Maruyama
    x1 += (v        - kappa*x1 - beta*x2) * dt + s * std::sqrt(dt) * R::rnorm(0.0, 1.0);
    x2 += ((1.0 - v) - kappa*x2 - beta*x1) * dt + s * std::sqrt(dt) * R::rnorm(0.0, 1.0);

    if (x1 < 0.0) x1 = 0.0;
    if (x2 < 0.0) x2 = 0.0;

    if (downsample_every <= 1 || (k % downsample_every) == 0) {
      t.push_back(time); x1_v.push_back(x1); x2_v.push_back(x2);
    }

    if (x1 >= b || x2 >= b) {
      int winner = (x1 >= b && x2 >= b) ? (x1 >= x2 ? 1 : 2) : (x1 >= b ? 1 : 2);
      return Rcpp::List::create(
        Rcpp::Named("hit")    = true,
        Rcpp::Named("rt")     = time,
        Rcpp::Named("winner") = winner,
        Rcpp::Named("t")      = t,
        Rcpp::Named("x1")     = x1_v,
        Rcpp::Named("x2")     = x2_v
      );
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("hit")    = false,
    Rcpp::Named("rt")     = NA_REAL,
    Rcpp::Named("winner") = NA_INTEGER,
    Rcpp::Named("t")      = t,
    Rcpp::Named("x1")     = x1_v,
    Rcpp::Named("x2")     = x2_v
  );
}

// [[Rcpp::export]]
Rcpp::DataFrame simulate_lca_cpp(int n, double v, double s, double b, double a,
                                 double kappa, double beta, double dt, double max_t) {
  Rcpp::NumericVector RT(n);
  Rcpp::LogicalVector correct(n);
  Rcpp::NumericVector start(n);
  int max_steps = (int)std::ceil(max_t / dt);

  for (int i = 0; i < n; ++i) {
    double x1 = R::runif(0.0, a);
    double x2 = R::runif(0.0, a);
    double t = 0.0; bool hit = false; int winner = NA_INTEGER;

    for (int k = 0; k < max_steps; ++k) {
      t += dt;
      x1 += (v        - kappa*x1 - beta*x2) * dt + s * std::sqrt(dt) * R::rnorm(0.0,1.0);
      x2 += ((1.0 - v) - kappa*x2 - beta*x1) * dt + s * std::sqrt(dt) * R::rnorm(0.0,1.0);
      if (x1 < 0.0) x1 = 0.0;
      if (x2 < 0.0) x2 = 0.0;
      if (x1 >= b || x2 >= b) {
        hit = true;
        winner = (x1 >= b && x2 >= b) ? (x1 >= x2 ? 1 : 2) : (x1 >= b ? 1 : 2);
        RT[i] = t;
        correct[i] = (winner == 1);
        start[i] = (winner == 1 ? x1 : x2);
        break;
      }
    }
    if (!hit) { RT[i] = NA_REAL; correct[i] = NA_LOGICAL; start[i] = NA_REAL; }
  }

  return Rcpp::DataFrame::create(Rcpp::Named("RT", RT),
                                 Rcpp::Named("correct", correct),
                                 Rcpp::Named("start", start));
}

// -------------------------------------------------------------------
// LDLIV (dual race) + path
// -------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::DataFrame simulate_ldliv_cpp(int n,
                                   double b,
                                   double qC, double qE,
                                   double k, double sigma,
                                   double zprop, double t0,
                                   double dt, double max_t) {

  b      = std::max(1e-9, b);
  k      = std::max(0.0, k);
  sigma  = std::max(0.0, sigma);
  zprop  = std::min(std::max(zprop, 0.0), 0.99);

  int max_steps = (int)std::ceil(max_t / dt);
  double sqrt_dt = std::sqrt(dt);
  double z0 = b * zprop;
  double two_sigma = 2.0 * sigma;

  Rcpp::NumericVector RT(n);
  Rcpp::LogicalVector correct(n);
  Rcpp::NumericVector start(n);

  for (int i = 0; i < n; ++i) {
    double xC = z0, xE = z0, t = 0.0;
    bool hit = false;

    for (int step = 0; step < max_steps; ++step) {
      t += dt;

      // noise sd uses *previous* state
      double sdC = std::sqrt(std::max(0.0, two_sigma * std::max(0.0, xC))) * sqrt_dt;
      double sdE = std::sqrt(std::max(0.0, two_sigma * std::max(0.0, xE))) * sqrt_dt;

      xC += (qC - k * xC) * dt + sdC * R::rnorm(0.0, 1.0);
      xE += (qE - k * xE) * dt + sdE * R::rnorm(0.0, 1.0);

      if (xC < 0.0) xC = 0.0;
      if (xE < 0.0) xE = 0.0;

      if (xC >= b || xE >= b) {
        // tie-break: larger state wins
        bool c_wins;
        if (xC >= b && xE >= b) c_wins = (xC >= xE);
        else                    c_wins = (xC >= b);

        RT[i]      = t0 + t;
        correct[i] = c_wins;
        start[i]   = z0;
        hit = true;
        break;
      }
    }

    if (!hit) {
      RT[i] = NA_REAL;
      correct[i] = NA_LOGICAL;
      start[i]   = z0;
    }
  }

  return Rcpp::DataFrame::create(Rcpp::Named("RT", RT),
                                 Rcpp::Named("correct", correct),
                                 Rcpp::Named("start", start));
}

// [[Rcpp::export]]
Rcpp::List ldliv_path_cpp(double rt, double b,
                          double q, double k, double sigma,
                          double z_abs, double dt,
                          int downsample_every) {
  if (!std::isfinite(rt) || rt <= 0.0) {
    return Rcpp::List::create(Rcpp::Named("t") = Rcpp::NumericVector(0),
                              Rcpp::Named("x") = Rcpp::NumericVector(0));
  }

  int n = (int)std::ceil(rt / dt);
  double sqrt_dt = std::sqrt(dt);
  double two_sigma = 2.0 * sigma;

  std::vector<double> T;  T.reserve(n+1);
  std::vector<double> X;  X.reserve(n+1);

  double x = z_abs;
  T.push_back(0.0); X.push_back(x);

  for (int i = 1; i <= n; ++i) {
    double sd = std::sqrt(std::max(0.0, two_sigma * std::max(0.0, x))) * sqrt_dt;
    x += (q - k * x) * dt + sd * R::rnorm(0.0, 1.0);
    if (x < 0.0) x = 0.0;
    if (x > b)   x = b;

    if (downsample_every <= 1 || (i % downsample_every) == 0) {
      T.push_back(i * dt);
      X.push_back(x);
    }
  }

  return Rcpp::List::create(Rcpp::Named("t") = T,
                            Rcpp::Named("x") = X);
}

// -------------------------------------------------------------------
// MDFT (simulate + one-trial path)
// -------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::DataFrame simulate_mdft_cpp(int n, int M,
                                  double b, double phi, double beta,
                                  double rho, double sigma_s,
                                  double v1, double vC, double sv,
                                  double s, double dt, double t0,
                                  double max_t, double dpos) {
  // guards
  M        = std::max(2, M);
  b        = std::max(1e-9, b);
  phi      = std::max(0.0, phi);
  beta     = std::max(0.0, beta);
  sigma_s  = std::max(1e-12, sigma_s);
  s        = std::max(0.0, s);
  dt       = std::max(1e-6, dt);
  max_t    = std::max(1e-3, max_t);
  dpos     = std::max(1e-12, dpos);

  int max_steps = (int)std::ceil(max_t / dt);
  double sqrt_dt = std::sqrt(dt);

  // positions and S_off (constant across trials)
  std::vector<double> pos(M);
  for (int i = 0; i < M; ++i) pos[i] = i * dpos;

  std::vector<double> S_off(M * M, 0.0);
  double denom = 2.0 * sigma_s * sigma_s;
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      if (i == j) continue;
      double d = std::abs(pos[i] - pos[j]);
      S_off[i*M + j] = std::exp(-(d*d) / denom);
    }
  }

  // outputs
  Rcpp::NumericVector RT(n);
  Rcpp::LogicalVector correct(n);

  // main trial loop
  for (int tr = 0; tr < n; ++tr) {
    // valences v: target=1 gets v1; others vC; sv noise; mean-center
    std::vector<double> v(M);
    for (int i = 0; i < M; ++i) {
      double base = (i == 0 ? v1 : vC);
      v[i] = base + (sv > 0.0 ? R::rnorm(0.0, sv) : 0.0);
    }
    double meanv = 0.0;
    for (int i = 0; i < M; ++i) meanv += v[i];
    meanv /= (double)M;
    for (int i = 0; i < M; ++i) v[i] -= meanv;

    // U = (I - rho*S_off) v  (diag of S_off is 0)
    std::vector<double> U(M, 0.0);
    for (int i = 0; i < M; ++i) {
      double acc = v[i];
      if (rho != 0.0) {
        double sdot = 0.0;
        for (int j = 0; j < M; ++j) if (j != i) sdot += S_off[i*M + j] * v[j];
        acc -= rho * sdot;
      }
      U[i] = acc;
    }

    // dynamics
    std::vector<double> x(M, 0.0), x_new(M, 0.0);
    double t = 0.0;
    bool hit = false;
    int winner = -1;

    for (int k = 0; k < max_steps; ++k) {
      t += dt;

      // Lx = (sum x)1 - x => drift_i = (-(phi - beta) x_i - beta*sumx + U_i)
      double sumx = 0.0;
      for (int i = 0; i < M; ++i) sumx += x[i];

      for (int i = 0; i < M; ++i) {
        double drift = (-(phi - beta) * x[i]) - beta * sumx + U[i];
        double noise = (s > 0.0 ? s * sqrt_dt * R::rnorm(0.0, 1.0) : 0.0);
        x_new[i] = x[i] + drift * dt + noise;
      }
      x.swap(x_new);

      // bound check
      int imax = 0;
      double xmax = x[0];
      for (int i = 1; i < M; ++i) {
        if (x[i] > xmax) { xmax = x[i]; imax = i; }
      }
      if (xmax >= b) {
        hit = true;
        winner = imax; // 0-based
        RT[tr] = t0 + t;
        correct[tr] = (winner == 0); // option 1 is target
        break;
      }
    }

    if (!hit) {
      RT[tr] = NA_REAL;
      correct[tr] = NA_LOGICAL;
    }
  }

  return Rcpp::DataFrame::create(Rcpp::Named("RT", RT),
                                 Rcpp::Named("correct", correct));
}

// [[Rcpp::export]]
Rcpp::List mdft_one_trial_path_cpp(int M,
                                   double b, double phi, double beta,
                                   double rho, double sigma_s,
                                   double v1, double vC, double sv,
                                   double s, double dt, double max_t,
                                   double dpos, int downsample_every) {
  // guards
  M        = std::max(2, M);
  b        = std::max(1e-9, b);
  phi      = std::max(0.0, phi);
  beta     = std::max(0.0, beta);
  sigma_s  = std::max(1e-12, sigma_s);
  s        = std::max(0.0, s);
  dt       = std::max(1e-6, dt);
  max_t    = std::max(1e-3, max_t);
  dpos     = std::max(1e-12, dpos);
  int ds   = std::max(1, downsample_every);

  int max_steps = (int)std::ceil(max_t / dt);
  int max_rows  = (max_steps / ds) + 2;      // + start row + safety
  double sqrt_dt = std::sqrt(dt);

  // positions & S_off (off-diagonal similarities)
  std::vector<double> pos(M);
  for (int i = 0; i < M; ++i) pos[i] = i * dpos;

  std::vector<double> S_off(M * M, 0.0);
  double denom = 2.0 * sigma_s * sigma_s;
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      if (i == j) continue;
      double d = std::abs(pos[i] - pos[j]);
      S_off[i*M + j] = std::exp(-(d*d) / denom);
    }
  }

  // trial-specific valences v, mean-centered
  std::vector<double> v(M);
  for (int i = 0; i < M; ++i) {
    double base = (i == 0 ? v1 : vC);
    v[i] = base + (sv > 0.0 ? R::rnorm(0.0, sv) : 0.0);
  }
  double meanv = 0.0;
  for (int i = 0; i < M; ++i) meanv += v[i];
  meanv /= (double)M;
  for (int i = 0; i < M; ++i) v[i] -= meanv;

  // U = (I - rho*S_off) v  (diag(S_off)=0)
  std::vector<double> U(M, 0.0);
  if (rho == 0.0) {
    for (int i = 0; i < M; ++i) U[i] = v[i];
  } else {
    for (int i = 0; i < M; ++i) {
      double sdot = 0.0;
      for (int j = 0; j < M; ++j) if (j != i) sdot += S_off[i*M + j] * v[j];
      U[i] = v[i] - rho * sdot;
    }
  }

  // preallocate storage for downsampled path
  Rcpp::NumericMatrix X(max_rows, M);
  Rcpp::NumericVector T(max_rows);
  int r = 0;
  auto push_row = [&](const std::vector<double>& x, double tt){
    if (r >= max_rows) return;
    for (int j = 0; j < M; ++j) X(r, j) = x[j];
    T[r] = tt;
    ++r;
  };

  // dynamics
  std::vector<double> x(M, 0.0), x_new(M, 0.0);
  push_row(x, 0.0);

  double t = 0.0;
  bool hit = false;
  int winner = -1;

  for (int k = 0; k < max_steps; ++k) {
    t += dt;

    double sumx = 0.0;
    for (int i = 0; i < M; ++i) sumx += x[i];

    for (int i = 0; i < M; ++i) {
      // drift: -phi*x - beta*((sum x) - x) + U
      double drift = (-(phi - beta) * x[i]) - beta * sumx + U[i];
      double noise = (s > 0.0 ? s * sqrt_dt * R::rnorm(0.0, 1.0) : 0.0);
      x_new[i] = x[i] + drift * dt + noise;
    }
    x.swap(x_new);

    if ((k % ds) == 0) push_row(x, t);

    // bound check
    int imax = 0; double xmax = x[0];
    for (int i = 1; i < M; ++i) if (x[i] > xmax) { xmax = x[i]; imax = i; }
    if (xmax >= b) { hit = true; winner = imax; break; }
  }

  // shrink to filled rows
  Rcpp::NumericMatrix Xout(r, M);
  for (int i = 0; i < r; ++i) for (int j = 0; j < M; ++j) Xout(i, j) = X(i, j);
  Rcpp::NumericVector Tout(r);
  for (int i = 0; i < r; ++i) Tout[i] = T[i];

  return Rcpp::List::create(
    Rcpp::Named("hit")    = hit,
    Rcpp::Named("rt")     = t,
    Rcpp::Named("winner") = (winner >= 0 ? winner + 1 : NA_INTEGER), // 1-based for R
    Rcpp::Named("t")      = Tout,
    Rcpp::Named("X")      = Xout
  );
}
', cacheDir = cache_dir)


# ────────────────────────────────────────────────────────────────────────────
#  Helpers: RNGs & simulators
rtruncnorm <- function(n, mean, sd, a = 0) {
  x <- rnorm(n, mean, sd)
  while (any(bad <- x <= a))
    x[bad] <- rnorm(sum(bad), mean, sd)
  x
}

lca_fixed_point_check <- function(v, kappa, beta, b, margin = 0.95) {
  det <- kappa^2 - beta^2
  if (det <= 0) {
    return(list(
      xstar = c(NA_real_, NA_real_), det = det,
      verdict = "inconclusive",
      msg = "Matrix not strictly contracting (kappa ≤ beta); cannot conclude from fixed point."
    ))
  }
  x1 <- ( kappa * v     - beta * (1 - v)) / det
  x2 <- ( kappa * (1-v) - beta * v      ) / det
  # Clamp negatives to 0 (due to rectification in your simulator)
  xstar <- pmax(c(x1, x2), 0)
  mx <- max(xstar)
  if (!is.finite(mx)) {
    return(list(xstar = c(NA_real_, NA_real_), det = det,
                verdict = "unknown", msg = "Numerical issue estimating fixed point."))
  }
  if (mx < margin * b) {
    return(list(
      xstar = xstar, det = det,
      verdict = "no-hit-likely",
      msg = sprintf("Fixed point below bound: max(x*)=%.3g < %.0f%% of b=%.3g.",
                    mx, 100*margin, b)
    ))
  }
  list(xstar = xstar, det = det,
       verdict = "inconclusive",
       msg = "Fixed point not clearly below bound; hits may still occur (esp. with noise).")
}

safe_scalar <- function(x, default) {
  if (length(x) != 1 || is.na(x) || !is.finite(x)) default else x
}

simulate_rw <- function(n, p, h, a, zprop, dt){
  zprop <- min(max(zprop, 1e-6), 1 - 1e-6)
  RT <- numeric(n); correct <- logical(n); start <- numeric(n); steps <- integer(n)
  z <- a * zprop
  max_steps <- max(1L, ceiling(20 / dt))
  for (i in seq_len(n)) {
    x <- z; t <- 0L
    for (k in seq_len(max_steps)) {
      t <- t + 1L
      x <- x + if (runif(1) < p) h else -h
      if (x >= a) { RT[i] <- t*dt; correct[i] <- TRUE;  start[i] <- z; steps[i] <- t; break }
      if (x <= 0) { RT[i] <- t*dt; correct[i] <- FALSE; start[i] <- z; steps[i] <- t; break }
    }
    if (steps[i] == 0L) { RT[i] <- NA_real_; correct[i] <- NA; start[i] <- z; steps[i] <- NA_integer_ }
  }
  data.frame(RT=RT, correct=correct, start=start, steps=steps)
}

simulate_lba <- function(n, v, s, b, a) {
  start_c <- runif(n, 0, a)
  start_e <- runif(n, 0, a)
  drift_c <- rtruncnorm(n, v,   s)
  drift_e <- rtruncnorm(n, 1-v, s)
  t_c <- (b - start_c) / drift_c
  t_e <- (b - start_e) / drift_e
  winC <- t_c < t_e
  data.frame(RT = ifelse(winC, t_c, t_e), correct = winC,
             start = ifelse(winC, start_c, start_e))
}

# LBA with collapsing upper bound: b(t) = max(b_min, b0 - k t)
simulate_lba_cb <- function(n, v, s, b0, a, k, bmin){
  a_eff   <- min(a, b0)
  start_c <- runif(n, 0, a_eff)
  start_e <- runif(n, 0, a_eff)
  drift_c <- rtruncnorm(n, v,   s)
  drift_e <- rtruncnorm(n, 1-v, s)
  tf <- if (k > 0) (b0 - bmin)/k else Inf
  hit_time_vec <- function(start, d){
    if (k == 0) return((b0 - start)/d)
    tc <- (b0 - start) / (d + k)
    use_sloped <- tc <= tf
    tstar <- (bmin - start) / d
    t <- ifelse(use_sloped, tc, pmax(tf, tstar))
    pmax(t, 0)
  }
  t_c <- hit_time_vec(start_c, drift_c)
  t_e <- hit_time_vec(start_e, drift_e)
  winC <- t_c < t_e
  data.frame(RT = ifelse(winC, t_c, t_e),
             correct = winC,
             start = ifelse(winC, start_c, start_e))
}

simulate_ddm <- function(n, v, s, a, zprop) {
  raw <- rdiffusion(n, a = a, v = v, s = s, z = a*zprop, t0 = 0)
  data.frame(RT = raw$rt, correct = raw$response == "upper",
             start = a * zprop)
}

# DDM with hyperbolic collapsing bound(s)
simulate_ddm_hcb <- function(n, v, s, a0, zprop, k, amin, dt, two_sided = TRUE){
  RT <- numeric(n); correct <- logical(n); start <- numeric(n)
  a0   <- max(1e-6, a0)
  amin <- max(0, min(amin, a0))
  max_steps <- max(1L, ceiling(20 / dt))
  for (i in seq_len(n)) {
    t <- 0.0; hit <- FALSE
    if (two_sided) {
      x0 <- (2*zprop - 1) * a0
      x  <- x0
      start[i] <- x0
      for (step in seq_len(max_steps)) {
        t   <- t + dt
        a_t <- max(amin, a0 / (1 + k * t))
        x   <- x + v*dt + s*sqrt(dt)*rnorm(1)
        if (x >=  a_t) { RT[i] <- t; correct[i] <- TRUE;  hit <- TRUE; break }
        if (x <= -a_t) { RT[i] <- t; correct[i] <- FALSE; hit <- TRUE; break }
      }
    } else {
      x0 <- a0 * zprop
      x  <- x0
      start[i] <- x0
      for (step in seq_len(max_steps)) {
        t   <- t + dt
        a_t <- max(amin, a0 / (1 + k * t))
        x   <- x + v*dt + s*sqrt(dt)*rnorm(1)
        if (x >= a_t) { RT[i] <- t; correct[i] <- TRUE;  hit <- TRUE; break }
        if (x <= 0)   { RT[i] <- t; correct[i] <- FALSE; hit <- TRUE; break }
      }
    }
    if (!hit) { RT[i] <- NA_real_; correct[i] <- NA }
  }
  data.frame(RT = RT, correct = correct, start = start)
}

simulate_cdm <- function(n, v, s, R, theta_deg) {
  R         <- safe_scalar(R, 1)
  theta_deg <- safe_scalar(theta_deg, 0)
  s         <- safe_scalar(s, 0.15)
  v         <- safe_scalar(v, 0.25)
  th <- theta_deg*pi/180
  as.data.frame(sim_cdm_cpp(n, v, s, R, 0.002, cos(th), sin(th)))
}

# LCA (2 accumulators)
simulate_lca <- function(n, v, s, b, a, k, beta, dt){
  RT <- numeric(n); correct <- logical(n); start <- numeric(n)
  max_steps <- max(1L, ceiling(20/dt))
  for (i in seq_len(n)){
    x1 <- runif(1, 0, a); x2 <- runif(1, 0, a)
    x1_0 <- x1; x2_0 <- x2
    t <- 0L; hit <- FALSE
    for (kstep in seq_len(max_steps)){
      t <- t + 1L
      x1 <- x1 + (v     - k*x1 - beta*x2)*dt + s*sqrt(dt)*rnorm(1)
      x2 <- x2 + ((1-v) - k*x2 - beta*x1)*dt + s*sqrt(dt)*rnorm(1)
      if (x1 < 0) x1 <- 0
      if (x2 < 0) x2 <- 0
      if (x1 >= b || x2 >= b){
        RT[i] <- t*dt
        if (x1 >= b && x2 >= b){
          w1 <- x1 >= x2
          correct[i] <- w1
          start[i] <- if (w1) x1_0 else x2_0
        } else if (x1 >= b){
          correct[i] <- TRUE;  start[i] <- x1_0
        } else {
          correct[i] <- FALSE; start[i] <- x2_0
        }
        hit <- TRUE
        break
      }
    }
    if (!hit){
      RT[i] <- NA_real_; correct[i] <- NA; start[i] <- if (runif(1)<0.5) x1_0 else x2_0
    }
  }
  data.frame(RT=RT, correct=correct, start=start)
}

# Wiener (across-trial variability)
simulate_wiener <- function(n, a, v, s, zprop, t0, sz = 0, sv = 0, st0 = 0) {
  zprop_i <- if (sz > 0) pmin(pmax(zprop + runif(n, -sz/2, sz/2), 1e-6), 1 - 1e-6) else rep(zprop, n)
  v_i     <- if (sv > 0) v + rnorm(n, 0, sv) else rep(v, n)
  t0_i    <- if (st0 > 0) t0 + runif(n, 0, st0) else rep(t0, n)
  z_abs   <- a * zprop_i
  raw <- rdiffusion(n, a = a, v = v_i, s = s, z = z_abs, t0 = t0_i)
  data.frame(RT = raw$rt, correct = raw$response == "upper", start = z_abs)
}

# ────────────────────────────────────────────────────────────────────────────
#  NEW: SCDM (sinusoidal field)
#  - Discretize circle into M angles; template f_kappa centered at theta
#  - Update y(θ,t) += v * f(θ) * dt + s * sqrt(dt) * N(0,1)
#  - Stop when max(y) >= b; response angle = argmax(y)
simulate_scdm <- function(n, v, s, b, kappa, M, theta_deg, dt){
  M <- max(8L, as.integer(round(M)))
  dt <- max(1e-4, dt)
  th0 <- theta_deg * pi/180
  # grid (exclude endpoint)
  ang <- seq(0, 2*pi, length.out = M + 1L); ang <- ang[-length(ang)]
  # von Mises template (fallback to cosine if kappa==0)
  if (kappa > 0) {
    vm <- exp(kappa * cos((ang - th0)))
    f  <- vm - mean(vm)
    f  <- f / max(f)  # crest = 1
  } else {
    f  <- cos(ang - th0)   # zero-mean already; crest = 1
  }
  max_steps <- max(1L, ceiling(20/dt))
  RT <- numeric(n); angle <- numeric(n); correct <- logical(n)
  for (i in seq_len(n)) {
    y <- rep(0.0, M); t <- 0.0; hit <- FALSE
    for (k in seq_len(max_steps)) {
      t <- t + dt
      y <- y + v * f * dt + s * sqrt(dt) * rnorm(M)
      if (max(y) >= b) {
        jstar    <- which.max(y)
        RT[i]    <- t
        angle[i] <- ang[jstar]
        # pseudo "correct": within ±30°
        err <- ((angle[i] - th0 + pi) %% (2*pi)) - pi
        correct[i] <- abs(err) <= (pi/6)
        hit <- TRUE
        break
      }
    }
    if (!hit) { RT[i] <- NA_real_; angle[i] <- NA_real_; correct[i] <- NA }
  }
  data.frame(RT = RT, angle = angle, correct = correct)
}

# Bernoulli–Weibull (descriptive)
#  RT ~ t0 + Weibull(k, λ);  Accuracy ~ Bernoulli(p)
simulate_bern_weibull <- function(n, k, lambda, p, t0 = 0) {
  rt  <- t0 + stats::rweibull(n, shape = k, scale = lambda)
  cor <- stats::rbinom(n, 1, prob = p) == 1
  data.frame(RT = rt, correct = cor)
}

# Poisson Counter (2-choice race of Poisson counters)
#  Each counter waits for N events at rate λ ⇒ T ~ Gamma(shape=N, rate=λ)
#  Decision = first counter to finish. Optional non-decision time t0.
simulate_poisson_counter <- function(n, N, lambda_c, lambda_e, t0 = 0) {
  N <- max(1L, as.integer(round(N)))
  lambda_c <- max(1e-9, lambda_c)
  lambda_e <- max(1e-9, lambda_e)
  t_c <- stats::rgamma(n, shape = N, rate = lambda_c)
  t_e <- stats::rgamma(n, shape = N, rate = lambda_e)
  winC <- t_c < t_e
  data.frame(RT = t0 + pmin(t_c, t_e), correct = winC)
}

# LDLIV (Linear Drift, Linear Infinitesimal Variance) dual-race
#  SDE per accumulator: dX = (q - k X) dt + sqrt(2 σ X) dW,   X ≥ 0 (0 is natural boundary)
#  Two independent accumulators (correct vs error) race to a common bound b.
simulate_ldliv <- function(n, b, q_c, q_e, k, sigma, zprop = 0, t0 = 0, dt = 0.001, max_t = 20) {
  b     <- max(1e-6, b)
  k     <- max(0, k)
  sigma <- max(0, sigma)
  zprop <- min(max(zprop, 0), 0.99)
  z0    <- b * zprop
  max_steps <- max(1L, ceiling(max_t / dt))
  
  RT <- numeric(n); correct <- logical(n); start <- numeric(n)
  for (i in seq_len(n)) {
    xC <- z0; xE <- z0; t <- 0.0; hit <- FALSE
    for (step in seq_len(max_steps)) {
      t <- t + dt
      # Euler–Maruyama increments with X clamped to R+
      sqrtC <- sqrt(pmax(xC, 0))
      sqrtE <- sqrt(pmax(xE, 0))
      xC <- xC + (q_c - k * xC) * dt + sqrt(2 * sigma * sqrtC^2) * sqrt(dt) * rnorm(1)
      xE <- xE + (q_e - k * xE) * dt + sqrt(2 * sigma * sqrtE^2) * sqrt(dt) * rnorm(1)
      if (xC < 0) xC <- 0
      if (xE < 0) xE <- 0
      if (xC >= b || xE >= b) {
        RT[i]     <- t0 + t
        correct[i] <- (xC >= b) && !(xE >= b && xE > xC)  # tie-break to winner with larger state
        start[i]   <- z0
        hit <- TRUE
        break
      }
    }
    if (!hit) { RT[i] <- NA_real_; correct[i] <- NA; start[i] <- z0 }
  }
  data.frame(RT = RT, correct = correct, start = start)
}

# Ornstein–Uhlenbeck (fixed 0..a bounds)
simulate_ou <- function(n, theta, mu_prop, s, a, zprop, t0 = 0, dt = 0.001, max_t = 20) {
  a     <- max(1e-6, a)
  theta <- max(0, theta)
  s     <- max(0, s)
  zprop <- min(max(zprop, 1e-6), 1 - 1e-6)
  mu    <- a * min(max(mu_prop, 0), 1)
  max_steps <- max(1L, ceiling(max_t / dt))
  
  RT <- numeric(n); correct <- logical(n); start <- numeric(n)
  for (i in seq_len(n)) {
    x <- a * zprop; t <- 0.0; hit <- FALSE
    for (step in seq_len(max_steps)) {
      t <- t + dt
      x <- x + theta * (mu - x) * dt + s * sqrt(dt) * rnorm(1)
      if (x >= a) { RT[i] <- t0 + t; correct[i] <- TRUE;  start[i] <- a * zprop; hit <- TRUE; break }
      if (x <= 0) { RT[i] <- t0 + t; correct[i] <- FALSE; start[i] <- a * zprop; hit <- TRUE; break }
    }
    if (!hit) { RT[i] <- NA_real_; correct[i] <- NA; start[i] <- a * zprop }
  }
  data.frame(RT = RT, correct = correct, start = start)
}

# ────────────────────────────────────────────────────────────────────────────
#  UI
ui <- fluidPage(
  titlePanel("Accumulator-Model Simulator"),
  add_busy_spinner(spin = "fading-circle", position = "top-right", timeout = 0),
  sidebarLayout(
    sidebarPanel(
      selectInput("model", "Choose model:",
                  c("Bernoulli–Weibull"                = "bern_weibull",
                    "Poisson Counter (Gamma race)"     = "poisson",
                    "Random Walk"                      = "rw",
                    "Linear Ballistic Accumulator"     = "lba",
                    "LBA (Collapsing bounds)"          = "lba_cb",
                    "Leaky Competing Accumulator"      = "lca",
                    "Diffusion Decision Model"         = "ddm",
                    "DDM (Hyperbolic collapsing)"      = "ddm_hcb",
                    "Wiener DDM (across-trial var)"    = "wiener",
                    "Ornstein–Uhlenbeck (fixed bound)" = "ou",
                    "LDLIV DDM"                        = "ldliv",
                    "Circular Diffusion Model"         = "cdm",
                    "SCDM (sinusoidal field)"          = "scdm", 
                    "MDFT (Decision Field Theory)"     = "mdft")),
      hr(),
      numericInput("n", "N-trials", 1000, 1, 10000, 1),
      uiOutput("paramUI"),
      actionButton("run", "Run simulation", class = "btn-primary"),
      actionButton("cancel", "Kill (stop now)", class = "btn-danger"),
      downloadButton("savecsv", "Save CSV", class = "btn-success")
    ),
    mainPanel(
      plotOutput("histPlot",   height = 300),
      plotOutput("extraPlot1", height = 300),
      plotOutput("extraPlot2", height = 300),
      plotOutput("extraPlot3", height = 400)
    )
  )
)

server <- function(input, output, session) {
  lcaCache <- reactiveVal(list(key = NULL, correct = NULL, error = NULL))
  key_lca <- function(p) paste(p$v, p$s, p$b, p$a, p$kappa, p$beta, p$dt, p$maxt, sep = "|")
  
  # dynamic parameter panel
  output$paramUI <- renderUI({
    switch(input$model,
           lba = tagList(
             sliderInput("v", "Drift rate  v", 0.01, 1, 0.5, 0.01),
             sliderInput("s", "Drift SD      s", 0.01, 1, 0.1, 0.01),
             sliderInput("b", "Threshold     b",   1, 10, 3, 0.1),
             sliderInput("a", "Start max     a",   0, 10, 1, 0.1)
           ),
           lba_cb = tagList(
             sliderInput("v",    "Drift rate       v", 0.01, 1, 0.5, 0.01),
             sliderInput("s",    "Drift SD         s", 0.01, 1, 0.1, 0.01),
             sliderInput("b",    "Initial bound   b0", 0.2, 10, 3, 0.1),
             sliderInput("a",    "Start max        a", 0,   10, 1, 0.1),
             sliderInput("kc",   "Collapse rate     k", 0,   5,  0.5, 0.01),
             sliderInput("bmin", "Min bound   b_min",  0,   5,  0.5, 0.05)
           ),
           lca = tagList(
             sliderInput("v",     "Drift (correct)    v", 0.01, 1, 0.5, 0.01),
             sliderInput("s",     "Noise SD           s", 0.01, 1, 0.1, 0.01),
             sliderInput("b",     "Threshold          b", 0.2, 10, 3, 0.1),
             sliderInput("a",     "Start max          a", 0, 5, 0.5, 0.05),
             sliderInput("kappa", "Leak               k", 0, 2, 0.1, 0.01),
             sliderInput("beta",  "Lateral Inhibition beta", 0, 2, 0.1, 0.01),
             sliderInput("dt",    "Time step          dt", 0.001, 0.02, 0.002, 0.001),
             sliderInput("lca_maxt", "Max time (s)", 0.5, 30, 30, 0.5)
           ),
           ddm = tagList(
             sliderInput("v", "Drift rate  v", -3, 3, 0.5, 0.05),
             sliderInput("s", "Noise σ      s", 0.1, 2, 1, 0.05),
             sliderInput("a", "Threshold   a", 0.1, 3, 1, 0.05),
             sliderInput("z", "Start prop  z/a", 0, 1, 0.5, 0.05)
           ),
           ddm_hcb = tagList(
             sliderInput("v",    "Drift rate       v",  -3, 3, 0.5, 0.05),
             sliderInput("s",    "Noise σ          s",   0.1, 2, 1, 0.05),
             sliderInput("a0",   "Initial bound   a0",   0.2, 3, 1, 0.05),
             sliderInput("z",    "Start prop    z/a0",   0, 1, 0.5, 0.05),
             sliderInput("kh",   "Collapse rate     k",  0, 5, 0.5, 0.01),
             sliderInput("amin", "Min bound  a_min",     0, 3, 0.3, 0.05),
             sliderInput("dt",   "Time step        dt",  0.0005, 0.01, 0.001, 0.0005),
             checkboxInput("hcb_sym", "Collapse both bounds (±a(t))", TRUE)
           ),
           wiener = tagList(
             sliderInput("v",   "Mean drift v",     -3, 3, 0.5, 0.05),
             sliderInput("s",   "Noise σ     s",     0.1, 2, 1, 0.05),
             sliderInput("a",   "Boundary a",        0.1, 3, 1, 0.05),
             sliderInput("z",   "Start prop z/a",    0, 1, 0.5, 0.05),
             sliderInput("t0",  "Non-decision t0",   0, 0.6, 0.2, 0.01),
             sliderInput("sz",  "Across-trial start range sz", 0, 1, 0.1, 0.01),
             sliderInput("sv",  "Across-trial drift SD sv",    0, 2, 0.3, 0.05),
             sliderInput("st0", "Across-trial t0 range st0",   0, 0.5, 0.05, 0.01)
           ),
           cdm = tagList(
             sliderInput("v", "Drift speed v", 0.01, 1, 0.25, 0.01),
             sliderInput("s", "Noise SD    s", 0.01, 1, 0.15, 0.01),
             sliderInput("R", "Radius      R", 0.5, 2, 1, 0.05),
             sliderInput("theta", "Target angle (°)", 0, 360, 0, 5)
           ),
           scdm = tagList(
             sliderInput("v",     "Drift amplitude v", 0.01, 1, 0.20, 0.01),
             sliderInput("s",     "Noise SD        s", 0.01, 1, 0.15, 0.01),
             sliderInput("b",     "Threshold       b", 0.1,  3, 1.00, 0.05),
             sliderInput("kappa", "Tuning (κ)",        0.0, 10, 2.00, 0.1),
             sliderInput("M",     "Angular bins   M",  16, 360, 72, 1),
             sliderInput("theta", "Target angle (°)",  0,  360, 0, 5),
             sliderInput("dt",    "Time step      dt", 0.001, 0.02, 0.005, 0.001)
           ),
           rw = tagList(
             sliderInput("p",  "Step up prob  p", 0, 1, 0.55, 0.01),
             sliderInput("h",  "Step size     h", 0.01, 0.2, 0.05, 0.01),
             sliderInput("a",  "Threshold     a", 0.1, 3, 1, 0.05),
             sliderInput("z",  "Start prop   z/a", 0, 1, 0.5, 0.05),
             sliderInput("dt", "Step time     dt", 0.001, 0.02, 0.005, 0.001)
           ),
           bern_weibull = tagList(
             sliderInput("k",      "Weibull shape  k",   0.2, 5, 1.8, 0.05),
             sliderInput("lambda", "Weibull scale  λ",   0.05, 2, 0.45, 0.01),
             sliderInput("p_acc",  "Accuracy  p",        0, 1, 0.75, 0.01),
             sliderInput("t0bw",   "Non-decision t0",    0, 0.6, 0.15, 0.01)
           ),
           poisson = tagList(
             sliderInput("Nthr",   "Threshold counts  N", 1, 50, 10, 1),
             sliderInput("lamC",   "Rate (correct)  λ_c", 0.5, 30, 8, 0.1),
             sliderInput("lamE",   "Rate (error)    λ_e", 0.5, 30, 5, 0.1),
             sliderInput("t0pois", "Non-decision t0",     0, 0.6, 0.15, 0.01)
           ),
           ldliv = tagList(
             sliderInput("b_ld",  "Threshold          b",   0.2, 10, 3.0, 0.1),
             sliderInput("qc",    "Input (correct)    q_c", 0.0, 10, 1.5, 0.05),
             sliderInput("qe",    "Input (error)      q_e", 0.0, 10, 1.0, 0.05),
             sliderInput("kld",   "Leak/decay         k",   0.0,  5, 0.5, 0.01),
             sliderInput("sigld", "Noise scale      σ",     0.0,  5, 1.0, 0.05),
             sliderInput("zld",   "Start prop        z/b",  0.0,  0.99, 0.10, 0.01),
             sliderInput("t0ld",  "Non-decision t0",        0.0,  0.6, 0.15, 0.01),
             sliderInput("dtld",  "Time step          dt",  0.0005, 0.01, 0.001, 0.0005)
           ),
           mdft = tagList(
             sliderInput("M_mdft",     "Number of options m",      2,   5,    3,   1),
             sliderInput("b_mdft",     "Threshold b",              0.2, 5.0,  1.5, 0.05),
             sliderInput("leak_mdft",  "Leak \u03C6",              0.0, 2.0,  0.30,0.01),
             sliderInput("beta_mdft",  "Lateral inhibition \u03B2",0.0, 2.0,  0.50,0.01),
             sliderInput("rho_mdft",   "Similarity strength \u03C1",0.0,1.0,  0.50,0.01),
             sliderInput("sigS_mdft",  "Similarity width \u03C3_s",0.1, 5.0,  1.00,0.1),
             sliderInput("v1_mdft",    "Valence (option 1)",     -2.0, 2.0,  0.80,0.05),
             sliderInput("vC_mdft",    "Valence (competitors)",  -2.0, 2.0,  0.30,0.05),
             sliderInput("sv_mdft",    "Across-trial valence SD", 0.0, 1.0,  0.05,0.01),
             sliderInput("s_mdft",     "Noise SD  s",             0.01,1.5,  0.20,0.01),
             sliderInput("dt_mdft",    "Time step  dt",           0.0005,0.02,0.002,0.0005),
             sliderInput("maxt_mdft",  "Max time (s)",            0.5, 30.0, 10.0,0.5),
             sliderInput("dpos_mdft",  "Option spacing \u0394pos",0.1,  5.0,  1.50,0.1),
             sliderInput("t0_mdft",    "Non-decision  t0",        0.0,  0.6,  0.15,0.01)
           ),
           ou = tagList(
             sliderInput("theta_ou", "Mean reversion  θ", 0.0, 10.0, 1.0, 0.05),
             sliderInput("mu_ou",    "Set point     μ/a", 0.0,  1.0, 0.6, 0.01),
             sliderInput("s_ou",     "Noise σ        s",  0.05, 2.0, 1.0, 0.05),
             sliderInput("a_ou",     "Boundary       a",  0.1,  3.0, 1.0, 0.05),
             sliderInput("z_ou",     "Start prop   z/a",  0.0,  1.0, 0.5, 0.05),
             sliderInput("t0_ou",    "Non-decision  t0",  0.0,  0.6, 0.15, 0.01),
             sliderInput("dt_ou",    "Time step      dt", 0.0005, 0.01, 0.001, 0.0005)
           )
    )
  })
  
  # state & kill handling
  simState  <- reactiveVal(NULL)
  running   <- reactiveVal(FALSE)
  cancelled <- FALSE
  
  observeEvent(input$model, {
    cancelled <<- TRUE
    running(FALSE)
    updateActionButton(session, "run", label = "Run simulation")
    simState(NULL)
    if (identical(input$model, "rw")) updateSliderInput(session, "z", value = 0.5)
  }, ignoreInit = TRUE)
  
  observeEvent(input$cancel, { if (running()) cancelled <<- TRUE })
  session$onSessionEnded(function(){ cancelled <<- TRUE })
  
  # Run (chunked + killable)
  observeEvent(input$run, {
    req(input$model)
    if (running()) return(invisible(NULL))
    running(TRUE); cancelled <<- FALSE
    updateActionButton(session, "run", label = "Running…")
    
    mdl <- isolate(input$model)
    n   <- isolate(input$n)
    
    if (mdl == "lba") {
      p <- list(v=input$v, s=input$s, b=input$b, a=input$a, n=n)
      sim_one_chunk <- function(m) simulate_lba(m, p$v, p$s, p$b, p$a)
      
    } else if (mdl == "rw") {
      z0 <- safe_scalar(input$z, 0.5); z0 <- min(max(z0, 1e-6), 1 - 1e-6)
      p <- list(p=input$p, h=input$h, a=input$a, z=z0, dt=input$dt, n=n)
      sim_one_chunk <- function(m) simulate_rw(m, p$p, p$h, p$a, p$z, p$dt)
      
    } else if (mdl == "lba_cb") {
      p <- list(v=input$v, s=input$s, b0=input$b, a=input$a, kc=input$kc, bmin=input$bmin, n=n)
      sim_one_chunk <- function(m) simulate_lba_cb(m, p$v, p$s, p$b0, p$a, p$kc, p$bmin)
      
    } else if (mdl == "lca") {
      p <- list(
        v=input$v, s=input$s, b=input$b, a=input$a,
        kappa=input$kappa, beta=input$beta, dt=input$dt, n=n,
        maxt = if (!is.null(input$lca_maxt)) input$lca_maxt else 5,
        seed_demo = as.integer(sample.int(.Machine$integer.max, 1))
      )
      #sim_one_chunk <- function(m) simulate_lca(m, p$v, p$s, p$b, p$a, p$kappa, p$beta, p$dt)
      sim_one_chunk <- function(m) {
        simulate_lca_cpp(m, p$v, p$s, p$b, p$a, p$kappa, p$beta, p$dt, p$maxt)
      }
    } else if (mdl == "ddm") {
      p <- list(v=input$v, s=input$s, a=input$a, z=input$z, n=n)
      sim_one_chunk <- function(m) simulate_ddm(m, p$v, p$s, p$a, p$z)
      
    } else if (mdl == "ddm_hcb") {
      p <- list(v=input$v, s=input$s, a0=input$a0, z=input$z, kh=input$kh, amin=input$amin,
                dt=input$dt, two_sided=input$hcb_sym, n=n)
      sim_one_chunk <- function(m) simulate_ddm_hcb(m, p$v, p$s, p$a0, p$z, p$kh, p$amin, p$dt, two_sided = p$two_sided)
      
    } else if (mdl == "wiener") {
      p <- list(v=input$v, s=input$s, a=input$a, z=input$z, t0=input$t0, sz=input$sz, sv=input$sv, st0=input$st0, n=n)
      sim_one_chunk <- function(m) simulate_wiener(m, p$a, p$v, p$s, p$z, p$t0, p$sz, p$sv, p$st0)
      
    } else if (mdl == "cdm") {
      p <- list(v=input$v, s=input$s, R=safe_scalar(input$R,1), theta=safe_scalar(input$theta,0), n=n)
      sim_one_chunk <- function(m) {
        dat <- simulate_cdm(m, p$v, p$s, p$R, p$theta)
        names(dat) <- c("RT","angle")
        target <- p$theta*pi/180
        err <- ((dat$angle - target + pi) %% (2*pi)) - pi
        dat$correct <- abs(err) <= pi/6
        dat
      }
      
    } else if (mdl == "scdm") {
      p <- list(v=input$v, s=input$s, b=input$b, kappa=input$kappa, M=input$M,
                theta=input$theta, dt=input$dt, n=n)
      sim_one_chunk <- function(m) simulate_scdm(m, p$v, p$s, p$b, p$kappa, p$M, p$theta, p$dt)
      
    } else if (mdl == "bern_weibull") {
      p <- list(k = input$k, lambda = input$lambda, p = input$p_acc, t0 = input$t0bw, n = n)
      sim_one_chunk <- function(m) simulate_bern_weibull(m, p$k, p$lambda, p$p, p$t0)
      
    } else if (mdl == "poisson") {
      p <- list(N = input$Nthr, lamC = input$lamC, lamE = input$lamE, t0 = input$t0pois, n = n)
      sim_one_chunk <- function(m) simulate_poisson_counter(m, p$N, p$lamC, p$lamE, p$t0)
      
    }  else if (mdl == "ldliv") {
      p <- list(b = input$b_ld, qC = input$qc, qE = input$qe,
                k = input$kld, sigma = input$sigld, z = input$zld,
                t0 = input$t0ld, dt = input$dtld, n = n)
      #sim_one_chunk <- function(m) simulate_ldliv(m, p$b, p$qC, p$qE, p$k, p$sigma, p$z, p$t0, p$dt)
      sim_one_chunk <- function(m) simulate_ldliv_cpp(m, p$b, p$qC, p$qE, p$k, p$sigma, p$z, p$t0, p$dt, max_t = 20)
    
    } else if (mdl == "ou") {
      p <- list(theta = input$theta_ou, mu = input$mu_ou, s = input$s_ou,
                a = input$a_ou, z = input$z_ou, t0 = input$t0_ou, dt = input$dt_ou, n = n)
      sim_one_chunk <- function(m) simulate_ou(m, p$theta, p$mu, p$s, p$a, p$z, p$t0, p$dt)
    
    } else if (mdl == "mdft") {
      p <- list(
        M     = input$M_mdft, b = input$b_mdft,
        phi   = input$leak_mdft, beta = input$beta_mdft,
        rho   = input$rho_mdft, sigma_s = input$sigS_mdft,
        v1    = input$v1_mdft,  vC = input$vC_mdft,
        sv    = input$sv_mdft,  s  = input$s_mdft,
        dt    = input$dt_mdft,  maxt = input$maxt_mdft,
        dpos  = input$dpos_mdft, t0 = input$t0_mdft,
        n     = n
      )
      sim_one_chunk <- function(m) simulate_mdft_cpp(
        n = m, M = p$M, b = p$b, phi = p$phi, beta = p$beta,
        rho = p$rho, sigma_s = p$sigma_s,
        v1 = p$v1, vC = p$vC, sv = p$sv,
        s = p$s, dt = p$dt, t0 = p$t0,
        max_t = p$maxt, dpos = p$dpos
      )
    } else {
      running(FALSE); updateActionButton(session, "run", label = "Run simulation")
      showNotification("Unknown model id; please reselect and run again.", type = "error")
      return(invisible())
    }
    
    pb <- shiny::Progress$new(session, min = 0, max = n)
    pb$set(message = paste("Simulating", toupper(mdl)), value = 0)
    
    buf  <- NULL
    done <- 0L
    chunk_size <- if (mdl %in% c("scdm","mdft")) {
      max(10L, min(300L, floor(n/60)))
    } else {
      max(20L, min(1000L, floor(n/30)))
    }
    
    step <- function(){
      if (cancelled) {
        pb$close(); running(FALSE)
        updateActionButton(session, "run", label = "Run simulation")
        if (is.null(buf)) simState(NULL) else simState(list(model = mdl, params = p, df = buf))
        return(invisible())
      }
      m <- min(chunk_size, n - done)
      if (m <= 0L) {
        pb$close(); running(FALSE)
        updateActionButton(session, "run", label = "Run simulation")
        simState(list(model = mdl, params = p, df = buf))
        return(invisible())
      }
      dfm <- sim_one_chunk(m)
      if (is.null(buf)) buf <<- dfm else buf <<- rbind(buf, dfm)
      done <<- done + nrow(dfm)
      pb$set(value = done)
      later::later(step, 0)
    }
    later::later(step, 0)
  })
  
  # helper: accuracy labels
  acc_label <- function(correct){
    acc <- mean(correct, na.rm = TRUE); N <- sum(correct, na.rm = TRUE)
    if (!is.finite(acc)) "Acc=NA" else sprintf("(%.1f%%, N=%i)", 100*acc, N)
  }
  err_label <- function(correct){
    acc <- mean(!correct, na.rm = TRUE); N <- sum(!correct, na.rm = TRUE)
    if (!is.finite(acc)) "Err=NA" else sprintf("(%.1f%%, N=%i)", 100*acc, N)
  }
  
  # ── Histogram / density routing -----------------------------------
  output$histPlot <- renderPlot({
    st <- simState()
    if (is.null(st)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    df <- st$df
    if (st$model %in% c("cdm","scdm")) {
      par(mfrow = c(1,1), mar = c(4,4,3,1))
      hist(df$RT, breaks = 40, col = "grey85", border = "white",
           main = "RT distribution", xlab = "RT (s)")
      abline(v = mean(df$RT, na.rm = TRUE), lwd = 3)
      legend("topright", legend = "mean RT", lwd = 3, bty = "n")
    } else {
      par(mfrow = c(1,2), mar = c(4,4,3,1))
      lab <- acc_label(df$correct); labe <- err_label(df$correct)
      if (any(df$correct, na.rm = TRUE)) {
        hist(df$RT[df$correct], breaks = 30, col = "steelblue", border = "white",
             main = paste("Correct RTs ", lab), xlab = "RT (s)")
      } else {
        plot.new(); title(paste("Correct RTs\n(no trials) ", lab), col.main = "grey50")
      }
      if (any(!df$correct, na.rm = TRUE)) {
        hist(df$RT[!df$correct], breaks = 30, col = "salmon", border = "white",
             main = paste("Error RTs ", labe), xlab = "RT (s)")
      } else {
        plot.new(); title(paste("Error RTs\n(no trials) ", labe), col.main = "grey50")
      }
    }
  })
  
  output$extraPlot1 <- renderPlot({
    st <- simState()
    if (is.null(st)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    if (st$model == "mdft") {
      mdft_path_example_plot(st$params)
      return(invisible())
    }
    if (st$model %in% c("lba","lba_cb","ddm","ddm_hcb","wiener","rw","lca","ldliv","ou")) {
      pathPlot(st$df, st$model, st$params)
    } else if (st$model %in% c("cdm","scdm")) {
      polarErrPlot(st$df, st$params)
    } else {  # bern_weibull, poisson
      densPlot(st$df, st$model, st$params)
    }
  })
  
  output$extraPlot2 <- renderPlot({
    st <- simState()
    if (is.null(st)) { plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); return() }
    if (st$model %in% c("lba","lba_cb","ddm","ddm_hcb","wiener","rw","lca","ldliv","ou")) {
      densPlot(st$df, st$model, st$params)
    } else if (st$model %in% c("bern_weibull","poisson")){ 
    } else {
      rtPolarPlot(st$df, st$params)
    }
  })
  
  output$extraPlot3 <- renderPlot({
    st <- simState()
    if (is.null(st)) { 
      plot.new(); title('Model changed. Click "Run simulation".', col.main='grey50'); 
      return()
    }
    
    if (st$model == "cdm") {
      cdmPathPlot(st$df, st$params)
      
    } else if (st$model == "scdm") {
      scdmWavePlot(st$params)
      
    } else if (st$model == "lca") {
      # reproducible RNG per run
      old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
      on.exit({ if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv) }, add = TRUE)
      if (!is.null(st$params$seed_demo)) set.seed(st$params$seed_demo)
      
      fp <- lca_fixed_point_check(st$params$v, st$params$kappa, st$params$beta, st$params$b)
      if (identical(fp$verdict, "no-hit-likely")) {
        plot.new()
        title("LCA dynamic — no-hit likely from fixed-point analysis", cex.main=.95)
        mtext(sprintf("x* ≈ (%.3g, %.3g); %s", fp$xstar[1], fp$xstar[2], fp$msg),
              side=3, line=0.5, col="grey40")
        return(invisible())
      }
      
      maxT <- if (is.null(st$params$maxt)) 5 else st$params$maxt
      tr <- lca_draw_one(st$params, want = "any", max_tries = 1000L, max_t = maxT)
      if (is.null(tr)) { plot.new(); title("No bound hits found.", col.main="grey40"); return(invisible()) }
      
      plot(NA, xlim = c(0, max(tr$t)), ylim = c(0, st$params$b),
           xlab = "Time (s)", ylab = "Activation",
           main = paste("Single Trial LCA dynamic —", if (tr$winner == 1L) "correct" else "error"))
      abline(h = st$params$b, lty = 3)
      lines(tr$t, tr$x1, col = "steelblue", lwd = 2)
      lines(tr$t, tr$x2, col = "salmon",   lwd = 2)
      legend("topleft", c("Accumulator 1", "Accumulator 2"),
             col = c("steelblue","salmon"), lwd = 2, bty = "n")
      
    } else {
      plot.new(); title("")
    }
  })
  
  
  
  # ── Save CSV -------------------------------------------------------
  output$savecsv <- downloadHandler(
    filename = function(){
      st <- simState()
      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
      if (is.null(st)) return(paste0("accumulator_no-data_", ts, ".csv"))
      mdl <- st$model; p <- st$params
      num <- function(x) format(round(x, 3), trim = TRUE, scientific = FALSE)
      key <- switch(mdl,
                    lba    = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s), "_b", num(p$b), "_a", num(p$a)),
                    lba_cb = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s),
                                    "_b0", num(p$b0), "_a", num(p$a), "_k", num(p$kc), "_bmin", num(p$bmin)),
                    lca    = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s),
                                    "_b", num(p$b), "_a", num(p$a),
                                    "_k", num(p$kappa), "_beta", num(p$beta), "_dt", num(p$dt)),
                    ddm    = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s), "_a", num(p$a), "_z", num(p$z)),
                    ddm_hcb= paste0("n", p$n, "_v", num(p$v), "_s", num(p$s),
                                    "_a0", num(p$a0), "_z", num(p$z),
                                    "_k", num(p$kh), "_amin", num(p$amin),
                                    "_dt", num(p$dt), "_sym", as.integer(isTRUE(p$two_sided))),
                    wiener = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s),
                                    "_a", num(p$a), "_z", num(p$z),
                                    "_t0", num(p$t0), "_sz", num(p$sz),
                                    "_sv", num(p$sv), "_st0", num(p$st0)),
                    rw     = paste0("n", p$n, "_p", num(p$p), "_h", num(p$h),
                                    "_a", num(p$a), "_z", num(p$z), "_dt", num(p$dt)),
                    cdm    = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s), "_R", num(p$R), "_theta", num(p$theta)),
                    scdm   = paste0("n", p$n, "_v", num(p$v), "_s", num(p$s), "_b", num(p$b),
                                    "_kappa", num(p$kappa), "_M", as.integer(p$M),
                                    "_theta", num(p$theta), "_dt", num(p$dt)),
                    bern_weibull = paste0("n", p$n, "_k", num(p$k), "_lam", num(p$lambda), "_p", num(p$p), "_t0", num(p$t0)),
                    poisson      = paste0("n", p$n, "_N", as.integer(p$N), "_lamC", num(p$lamC), "_lamE", num(p$lamE), "_t0", num(p$t0)),
                    ldliv        = paste0("n", p$n, "_b", num(p$b), "_qC", num(p$qC), "_qE", num(p$qE),
                                          "_k", num(p$k), "_sig", num(p$sigma), "_z", num(p$z), "_t0", num(p$t0), "_dt", num(p$dt)),
                    ou = paste0("n", p$n, "_th", num(p$theta), "_mu", num(p$mu), "_s", num(p$s),
                                "_a", num(p$a), "_z", num(p$z), "_t0", num(p$t0), "_dt", num(p$dt)),
                    mdft = paste0("n", p$n, "_m", as.integer(p$M),
                                  "_b", num(p$b), "_phi", num(p$phi), "_beta", num(p$beta),
                                  "_rho", num(p$rho), "_sigS", num(p$sigma_s),
                                  "_v1", num(p$v1), "_vC", num(p$vC),
                                  "_sv", num(p$sv), "_s", num(p$s), "_dt", num(p$dt),
                                  "_t0", num(p$t0), "_dpos", num(p$dpos)),
                    paste0("n", p$n)
      )
      paste0("accumulator_", mdl, "_", key, "_", ts, ".csv")
    },
    content = function(file){
      st <- simState()
      if (is.null(st)) {
        write.csv(data.frame(note = "No data yet. Click 'Run simulation' first."), file, row.names = FALSE)
      } else {
        write.csv(st$df, file, row.names = FALSE)
      }
    }
  )
  
  # ---------- helper plotting functions ------------------------------
  pathPlot <- function(df, model, params) {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    if (model == "lca") {
      par(mfrow = c(1, 2), mar = c(4,4,3,1))
      
      need <- key_lca(params)
      cache <- lcaCache()
      if (is.null(cache$key) || cache$key != need || is.null(cache$correct) || is.null(cache$error)) {
        # compute ONCE, store; aim for 1 path each to keep it snappy
        trC <- lca_one_trial_path_fast(params$v, params$s, params$b, params$a,
                                       params$kappa, params$beta, params$dt, params$maxt,
                                       target = "correct")
        trE <- lca_one_trial_path_fast(params$v, params$s, params$b, params$a,
                                       params$kappa, params$beta, params$dt, params$maxt,
                                       target = "error")
        lcaCache(list(key = need, correct = trC, error = trE))
        cache <- lcaCache()
      }
      
      draw_one <- function(tr, title_txt) {
        if (is.null(tr)) { plot.new(); title(paste(title_txt, "\n(no example)"), col.main="grey50"); return() }
        plot(NA, xlim = c(0, max(tr$t)), ylim = c(0, params$b),
             xlab = "Time (s)", ylab = "Activation", main = title_txt)
        abline(h = params$b, lty = 3)
        lines(tr$t, tr$x1, col = "steelblue", lwd = 2.0)
        lines(tr$t, tr$x2, col = "salmon",    lwd = 1.8)
        legend("topleft", c("Accumulator 1 (target)","Accumulator 2"),
               col = c("steelblue","salmon"), lwd = 2, bty = "n", cex = .9)
        abline(v = max(tr$t), lty = 2)
      }
      
      draw_one(cache$correct, "LCA dynamics — Correct example")
      draw_one(cache$error,   "LCA dynamics — Error example")
      return(invisible())
    }
    ok <- is.finite(df$RT) & df$RT > 0
    if (!any(ok)) {
      par(mfrow = c(1,1), mar = c(4,4,3,1))
      plot.new()
      title("No finite RTs.\nTry lower thresholds / more drift / higher noise.", col.main = "grey50")
      return(invisible())
    }
    xmax <- max(df$RT, na.rm = TRUE)
    
    if (model %in% c("ddm","wiener","ou")) {
      for (flag in c(TRUE, FALSE)) {
        plot(NA, xlim = c(0, xmax), ylim = c(0, params$a),
             xlab = "Time", ylab = "Evidence",
             main = paste(ifelse(flag, "Eg. Correct", "Eg. Error"), "paths"))
        abline(h = params$a, lty = 3)
        idx   <- which(df$correct == flag)
        nshow <- min(100L, length(idx))
        if (nshow == 0L) { text(xmax/2, params$a/2, "No trials", col = "grey50"); next }
        cols <- adjustcolor(ifelse(flag, "steelblue", "salmon"), 0.18)
        for (j in idx[seq_len(nshow)]) {
          if (model == "ou") {
            p <- ou_path(
              rt       = df$RT[j],
              a        = params$a,
              z        = df$start[j],
              to_upper = df$correct[j],
              theta    = params$theta,
              mu_prop  = params$mu,
              s        = params$s,
              dt       = max(0.0005, min(0.01, params$dt))
            )
          } else {
            p <- ddm_path(
              rt       = df$RT[j],
              a        = params$a,
              z        = df$start[j],
              to_upper = df$correct[j],
              v        = params$v,
              s        = params$s,
              dt       = 0.002
            )
          }
          if (!is.null(p)) lines(p$t, p$x, col = cols)
        }
      }
      
    } else if (model == "ddm_hcb") {
      a0 <- params$a0; amin <- params$amin; kh <- params$kh
      two <- isTRUE(params$two_sided)
      thr_fun <- function(t) pmax(amin, a0/(1 + kh*t))
      tline <- seq(0, xmax, length.out = 200)
      ylim <- if (two) c(-a0, a0) else c(0, a0)
      
      for (flag in c(TRUE, FALSE)) {
        plot(NA, xlim = c(0, xmax), ylim = ylim,
             xlab = "Time", ylab = "Evidence",
             main = paste(ifelse(flag, "Eg. Correct", "Eg. Error"),
                          "paths (hyperbolic Brownian)"))
        lines(tline, thr_fun(tline), lty = 3)
        if (two) lines(tline, -thr_fun(tline), lty = 3) else abline(h = 0, lty = 3)
        
        idx   <- which(df$correct == flag)
        nshow <- min(100L, length(idx))
        if (nshow > 0L) {
          cols <- adjustcolor(ifelse(flag, "steelblue", "salmon"), 0.20)
          for (j in idx[seq_len(nshow)]) {
            pth <- ddm_hcb_path(
              rt       = df$RT[j],
              a0       = a0,
              amin     = amin,
              k        = kh,
              z        = params$z,
              to_upper = isTRUE(df$correct[j]),
              v        = params$v,
              s        = params$s,
              dt       = params$dt,
              two_sided= two
            )
            if (!is.null(pth)) lines(pth$t, pth$x, col = cols)
          }
        } else {
          text(xmax/2, if (two) 0 else a0/2, "No trials", col = "grey50")
        }
      }
      
    } else if (model == "rw") {
      for (flag in c(TRUE, FALSE)) {
        plot(NA, xlim = c(0, xmax), ylim = c(0, params$a),
             xlab = "Time", ylab = "Evidence",
             main = paste(ifelse(flag, "Eg. Correct", "Eg. Error"), "paths (RW)"))
        abline(h = params$a, lty = 3)
        idx   <- which(df$correct == flag)
        nshow <- min(100L, length(idx))
        if (nshow == 0L) { text(xmax/2, params$a/2, "No trials", col = "grey50"); next }
        cols <- grDevices::adjustcolor(ifelse(flag, "steelblue", "salmon"), 0.18)
        for (j in idx[seq_len(nshow)]) {
          stps <- df$steps[j]
          if (is.null(stps) || is.na(stps)) stps <- max(1L, round(df$RT[j] / params$dt))
          p <- rw_path(
            rt       = df$RT[j],
            a        = params$a,
            z        = df$start[j],
            to_upper = isTRUE(df$correct[j]),
            h        = params$h,
            dt       = params$dt,
            steps    = stps
          )
          if (!is.null(p)) lines(p$t, p$x, col = cols)
        }
      }
      
    } else if (model %in% c("lba","lba_cb","ldliv")) {
      thr_fun <- NULL; thr <- NULL
      if (model %in% c("lba","lca","ldliv")) {
        thr <- if (model == "ldliv") params$b else params$b
      } else if (model == "lba_cb") {
        b0 <- params$b0; bmin <- params$bmin; kc <- params$kc
        thr_fun <- function(t) pmax(bmin, b0 - kc*t)
      }
      for (flag in c(TRUE, FALSE)) {
        if (is.null(thr_fun)) {
          plot(NA, xlim = c(0, xmax), ylim = c(0, thr),
               xlab = "Time", ylab = "Accumulated Evidence",
               main = paste(ifelse(flag, "Eg. Correct", "Eg. Error"), "paths"))
          abline(h = thr, lty = 3)
        } else {
          plot(NA, xlim = c(0, xmax), ylim = c(0, params$b0),
               xlab = "Time", ylab = "Accumulated Evidence",
               main = paste(ifelse(flag, "Eg. Correct", "Eg. Error"), "paths (collapsing)"))
          tline <- seq(0, xmax, length.out = 200)
          lines(tline, thr_fun(tline), lty = 3)
        }
        idx   <- which(df$correct == flag)
        nshow <- min(100L, length(idx))
        if (nshow > 0L) {
          for (j in idx[seq_len(nshow)]) {
            if (model == "ldliv") {
              # draw one simulated LDLIV path conditioned on RT
              pth <- ldliv_path_fast(
                rt  = df$RT[j],                  # t0 already included in df$RT? If so, pass (df$RT[j] - p$t0)
                b   = params$b,
                q   = if (flag) params$qC else params$qE,
                k   = params$k,
                sigma = params$sigma,
                z   = params$z * params$b,
                dt  = params$dt
              )
              if (!is.null(pth) && length(pth$t)) {
                lines(pth$t, pth$x, col = adjustcolor(ifelse(flag, "steelblue","salmon"), 0.15))
              }
            } else {
              rt <- df$RT[j]
              y1 <- if (is.null(thr_fun)) thr else thr_fun(rt)
              lines(c(0, rt), c(df$start[j], y1),
                    col = adjustcolor(ifelse(flag, "steelblue", "salmon"), 0.15))
            }
          }
        } else {
          ymid <- if (is.null(thr_fun)) if (is.null(thr)) 1 else thr/2 else max(thr_fun(xmax), 0)/2
          text(xmax/2, ymid, "No trials", col = "grey50")
        }
      }
    }
  }
  
  # Brownian bridge for fixed-bound DDM (paths only)
  ddm_path <- function(rt, a, z, to_upper, v, s, dt = 0.002) {
    if (is.na(rt) || rt <= 0) return(NULL)
    n  <- max(1L, ceiling(rt / dt))
    t  <- seq(0, rt, length.out = n + 1)
    drift <- v * t
    W     <- c(0, cumsum(rnorm(n, mean = 0, sd = s * sqrt(dt))))
    x     <- z + drift + W
    target <- if (to_upper) a else 0
    corr   <- (target - x[length(x)]) / rt
    x      <- x + corr * t
    x <- pmin(pmax(x, 0), a)
    data.frame(t = t, x = x)
  }
  
  # Brownian-like sample path under hyperbolic collapsing bounds (plotting only)
  ddm_hcb_path <- function(rt, a0, amin, k, z, to_upper, v, s,
                           dt = 0.001, two_sided = TRUE) {
    if (is.na(rt) || rt <= 0) return(NULL)
    n  <- max(1L, ceiling(rt / dt))
    t  <- seq(0, rt, length.out = n + 1)
    a_t <- pmax(amin, a0 / (1 + k * t))           # hyperbolic bound over time
    
    # starting point in absolute units
    y0 <- if (two_sided) (2 * z - 1) * a0 else a0 * z
    
    # endpoint pinned to the hit boundary at rt
    a_rt  <- a_t[length(a_t)]
    y_end <- if (two_sided) { if (isTRUE(to_upper))  a_rt else -a_rt } else { if (isTRUE(to_upper)) a_rt else 0 }
    
    # raw Brownian motion with drift
    W <- c(0, cumsum(rnorm(n, mean = 0, sd = s * sqrt(dt))))
    x <- y0 + v * t + W
    
    # bridge: pin to the known endpoint
    corr <- (y_end - x[length(x)]) / rt
    x <- x + corr * t
    
    # clip to the (time-varying) bounds to avoid fake pre-hits in the drawing
    if (two_sided) {
      x <- pmin(pmax(x, -a_t), a_t)
    } else {
      x <- pmin(pmax(x, 0), a_t)
    }
    
    data.frame(t = t, x = x, up = a_t, lo = if (two_sided) -a_t else rep(0, length(a_t)))
  }
  
  # Densities with optional parametric overlays for Bernoulli–Weibull & Poisson
  densPlot <- function(df, model = NULL, params = NULL){
    par(mfrow = c(1,2), mar = c(4,4,3,1))
    plot_one <- function(flag, col){
      ttl <- paste(ifelse(flag, "Correct", "Error"), "RT Density")
      idx <- (!is.na(df$correct)) & (df$correct == flag) & is.finite(df$RT)
      x   <- df$RT[idx]
      if (length(x) < 2L) {
        plot.new(); title(paste0(ttl, " (not enough data)"), col.main = "grey50"); return(invisible())
      }
      xr <- range(x); pad <- diff(xr) * 0.05
      if (!is.finite(sd(x)) || diff(xr) == 0) {
        pad <- max(1e-6, abs(xr[1]) * 1e-3)
        plot(NA, xlim = xr + c(-pad, pad), ylim = c(0, 1),
             xlab = "RT (s)", ylab = "Density", main = paste0(ttl, " (degenerate)"))
        abline(v = xr[1], col = col, lwd = 2); rug(x, col = col); return(invisible())
      }
      d <- tryCatch(stats::density(x), warning = function(w) NULL, error = function(e) NULL)
      if (is.null(d) || any(!is.finite(d$y))) {
        jit_amt <- max(1e-6, 0.001 * sd(x))
        d <- stats::density(jitter(x, factor = 1, amount = jit_amt))
      }
      plot(d, col = col, lwd = 2, main = ttl, xlab = "RT (s)")
      
      # ── Parametric overlays
      if (!is.null(model) && !is.null(params)) {
        xx <- seq(max(min(x) - 3*sd(x), 1e-6), max(x) + 3*sd(x), length.out = 400)
        if (model == "bern_weibull") {
          k  <- params$k; lam <- params$lambda; t0 <- params$t0
          # Shifted Weibull: f(t) = 0 for t<=t0; else f(t-t0; k, λ)
          ww <- ifelse(xx <= t0, 0, dweibull(xx - t0, shape = k, scale = lam))
          lines(xx, ww, lwd = 2, lty = 2)
          legend("topright", legend = c("Kernel density", "Weibull PDF"),
                 col = c(col, "black"), lwd = c(2,2), lty = c(1,2), bty = "n", cex = .9)
        } else if (model == "poisson") {
          N   <- max(1L, as.integer(round(params$N)))
          lc  <- max(1e-9, params$lamC)
          le  <- max(1e-9, params$lamE)
          t0  <- params$t0
          t   <- pmax(xx - t0, 0)
          # f_correct(t) = f_c(t) * S_e(t); f_error(t) = f_e(t) * S_c(t)
          f_c <- dgamma(t, shape = N, rate = lc) * (1 - pgamma(t, shape = N, rate = le))
          f_e <- dgamma(t, shape = N, rate = le) * (1 - pgamma(t, shape = N, rate = lc))
          if (isTRUE(flag)) {
            lines(xx, f_c, lwd = 2, lty = 2)
            legend("topright", legend = c("Kernel density", "Gamma-race PDF (correct)"),
                   col = c(col, "black"), lwd = c(2,2), lty = c(1,2), bty = "n", cex = .9)
          } else {
            lines(xx, f_e, lwd = 2, lty = 2)
            legend("topright", legend = c("Kernel density", "Gamma-race PDF (error)"),
                   col = c(col, "black"), lwd = c(2,2), lty = c(1,2), bty = "n", cex = .9)
          }
        }
      }
    }
    plot_one(TRUE,  "steelblue")
    plot_one(FALSE, "salmon")
  }
  
  polarErrPlot <- function(df, params){
    # works for cdm & scdm (both provide angle)
    if (!("angle" %in% names(df))) { plot.new(); title("No angles to plot", col.main="grey50"); return(invisible()) }
    target <- safe_scalar(params$theta, 0) * pi/180
    ang  <- circular(df$angle, units="radians")
    dens <- density.circular(ang, bw = 20)
    maxR <- max(1,dens$y)*1.2
    plot(dens, shrink=1, xlim=c(-maxR,maxR), ylim=c(-maxR,maxR),
         main="Angular-error density", col="darkorange", lwd=2)
    arrows(0,0, cos(target), sin(target), length=.08, lwd=2)
  }
  
  rtPolarPlot <- function(df, params){
    if (!("angle" %in% names(df))) { plot.new(); title("No angles to plot"); return(invisible()) }
    if (nrow(df) < 2) { plot.new(); title("Not enough data"); return(invisible()) }
    ang <- (df$angle + 2*pi) %% (2*pi)
    sector <- factor(cut(ang, breaks = seq(0, 2*pi, length.out = 13),
                         include.lowest = TRUE, labels = FALSE), levels = 1:12)
    mu <- tapply(df$RT, sector, mean, na.rm = TRUE)
    edges   <- seq(0, 2*pi, length.out = 13)
    centers <- edges[-length(edges)] + diff(edges)[1]/2
    theta   <- c(centers, centers[1]); r <- c(mu, mu[1])
    if (all(is.na(r))) { plot.new(); title("No RTs in any sector", col.main="grey50"); return(invisible()) }
    maxR <- max(r, na.rm = TRUE) * 1.15
    plot(0, 0, type = "n", xlim = c(-maxR, maxR), ylim = c(-maxR, maxR),
         asp = 1, axes = FALSE, main = "Mean RT by angle (12 bins)")
    rgrid <- pretty(c(0, maxR), 4); rgrid <- rgrid[rgrid > 0]
    if (length(rgrid)) {
      symbols(rep(0, length(rgrid)), rep(0, length(rgrid)),
              circles = rgrid, inches = FALSE, add = TRUE, fg = "grey90")
      text(rgrid, 0, labels = round(rgrid, 2), col = "grey60", pos = 4, cex = 0.7)
    }
    xs <- r * cos(theta); ys <- r * sin(theta)
    lines(xs, ys, col = "purple", lwd = 2); points(xs, ys, pch = 16, col = "purple")
    tar <- safe_scalar(params$theta, 0) * pi/180
    arrows(0, 0, maxR*cos(tar), maxR*sin(tar), length = 0.06, lwd = 2)
  }
  
  cdmPathPlot <- function(df, params){
    get_num <- function(x, default) {
      y <- suppressWarnings(as.numeric(x))
      if (length(y) != 1 || is.na(y) || !is.finite(y)) default else y
    }
    R     <- get_num(params[["R"]], 1)
    theta <- get_num(params[["theta"]], 0) * pi/180
    v     <- get_num(params[["v"]], 0.25)
    s     <- max(0, get_num(params[["s"]], 0.15))
    if (!is.finite(R) || R <= 0 || !is.finite(v) || !is.finite(s)) {
      plot.new(); title("Invalid CDM params. Click Run again.", col.main = "grey50"); return(invisible())
    }
    nshow  <- 50
    par(mar = c(1,1,3,1))
    plot(NA, xlim=c(-R*1.1, R*1.1), ylim=c(-R*1.1, R*1.1),
         asp=1, axes=FALSE, main=paste("Vector diffusion paths"))
    symbols(0,0, circles=R, inches=FALSE, add=TRUE, fg="grey60")
    ramp <- grDevices::colorRampPalette(c("steelblue","grey80"))
    dt <- 0.02
    for (i in seq_len(min(nshow, nrow(df)))) {
      x <- y <- 0; xs <- ys <- numeric(0)
      iter <- 0L; iter_max <- 1e6L
      while ((x^2 + y^2) < (R^2)) {
        if (iter >= iter_max) break
        iter <- iter + 1L
        xs <- c(xs, x); ys <- c(ys, y)
        x  <- x + v*dt*cos(theta) + rnorm(1, 0, s*sqrt(dt))
        y  <- y + v*dt*sin(theta) + rnorm(1, 0, s*sqrt(dt))
      }
      ang  <- atan2(y, x)
      delta<- ((ang - theta + pi) %% (2*pi)) - pi
      errA <- abs(delta)
      idx  <- min(101, max(1, floor(errA/pi*100)+1))
      col  <- grDevices::adjustcolor(ramp(101)[idx], 0.35)
      if (length(xs) > 1) lines(xs, ys, col = col)
    }
    points(0,0,pch=3)
  }
  
  # NEW: SCDM wave snapshots (demo of one trial)
  scdmWavePlot <- function(params){
    v <- params$v; s <- params$s; b <- params$b
    kappa <- params$kappa; M <- max(8L, as.integer(round(params$M)))
    theta <- params$theta * pi/180; dt <- params$dt
    ang <- seq(0, 2*pi, length.out = M + 1L); ang <- ang[-length(ang)]
    if (kappa > 0) {
      vm <- exp(kappa * cos((ang - theta))); f <- vm - mean(vm); f <- f / max(f)
    } else {
      f <- cos(ang - theta)
    }
    max_steps <- max(1L, ceiling(5/dt))  # just a short illustrative run
    y <- rep(0.0, M); t <- 0.0
    snaps <- list(); times <- c()
    # take ~5 snapshots until hit
    snap_idx <- unique(round(seq(1, max_steps, length.out = 5)))
    hit <- FALSE
    for (k in seq_len(max_steps)) {
      t <- t + dt
      y <- y + v * f * dt + s * sqrt(dt) * rnorm(M)
      if (k %in% snap_idx || max(y) >= b) {
        snaps[[length(snaps)+1]] <- pmin(pmax(y, 0), b)  # clamp for plotting
        times <- c(times, t)
      }
      if (max(y) >= b) { hit <- TRUE; break }
    }
    rmax <- b * 1.1
    par(mar=c(1,1,3,1))
    plot(NA, xlim=c(-rmax,rmax), ylim=c(-rmax,rmax), asp=1, axes=FALSE,
         main="SCDM: wave snapshots (crest→bound)")
    symbols(0,0, circles=b, inches=FALSE, add=TRUE, fg="grey80")
    cols <- grDevices::adjustcolor(colorRampPalette(c("grey70","steelblue"))(length(snaps)), 0.9)
    for (i in seq_along(snaps)) {
      r <- snaps[[i]]
      xs <- r * cos(ang); ys <- r * sin(ang)
      lines(xs, ys, lwd = 2, col = cols[i])
    }
    # target direction
    arrows(0,0, rmax*cos(theta), rmax*sin(theta), length=0.06, lwd=2)
    legend("bottomleft", legend = paste("t=", round(times,3), "s"), lwd = 2, col = cols, bty="n", cex=.8)
  }
  
  rw_path <- function(rt, a, z, to_upper, h, dt, steps = NULL) {
    if (is.na(rt) || rt <= 0) return(NULL)
    n  <- if (is.null(steps)) max(1L, round(rt / dt)) else as.integer(steps)
    t  <- (0:n) * dt
    target <- if (isTRUE(to_upper)) a else 0
    net    <- target - z
    k      <- round(net / h)
    U      <- (n + k) / 2
    D      <- n - U
    if (U < 0 || D < 0 || abs(U - round(U)) > 1e-8) {
      x <- z + (net / (n*dt)) * t
      x <- pmin(pmax(x, 0), a)
      return(data.frame(t = t, x = x))
    }
    U <- as.integer(round(U)); D <- as.integer(round(D))
    stepsigns <- sample(c(rep(1L, U), rep(-1L, D)))
    x <- numeric(n + 1); x[1] <- z
    for (i in 1:n) x[i + 1] <- x[i] + h * stepsigns[i]
    x <- pmin(pmax(x, 0), a)
    data.frame(t = t, x = x)
  }
  
  # Sample path generator for LDLIV (for plotting only; not used in simulation core)
  ldliv_path <- function(rt, b, q, k, sigma, z, dt = 0.001) {
    if (is.na(rt) || rt <= 0) return(NULL)
    n <- max(1L, ceiling(rt / dt))
    t <- seq(0, rt, length.out = n + 1L)
    x <- numeric(n + 1L); x[1] <- z
    for (i in 1:n) {
      xi <- max(0, x[i])
      x[i+1] <- xi + (q - k * xi) * dt + sqrt(2 * sigma * xi) * sqrt(dt) * rnorm(1)
      if (x[i+1] < 0) x[i+1] <- 0
      if (x[i+1] > b) x[i+1] <- b
    }
    data.frame(t = t, x = x)
  }
  
  # OU sample path (for plotting only; lightly "bridged" to the endpoint)
  ou_path <- function(rt, a, z, to_upper, theta, mu_prop, s, dt = 0.002) {
    if (is.na(rt) || rt <= 0) return(NULL)
    n  <- max(1L, ceiling(rt / dt))
    t  <- seq(0, rt, length.out = n + 1)
    mu <- a * min(max(mu_prop, 0), 1)
    
    x <- numeric(n + 1); x[1] <- z
    for (i in 1:n) {
      x[i+1] <- x[i] + theta * (mu - x[i]) * dt + s * sqrt(dt) * rnorm(1)
    }
    
    target <- if (isTRUE(to_upper)) a else 0
    corr   <- (target - x[length(x)]) / rt
    x      <- x + corr * t
    x      <- pmin(pmax(x, 0), a)
    data.frame(t = t, x = x)
  }
  
  # --- LCA single-trial path (for visualization only; not used in main sim) ---
  lca_one_trial_path <- function(v, s, b, a, kappa, beta, dt,
                                 max_t = 10, stall_eps = 1e-7, stall_run = 500L) {
    max_steps <- max(1L, ceiling(max_t / dt))
    
    x1 <- runif(1, 0, a); x2 <- runif(1, 0, a)
    X1 <- numeric(max_steps + 1L); X2 <- numeric(max_steps + 1L); TT <- numeric(max_steps + 1L)
    X1[1] <- x1; X2[1] <- x2; TT[1] <- 0
    
    stall <- 0L
    
    for (st in 2:(max_steps + 1L)) {
      # Kill switch: bail immediately if user pressed "Kill"
      if (isTRUE(cancelled)) {
        last <- st - 1L
        return(list(t = TT[1:last], x1 = X1[1:last], x2 = X2[1:last],
                    winner = NA_integer_, rt = NA_real_, hit = FALSE, cancelled = TRUE))
      }
      
      x1_prev <- x1; x2_prev <- x2
      x1 <- x1 + (v     - kappa * x1 - beta * x2) * dt + s * sqrt(dt) * rnorm(1)
      x2 <- x2 + ((1-v) - kappa * x2 - beta * x1) * dt + s * sqrt(dt) * rnorm(1)
      if (x1 < 0) x1 <- 0
      if (x2 < 0) x2 <- 0
      
      TT[st] <- TT[st-1] + dt
      X1[st] <- x1; X2[st] <- x2
      
      if (x1 >= b || x2 >= b) {
        winner <- if (x1 >= b && x2 >= b) { if (x1 >= x2) 1L else 2L } else if (x1 >= b) 1L else 2L
        last <- st
        return(list(t = TT[1:last], x1 = X1[1:last], x2 = X2[1:last],
                    winner = winner, rt = TT[last], hit = TRUE, cancelled = FALSE))
      }
      
      if (abs(x1 - x1_prev) + abs(x2 - x2_prev) < stall_eps) {
        stall <- stall + 1L
        if (stall >= stall_run) break
      } else {
        stall <- 0L
      }
    }
    
    last <- which.max(TT)  # end of simulated time
    list(t = TT[1:last], x1 = X1[1:last], x2 = X2[1:last],
         winner = NA_integer_, rt = NA_real_, hit = FALSE, cancelled = FALSE)
  }
  
  
  # --- LCA dynamics panels: correct vs error ---
  lca_traj_panels <- function(params, n_each = 6, max_t = NULL) {
    # Optional fixed-point precheck: skip drawing if clearly sub-threshold
    fp <- lca_fixed_point_check(params$v, params$kappa, params$beta, params$b)
    if (identical(fp$verdict, "no-hit-likely")) {
      plot.new()
      title("LCA dynamics — no-hit likely from fixed-point analysis", cex.main=.95)
      mtext(sprintf("x* ≈ (%.3g, %.3g); %s", fp$xstar[1], fp$xstar[2], fp$msg),
            side=3, line=0.5, col="grey40")
      return(invisible())
    }
    
    # Reproducible RNG per run (same trick you use in extraPlot3)
    old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
    on.exit({ if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv) }, add = TRUE)
    if (!is.null(params$seed_demo)) set.seed(params$seed_demo)
    
    ex_correct <- list(); ex_error <- list()
    tries <- 0L; cap <- 2000L
    maxt <- if (is.null(max_t)) if (is.null(params$maxt)) 5 else params$maxt else max_t
    
    while ((length(ex_correct) < n_each || length(ex_error) < n_each) && tries < cap) {
      tries <- tries + 1L
      tr <- lca_one_trial_path(v = params$v, s = params$s, b = params$b, a = params$a,
                               kappa = params$kappa, beta = params$beta, dt = params$dt,
                               max_t = maxt)
      if (is.null(tr) || !isTRUE(tr$hit)) next
      if (tr$winner == 1L && length(ex_correct) < n_each) {
        ex_correct[[length(ex_correct) + 1L]] <- tr
      } else if (tr$winner == 2L && length(ex_error) < n_each) {
        ex_error[[length(ex_error) + 1L]] <- tr
      }
    }
    
    par(mfrow = c(1, 2), mar = c(4,4,3,1))
    draw_panel <- function(ex_list, title_txt) {
      if (!length(ex_list)) { plot.new(); title(paste(title_txt, "\n(no examples)"), col.main="grey50"); return(invisible()) }
      xmax <- max(vapply(ex_list, function(z) max(z$t), numeric(1)))
      plot(NA, xlim = c(0, xmax), ylim = c(0, params$b),
           xlab = "Time (s)", ylab = "Activation", main = title_txt)
      abline(h = params$b, lty = 3)
      for (e in ex_list) {
        lines(e$t, e$x1, col = adjustcolor("steelblue", 0.85), lwd = 1.5)
        lines(e$t, e$x2, col = adjustcolor("salmon",    0.85), lwd = 1.2)
      }
      legend("topleft", c("Accumulator 1 (target)", "Accumulator 2 (competitor)"),
             col = c("steelblue", "salmon"), lwd = c(2,2), bty = "n", cex = .9)
    }
    draw_panel(ex_correct[seq_len(min(n_each, length(ex_correct)))], "LCA dynamics — Correct trials")
    draw_panel(ex_error[seq_len(min(n_each, length(ex_error)))],     "LCA dynamics — Error trials")
  }
  
  
  lca_draw_one <- function(params, want = c("any","correct","error"),
                           max_tries = 2000L, max_t = 10) {
    want <- match.arg(want)
    for (i in seq_len(max_tries)) {
      if (isTRUE(cancelled)) return(list(cancelled = TRUE))
      tr <- lca_one_trial_path(params$v, params$s, params$b, params$a,
                               params$kappa, params$beta, params$dt,
                               max_t = max_t)
      if (isTRUE(tr$cancelled)) return(tr)
      if (!tr$hit) next
      if (want == "any" ||
          (want == "correct" && tr$winner == 1L) ||
          (want == "error"   && tr$winner == 2L)) return(tr)
    }
    NULL
  }
  
  lca_one_trial_path_fast <- function(v, s, b, a, kappa, beta, dt, max_t,
                                      target = c("any","correct","error")) {
    target <- match.arg(target)
    # ~200 points per second is plenty smooth for teaching plots
    downsample_every <- max(1L, floor(0.005 / dt))
    
    for (i in 1:2000) {
      tr <- lca_one_trial_path_cpp(v, s, b, a, kappa, beta, dt, max_t, downsample_every)
      if (!isTRUE(tr$hit)) next
      if (target == "any" ||
          (target == "correct" && tr$winner == 1L) ||
          (target == "error"   && tr$winner == 2L)) {
        tr$t  <- as.numeric(tr$t)
        tr$x1 <- as.numeric(tr$x1)
        tr$x2 <- as.numeric(tr$x2)
        return(tr)
      }
    }
    NULL
  }
  
  
  # ────────────────────────────────────────────────────────────────────────────
  # MDFT helpers & simulator
  # Dynamics (per step): 
  #   x <- x + ( -phi * x  - beta * L %*% x  + U ) * dt  + s * sqrt(dt) * N(0, I)
  #   where U = (I - rho * S_off) %*% v  and  S_off[i!=j] = exp( -d_ij^2 / (2*sigma_s^2) )
  # Positions are on a 1D line with spacing Δpos; option 1 is the "target".
  mdft_build_mats <- function(M, dpos, rho, sigma_s){
    M <- max(2L, as.integer(M))
    pos <- seq(0, (M-1)*dpos, length.out = M)
    D   <- abs(outer(pos, pos, "-"))
    S_off <- exp(-(D^2) / (2 * sigma_s^2))
    diag(S_off) <- 0  # off-diagonal similarities only
    L <- matrix(1, M, M); diag(L) <- 0  # pure lateral structure (sum competitors)
    list(S_off = S_off, L = L)
  }
  
  mdft_one_trial_path <- function(M, b, phi, beta, rho, sigma_s,
                                  v1, vC, sv, s, dt, max_t, dpos){
    mats <- mdft_build_mats(M, dpos, rho, sigma_s)
    S_off <- mats$S_off; L <- mats$L
    
    # trial-specific valences (centered)
    v <- rep(vC, M); v[1] <- v1
    if (sv > 0) v <- v + rnorm(M, 0, sv)
    v <- v - mean(v)
    U <- (diag(M) - rho * S_off) %*% v
    
    nmax <- max(1L, ceiling(max_t / dt))
    X <- matrix(NA_real_, nrow = nmax + 1L, ncol = M)
    T <- numeric(nmax + 1L)
    x <- rep(0, M)     # preferences start at 0 (standard in DFT)
    X[1,] <- x; T[1] <- 0
    hit <- FALSE; winner <- NA_integer_; t <- 0
    
    for (k in 1:nmax){
      t <- t + dt
      x <- x + (-phi * x - beta * (L %*% x) + as.vector(U)) * dt + s * sqrt(dt) * rnorm(M)
      X[k+1,] <- x; T[k+1] <- t
      i <- which.max(x)
      if (x[i] >= b){
        hit <- TRUE; winner <- i
        X <- X[1:(k+1), , drop=FALSE]; T <- T[1:(k+1)]
        break
      }
    }
    list(hit = hit, rt = t, winner = winner, t = T, X = X)
  }
  
  simulate_mdft <- function(n, M, b, phi, beta, rho, sigma_s,
                            v1, vC, sv, s, dt, t0, max_t, dpos){
    M     <- max(2L, as.integer(M))
    b     <- max(1e-6, b)
    phi   <- max(0, phi)
    beta  <- max(0, beta)
    s     <- max(0, s)
    dt    <- max(1e-5, dt)
    max_t <- max(1e-3, max_t)
    dpos  <- max(1e-6, dpos)
    
    RT <- numeric(n); correct <- logical(n)
    for (i in seq_len(n)){
      tr <- mdft_one_trial_path(M, b, phi, beta, rho, sigma_s, v1, vC, sv, s, dt, max_t, dpos)
      if (isTRUE(tr$hit)){
        RT[i] <- t0 + tr$rt
        correct[i] <- (tr$winner == 1L)
      } else {
        RT[i] <- NA_real_; correct[i] <- NA
      }
    }
    data.frame(RT = RT, correct = correct)
  }
  
  # Single illustrative MDFT path (for students)
  mdft_path_example_plot <- function(params){
    M <- params$M; b <- params$b; phi <- params$phi; beta <- params$beta
    rho <- params$rho; sigma_s <- params$sigma_s; v1 <- params$v1
    vC <- params$vC; sv <- params$sv; s <- params$s; dt <- params$dt
    maxt <- params$maxt; dpos <- params$dpos
    
    # ~200 points/sec is smooth and light to draw
    ds <- max(1L, floor(0.005 / dt))
    tr <- mdft_one_trial_path_cpp(M, b, phi, beta, rho, sigma_s,
                                  v1, vC, sv, s, dt, maxt, dpos, ds)
    T  <- as.numeric(tr$t)
    X  <- tr$X
    if (length(T) == 0 || nrow(X) == 0) {
      plot.new(); title("No path samples (try larger max_t or smaller dt)", col.main="grey50")
      return(invisible())
    }
    M  <- ncol(X)
    
    cols <- c("steelblue", rep("grey55", max(0, M-1)))
    lwds <- c(2.4,        rep(1.6,       max(0, M-1)))
    
    ymin <- min(0, min(X))
    plot(NA, xlim = c(0, ifelse(length(T), max(T), maxt)), ylim = c(ymin, b),
         xlab = "Time (s)", ylab = "Preference",
         main = "MDFT: Single-trial preference dynamics (downsampled)")
    abline(h = b, lty = 3)
    
    for (j in seq_len(M)) lines(T, X[, j], col = cols[j], lwd = lwds[j])
    
    labs <- paste("Option", seq_len(M)); if (M >= 1) labs[1] <- "Option 1 (target)"
    legend("topleft", labs, col = cols, lwd = lwds, bty = "n", cex = .9)
    
    if (isTRUE(tr$hit)) {
      abline(v = T[length(T)], lty = 2)
      mtext(sprintf("Winner: %s   RT: %.3fs",
                    labs[suppressWarnings(as.integer(tr$winner))], tr$rt),
            side = 3, line = 0.2, cex = .9)
    } else {
      mtext("No bound hit (adjust v/φ/β/ρ/σs/b).", side = 3, line = 0.2, cex = .9, col = "grey40")
    }
  }
  
  
  ldliv_path_fast <- function(rt, b, q, k, sigma, z, dt) {
    # ~200 points per second is visually smooth
    ds <- max(1L, floor(0.005 / dt))
    ldliv_path_cpp(rt, b, q, k, sigma, z, dt, ds)
  }
  
  
}

shinyApp(ui, server)
