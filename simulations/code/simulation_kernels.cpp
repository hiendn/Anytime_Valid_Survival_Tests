#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::plugins(cpp17)]]

namespace {

inline double log_add_exp(const double a, const double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  if (a > b) return a + std::log1p(std::exp(b - a));
  return b + std::log1p(std::exp(a - b));
}

inline double log_choose(const int n, const int k) {
  if (k < 0 || k > n) return R_NegInf;
  return std::lgamma(static_cast<double>(n) + 1.0) -
    std::lgamma(static_cast<double>(k) + 1.0) -
    std::lgamma(static_cast<double>(n - k) + 1.0);
}

inline double single_table_loge(
    const int n1,
    const int n0,
    const int d,
    const int x,
    const double u,
    const double tol) {
  const int n = n1 + n0;
  if (d <= 0 || d >= n || n1 <= 0 || n0 <= 0 ||
      !std::isfinite(u) || std::abs(u) <= tol) {
    return 0.0;
  }

  const int lo = std::max(0, d - n0);
  const int hi = std::min(n1, d);
  int r = lo;
  double log_term = log_choose(n1, r) + log_choose(n0, d - r) + u * r;
  double log_z = log_term;
  while (r < hi) {
    log_term +=
      std::log(static_cast<double>(n1 - r)) -
      std::log(static_cast<double>(r + 1)) +
      std::log(static_cast<double>(d - r)) -
      std::log(static_cast<double>(n0 - d + r + 1)) + u;
    log_z = log_add_exp(log_z, log_term);
    ++r;
  }
  return u * x + log_choose(n, d) - log_z;
}

}  // namespace


// [[Rcpp::export]]
Rcpp::NumericVector table_loge_cpp(
    const Rcpp::IntegerVector& n1,
    const Rcpp::IntegerVector& n0,
    const Rcpp::IntegerVector& d,
    const Rcpp::IntegerVector& x,
    const Rcpp::NumericVector& u,
    const double tol = 1e-15) {
  const R_xlen_t size = n1.size();
  if (n0.size() != size || d.size() != size || x.size() != size || u.size() != size) {
    Rcpp::stop("n1, n0, d, x, and u must have identical lengths");
  }
  Rcpp::NumericVector out(size);
  for (R_xlen_t i = 0; i < size; ++i) {
    out[i] = single_table_loge(n1[i], n0[i], d[i], x[i], u[i], tol);
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericVector calibration_curve_cpp(
    const Rcpp::IntegerVector& n1,
    const Rcpp::IntegerVector& n0,
    const Rcpp::IntegerVector& d,
    const Rcpp::IntegerVector& x,
    const Rcpp::NumericVector& weight,
    const Rcpp::NumericVector& eta_grid,
    const int reps,
    const double tol = 1e-15) {
  const R_xlen_t size = n1.size();
  if (n0.size() != size || d.size() != size || x.size() != size || weight.size() != size) {
    Rcpp::stop("calibration inputs must have identical lengths");
  }
  if (reps <= 0) Rcpp::stop("reps must be positive");

  Rcpp::NumericVector out(eta_grid.size());
  for (R_xlen_t j = 0; j < eta_grid.size(); ++j) {
    const double eta = eta_grid[j];
    long double total = 0.0L;
    for (R_xlen_t i = 0; i < size; ++i) {
      total += single_table_loge(
        n1[i], n0[i], d[i], x[i], eta * weight[i], tol
      );
    }
    out[j] = static_cast<double>(total / static_cast<long double>(reps));
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix row_cumsums_cpp(const Rcpp::NumericMatrix& input) {
  const int nr = input.nrow();
  const int nc = input.ncol();
  Rcpp::NumericMatrix out(nr, nc);
  for (int r = 0; r < nr; ++r) {
    double running = 0.0;
    for (int c = 0; c < nc; ++c) {
      running += input(r, c);
      out(r, c) = running;
    }
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::List build_updates_cpp(
    const Rcpp::IntegerMatrix& group,
    const Rcpp::NumericMatrix& event_time,
    const Rcpp::IntegerMatrix& last_visit,
    const Rcpp::NumericMatrix& regime_uniform,
    const Rcpp::IntegerMatrix& exit_order,
    const double rho,
    const int horizon,
    const int method_code,
    const double plan_cure,
    const double plan_rate) {
  const int reps = group.nrow();
  const int n = group.ncol();
  if (event_time.nrow() != reps || event_time.ncol() != n ||
      last_visit.nrow() != reps || last_visit.ncol() != n ||
      regime_uniform.nrow() != reps || regime_uniform.ncol() != n ||
      exit_order.nrow() != reps || exit_order.ncol() != n) {
    Rcpp::stop("latent-data matrices must have identical dimensions");
  }
  if (horizon <= 0) Rcpp::stop("horizon must be positive");
  if (method_code < 0 || method_code > 3) Rcpp::stop("invalid method code");

  const int max_updates = n + horizon;
  Rcpp::IntegerMatrix n1_out(reps, max_updates);
  Rcpp::IntegerMatrix n0_out(reps, max_updates);
  Rcpp::IntegerMatrix d_out(reps, max_updates);
  Rcpp::IntegerMatrix x_out(reps, max_updates);
  Rcpp::NumericMatrix time_out(reps, max_updates);
  Rcpp::NumericMatrix weight_out(reps, max_updates);
  Rcpp::IntegerVector count_out(reps);
  std::fill(time_out.begin(), time_out.end(), R_PosInf);
  std::fill(weight_out.begin(), weight_out.end(), 1.0);

  for (int rep = 0; rep < reps; ++rep) {
    std::vector<int> exact_n1(n, 0), exact_n0(n, 0), exact_d(n, 0), exact_x(n, 0);
    std::vector<double> exact_time(n, 0.0), exact_weight(n, 1.0);
    std::vector<int> panel_n1(horizon, 0), panel_n0(horizon, 0),
      panel_d(horizon, 0), panel_x(horizon, 0);
    std::vector<double> panel_time(horizon, 0.0), panel_weight(horizon, 1.0);
    std::vector<unsigned char> include_exact(n, 0), include_panel(n, 0);
    int n_exact_updates = 0;
    int n_panel_updates = 0;

    for (int i = 0; i < n; ++i) {
      const bool is_exact = regime_uniform(rep, i) < rho;
      if (method_code == 1) {
        include_panel[i] = 1;
      } else if (method_code == 2) {
        if (is_exact) include_exact[i] = 1;
      } else if (method_code == 3) {
        include_exact[i] = 1;
      } else if (is_exact) {
        include_exact[i] = 1;
      } else {
        include_panel[i] = 1;
      }
    }

    int risk1 = 0;
    int risk0 = 0;
    for (int i = 0; i < n; ++i) {
      if (!include_exact[i]) continue;
      if (group(rep, i) == 1) ++risk1;
      else ++risk0;
    }

    if (risk1 + risk0 > 0) {
      for (int pos = 0; pos < n; ++pos) {
        const int i = exit_order(rep, pos) - 1;  // R indices are one-based.
        if (i < 0 || i >= n) Rcpp::stop("exit_order contains an invalid index");
        if (!include_exact[i]) continue;
        const double ti = event_time(rep, i);
        const double ci = static_cast<double>(last_visit(rep, i));
        if (ti <= ci) {
          exact_n1[n_exact_updates] = risk1;
          exact_n0[n_exact_updates] = risk0;
          exact_d[n_exact_updates] = 1;
          exact_x[n_exact_updates] = group(rep, i) == 1 ? 1 : 0;
          exact_time[n_exact_updates] = ti;
          const double g_plan = plan_cure + (1.0 - plan_cure) * std::exp(-plan_rate * ti);
          exact_weight[n_exact_updates] = 1.0 / std::max(g_plan, 0.10);
          ++n_exact_updates;
        }
        if (group(rep, i) == 1) --risk1;
        else --risk0;
      }
    }

    for (int k = 1; k <= horizon; ++k) {
      int r1 = 0;
      int r0 = 0;
      int events = 0;
      int treatment_events = 0;
      const double left = static_cast<double>(k - 1);
      const double right = static_cast<double>(k);
      for (int i = 0; i < n; ++i) {
        if (!include_panel[i]) continue;
        if (last_visit(rep, i) >= k && event_time(rep, i) > left) {
          if (group(rep, i) == 1) ++r1;
          else ++r0;
          if (event_time(rep, i) <= right) {
            ++events;
            if (group(rep, i) == 1) ++treatment_events;
          }
        }
      }
      if (events > 0 && r1 + r0 > 0) {
        panel_n1[n_panel_updates] = r1;
        panel_n0[n_panel_updates] = r0;
        panel_d[n_panel_updates] = events;
        panel_x[n_panel_updates] = treatment_events;
        panel_time[n_panel_updates] = right;
        const double g_plan = plan_cure + (1.0 - plan_cure) * std::exp(-plan_rate * right);
        panel_weight[n_panel_updates] = 1.0 / std::max(g_plan, 0.10);
        ++n_panel_updates;
      }
    }

    int a = 0;
    int b = 0;
    int h = 0;
    while (a < n_exact_updates || b < n_panel_updates) {
      const bool take_exact =
        b >= n_panel_updates ||
        (a < n_exact_updates && exact_time[a] < panel_time[b]);
      if (take_exact) {
        n1_out(rep, h) = exact_n1[a];
        n0_out(rep, h) = exact_n0[a];
        d_out(rep, h) = exact_d[a];
        x_out(rep, h) = exact_x[a];
        time_out(rep, h) = exact_time[a];
        weight_out(rep, h) = exact_weight[a];
        ++a;
      } else {
        n1_out(rep, h) = panel_n1[b];
        n0_out(rep, h) = panel_n0[b];
        d_out(rep, h) = panel_d[b];
        x_out(rep, h) = panel_x[b];
        time_out(rep, h) = panel_time[b];
        weight_out(rep, h) = panel_weight[b];
        ++b;
      }
      ++h;
    }
    count_out[rep] = h;
  }

  return Rcpp::List::create(
    Rcpp::Named("n1") = n1_out,
    Rcpp::Named("n0") = n0_out,
    Rcpp::Named("d") = d_out,
    Rcpp::Named("x") = x_out,
    Rcpp::Named("time") = time_out,
    Rcpp::Named("gt_weight") = weight_out,
    Rcpp::Named("count") = count_out
  );
}
