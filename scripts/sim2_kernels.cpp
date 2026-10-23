#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
// [[Rcpp::plugins(cpp17)]]

namespace sim2 {
struct Moments { double mean, variance, logz; };
inline double logadd(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  return a > b ? a + std::log1p(std::exp(b-a)) : b + std::log1p(std::exp(a-b));
}
inline double logistic(double a) {
  if (a >= 0) return 1.0/(1.0+std::exp(-a));
  const double e = std::exp(a); return e/(1.0+e);
}
struct Table {
  int n1, n0, d, x, lo, hi;
  double logz0, logn1, logn0, logratio;
  std::vector<double> logw;
  Table(int b, int a, int events, int treatment) : n1(b), n0(a), d(events), x(treatment),
      lo(std::max(0, events-a)), hi(std::min(b, events)), logz0(0.0),
      logn1(b>0?std::log(double(b)):R_NegInf),
      logn0(a>0?std::log(double(a)):R_NegInf),logratio(logn1-logn0) {
    if (b < 0 || a < 0 || events < 0 || events > b+a || treatment < lo || treatment > hi)
      Rcpp::stop("Invalid fixed-margin allocation");
    if (hi > lo && d > 1) {
      logz0 = R::lchoose(b+a, events);
      logw.resize(hi-lo+1);
      logw[0] = R::lchoose(b,lo)+R::lchoose(a,events-lo);
      for (int y=lo; y<hi; ++y) logw[y-lo+1] = logw[y-lo] +
          std::log(double(b-y))-std::log(double(y+1)) +
          std::log(double(events-y))-std::log(double(a-events+y+1));
    }
  }
  bool informative() const { return hi > lo; }
  Moments moments(double eta, bool need_logz=false) const {
    if (!informative()) return {double(lo), 0.0, 0.0};
    if (d == 1) {
      double logodds=eta+logratio;
      double p=logistic(logodds);
      return {p, p*(1-p), need_logz?logadd(logn0,eta+logn1):0.0};
    }
    double maxterm=R_NegInf;
    for (int y=lo; y<=hi; ++y) maxterm=std::max(maxterm,logw[y-lo]+eta*y);
    long double s=0, m=0, v=0;
    // Center at lo; variance of the centered support has smaller cancellation.
    for (int y=lo; y<=hi; ++y) {
      const double p=std::exp(logw[y-lo]+eta*y-maxterm);
      const double z=y-lo; s+=p; m+=p*z; v+=p*z*z;
    }
    const double mm=double(m/s);
    return {lo+mm, std::max(0.0,double(v/s)-mm*mm), maxterm+std::log(double(s))};
  }
  double loge(double eta) const {
    if (!informative() || eta == 0) return 0.0;
    const Moments m=moments(eta,true);
    return eta*x+(d==1 ? std::log(double(n1+n0)) : logz0)-m.logz;
  }
};
struct Derivatives { double score, information; };
Derivatives prior(double eta, int a, int b) {
  double p1=logistic(eta+std::log(double(b+1))-std::log(double(a+1)));
  double p2=a==0 ? 1.0 : logistic(eta+std::log(double(b+1))-std::log(double(a)));
  return {1.0-p1-p2,p1*(1-p1)+p2*(1-p2)};
}
Derivatives derivatives(const std::vector<Table>& history, double eta, int a, int b) {
  Derivatives p=prior(eta,a,b);
  long double score=p.score, info=p.information;
  // In exact-event histories theta is shared by the entire likelihood scan.
  // Compute its exponential once; use log odds only for extreme eta values.
  const bool ordinary_eta=(eta>-500.0 && eta<500.0);
  const double theta=ordinary_eta?std::exp(eta):0.0;
  for (const auto& tab:history) {
    if (tab.d==1 && ordinary_eta) {
      const double q=theta*tab.n1/(tab.n0+theta*tab.n1);
      score+=tab.x-q;info+=q*(1-q);
    } else {
      const Moments m=tab.moments(eta); score+=tab.x-m.mean; info+=m.variance;
    }
  }
  return {double(score),double(info)};
}
// Strictly concave smoothed log likelihood; no parameter clipping. Damping
// bounds the Newton STEP, and is not a bound on the fitted log odds.
double fit(const std::vector<Table>& history, int a, int b, double eta, Derivatives& state) {
  double lower=R_NegInf, upper=R_PosInf;
  for (int iter=0; iter<100; ++iter) {
    if (std::abs(state.score)<1e-10) return eta;
    if (state.score>0) lower=eta; else upper=eta;
    if (!(state.information>0) || !std::isfinite(state.information)) Rcpp::stop("Prequential MLE has zero information");
    double step=std::max(-2.0,std::min(2.0,state.score/state.information));
    double next=eta+step;
    if (next<=lower || next>=upper) next=(lower+upper)/2;
    if (!std::isfinite(next)) Rcpp::stop("Nonfinite prequential MLE step");
    if (std::abs(next-eta)<1e-12) return eta;
    eta=next; state=derivatives(history,eta,a,b);
  }
  Rcpp::stop("Prequential smoothed MLE failed to converge");
  return eta;
}
struct Process {
  Rcpp::LogicalVector reject;
  Rcpp::NumericVector stop, terminal, maximum;
  Rcpp::NumericMatrix path, eta_path;
  Process(int nr,int nc,double horizon,bool keep): reject(nr,false),stop(nr,horizon),
    terminal(nr,0.0),maximum(nr,0.0),path(keep?nr:0,keep?nc:0),eta_path(keep?nr:0,keep?nc:0) {
    if (keep) std::fill(eta_path.begin(),eta_path.end(),NA_REAL);
  }
  void release(int row,int start,int end,double increment,double eta,double time,double boundary,bool keep) {
    terminal[row]+=increment;
    maximum[row]=std::max(maximum[row],terminal[row]);
    if (!reject[row] && terminal[row]>=boundary) { reject[row]=true;stop[row]=time; }
    if (keep) for (int col=start;col<end;++col) { path(row,col)=terminal[row];eta_path(row,col)=eta; }
  }
  Rcpp::List result(bool keep) {
    Rcpp::List out=Rcpp::List::create(Rcpp::_ ["reject"]=reject,Rcpp::_ ["stop_time"]=stop,
      Rcpp::_ ["terminal_loge"]=terminal,Rcpp::_ ["maximum_loge"]=maximum);
    if (keep) { out["loge_path"]=path;out["eta_path"]=eta_path; }
    return out;
  }
};
}

// [[Rcpp::export]]
Rcpp::List sim2_comparators_cpp(const Rcpp::IntegerMatrix& n1,const Rcpp::IntegerMatrix& n0,
    const Rcpp::IntegerMatrix& d,const Rcpp::IntegerMatrix& x,const Rcpp::NumericMatrix& time,
    const Rcpp::IntegerVector& count,const Rcpp::IntegerVector& initial_n1,
    const Rcpp::IntegerVector& initial_n0,double horizon,double fixed_eta,double alpha,bool keep_path=false) {
  const int nr=n1.nrow(),nc=n1.ncol(); const double boundary=-std::log(alpha);
  sim2::Process fixed(nr,nc,horizon,keep_path), preq(nr,nc,horizon,keep_path), protective(nr,nc,horizon,keep_path);
  for (int r=0;r<nr;++r) {
    if (r%16==0) Rcpp::checkUserInterrupt();
    if (count[r]<0 || count[r]>nc) Rcpp::stop("Invalid update count");
    std::vector<sim2::Table> history; history.reserve(count[r]);
    const int a=initial_n0[r],b=initial_n1[r];
    const bool one_group=(a==0 || b==0);
    double eta=one_group ? 0.0 : 0.5*(std::log(double(a))+std::log(double(a+1)))-std::log(double(b+1));
    sim2::Derivatives state=sim2::prior(eta,a,b);
    int col=0;
    while (col<count[r]) {
      const double release=time(r,col);
      if (!std::isfinite(release) || release<0 || release>horizon || (col>0 && release<time(r,col-1)))
        Rcpp::stop("Release times must be finite, ordered, and within horizon");
      int end=col+1; while (end<count[r] && time(r,end)==release) ++end;
      if (!one_group) eta=sim2::fit(history,a,b,eta,state);
      const double cap_eta=std::min(eta,0.0);
      double fi=0,pi=0,ci=0;
      std::vector<sim2::Table> current; current.reserve(end-col);
      for (int c=col;c<end;++c) {
        current.emplace_back(n1(r,c),n0(r,c),d(r,c),x(r,c));
        const auto& tab=current.back();
        fi+=tab.loge(fixed_eta); const double value=tab.loge(eta);pi+=value;
        if (eta<=0) ci+=value;
      }
      fixed.release(r,col,end,fi,fixed_eta,release,boundary,keep_path);
      preq.release(r,col,end,pi,eta,release,boundary,keep_path);
      protective.release(r,col,end,ci,cap_eta,release,boundary,keep_path);
      // Only AFTER all current predictions and decisions enter the allocation
      // cells into history. This prevents leakage across simultaneous cells.
      for (auto& tab:current) if (tab.informative()) {
        const sim2::Moments m=tab.moments(eta);
        state.score+=tab.x-m.mean;state.information+=m.variance;
        history.push_back(std::move(tab));
      }
      col=end;
    }
    if (keep_path) for (int c=count[r];c<nc;++c) {
      fixed.path(r,c)=fixed.terminal[r];preq.path(r,c)=preq.terminal[r];protective.path(r,c)=protective.terminal[r];
    }
  }
  return Rcpp::List::create(Rcpp::_ ["fixed"]=fixed.result(keep_path),
    Rcpp::_ ["prequential"]=preq.result(keep_path),Rcpp::_ ["prequential_protective"]=protective.result(keep_path));
}
