#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double ppexpV(double q, const arma::vec& x, const arma::vec& cuts) {

  int nc = cuts.n_elem;
  double ret;
  double mi = nc;
  arma::vec cp1 = ones<arma::vec>(1);
  arma::vec cp0 = zeros<arma::vec>(1);

  // Find cut point where q >= cuts
  double ind;
  arma::vec indcut(nc);

  for(int i = 0; i < nc; i++){
    indcut(i) = q >= cuts(i);
  }
  ind = sum(indcut) -1;

  // Compute cdf at cuts[q]
  ret = R::pexp(q - cuts(ind), 1/x(ind), TRUE, FALSE);

  if( mi > ind){
    mi = ind;
  }

  if (nc > 1) {
    arma::vec dt = cuts(find(cuts>0)) - cuts(find(cuts != cuts(mi)));

    arma::vec xc = x(find(cuts != cuts(mi)));

    int nc1 = nc-1;
    arma::vec pe(nc1);
    for(int i = 0; i < nc1; i++){
      pe(i) = R::pexp(dt(i), 1/xc(i), TRUE, FALSE);;
    }
    arma::vec cp = arma::cumprod(1 - pe);
    cp = arma::join_vert(cp1,cp);

    int ncp = cp.n_elem;
    arma::vec ret0 = arma::cumsum(cp.head(ncp-1)%pe);
    ret0 = arma::join_vert(cp0, ret0);
    double ret1 = ret0(ind) + ret*cp(ind);
    return(ret1);
  } else{
    return(ret);
  }
}


// [[Rcpp::export]]
arma::vec ppexpM(double q, const arma::mat& x, const arma::vec& cuts) {

  arma::mat y = x.t();

  int ny = y.n_cols;
  arma::vec ret(ny);
  for(int i = 0; i < ny; i++){
    ret(i) = ppexpV(q, y.col(i), cuts);
  }
  return(ret);
}


// [[Rcpp::export]]
arma::mat ppexpMV(const arma::vec& q, const arma::mat& x, const arma::vec& cuts) {

  // Transpose the hazard matrix once, then evaluate the cdf at every query
  // time. Returns a matrix with one row per query time and one column per
  // hazard draw (i.e. nrow(q) x nrow(x)).
  arma::mat y = x.t();

  int nq = q.n_elem;
  int ny = y.n_cols;
  arma::mat ret(nq, ny);
  for(int j = 0; j < ny; j++){
    arma::vec xj = y.col(j);
    for(int i = 0; i < nq; i++){
      ret(i, j) = ppexpV(q(i), xj, cuts);
    }
  }
  return(ret);
}
