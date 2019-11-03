#include <RcppArmadillo.h>
#include "saem_class_rcpp.hpp"


using namespace std;
using namespace arma;

extern "C" {
  typedef void (*ode_solver_c)(int *neq, double *theta, double *time, int *evid,
                               int *ntime, double *inits, double *dose,
                               double *ret, double *atol, double *rtol,
                               int *stiff, int *transit_abs, int *nlhs,
                               double *lhs, int *rc);

  void rx_c133fd76d01c28263281c68fa73d8fa6__ode_solver(int *neq, double *theta, double *time, int *evid, int *ntime,
                       double *inits, double *dose, double *ret, double *atol,
                       double *rtol, int *stiff, int *transit_abs, int *nlhs,
                       double *lhs, int *rc){
    static ode_solver_c fun=NULL;
    if (fun == NULL) fun = (ode_solver_c) R_GetCCallable("rx_c133fd76d01c28263281c68fa73d8fa6_","rx_c133fd76d01c28263281c68fa73d8fa6__ode_solver");
    fun(neq, theta, time, evid, ntime, inits, dose, ret, atol,
        rtol, stiff, transit_abs, nlhs, lhs, rc);
  }

}

Function ff("sd");

RObject mat2NumMat(const mat &m) {
	RObject x = wrap( m.memptr() , m.memptr() + m.n_elem ) ;
	x.attr( "dim" ) = Dimension( m.n_rows, m.n_cols ) ;
	return x;
}

vec Ruser_function(const mat &phi_, const mat &evt_, const List &opt) {
  RObject phi, evt;
  phi = mat2NumMat(phi_);
  evt = mat2NumMat(evt_);
  NumericVector g;
  g = ff(phi, evt);
  vec yp(g);

  return yp;
}


vec user_function(const mat &phi, const mat &evt, const List &opt) {
  uvec ix;
  vec id = evt.col(0);
  mat wm;
  vec wv;

  ix = find(evt.col(2) == 0);
  vec yp(ix.n_elem);
  double *p=yp.memptr();
  int N=id.max()+1;



  for (int i=0; i<N; i++) {


    wm = evt.rows( find(id == i) );

    vec time__;
    time__ = wm.col(1);
    int ntime = time__.n_elem;
    wv = wm.col(2);
    ivec evid(ntime);
    for (int k=0; k<ntime; ++k) evid(k) = wv(k);
    vec amt;
    amt = wm.col(3);
    amt = amt( find(evid > 0) );

    int neq=as<int>(opt["neq"]);
    vec inits(neq);
    inits.zeros(); //as<vec>(opt["inits"]);	//FIXME

  double
  ka,
  tka,
  cl,
  tcl,
  v2,
  tv2,
  v3,
  tv3,
  q,
  tq;

  tka = phi(i, 0);
  tcl = phi(i, 1);
  tv2 = phi(i, 2);
  tv3 = phi(i, 3);
  tq = phi(i, 4);

  ka = exp( tka);
  cl = exp( tcl);
  v2 = exp( tv2);
  v3 = exp( tv3);
  q = exp( tq);

    vec params(5);
    params(0) = ka;
    params(1) = cl;
    params(2) = v2;
    params(3) = q;
    params(4) = v3;



    int stiff=as<int>(opt["stiff"]);
    int transit_abs=as<int>(opt["transitAbs"]);
    int nlhs=as<int>(opt["nlhs"]);
    double atol=as<double>(opt["atol"]);
    double rtol=as<double>(opt["rtol"]);
    int rc=0;

    mat ret(neq, ntime);
    mat lhs(nlhs, ntime);
	rx_c133fd76d01c28263281c68fa73d8fa6__ode_solver(&neq, params.memptr(), time__.memptr(),
	    evid.memptr(), &ntime, inits.memptr(), amt.memptr(), ret.memptr(),
	    &atol, &rtol, &stiff, &transit_abs, &nlhs, lhs.memptr(), &rc);
	ret = join_cols(join_cols(time__.t(), ret), lhs).t();
	ret = ret.rows( find(evid == 0) );

vec time;
time=ret.col(0);
vec depot;
depot=ret.col(1);
vec center;
center=ret.col(2);
vec periph;
periph=ret.col(3);
vec cp;
cp=ret.col(4);
vec nlmixr_pred;
nlmixr_pred=ret.col(5);


vec g;
g = nlmixr_pred;
if (g.has_nan()) {
	Rcout <<"WARNING: NAN in pred" << endl;
	Rcout  << i << params << endl;
	Rcout << join_rows(time,g) << endl;
	Rcout << "Consider to relax atol & rtol" << endl;
	g.replace(datum::nan, 1.0e9);
}

    int no = g.n_elem;
    memcpy(p, g.memptr(), no*sizeof(double));
    p += no;
  }

  return yp;
}

// definition
RcppExport SEXP dopred( SEXP in_phi, SEXP in_evt, SEXP in_opt ) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< mat& >::type phi(in_phi);
    Rcpp::traits::input_parameter< mat& >::type evt(in_evt);
    List opt(in_opt);

    vec g = user_function(phi, evt, opt);
    return Rcpp::wrap(g);
END_RCPP
}

RcppExport SEXP saem_fit(SEXP xSEXP) {
BEGIN_RCPP
  List x(xSEXP);

  SAEM saem;
  saem.inits(x);

  if(x.containsElementNamed("Rfn")) {
    ff = as<Function>(x["Rfn"]);
    saem.set_fn(Ruser_function);
  } else {
    saem.set_fn(user_function);
  }

  saem.saem_fit();

  List out = List::create(
    Named("mpost_phi") = saem.get_mpost_phi(),
    Named("Gamma2_phi1") = saem.get_Gamma2_phi1(),
    Named("Plambda") = saem.get_Plambda(),
    Named("Ha") = saem.get_Ha(),
    Named("sig2") = saem.get_sig2(),
    Named("eta") = saem.get_eta(),
    Named("par_hist") = saem.get_par_hist()
  );
  out.attr("saem.cfg") = x;
  out.attr("class") = "saemFit";
  return out;
END_RCPP
}

