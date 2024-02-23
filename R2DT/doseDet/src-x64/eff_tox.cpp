#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
int efftox_add(NumericMatrix ep, NumericMatrix tp, double r, double pie, double pit,
             double ptt, double pet, double ptc, double pec, int h_treated) {
  int nrow = ep.nrow(), out=0;
  int max_dose_search = std::min(ep.ncol(), h_treated + 1);
  double u1= 1 - pow( pow( 1  / (1 - pie), r) + pow( 1 / pit, r), 1/r);
  double e_thres = nrow * (1 - pec);
  double t_thres = nrow * (1 - ptc);


  for (int j = 0; j < max_dose_search; j++) {
    double total = 0, m_pe=0, m_pt=0;
    int e_add=0, t_add=0; 
    for (int i = 0; i < nrow; i++) {

      m_pe += ep(i, j);
      m_pt += tp(i, j);
      e_add += ceil(pet - ep(i, j));
      t_add += ceil(tp(i, j) - ptt);
    }
    
    m_pe = m_pe / nrow;
    m_pt = m_pt / nrow;
    total = 1 - pow( pow( (1 - m_pe) / (1 - pie), r) + pow( m_pt / pit, r), 1/r);

     if ((total >= u1) && (e_add <= e_thres) && (t_add <= t_thres)) {
      out=j+1;
      u1=total;
    } 
    
  }
  return(out);
}


