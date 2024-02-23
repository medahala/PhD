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
int prospect_pow_addC(NumericMatrix ep, double er, double el, double eag, double eal,
            NumericMatrix tp, double tr, double tl, double tag, double tal, 
            double k1, double k2, double ptt, double pet, double ptc, double pec, int h_treated) {
  int nrow = ep.nrow(), ndose=ep.ncol(), out=0;
  double u1= 0;
  double e_thres = nrow * (1 - pec);
  double t_thres = nrow * (1 - ptc);
  double A = -el * pow(er, eal);
  double B = pow(1 - er,eag) - A;
  double C = -tl * pow(abs(tr - 1), tal);
  double D = pow(tr,tag) - C;
  double ue;
  double ut;

  for (int j = 0; j < ndose; j++) {
    double total = 0;
    int e_add=0, t_add=0; 
    for (int i = 0; i < nrow; i++) {

            if (ep(i, j) >= er) {
            ue = (pow((ep(i, j) - er),eag) - A)/B;
            } else {
            ue = ((-el * pow(abs(ep(i, j)-er),eal)) - A)/B;
            }
            if (tp(i, j) >= tr) {
            ut = ((-tl * pow(abs(tr - tp(i, j)),tal)) - C)/D;
            } else {
            ut = (pow((tr- tp(i, j)),tag) - C)/D;
            }

      total += k1*ue  +  k2*ut + (1 - k1 -k2) *ue *ut;
      e_add += ceil(pet - ep(i, j));
      t_add += ceil(tp(i, j) - ptt);
    }
    
    if ((total >= u1) && (e_add <= e_thres) && (t_add <= t_thres)) {
      out=j+1;
      u1=total;
    } 
    
  }
 return(std::min(out, h_treated + 1));
}

