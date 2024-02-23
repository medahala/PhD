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
NumericVector UdetC(NumericVector ep, double er, double el, double eag, double eal,
                    NumericVector tp, double tr, double tl, double tag, double tal, 
                    double k1, double k2) {
  int n = ep.size();
  NumericVector out(n);
  double A = -el * pow(er, eal);
  double B = pow(1 - er,eag) - A;
  double C = -tl * pow(abs(tr - 1), tal);
  double D = pow(tr,tag) - C;
  double ue;
  double ut;
  
  for(int i = 0; i < n; ++i) {
    
    if (ep[i] >= er) {
      ue = (pow((ep[i] - er),eag) - A)/B;
    } else {
      ue = ((-el * pow(abs(ep[i]-er),eal)) - A)/B;
    }
    if (tp[i] >= tr) {
      ut = ((-tl * pow(abs(tr - tp[i]),tal)) - C)/D;
    } else {
      ut = (pow((tr- tp[i]),tag) - C)/D;
    }
    
    out[i] = k1*ue  +  k2*ut + (1 - k1 -k2) *ue *ut;
  }
  return out;
}


