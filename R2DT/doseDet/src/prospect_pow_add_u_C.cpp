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
int prospect_pow_add_u_C(NumericMatrix ep, double er, double el, double eag, double eal,
            NumericMatrix tp, double tr, double tl, double tag, double tal, 
            double k1, double k2, double ustop, double pstop, int h_treated) {
  int nrow = ep.nrow(), ndose=ep.ncol(), out=0;
  double u1= 0;
  double u_thres = nrow * (1 - pstop);
  double A = -el * pow(er, eal);
  double B = pow(1 - er,eag) - A;
  double C = -tl * pow(abs(tr - 1), tal);
  double D = pow(tr,tag) - C;
  double ue;
  double ut;
  double u;

  for (int j = 0; j < ndose; j++) {
    double total = 0;
    int u_add=0; 
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
      u = k1*ue  +  k2*ut + (1 - k1 -k2) *ue *ut;
      total += u;
      u_add += ceil(ustop - u);
    }
    
    if ((total >= u1) && (u_add <= u_thres)) {
      out=j+1;
      u1=total;
    } 
    
  }
 return(std::min(out, h_treated + 1));
}


