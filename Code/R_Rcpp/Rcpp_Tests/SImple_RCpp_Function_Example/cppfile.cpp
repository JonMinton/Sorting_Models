#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;


//R calls f1 by giving it an integer value
//f1 calls f2, passing the integer value to it
//f2 returns double the value passed to it
//f1 returns the value passed to it by f2


double f2
  (
    double input
  )
  {
    double out;
    out = input * 2;
    return(out);
  }


// [[Rcpp::export]]
double f1 
  (
    double input
   )
  {
    
    double out;
    out = f2(input);
    
    return(out);
  }
  
  

  