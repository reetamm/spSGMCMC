#ifndef COVMATRIX_FUNS_H
#define COVMATRIX_FUNS_H

#include <RcppArmadillo.h>
#include <iostream>
#include <vector>

//cov matrix file
#include "cov_header.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]


void get_covfun(std::string covfun_name_string,  mat (*p_covfun[1])(arma::vec, arma::mat), cube (*p_d_covfun[1])(arma::vec, arma::mat)  )
{
  
  if( covfun_name_string.compare("matern_isotropic") == 0 )
  { 
    p_covfun[0] = matern_isotropic; 
    p_d_covfun[0] = d_matern_isotropic;
  } 
  else if( covfun_name_string.compare("exponential_isotropic") == 0 )
  { 
    p_covfun[0] = exponential_isotropic; 
    p_d_covfun[0] = d_exponential_isotropic;
  }
  else if( covfun_name_string.compare("matern15_isotropic") == 0 )
  { 
    p_covfun[0] = matern15_isotropic; 
    p_d_covfun[0] = d_matern15_isotropic;
  } 
  else { // stop the program
    Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
  }
  
}



#endif
