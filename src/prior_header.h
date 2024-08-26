#ifndef PRIOR_HEADER_H
#define PRIOR_HEADER_H

#include <Rcpp.h>


using namespace Rcpp;



//' log-prior for regression coefficients \eqn{\beta}
//' We consider a standard normal prior on all the components 
//' @param beta, the current value of the regression parameters
//' @return a scalar i.e the log-prior at beta
//' @export
// [[Rcpp::export]]
double beta_logpior( NumericVector& beta, double var = 1.0 ){
   return Rcpp::sum( Rcpp::dnorm(beta, 0.0, sqrt(var), true) );
 }
 
//' gradient of log_prior on \eqn{\beta}
//' @param beta, current value
//' @return vector with same size as \eqn{\beta}
//' @export
// [[Rcpp::export]]
NumericVector beta_grad_logprior(NumericVector& beta, double var = 1.0){
  NumericVector grad_lp = -beta/var;
  return grad_lp;
}
 
 
//' log-prior for matern cov parameters
//' @param covparms \eqn{\sigma^2, \alpha, \nu, \tau^2}
//' the nugget is \eqn{\sigma^2  \tau^2}
//' the priors are log-normal (1,1) for the smoothness \eqn{\nu}
//' Gamma(9,2) for the range \eqn{\alpha}
//' Gamma(.1,.1) for the spatial variance \eqn{\sigma^2}
//' and  Gamma(.1, .1) for the scaled nugget
//' @return a vector of log pdf
//' @export
// [[Rcpp::export]]
NumericVector matern_parms_logprior(NumericVector& covparms, NumericMatrix& prior_params){
  NumericVector covparms_lp = { R::dgamma( covparms(0), prior_params(0,0), 1.0/prior_params(0,1),  true),
                                R::dgamma(covparms(1), prior_params(1,0), 1.0/prior_params(1,1),  true),
                                R::dlnorm(covparms(2), prior_params(2,0), prior_params(2,1), true),
                                R::dgamma(covparms(3), prior_params(3,0), 1.0/prior_params(3,1), true)};
  return covparms_lp;
}
 
//' gradient of log_prior for matern cov_parameters
//' @param covparms \eqn{\sigma^2, \alpha, \nu, \tau^2}
//' the nugget is \eqn{\sigma^2  \tau^2}
//' the priors are log-normal (1,1) for the smoothness \eqn{\nu}
//' Gamma(9,2) for the range \eqn{\alpha}
//' Gamma(.1,.1) for the spatial variance \eqn{\sigma^2}
//' and  Gamma(.1, .1) for the scaled nugget
//' @return a vector same length as covparms.
//' @export
// [[Rcpp::export]]
NumericVector matern_parms_prior_grad(NumericVector& covparms, NumericMatrix& prior_params){
  NumericVector covparms_grad(covparms.length());
  covparms_grad(0) = ( prior_params(0,0) - 1.0)/covparms(0) - prior_params(0,1);
  covparms_grad(1) = ( prior_params(1,0) - 1.0)/covparms(1) - prior_params(1,1);
  covparms_grad(2) = -1.0*(1.0 + (log(covparms[2]) - prior_params(2,0))/pow( prior_params(2,1), 2.0 ) )/covparms(2); // 
  covparms_grad(3) = ( prior_params(3,0) - 1.0)/covparms(3) - prior_params(3,1);
  return covparms_grad;
}

//' Bijector for the cov parameters
//' all the parameters are sampled on the log scale
//' @export
// [[Rcpp::export]]
NumericVector parms_link( NumericVector& covparms){
  NumericVector logparms(covparms.length());
  logparms = Rcpp::log(covparms);
  return logparms;
}

//' Inverse Bijector for the cov parameters
//' all the parameters are sampled on the log scale
//' @export
// [[Rcpp::export]]
NumericVector parms_invlink( NumericVector& logparms){
  NumericVector covparms(logparms.length());
  covparms = Rcpp::exp(logparms);
  return covparms;
}


//' Bijector gradient for the cov parameters
//' @export
// [[Rcpp::export]]
NumericVector parms_link_grad(NumericVector& covparms){
  NumericVector dlink = 1.0/covparms;
  return dlink;
}

//' Inverse Bijector gradient for the cov parameters
//' @export
// [[Rcpp::export]]
NumericVector parms_invlink_grad(NumericVector& logparms){
  NumericVector dinvlink = Rcpp::exp(logparms);
  return dinvlink;
}

//' log prior of transformed matern parameters
//' @export
// [[Rcpp::export]]
NumericVector transformed_matern_parms_logprior(NumericVector& logparms, NumericMatrix& prior_params){
  NumericVector covparms = parms_invlink(logparms);
  NumericVector log_jacob = Rcpp::log( parms_invlink_grad(logparms) );
  NumericVector lp = matern_parms_logprior(covparms, prior_params);
  return lp + log_jacob;
}

//' log prior gradient of transformed matern parameters
//' @export
// [[Rcpp::export]]
NumericVector transformed_matern_parms_logprior_grad(NumericVector& logparms, NumericMatrix& prior_params){
  NumericVector covparms = parms_invlink(logparms);
  NumericVector jacob = parms_invlink_grad(logparms) ;
  NumericVector dlog_jacob = {1.0, 1.0, 1.0, 1.0};
  NumericVector cov_grad = matern_parms_prior_grad(covparms, prior_params);
  return cov_grad*jacob + dlog_jacob;
}

#endif
