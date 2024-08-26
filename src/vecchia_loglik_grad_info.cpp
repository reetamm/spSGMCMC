
#define BOOST_DISABLE_ASSERTS

#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "onepass.h"
#include "prior_header.h"
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]

using namespace std;
using namespace Rcpp;
using namespace arma;


//' @param X Design matrix of covariates. Row \code{i} of \code{X} contains
//' the covariates for the observation at row \code{i} of \code{locs}.
//' @return A list containing
//' \itemize{
//'     \item \code{loglik}: the loglikelihood
//'     \item \code{grad}: gradient with respect to covariance parameters
//'     \item \code{info}: Fisher information for covariance parameters
//'     \item \code{betahat}: profile likelihood estimate of regression coefs
//'     \item \code{betainfo}: information matrix for \code{betahat}.
//'     \item \code{grad_beta}: gradient of the loglikelihood with respect to mean parameters
//' }
//' The covariance matrix for \code{$betahat} is the inverse of \code{$betainfo}.
//' @export
// [[Rcpp::export]]
List vecchia_profbeta_loglik_grad_info(
    IntegerVector batch_id,
    NumericVector covparms,
    StringVector covfun_name,
    NumericVector& y,
    NumericMatrix& X,
    NumericVector current_beta,
    const NumericMatrix& locs,
    NumericMatrix& NNarray ){

    NumericVector ll(1);
    NumericVector grad_t( covparms.length() );
    NumericVector grad( covparms.length() );
    NumericVector grad_beta (X.ncol());
    NumericVector betahat( X.ncol() );
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );
    
    // data dimensions
    //int n = batch_id.length();//y.length();
    //int m = NNarray.ncol();
    //int p = X.ncol();
    //int nparms = covparms.length();
    //int dim = locs.ncol();
    
    // // likelihood objects
    // arma::mat XSX = arma::mat(p, p, fill::zeros);
    // arma::vec ySX = arma::vec(p, fill::zeros);
    // double ySy = 0.0;
    // double logdet = 0.0;
    // 
    // // gradient objects
    // arma::cube dXSX = arma::cube(p,p,nparms,fill::zeros);
    // arma::mat dySX = arma::mat(p, nparms, fill::zeros);
    // arma::vec dySy = arma::vec(nparms, fill::zeros);
    // arma::vec dlogdet = arma::vec(nparms, fill::zeros);
    // // fisher information
    // arma::mat ainfo = arma::mat(nparms, nparms, fill::zeros);
    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    my_synthesize(covparms, covfun_name, locs, batch_id, NNarray, y, X, current_beta,
        &ll, &betahat, &grad, &grad_t, &grad_beta, &info, &betainfo, true, true);//, XSX, ySX, ySy,
        //logdet, dXSX, dySX, dySy, dlogdet, ainfo);
    List pieces = List::create( Named("loglik") = ll,
                                Named("betahat") = betahat,
                                Named("grad") = grad,
                                Named("grad_t") = grad_t,
                                Named("grad_beta") = grad_beta,
                                Named("info") = info,
                                Named("betainfo") = betainfo );
    //Named("logdet") = logdet, Named("dlogdet") = dlogdet, Named("ySy") = ySy, Named("dySy") = dySy,Named("XSX") = XSX,Named("dXSX") = dXSX,Named("ySX") = ySX,Named("dySX") = dySX,
    return pieces;

}

//' Function to transform the gradients and fisher information
//' of the likelihood from the constrained to unconstrained space
//' @param \code{grad} gradient with respect to cov_params
//' @param \code{info} Fisher information for cov_params
//' @param \code{grad_phi} reference to vector for storing unconstrained gradients
//' @param \code{info_phi} reference to matrix for unconstrained Fisher info
//' @export
// [[Rcpp::export]]
void reparameterized_quantities(
    NumericVector& cov_params,
    NumericVector& grad,
    NumericMatrix& info,
    NumericVector& grad_phi,
    NumericVector& info_phi){
    
    arma::vec cov_params_c = arma::vec(cov_params.begin(), cov_params.length());
    arma::vec grad_c = arma::vec(grad.begin(), grad.length(), false);
    arma::mat info_c = arma::mat(info.begin(), info.nrow(), info.ncol(), false);
    
    for(int i=0; i<grad.length(); i++){
      grad_phi(i) = cov_params(i)*grad(i);
      info_phi(i,i) = cov_params(i)*cov_params(i)*info(i,i);
      for(int j=0; j<i; j++){
        info_phi(i,j) = cov_params(i)*cov_params(j)*info(i,j);
        info_phi(j,i) = info_phi(i,j);
      }
    }
}



//' Function to take an SGDRLD step
//' for the covariance parameters
//' This function returns \eqn{( \log(\sigma^2), \log(\alpha), \log(\nu), \log(\tau^2) )_{t+1}}
//' @param \eqn{\epsilon} the step size
//' @param \code{info} the fisher information matrix i.e. preconditionner for the cov parameters
//' @parm \code{cov_params} the current state of the cov parameters
//' @parm \code{grad} sample gradient with respect to \eqn{( \log(\sigma^2), \log(\alpha), \log(\nu), \log(\tau^2) )_{t+1}}
//' @export
// [[Rcpp::export]]
NumericVector SGRLD_step(
    double epsilon,
    NumericMatrix& info,
    NumericVector& cov_params,
    NumericVector& grad){
    
    arma::mat info_c = arma::mat(info.begin(), info.nrow(), info.ncol(), false);
    arma::mat grad_c = arma::vec(grad.begin(), grad.length(), false);
    arma::mat chol_info = chol(info_c, "upper");
    arma::vec noise = arma::solve( arma::trimatu(chol_info), arma::randn(cov_params.length(),
                                                 arma::distr_param(0., sqrt(epsilon))));
    arma::vec step = 0.5*epsilon*arma::solve(info_c, grad_c);
    arma::vec cov_params_c = arma::vec(cov_params.begin(), cov_params.length());
    cov_params_c = cov_params_c + step + noise;
    return Rcpp::NumericVector(cov_params_c.begin(), cov_params_c.end());
    
}

//' Function to run spSGMCMC for a number of iterations
//' @param \code{y}
//' @param \code{X}
//' @param \code{NNarray}
// [[Rcpp::export]]
List SGRLD_loop(
    NumericVector& y,
    NumericMatrix& X,
    NumericMatrix& NNarray,
    StringVector covfun_name,
    NumericMatrix& locs,
    NumericVector beta_0,
    NumericVector covparams0,
    NumericMatrix prior_params,
    IntegerVector& indexes,
    int n_epochs,
    int n_batch,
    int n_burn,
    double lr,
    int thin ){
    //Start of the main function
    int N = y.length();
    int p = X.ncol();
    int d = covparams0.length();
    int n_iter = N*n_epochs/n_batch;
    int n_samples = (n_iter - n_burn)/thin;
    double grad_factor =  N/(double)n_batch;
    int n_batches = (int)(grad_factor+0.5);
    //matrices to save samples
    arma::mat covparams_samples(n_samples, d, arma::fill::zeros );
    arma::mat beta_samples(n_samples, p, arma::fill::zeros);
    arma::vec llk_trace(n_samples, arma::fill::zeros);
    NumericVector current_covparams = covparams0;//arma::vec(covparams0.begin(), covparams0.length());
    arma::vec current_covparams_c = arma::vec(current_covparams.begin(), current_covparams.length());
    NumericVector current_beta = beta_0;//arma::vec(beta_0.begin(), beta_0.length());
    NumericVector grad_phi(covparams0.length());
    NumericMatrix info_phi(d,d);
    NumericVector grad_theta(covparams0.length());
    //NumericMatrix info_theta(d,d);
    int iter = 1;
    //Now start the loop
    for( int epoch=0; epoch<n_epochs; epoch++){
      //shuffle the dataset
      IntegerVector Shuffle_id = Rcpp::sample(indexes, N);
      //number of minibatches per epoch
      int start = 0;
      for(int i = 0; i<n_batches; i++){
        //Now one step spSGMCMC
        int end = std::min(n_batch*(i+1), N) - 1;
        IntegerVector batch_id = Shuffle_id[Rcpp::seq(start, end)];
        std::sort(batch_id.begin(), batch_id.end());
        start = n_batch*(i+1);
        //Now we have the ordered_batch_id
        List pieces = vecchia_profbeta_loglik_grad_info(batch_id, current_covparams, covfun_name, y, X, current_beta, locs, NNarray);
        //Now we need to reparametereize the quantities
        double pieces3 = as<NumericVector>(pieces[3])[0];
        grad_theta = grad_factor*pieces3;
        NumericMatrix info_theta = pieces['info'];
        reparameterized_quantities(current_covparams, grad_theta, info_theta, grad_phi, info_phi );
        //info_phi.diag()+=0.0001;
        current_beta = pieces[1];
        Rcpp::NumericVector phi = Rcpp::log(current_covparams);
        phi = SGRLD_step(lr, info_phi, phi , grad_phi);
        current_covparams = Rcpp::exp(phi);
        if(iter > n_burn){
          int j = iter - n_burn;
          int r = j%thin;
          if(r==0){
            int k = j/thin;
            beta_samples.row(k) = arma::vec(current_beta.begin(), current_beta.length());
            covparams_samples.row(k) = arma::vec(current_covparams.begin(), current_covparams.length());
            llk_trace(k) = pieces[0];
            }
          }
        iter+=1;
        }
      if((epoch%20==0) && (epoch>0)){
        Rcpp::Rcout << "Iteration " << iter << ", and Ecpoch " << epoch+1 <<"\n";
        }
      if((epoch%10==0) && (epoch>0) && (lr>1e-5) ){
        lr/=2.;
        }
      }
    Rcpp::List traces =  Rcpp::List::create( Named("beta_samples") = beta_samples,
                                             Named("theta_samples") = covparams_samples,
                                             Named("llk_trace") = llk_trace);
    return traces;
}


//[[Rcpp::export]]
List sample_pieces(
  IntegerVector batch_id,
  NumericMatrix NNarray,
  NumericVector covparms,
  NumericMatrix X,
  NumericVector y,
  NumericMatrix locs,
  StringVector covfun_name){
  
  arma::vec covparms_c = arma::vec(covparms.begin(),covparms.length());
  arma::mat locs_c = arma::mat(locs.begin(),locs.nrow(),locs.ncol());
  arma::mat NNarray_c = arma::mat(NNarray.begin(),NNarray.nrow(),NNarray.ncol());
  arma::vec y_c = arma::vec(y.begin(),y.length());
  arma::mat X_c = arma::mat(X.begin(),X.nrow(),X.ncol());
  int n = batch_id.length(); // y.n_elem;
  int m = NNarray_c.n_cols;
  int p = X_c.n_cols;
  int nparms = covparms_c.n_elem;
  int dim = locs_c.n_cols;
  
  // convert StringVector to std::string to use .compare() below
  std::string covfun_name_string;
  covfun_name_string = covfun_name[0];
  
  // assign covariance fun and derivative based on covfun_name_string
  
  /* p_covfun is an array of length 1. Its entry is a pointer to a function which takes
   in arma::vec and arma::mat and returns mat. p_d_covfun is analogous. This was a workaround for the solaris bug*/
  
  mat (*p_covfun[1])(arma::vec, arma::mat);
  cube (*p_d_covfun[1])(arma::vec, arma::mat);
  get_covfun(covfun_name_string, p_covfun, p_d_covfun);
  
  arma::mat l_XSX = arma::mat(p, p, fill::zeros);
  arma::vec l_ySX = arma::vec(p, fill::zeros);
  //double l_ySy = 0.0;
  //double l_logdet = 0.0;
  arma::cube l_dXSX = arma::cube(p,p, nparms, fill::zeros);
  arma::mat l_dySX = arma::mat(p, nparms, fill::zeros);
  arma::vec l_dySy = arma::vec(nparms, fill::zeros);
  arma::vec l_dlogdet = arma::vec(nparms, fill::zeros);
  arma::mat l_ainfo = arma::mat(nparms, nparms, fill::zeros);
  bool profbeta = true;
  bool grad_info = false;
  int bsize = std::min( batch_id[0],m);
  arma::mat locsub(bsize, dim);
  arma::vec ysub(bsize);
  arma::mat X0( bsize, p );
  for(int i=0; i<n; i++){
    //here we need to modify this because neighbour size depends on the batch_id
    // batch_id[i] is the order of the observation
    // so the neighbor set size is min(batch_id[i], m)
    
    
    //std::vector<std::chrono::steady_clock::time_point> tt;
    
    //tt.push_back( std::chrono::steady_clock::now() );
    
    // first, fill in ysub, locsub, and X0 in reverse order
    // arma::mat locsub(bsize, dim);
    // arma::vec ysub(bsize);
    // arma::mat X0( bsize, p );
    for(int j=bsize-1; j>=0; j--){
      ysub(bsize-1-j) = y( NNarray_c(batch_id[i]-1,j)-1 );
      for(int k=0;k<dim;k++){ locsub(bsize-1-j,k) = locs_c( NNarray(batch_id[i]-1,j)-1, k ); }
      if(profbeta){
        for(int k=0;k<p;k++){ X0(bsize-1-j,k) = X_c( NNarray(batch_id[i]-1,j)-1, k ); }
      }
    }
  }
    
    // compute covariance matrix and derivatives and take cholesky
    arma::mat covmat = p_covfun[0]( covparms, locsub );
    
    arma::cube dcovmat;
    if(grad_info){
      dcovmat = p_d_covfun[0]( covparms, locsub );
    }
    List results = List::create( Named("v") = ysub,
                                 Named("covmat") = covmat,
                                 Named("R") = X0,
                                 Named("locsub") = locsub
                                 );
    return results;
}

