#' Spatial SGMCMC using SGRLD
#' @description
#' Main function from the paper to draw MCMC samples using SGRLD. Initial values usually obtained
#' from GpGp. Order of covariance parameters same as GpGp. Thoroughly tested with isotropic_matern.
#' Might struggle with other covariance structures for the time being.
#' 
#' @param y Vector of responses
#' @param X Matrix of covariates - usually the first column is going to be 1's for the intercept
#' @param NNarray Nearest neighbor object output from GpGp's find_ordered_nn function
#' @param locs Matrix of locations with each row corresponding to a location
#' @param beta_0 Initial values for GP mean parameters
#' @param covparams0 Initial values for covariance parameters (same structure as GpGp)
#' @param covfun_name Supports "matern_isotropic" and "exponential_isotropic"
#' @param algorithm Stochastic gradient algorithm. Choose from 'ADAMSGLD', 'MSGLD', 'PSGLD', 'SGRLD'. Default is 'SGRLD'.
#' @param lr Learning rate; 1e-3 by default
#' @param lr_min Lower bound for learning rate; 2e-6 by default
#' @param n_epochs Number of epochs
#' @param n_batch Size of each batch
#' @param n_burn Burn-in period for MCMC
#' @param thin Thin posterior samples
#' @param covparams_prior_params Check papaer for the distributions; order same as GpGp
#' @param silent Binary operator. `FALSE` by default
#' @return A list with 4 components - draws of beta, draws of covariance paramters, trace of loglik, and time
#' @export
fit_model_sgmcmc <- function(y, X, NNarray, locs, beta_0, covparams0, covfun_name = "matern_isotropic",
                             algorithm = "SGRLD",
                       lr = 1e-3, lr_min = 2e-6, n_epochs=100, n_batch = 250, n_burn = 2000,
                       thin = 5, covparams_prior_params, silent = F){
  if(algorithm=="SGRLD"){
    output <- sgrld_mcmc(y=y, X = X, NNarray = NNarray, locs = locs, beta_0 = beta_0,
                         covparams0 = covparams0, covfun_name = covfun_name, lr = lr,
                         lr_min = lr_min, n_epochs = n_epochs, n_batch = n_batch, n_burn = n_burn,
                         thin = thin, covparams_prior_params = covparams_prior_params, silent = silent)
    return(output)
  }
  if(algorithm=="MSGLD"){
    output <- msgld_mcmc(y=y, X = X, NNarray = NNarray, locs = locs, beta_0 = beta_0,
                         covparams0 = covparams0, covfun_name = covfun_name, lr = lr,
                         lr_min = lr_min, n_epochs = n_epochs, n_batch = n_batch, n_burn = n_burn,
                         thin = thin, covparams_prior_params = covparams_prior_params, silent = silent)
    return(output)
  }
  if(algorithm=="PSGLD"){
    output <- psgld_mcmc(y=y, X = X, NNarray = NNarray, locs = locs, beta_0 = beta_0,
                         covparams0 = covparams0, covfun_name = covfun_name, lr = lr,
                         lr_min = lr_min, n_epochs = n_epochs, n_batch = n_batch, n_burn = n_burn,
                         thin = thin, covparams_prior_params = covparams_prior_params, silent = silent)
    return(output)
  }
  if(algorithm=="ADAMSGLD"){
    output <- adamsgld_mcmc(y=y, X = X, NNarray = NNarray, locs = locs, beta_0 = beta_0,
                         covparams0 = covparams0, covfun_name = covfun_name, lr = lr,
                         lr_min = lr_min, n_epochs = n_epochs, n_batch = n_batch, n_burn = n_burn,
                         thin = thin, covparams_prior_params = covparams_prior_params, silent = silent)
    return(output)
  }
  print("Please select an algorithm from one of the supported ones")
}