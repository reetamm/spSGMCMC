#' get default starting values of covariance parameters
#'
#' @param y response
#' @param X design matrix
#' @param locs locations
#' @param covfun_name string name of covariance function
get_start_parms <- function(y,X,locs,covfun_name){
  
  fitlm <- stats::lm(y ~ X - 1 )
  start_var <- summary(fitlm)$sigma^2
  start_smooth <- 0.8
  start_nug <- 0.1
  n <- length(y)
  
  randinds <- sample(1:n, min(n,200))
  dmat <- fields::rdist(locs[randinds,])
  
  if(covfun_name %in% c("exponential_isotropic","exponential_isotropic_fast")){
    start_range <- mean( dmat )/4
    start_parms <- c(start_var, start_range, start_nug)
  }
  if(covfun_name == "matern_isotropic"){
    start_range <- mean( dmat )/4
    start_parms <- c(start_var, start_range, start_smooth, start_nug)
  }
  if(covfun_name == "matern15_isotropic"){
    start_range <- mean( dmat )/4
    start_parms <- c(start_var, start_range, start_nug)
  }
  return( list( start_parms = start_parms ) )
}


