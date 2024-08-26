
context("Covariance Functions")

covfun_names <- c(
  "matern_isotropic",
  "exponential_isotropic",
  "matern15_isotropic"
  )

get_test_locs <- function(covfun_name,n){
  
  nside <- round(sqrt(n))
  longrid <- seq(1,360,length.out=nside)
  latgrid <- seq(-80,80,length.out=nside)
  lonlat <- as.matrix(expand.grid(longrid,latgrid))
  if(covfun_name=="exponential_isotropic"){
    locs <- matrix(runif(3*n),n,3)           
  } else if(covfun_name=="matern_isotropic"){
    locs <- matrix(runif(2*n),n,2)
  } else if(covfun_name=="matern15_isotropic"){
    locs <- matrix(runif(2*n),n,2)
  } else {
    stop("unrecognized covariance in testing function")
  }
  return(locs)
}

test_that("covariance functions return positive definite matrix", {
  
  
  n <- 100    
  for(j in 1:length(covfun_names)){
    locs <- get_test_locs(covfun_names[j],n)
    covparms <- get_start_parms(rnorm(n),rep(1,n),locs,covfun_names[j])
    covfun <- get( covfun_names[j] )
    covmat <- covfun( covparms$start_parms, locs )
    cholmat <- t(chol(covmat))
    logdet <- 2*sum(log(diag(cholmat)))
    expect_lt( logdet, sum(log(diag(covmat))) )
  }
  
})



test_that("covariance function derivatives match finite differencing", {
  
  n <- 100    
  
  for(j in 1:length(covfun_names)){
    
    locs <- get_test_locs(covfun_names[j],n)
    covparms <- get_start_parms(rnorm(n),rep(1,n),locs,covfun_names[j])
    covparms <- covparms$start_parms
    nparms <- length(covparms)
    covfun <- get( covfun_names[j] )
    dcovfun <- get(paste0("d_",covfun_names[j]))
    covmat <- covfun( covparms, locs )
    dcovmat <- dcovfun( covparms,locs )
    ddcov <- array(NA, c(n,n,nparms) )
    eps <- 1e-8
    for(k in 1:nparms){
      dcovparms <- covparms
      dcovparms[k] <- covparms[k] + eps
      cov <- covfun( dcovparms, locs )
      ddcov[,,k] <- (cov - covmat)/eps
    }
    denom <- covmat[1,1]
    expect_equal( dcovmat/denom, ddcov/denom, tolerance = 1e-4 )
    
  }
})
