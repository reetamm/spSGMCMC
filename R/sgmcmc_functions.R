
#' compute condition number of matrix
#'
#' @param info matrix
condition_number <- function(info){
  # assumes that information matrix has finite numbers in it
  if(max(diag(info))/min(diag(info)) > 1e6){
    return( max(diag(info))/min(diag(info)) )
  } else {
    ee <- eigen(info)
    return( max(ee$values)/min(ee$values) )
  }
}

#' Update preconditioner
#' @export
update_preconditioner <- function(V, grad, alpha=0.99, lambda = 1e-5, diag_V = T){
  
  if(diag_V){
    V <- alpha*V + (1-alpha)*grad*grad
    G <- 1/(lambda + sqrt(V))
  }else{
    V <- alpha*V + (1-alpha)*diag(grad*grad, nrow = length(grad))
    G <- chol2inv(chol(V + lambda*diag(nrow = length(grad))))
  }
  
  return(list(V = V, G = G))
}

#' update State
#' Not sure if used
updat_state <- function(phi, grad, G, lr, add_noise = T, diag_G = T){
  if(diag_G){
    new_phi <- phi + 0.5*lr*grad*G
    if(add_noise) new_phi = new_phi + rnorm(length(phi), mean = 0, sd = sqrt(lr*G))
  }else{
    new_phi <- phi + 0.5*lr*G%*%matrix(grad, ncol = 1)
    if(add_noise) new_phi = new_phi + t(chol(G))%*%rnorm(length(phi), mean = 0, sd = lr)
  }
    
  return(drop(new_phi))
}

#' SGRLD steps for covariance parameters
r_sgrld_step <- function(lr, info, logparms, grad){
  chol_info <- chol(info)
  G <- chol2inv(chol_info)
  if(length(lr) > 1){
    step <- 0.5*diag(lr)%*%G%*%grad
    noise <- backsolve(chol_info, rnorm(length(logparms), mean = 0, sd = sqrt(lr) ))
    return(drop(step) + drop(noise) + logparms)
  }else{
    step <- 0.5*lr*G%*%grad
    noise <- backsolve(chol_info, rnorm(length(logparms), mean = 0, sd = sqrt(lr) ))
    return(drop(step) + drop(noise) + logparms)
  }

}

#' SGRLD step for betas
#' Probably not being used anywhere
r_sgrld_step_beta <- function(lr, info, beta_c, beta_hat){
  chol_info <- chol(info)
  noise <- backsolve(chol_info, rnorm(length(beta_c), mean = 0, sd = sqrt(lr) ))
  return( (1-lr)*beta_c + lr*beta_hat + drop(noise) )
}

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
#' @param covparams_0 Initial values for covariance parameters (same structure as GpGp)
#' @param covfun_name Supports "matern_isotropic" and "exponential_isotropic"
#' @param lr Learning rate; 1e-3 by default
#' @param lr_min Lower bound for learning rate; 2e-6 by default
#' @param epochs Number of epochs
#' @param n_batch Size of each batch
#' @param n_burn Burn-in period for MCMC
#' @param thin Thin posterior samples
#' @param covparams_prior_params Check papaer for the distributions; order same as GpGp
#' @return A list with 4 components - draws of beta, draws of covariance paramters, trace of loglik, and time
#' @export
sgrld_mcmc <- function(y, X, NNarray, locs, beta_0, covparams0, covfun_name = "matern_isotropic",
                       lr = 1e-3, lr_min = 2e-6, n_epochs=100, n_batch = 250, n_burn = 2000,
                       thin = 5, covparams_prior_params, silent = F){
  n_iter <- n_epochs*length(y)%/%n_batch
  nmc_samples <- (n_iter - n_burn)%/%thin
  covparams_matrix <- matrix(0, nrow=nmc_samples, ncol = 4)
  beta_matrix <- matrix(0, nrow=nmc_samples, ncol = length(beta_0))
  llk_matrix <- rep(0, nmc_samples)
  n <- length(y)
  tic <- proc.time()
  iter_count <- 1
  restart_count <- 0
  while(iter_count <=n_iter && restart_count <=10){
    batch_ind <- sample(1:n, size = n_batch, replace = F)
    ordered_batch_ind <- sort(batch_ind)

    pass_list_all_batch_ordered <- vecchia_profbeta_loglik_grad_info(batch_id = ordered_batch_ind, covparms = covparams0,
                                                                        covfun_name = "matern_isotropic", y = y,
                                                                        X = X, current_beta = beta_0, locs = locs, NNarray = NNarray)
    grad_theta <- (n/n_batch)*pass_list_all_batch_ordered$grad_t
    tol <- 1e-10
    grad_phi <- covparams0*grad_theta + transformed_matern_parms_logprior_grad(log(covparams0), covparams_prior_params)
    if( any(is.infinite(grad_phi)) || any(is.na(grad_phi)) ){
      cat(paste("Iteration ", iter_count, " non-numeric results in gradient.\n"))
      cat(paste("Current covparams: ", covparams0, ".\n", sep=""))
      stop()
    }
    info_phi <- diag(covparams0)%*%pass_list_all_batch_ordered$info%*%diag(covparams0)
    if (condition_number(info_phi) > 1 / tol) {
      if (!silent) cat("Cond # of info matrix > 1/tol \n")
      #info <- 1.0*max(likobj0$info)*diag(nrow(likobj0$info))
      # regularize
      ee <- eigen(info_phi)
      ee_ratios <- ee$values/max(ee$values)
      ee_ratios[ ee_ratios < 1e-5 ] <- 1e-5
      ee$values <- max(ee$values)*ee_ratios
      info_phi <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)

      #diag(info) <- diag(info) + tol*max(diag(info))
    }
    #reparameterized_quantities(covparams0, pass_list_all_batch_ordered$grad_t, pass_list_all_batch_ordered$info, grad_phi, info_phi)
    covparams_new = r_sgrld_step(lr = lr, logparms = log(covparams0), info = info_phi, grad = grad_phi)
    counter_while <- 1

    while( any( is.infinite(exp(covparams_new))) || any(is.na(exp(covparams_new)))){
      cat("Need to reduce step size.\n")
      # which cov_params
      #if(length(lr) ==1) lr  <- rep(lr, length(covparams0))
      current_theta <- exp(covparams_new)
      inf_id <- which(is.infinite(current_theta))
      na_id <- which(is.na(current_theta))
      if(length(inf_id)>0){
        cat(paste("Indexes with Infinite values ", inf_id, ".\n", sep="" ))
        cat(paste("Current parameter values ", covparams0, ".\n", sep = "" ))
        cat(paste("Current gradient values ", grad_theta, ".\n", sep = ""))
        #lr[inf_id] <- lr[inf_id]/10
      }
      if(length(na_id)>0){
        cat(paste("Indexes with NaN values ", na_id, ".\n", sep="" ))
        cat(paste("Current parameter values ", covparams0, ".\n", sep = "" ))
        cat(paste("Current gradient values ", grad_theta, ".\n", sep = ""))
        #lr[na_id] <- lr[na_id]/10
      }
      cat("Restarting sampling and reduce learning rate.\n")
      if(covparams0[1]>10 || covparams0[1]<.15) covparams0[1] <- 1.
      if(covparams0[2]>3 || covparams0[2]<0.1) covparams0[2] <- 0.5
      if(covparams0[3]>8 || covparams0[3]<0.25) covparams0[3] <- 1.
      if(covparams0[4]>10 || covparams0[4]<0.05) covparams0[4] <- 0.5
      lr <- lr/2
      lr_min <- lr_min/10
      iter_count <- 1
      restart_count <- restart_count + 1
      if(restart_count > 11) break
      pass_list_all_batch_ordered <- vecchia_profbeta_loglik_grad_info(batch_id = ordered_batch_ind, covparms = covparams0,
                                                                          covfun_name = "matern_isotropic", y = y,
                                                                          X = X, current_beta = beta_0, locs = locs, NNarray = NNarray)
      grad_theta <- (n/n_batch)*pass_list_all_batch_ordered$grad_t
      tol <- 1e-10
      grad_phi <- covparams0*grad_theta + transformed_matern_parms_logprior_grad(log(covparams0), covparams_prior_params)
      if( any(is.infinite(grad_phi)) || any(is.na(grad_phi)) ){
        cat(paste("Iteration ", iter_count, " non-numeric results in gradient.\n"))
        cat(paste("Current covparams: ", covparams0, ".\n", sep=""))
        stop()
      }
      info_phi <- diag(covparams0)%*%pass_list_all_batch_ordered$info%*%diag(covparams0)
      if (condition_number(info_phi) > 1 / tol) {
        if (!silent) cat("Cond # of info matrix > 1/tol \n")
        #info <- 1.0*max(likobj0$info)*diag(nrow(likobj0$info))
        # regularize
        ee <- eigen(info_phi)
        ee_ratios <- ee$values/max(ee$values)
        ee_ratios[ ee_ratios < 1e-5 ] <- 1e-5
        ee$values <- max(ee$values)*ee_ratios
        info_phi <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
        
        #diag(info) <- diag(info) + tol*max(diag(info))
      }
      covparams_new = r_sgrld_step(lr = lr, logparms = log(covparams0), info = info_phi, grad = grad_phi)
      counter_while <- counter_while + 1
    }
    if(counter_while> 10) cat(paste("Reduced the step size by ", counter_while, ".\n"))
    covparams0 <- covparams_new
    covparams0 <- drop(exp(covparams0))
    beta_0 <- pass_list_all_batch_ordered$betahat
    if(iter_count>n_burn){
      #browser()
      j = iter_count - n_burn
      if(j%%thin==0){
        beta_matrix[j%/%thin,] <- beta_0
        covparams_matrix[j%/%thin,] <- covparams0
        llk_matrix[j%/%thin] <- pass_list_all_batch_ordered$loglik
      }
    }


    if(iter_count%%(n%/%n_batch)==0) {

      if(iter_count%%(10*n%/%n_batch)==0){
        epoch_c <- iter_count/(10*n%/%n_batch)
        if(lr >= lr_min) lr = lr/epoch_c

        cat(paste("Iteration ", iter_count, ", Epoch ", iter_count%/%(n%/%n_batch), sep = ""), "\n")
      }
    }
    iter_count <- iter_count + 1
  }
  toc <- proc.time()
  return( list(beta_samples = beta_matrix, theta_samples = covparams_matrix, llk_trace = llk_matrix, elapsed_time = toc["elapsed"] - tic["elapsed"]))
}



# ###
# ### Function for diagonal adaptive RMSPROP sgld
# ###
# psgld_mcmc <- function(y, X, NNarray, locs, beta_0, covparams0, covfun_name = "matern_isotropic",
#                          lr = 1e-3, lr_min = 2e-6, n_epochs=100, n_batch = 250, n_burn = 2000,
#                          thin = 5, alpha = 0.99, lambda=1e-5, covparams_prior_params, silent = F, 
#                          initial_V = NULL, diag_V = T, diag_G = T){
#   if(is.null(initial_V)){
#     V <-0
#     diag_V <- T
#     diag_G <-T
#   }else{
#     V <- initial_V
#   }
#   #V = 0
#   add_noise <- T
#   n_iter <- n_epochs*length(y)/n_batch
#   nmc_samples <- (n_iter - n_burn)/thin
#   covparams_matrix <- matrix(0, nrow=nmc_samples, ncol = 4)
#   beta_matrix <- matrix(0, nrow=nmc_samples, ncol = 2)
#   llk_matrix <- rep(0, nmc_samples)
#   n <- length(y)
#   tic <- proc.time()
#   iter_count <- 1
#   restart_count <- 0
#   while(iter_count <= n_iter && restart_count <=10 ){
#     batch_ind <- sample(1:n, size = n_batch, replace = F)
#     ordered_batch_ind <- sort(batch_ind)
# 
#     pass_list_all_batch_ordered <- vecchia_profbeta_loglik_grad_info(batch_id = ordered_batch_ind, covparms = covparams0,
#                                                                         covfun_name = "matern_isotropic", y = y,
#                                                                         X = X, current_beta = beta_0, locs = locs, NNarray = NNarray)
#     #browser()
#     grad_theta <- (n/n_batch)*pass_list_all_batch_ordered$grad_t
#     grad_theta <- grad_theta/n
#     grad_phi <- covparams0*grad_theta
#     if( any(is.infinite(grad_phi)) || any(is.na(grad_phi)) ){
#       cat(paste("Iteration ", iter_count, " non-numeric results in gradient.\n"))
#       cat(paste("Current covparams: ", covparams0, ".\n", sep=""))
#       stop()
#     }
#     prior_grad <- transformed_matern_parms_logprior_grad(log(covparams0), covparams_prior_params)
#     new_cond <- update_preconditioner(V, grad_phi, alpha, lambda, diag_V = diag_V)
#     V <- new_cond$V
#     G <- new_cond$G
#     new_phi <- updat_state(log(covparams0), grad = n*grad_phi + prior_grad, G = G, lr = lr, add_noise = add_noise, diag_G = diag_G)
#     counter_while <- 1
#     while( any( is.infinite(exp(new_phi))) || any(is.na(exp(new_phi)))){
#       cat("Need to reduce step size.\n")
#       # which cov_params
#       #if(length(lr) ==1) lr  <- rep(lr, length(covparams0))
#       current_theta <- exp(new_phi)
#       inf_id <- which(is.infinite(current_theta))
#       na_id <- which(is.na(current_theta))
#       if(length(inf_id)>0){
#         cat(paste("Indexes with Infinite values ", inf_id, ".\n", sep="" ))
#         cat(paste("Current parameter values ", covparams0, ".\n", sep = "" ))
#         cat(paste("Current gradient values ", grad_theta, ".\n", sep = ""))
#         #lr[inf_id] <- lr[inf_id]/2
#       }
#       if(length(na_id)>0){
#         cat(paste("Indexes with NaN values ", na_id, ".\n", sep="" ))
#         cat(paste("Current parameter values ", covparams0, ".\n", sep = "" ))
#         cat(paste("Current gradient values ", grad_theta, ".\n", sep = ""))
#         #lr[na_id] <- lr[na_id]/2
#       }
#       cat("Restarting sampling and reduce learning rate.\n")
#       if(covparams0[1]>10 || covparams0[1]<.15) covparams0[1] <- 1.
#       if(covparams0[2]>3 || covparams0[2]<0.1) covparams0[2] <- 0.5
#       if(covparams0[3]>8 || covparams0[3]<0.25) covparams0[3] <- 1.
#       if(covparams0[4]>10 || covparams0[4]<0.05) covparams0[4] <- 0.5
#       lr = lr/2
#       lr_min = lr/10
#       iter_count <- 1
#       restart_count <- restart_count + 1
#       if(restart_count > 10) break
#       V <- initial_V
#       pass_list_all_batch_ordered <- vecchia_profbeta_loglik_grad_info(batch_id = ordered_batch_ind, covparms = covparams0,
#                                                                           covfun_name = "matern_isotropic", y = y,
#                                                                           X = X, current_beta = beta_0, locs = locs, NNarray = NNarray)
#       #browser()
#       grad_theta <- (n/n_batch)*pass_list_all_batch_ordered$grad_t
#       grad_theta <- grad_theta/n
#       grad_phi <- covparams0*grad_theta
#       if( any(is.infinite(grad_phi)) || any(is.na(grad_phi)) ){
#         cat(paste("Iteration ", iter_count, " non-numeric results in gradient.\n"))
#         cat(paste("Current covparams: ", covparams0, ".\n", sep=""))
#         stop()
#       }
#       prior_grad <- transformed_matern_parms_logprior_grad(log(covparams0), covparams_prior_params)
#       new_cond <- update_preconditioner(V, grad_phi, alpha, lambda, diag_V = diag_V)
#       V <- new_cond$V
#       G <- new_cond$G
#       new_phi <- updat_state(log(covparams0), grad = n*grad_phi + prior_grad, G = G, lr = lr, add_noise = add_noise, diag_G = diag_G)
#       counter_while <- counter_while + 1
#     }
#     if(counter_while> 10) cat(paste("Reduced the step size by ", counter_while, ".\n"))
#     covparams0 <- new_phi
#     covparams0 <- drop(exp(covparams0))
#     beta_0 <- pass_list_all_batch_ordered$betahat
#     if(iter_count>n_burn){
#       j = iter_count - n_burn
#       if(j%%thin==0){
#         beta_matrix[j%/%thin,] <- beta_0
#         covparams_matrix[j%/%thin,] <- covparams0
#         llk_matrix[j%/%thin] <- pass_list_all_batch_ordered$loglik
#       }
#     }
# 
# 
#     if(iter_count%%(n%/%n_batch)==0) {
# 
#       if(iter_count%%(10*n%/%n_batch)==0){
#         epoch_c <- iter_count/(10*n%/%n_batch)
#         if(lr >= lr_min) lr = lr/epoch_c
# 
#         cat(paste("Iteration ", iter_count, ", Epoch ", iter_count%/%(n%/%n_batch), sep = ""), "\n")
#       }
#     }
#     iter_count <- iter_count + 1
#   }
#   toc <- proc.time()
#   return( list(beta_samples = beta_matrix, theta_samples = covparams_matrix, llk_trace = llk_matrix, elapsed_time = toc["elapsed"] - tic["elapsed"]))
# }
