#' ADAM SGLD update step
#' 
adam_sgld_update <- function(phi, grad_t, lr,
                             momentum_t, V_t,
                             beta_1 = 0.9, 
                             beta_2 = 0.9, bias_factor=4
                             ){
  ### first update momentum
  new_momentum <- beta_1*momentum_t + (1-beta_1)*(grad_t)
  V <- beta_2*V_t + (1-beta_2)*crossprod(matrix(grad_t, nrow = 1))
  R <- chol(V + 1e-5*diag(nrow = length(grad_t)))
  A <- forwardsolve(t(R), new_momentum)
  new_phi <- phi + lr*(grad_t + bias_factor*A) + sqrt(2*lr)*rnorm(n = length(phi))
  return(list(new_phi = new_phi, momentum = new_momentum, V = V))
}

#' Main ADAM SGLD function
#' @export
adamsgld_mcmc <- function(y, X, NNarray, locs, beta_0, covparams0, covfun_name = "matern_isotropic",
                      lr = 1e-3, lr_min = 2e-6, n_epochs=100, n_batch = 250, n_burn = 2000,
                      thin = 5, beta1_factor = 0.9, beta2_factor = 0.99, 
                      bias_factor = 2, lambda=1e-5, covparams_prior_params, silent = F, initial_V = NULL){
  momentum <- 0
  if(is.null(initial_V)){
    V <- 0
  }else{
    V <- initial_V
  }
  add_noise <- T
  n_iter <- n_epochs*length(y)%/%n_batch
  nmc_samples <- (n_iter - n_burn)%/%thin
  covparams_matrix <- matrix(0, nrow=nmc_samples, ncol = 4)
  beta_matrix <- matrix(0, nrow=nmc_samples, ncol = length(beta_0))
  llk_matrix <- rep(0, nmc_samples)
  n <- length(y)
  tic <- proc.time()
  iter_count <- 1
  restart_count <- 0
  while(iter_count <= n_iter && restart_count <=3 ){
    batch_ind <- sample(1:n, size = n_batch, replace = F)
    ordered_batch_ind <- sort(batch_ind)
    
    pass_list_all_batch_ordered <- vecchia_profbeta_loglik_grad_info(batch_id = ordered_batch_ind, covparms = covparams0,
                                                                        covfun_name = "matern_isotropic", y = y,
                                                                        X = X, current_beta = beta_0, locs = locs, NNarray = NNarray)
    #browser()
    grad_theta <- (n/n_batch)*pass_list_all_batch_ordered$grad_t
    grad_theta <- grad_theta/n
    grad_phi <- covparams0*grad_theta
    if( any(is.infinite(grad_phi)) || any(is.na(grad_phi)) ){
      cat(paste("Iteration ", iter_count, " non-numeric results in gradient.\n"))
      cat(paste("Current covparams: ", covparams0, ".\n", sep=""))
      stop()
    }
    prior_grad <- transformed_matern_parms_logprior_grad(log(covparams0), covparams_prior_params)

    new_results <- adam_sgld_update(phi = log(covparams0), grad_t = n*grad_phi + prior_grad,
                                    momentum_t = momentum, V_t = V,  beta_1 = beta1_factor,
                                    beta_2 = beta2_factor, lr = lr,
                                    bias_factor = bias_factor)
    new_phi <- new_results$new_phi
    new_momentum <- new_results$momentum
    new_V <- new_results$V
    counter_while <- 1

    while( any( is.infinite(exp(new_phi))) || any(is.na(exp(new_phi)))){
      # which cov_params
      #if(length(lr) ==1) lr  <- rep(lr, length(covparams0))
      current_theta <- exp(new_phi)
      inf_id <- which(is.infinite(current_theta))
      na_id <- which(is.na(current_theta))
      if(length(inf_id)>0){
        cat(paste("Indexes with Infinite values ", inf_id, ".\n", sep="" ))
        cat(paste("Current parameter values ", covparams0, ".\n", sep = "" ))
        cat(paste("Current gradient values ", grad_theta, ".\n", sep = ""))
        #lr[inf_id] <- lr[inf_id]/2
      }
      if(length(na_id)>0){
        cat(paste("Indexes with NaN values ", na_id, ".\n", sep="" ))
        cat(paste("Current parameter values ", covparams0, ".\n", sep = "" ))
        cat(paste("Current gradient values ", grad_theta, ".\n", sep = ""))
        #lr[na_id] <- lr[na_id]/2
      }
      cat("Restarting sampling and reduce learning rate.\n")
      if(covparams0[1]>10 || covparams0[1]<.15) covparams0[1] <- 1.
      if(covparams0[2]>3 || covparams0[2]<0.1) covparams0[2] <- 0.5
      if(covparams0[3]>8 || covparams0[3]<0.25) covparams0[3] <- 1.
      if(covparams0[4]>10 || covparams0[4]<0.05) covparams0[4] <- 0.5
      lr = lr/1.5
      lr_min = lr/1.5
      iter_count <- 1
      restart_count <- restart_count + 1
      if(restart_count > 10) break
      momentum <- 0
      V <- 0
      pass_list_all_batch_ordered <- vecchia_profbeta_loglik_grad_info(batch_id = ordered_batch_ind, covparms = covparams0,
                                                                          covfun_name = "matern_isotropic", y = y,
                                                                          X = X, current_beta = beta_0, locs = locs, NNarray = NNarray)
      #browser()
      grad_theta <- (n/n_batch)*pass_list_all_batch_ordered$grad_t
      grad_theta <- grad_theta/n
      grad_phi <- covparams0*grad_theta
      if( any(is.infinite(grad_phi)) || any(is.na(grad_phi)) ){
        cat(paste("Iteration ", iter_count, " non-numeric results in gradient.\n"))
        cat(paste("Current covparams: ", covparams0, ".\n", sep=""))
        stop()
      }
      prior_grad <- transformed_matern_parms_logprior_grad(log(covparams0), covparams_prior_params)
      new_results <- adam_sgld_update(phi = log(covparams0), grad_t = n*grad_phi + prior_grad,
                                      momentum_t = momentum, V_t = V,  beta_1 = beta1_factor,
                                      beta_2 = beta2_factor, lr = lr,
                                      bias_factor = bias_factor)
      new_phi <- new_results$new_phi
      new_momentum <- new_results$momentum
      new_V <- new_results$V
      counter_while <- counter_while + 1
    }
    if(counter_while> 10) cat(paste("Reduced the step size by ", counter_while, ".\n"))
    covparams0 <- drop(exp(new_phi))
    V <- new_V
    momentum <- new_momentum
    beta_0 <- pass_list_all_batch_ordered$betahat
    if(iter_count>n_burn){
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

