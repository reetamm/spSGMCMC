#library(testsgld)
####
#### Function to run different flavors of preconditionned SGRLD
####
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

r_sgrld_step_beta <- function(lr, info, beta_c, beta_hat){
  chol_info <- chol(info)
  noise <- backsolve(chol_info, rnorm(length(beta_c), mean = 0, sd = sqrt(lr) ))
  return( (1-lr)*beta_c + lr*beta_hat + drop(noise) )
}

###
### Funtion to use fisher Info as preconditionner
###
sgrld_no_gamma <- function(y, X, NNarray, locs, beta_0, covparams0, covfun_name = "matern_isotropic",
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

    # while(any( is.infinite(exp(covparams_new))) || any(is.na(exp(covparams_new))) ){
    #   if (counter_while==1) {
    #     cat("Need to reduce step size.\n")
    #     # which cov_params
    #     if(length(lr) ==1) lr  <- rep(lr, length(covparams0))
    #     current_theta <- exp(covparams_new)
    #     inf_id <- which(is.infinite(current_theta))
    #     na_id <- which(is.na(current_theta))
    #     if(length(inf_id)>0){
    #       cat(paste("Indexes with Infinite values ", inf_id, ".\n", sep="" ))
    #       lr[inf_id] <- lr[inf_id]/10
    #     }
    #     if(length(na_id)>0){
    #       cat(paste("Indexes with NaN values ", na_id, ".\n", sep="" ))
    #       lr[na_id] <- lr[na_id]/10
    #     }
    # 
    #   }else{
    #     current_theta <- exp(covparams_new)
    #     inf_id <- which(is.infinite(current_theta))
    #     na_id <- which(is.na(current_theta))
    #     if(length(inf_id)>0){
    #       #cat(paste("Indexes with Infinite values ", inf_id, ".\n", sep="" ))
    #       lr[inf_id] <- lr[inf_id]/10
    #     }
    #     if(length(na_id)>0){
    #       #cat(paste("Indexes with Infinite values ", na_id, ".\n", sep="" ))
    #       lr[na_id] <- lr[na_id]/10
    #     }
    #   }
    # 
    #   covparams_new = r_sgrld_step(lr = lr, logparms = log(covparams0), info = info_phi, grad = grad_phi)
    #   if(counter_while==200){
    #     cat(paste("Learning rate might be too large!\n"))
    #     break()
    #   }
    #   counter_while <- counter_while + 1
    # }
    # if(counter_while> 10) cat(paste("Reduced the step size by ", counter_while, ".\n"))
    # covparams0 <- covparams_new
    # covparams0 <- drop(exp(covparams0))
    # if(counter_while>100){
    #   cat("Restarting sampling.\n")
    #   if(covparams0[1]>8 || covparams0[1]<.5) covparams0[1] <- 1.
    #   if(covparams0[2]>3 || covparams0[2]<0.1) covparams0[2] <- 0.5
    #   if(covparams0[3]>8 || covparams0[3]<0.5) covparams0[3] <- 1.
    #   if(covparams0[4]>10 || covparams0[4]<0.05) covparams0[4] <- 0.5
    #   #browser()
    # }
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



###
### Function for diagonal adaptive RMSPROP sgld
###
sgld_rmsprop <- function(y, X, NNarray, locs, beta_0, covparams0, covfun_name = "matern_isotropic",
                         lr = 1e-3, lr_min = 2e-6, n_epochs=100, n_batch = 250, n_burn = 2000,
                         thin = 5, alpha = 0.99, lambda=1e-5, covparams_prior_params, silent = F, 
                         initial_V = NULL, diag_V = T, diag_G = T){
  if(is.null(initial_V)){
    V <-0
    diag_V <- T
    diag_G <-T
  }else{
    V <- initial_V
  }
  #V = 0
  add_noise <- T
  n_iter <- n_epochs*length(y)/n_batch
  nmc_samples <- (n_iter - n_burn)/thin
  covparams_matrix <- matrix(0, nrow=nmc_samples, ncol = 4)
  beta_matrix <- matrix(0, nrow=nmc_samples, ncol = 2)
  llk_matrix <- rep(0, nmc_samples)
  n <- length(y)
  tic <- proc.time()
  iter_count <- 1
  restart_count <- 0
  while(iter_count <= n_iter && restart_count <=10 ){
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
    new_cond <- update_preconditioner(V, grad_phi, alpha, lambda, diag_V = diag_V)
    V <- new_cond$V
    G <- new_cond$G
    new_phi <- updat_state(log(covparams0), grad = n*grad_phi + prior_grad, G = G, lr = lr, add_noise = add_noise, diag_G = diag_G)
    counter_while <- 1
    ##
    # while(any( is.infinite(exp(new_phi))) || any(is.na(exp(new_phi))) ){
    #   if (counter_while==1) {
    #     cat("Need to reduce step size.\n")
    #     # which cov_params
    #     if(length(lr) ==1) lr  <- rep(lr, length(covparams0))
    #     current_theta <- exp(new_phi)
    #     inf_id <- which(is.infinite(current_theta))
    #     na_id <- which(is.na(current_theta))
    #     if(length(inf_id)>0){
    #       cat(paste("Indexes with Infinite values ", inf_id, ".\n", sep="" ))
    #       lr[inf_id] <- lr[inf_id]/2
    #     }
    #     if(length(na_id)>0){
    #       cat(paste("Indexes with NaN values ", na_id, ".\n", sep="" ))
    #       lr[na_id] <- lr[na_id]/2
    #     }
    # 
    #   }else{
    #     current_theta <- exp(new_phi)
    #     inf_id <- which(is.infinite(current_theta))
    #     na_id <- which(is.na(current_theta))
    #     if(length(inf_id)>0){
    #       #cat(paste("Indexes with Infinite values ", inf_id, ".\n", sep="" ))
    #       lr[inf_id] <- lr[inf_id]/2
    #     }
    #     if(length(na_id)>0){
    #       #cat(paste("Indexes with Infinite values ", na_id, ".\n", sep="" ))
    #       lr[na_id] <- lr[na_id]/2
    #     }
    #   }
    #   new_phi <- updat_state(log(covparams0), grad = grad_phi, G = G, lr = lr, add_noise = add_noise)
    #   if(counter_while==100){
    #     cat(paste("Learning rate might be too large!\n"))
    #     #browser()
    #   }
    #   counter_while <- counter_while + 1
    # }
    # if(counter_while> 10) cat(paste("Reduced the step size by ", counter_while, ".\n"))
    # #covparams0 <- new_phi
    # #covparams0 <- drop(exp(covparams0))
    # if(counter_while>100){
    #   cat("Restarting sampling.\n")
    #   if(covparams0[1]>8 || covparams0[1]<.5) covparams0[1] <- 1.
    #   if(covparams0[2]>3 || covparams0[2]<0.1) covparams0[2] <- 0.5
    #   if(covparams0[3]>8 || covparams0[3]<0.5) covparams0[3] <- 1.
    #   if(covparams0[4]>10 || covparams0[4]<0.05) covparams0[4] <- 0.4
    # }
    while( any( is.infinite(exp(new_phi))) || any(is.na(exp(new_phi)))){
      cat("Need to reduce step size.\n")
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
      lr = lr/2
      lr_min = lr/10
      iter_count <- 1
      restart_count <- restart_count + 1
      if(restart_count > 10) break
      V <- initial_V
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
      new_cond <- update_preconditioner(V, grad_phi, alpha, lambda, diag_V = diag_V)
      V <- new_cond$V
      G <- new_cond$G
      new_phi <- updat_state(log(covparams0), grad = n*grad_phi + prior_grad, G = G, lr = lr, add_noise = add_noise, diag_G = diag_G)
      counter_while <- counter_while + 1
    }
    if(counter_while> 10) cat(paste("Reduced the step size by ", counter_while, ".\n"))
    covparams0 <- new_phi
    covparams0 <- drop(exp(covparams0))
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


###
### Funtion to use fisher Info as preconditionner
###
sgrld_fixed_G <- function(y, X, NNarray, locs, beta_0, covparams0, covfun_name = "matern_isotropic",
                           lr = 1e-3, lr_min = 2e-6, n_epochs=100, n_batch = 250, n_burn = 2000,
                           thin = 5, covparams_prior_params, silent = F, gpgp_subsample_size = 5000){

  ###
  ### First run GpGp fit model
  ###

  n_iter <- n_epochs*length(y)/n_batch
  nmc_samples <- (n_iter - n_burn)/thin
  covparams_matrix <- matrix(0, nrow=nmc_samples, ncol = 4)
  beta_matrix <- matrix(0, nrow=nmc_samples, ncol = 2)
  llk_matrix <- rep(0, nmc_samples)
  n <- length(y)
  tic <- proc.time()
  sub_sample <- sample(1:length(y), size = gpgp_subsample_size, replace = F)
  gpgp_fit <- GpGp::fit_model(y[sub_sample], locs[sub_sample,], X[sub_sample,], covfun_name, reorder = T, m_seq = c(30), max_iter = 20,silent = F)
  beta_0 <- gpgp_fit$betahat; covparams0 <- gpgp_fit$covparms; info_theta <- gpgp_fit$info
  for(i in 1:n_iter){
    batch_ind <- sample(1:n, size = n_batch, replace = F)
    ordered_batch_ind <- sort(batch_ind)

    pass_list_all_batch_ordered <- vecchia_profbeta_loglik_grad_info(batch_id = ordered_batch_ind, covparms = covparams0,
                                                                        covfun_name = "matern_isotropic", y = y,
                                                                        X = X, current_beta = beta_0, locs = locs, NNarray = NNarray)
    grad_theta <- (n/n_batch)*pass_list_all_batch_ordered$grad
    tol <- 1e-10
    grad_phi <- covparams0*grad_theta + transformed_matern_parms_logprior_grad(log(covparams0), covparams_prior_params)
    if( any(is.infinite(grad_phi)) || any(is.na(grad_phi)) ){
      cat(paste("Iteration ", i, " non-numeric results in gradient.\n"))
      cat(paste("Current covparams: ", covparams0, ".\n", sep=""))
      stop()
    }
    info_phi <- diag(covparams0)%*%info_theta%*%diag(covparams0)
    if (condition_number(info_phi) > 1 / tol) {
      #if (!silent) cat("Cond # of info matrix > 1/tol \n")
      #info <- 1.0*max(likobj0$info)*diag(nrow(likobj0$info))
      # regularize
      ee <- eigen(info_phi)
      ee_ratios <- ee$values/max(ee$values)
      ee_ratios[ ee_ratios < 1e-5 ] <- 1e-5
      ee$values <- max(ee$values)*ee_ratios
      info_phi <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)

      #diag(info) <- diag(info) + tol*max(diag(info))
    }
    #reparameterized_quantities(covparams0, pass_list_all_batch_ordered$grad_t, pass_list_all_batch_ordered$info, grad_phi, info_ph
    covparams_new = r_sgrld_step(lr = lr, logparms = log(covparams0), info = info_phi, grad = grad_phi)
    counter_while <- 1
    while(any( is.infinite(exp(covparams_new))) || any(is.na(exp(covparams_new))) ){
      if (counter_while==1) {
        cat("Need to reduce step size.\n")
        # which cov_params
        if(length(lr) ==1) lr  <- rep(lr, length(covparams0))
        current_theta <- exp(covparams_new)
        inf_id <- which(is.infinite(current_theta))
        na_id <- which(is.na(current_theta))
        if(length(inf_id)>0){
          cat(paste("Indexes with Infinite values ", inf_id, ".\n", sep="" ))
          lr[inf_id] <- lr[inf_id]/10
        }
        if(length(na_id)>0){
          cat(paste("Indexes with NaN values ", na_id, ".\n", sep="" ))
          lr[na_id] <- lr[na_id]/10
        }

      }else{
        current_theta <- exp(covparams_new)
        inf_id <- which(is.infinite(current_theta))
        na_id <- which(is.na(current_theta))
        if(length(inf_id)>0){
          #cat(paste("Indexes with Infinite values ", inf_id, ".\n", sep="" ))
          lr[inf_id] <- lr[inf_id]/10
        }
        if(length(na_id)>0){
          #cat(paste("Indexes with Infinite values ", na_id, ".\n", sep="" ))
          lr[na_id] <- lr[na_id]/10
        }
      }
      covparams_new = r_sgrld_step(lr = lr, logparms = log(covparams0), info = info_phi, grad = grad_phi)
      if(counter_while==100){
        cat(paste("Learning rate might be too large!\n"))
        #browser()
      	#stop()
      }
      counter_while <- counter_while + 1
    }

    if(counter_while> 10) cat(paste("Reduced the step size by ", counter_while, ".\n"))
    covparams0 <- covparams_new
    covparams0 <- drop(exp(covparams0))
    if(counter_while>200){
      cat("Restarting sampling.\n")
      if(covparams0[1]>8 || covparams0[1]<1) covparams0[1] <- 2.
      if(covparams0[2]>3 || covparams0[2]<0.3) covparams0[2] <- 1.9
      if(covparams0[3]>5 || covparams0[3]<0.5) covparams0[3] <- 1.
      if(covparams0[4]>2 || covparams0[4]<0.05) covparams0[4] <- 0.4
    }
    beta_0 <- pass_list_all_batch_ordered$betahat
    if(i>n_burn){
      j = i - n_burn
      if(j%%thin==0){
        beta_matrix[j%/%thin,] <- beta_0
        covparams_matrix[j%/%thin,] <- covparams0
        llk_matrix[j%/%thin] <- pass_list_all_batch_ordered$loglik
      }
    }


    if(i%%(n%/%n_batch)==0) {

      if(i%%(10*n%/%n_batch)==0){
        epoch_c <- i/(10*n%/%n_batch)
        if(lr >= lr_min) lr = lr/epoch_c

        cat(paste("Iteration ", i, ", Epoch ", i%/%(n%/%n_batch), sep = ""), "\n")
      }
    }
  }
  toc <- proc.time()
  return( list(beta_samples = beta_matrix, theta_samples = covparams_matrix, llk_trace = llk_matrix, elapsed_time = toc["elapsed"] - tic["elapsed"]))
}





test_funs <- F
if (test_funs){
  set.seed(12345)

  n1 <- 100
  n2 <- 100
  n <- n1*n2
  locs <- as.matrix( expand.grid( 1:n1, 1:n2 ) )
  ord <- order_maxmin(locs)
  locsord <- locs[ord,]
  m <- 30
  NNarray <- find_ordered_nn(locsord,m=m)
  Xord <- as.matrix(cbind(rep(1,n), runif(n)))
  beta_true <- c(3, -5)
  mu <- Xord%*%beta_true


  covparms <- c(5,1,1.5, 0.2)
  spatial_effect <- fast_Gp_sim(covparms,"matern_isotropic",locsord)

  y_resp <- drop(mu) + spatial_effect #+ rnorm(length(spatial_effect))



  ###
  ### Now start the algorithm
  ###
  covparams0 <- GpGp::get_start_parms(y_resp, Xord, locsord, covfun_name = "matern_isotropic")$start_parms
  sub_sample <- sort(sample(1:n,500))
  all_active <- unique(as.vector(NNarray[sub_sample,]))
  non_nas <- which(!is.na(all_active))
  all_active <- all_active[non_nas]
  all_active <- sort(all_active)
  gpgp_fit <- GpGp::fit_model(y_resp[all_active], locsord[all_active,], Xord[all_active,], covfun_name, reorder = F, m_seq = c(15), max_iter = 20,silent = F, convtol = 1e-8)
  beta_c <- gpgp_fit$betahat; covparams0 <- gpgp_fit$covparms; initial_V <- gpgp_fit$info
  if(covparams0[1]>8 || covparams0[1]<1) covparams0[1] <- 2.
  if(covparams0[2]>3 || covparams0[2]<0.1) covparams0[2] <- 1.9
  if(covparams0[3]>8 || covparams0[3]<0.5) covparams0[3] <- 1.
  if(covparams0[4]>20 || covparams0[4]<0.05) covparams0[4] <- 0.5

  beta_c <- c(1., -3.)
  covparams_prior_params <- cbind(c(.1, 9, 1, .1), c(.1, 2, 1, .1) )
  test_no_gamma <- sgrld_no_gamma(y=y_resp, X = Xord, NNarray = NNarray, locs = locsord, beta_0 = beta_c,
                                  covparams0 = covparams0, covfun_name ="matern_isotropic", lr = lr_min_rmsprop,
                                  lr_min = lr_min_rmsprop, n_epochs = 50, n_batch = n_batch, n_burn = n_burn,
                                  thin = thin, covparams_prior_params = covparams_prior_params, silent = F)
  op <- par(no.readonly = T)
  library(latex2exp)
  par(mfrow = c(2,2))
  plot(test_no_gamma$theta_samples[,1], xlab = "", ylab = TeX("$\\sigma^2$"), main = "", type = "l",  )
  abline(h=covparms[1], col="red")
  plot(test_no_gamma$theta_samples[,2], xlab = "", ylab = TeX("$\\rho$"), main = "", type = "l",  )
  abline(h=covparms[2], col="red")
  plot(test_no_gamma$theta_samples[,3], xlab = "", ylab = TeX("$\\nu$"), main = "", type = "l",  )
  abline(h=covparms[3], col="red")
  plot(test_no_gamma$theta_samples[,4]*test_no_gamma$theta_samples[,1], xlab = "", ylab = TeX("$\\tau^2$"), main = "", type = "l",)
  abline(h=covparms[4]*covparms[1], col="red")
  par(mfrow = c(2,1))
  plot(test_no_gamma$beta_samples[,1], xlab = "", ylab = TeX("$\\beta_1$"), main = "", type = "l" )
  abline(h = beta_true[1], col = "red")
  plot(test_no_gamma$beta_samples[,2], xlab = "", ylab = TeX("$\\beta_2$"), main = "", type = "l" )
  abline(h=beta_true[2], col = "red")
  
  test_rmsprop <- sgld_rmsprop(y=y_resp, X = Xord, NNarray = NNarray, locs = locsord, beta_0 = beta_c,
                              covparams0 = covparams0, covfun_name ="matern_isotropic", lr = lr_min_rmsprop,
                              lr_min = lr_min_rmsprop, n_epochs = 50, n_batch = n_batch, n_burn = n_burn,
                              thin = thin, covparams_prior_params = covparams_prior_params, alpha =0.9,
                              lambda = 1e-5, silent = F,initial_V = diag(initial_V), diag_V = T, diag_G = T)
  
  par(mfrow = c(2,2))
  plot(test_rmsprop$theta_samples[,1], xlab = "", ylab = TeX("$\\sigma^2$"), main = "", type = "l",  )
  abline(h=covparms[1], col="red")
  plot(test_rmsprop$theta_samples[,2], xlab = "", ylab = TeX("$\\rho$"), main = "", type = "l",  )
  abline(h=covparms[2], col="red")
  plot(test_rmsprop$theta_samples[,3], xlab = "", ylab = TeX("$\\nu$"), main = "", type = "l",  )
  abline(h=covparms[3], col="red")
  plot(test_rmsprop$theta_samples[,4]*test_rmsprop$theta_samples[,1], xlab = "", ylab = TeX("$\\tau^2$"), main = "", type = "l",)
  abline(h=covparms[4]*covparms[1], col="red")
  par(mfrow = c(2,1))
  plot(test_rmsprop$beta_samples[,1], xlab = "", ylab = TeX("$\\beta_1$"), main = "", type = "l" )
  abline(h = beta_true[1], col = "red")
  plot(test_rmsprop$beta_samples[,2], xlab = "", ylab = TeX("$\\beta_2$"), main = "", type = "l" )
  abline(h=beta_true[2], col = "red")
  
  test_fixed_G <- sgrld_fixed_G(y=y_resp, X = Xord, NNarray = NNarray, locs = locsord, beta_0 = beta_c,
                                covparams0 = covparams0, covfun_name ="matern_isotropic", lr = lr,
                                lr_min = lr_min, n_epochs = n_epoch, n_batch = n_batch, n_burn = n_burn,
                                thin = thin, covparams_prior_params = covparams_prior_params, silent = F,
                                gpgp_subsample_size = 1000)



  op <- par(no.readonly = T)
  library(latex2exp)
  par(mfrow = c(2,2))
  plot(test_rmsprop$theta_samples[,1], xlab = "", ylab = TeX("$\\sigma^2$"), main = "", type = "l",  )
  abline(h=covparms[1], col="red")
  plot(test_rmsprop$theta_samples[,2], xlab = "", ylab = TeX("$\\rho$"), main = "", type = "l",  )
  abline(h=covparms[2], col="red")
  plot(test_rmsprop$theta_samples[,3], xlab = "", ylab = TeX("$\\nu$"), main = "", type = "l",  )
  abline(h=covparms[3], col="red")
  plot(test_rmsprop$theta_samples[,4]*test_rmsprop$theta_samples[,1], xlab = "", ylab = TeX("$\\tau^2$"), main = "", type = "l",)
  abline(h=covparms[4]*covparms[1], col="red")
  par(mfrow = c(2,1))
  plot(test_rmsprop$beta_samples[,1], xlab = "", ylab = TeX("$\\beta_1$"), main = "", type = "l" )
  abline(h = beta_true[1], col = "red")
  plot(test_rmsprop$beta_samples[,2], xlab = "", ylab = TeX("$\\beta_2$"), main = "", type = "l" )
  abline(h=beta_true[2], col = "red")
  par(mfrow=c(1,1))
  plot(test_rmsprop$llk_trace, xlab = "", ylab = TeX("$\\ell(\\theta)$"), main = "", type = "l")


  matern_vals <- fields::Matern(d = seq(0, 100, by = 0.5), range = 2, smoothness = 1.5 )
  par(mfrow = c(1,1))
  plot(seq(0, 100, by = 0.5), matern_vals, type = "l")
}

# i <- 1
# counter <- 1
# while( i <= 10 && counter <= 20){
#   if(i>7){
#     i <-1
#   }
#   cat(paste("i = ", i, ".\n", sep=""))
#   i <- i+1
#   counter <- counter + 1
# }
