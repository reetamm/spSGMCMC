rm(list=ls())
#library(GpGp)
library(ggplot2)
library(tidyverse)
library(viridis)
library(spSGMCMC)
####
#### Plot data
####
argo2016 <- GpGp::argo2016
argo2016$lon <- ifelse(argo2016$lon>180,argo2016$lon-360,argo2016$lon)
world <- map_data("world")
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "gray", fill = "white", size = 0.01
  ) +
  geom_point(
    data = argo2016,
    aes(lon, lat, color = temp100),
    alpha = 0.7,size=0.1
  ) + coord_fixed()+theme_bw() +
  scale_color_viridis_c(option = "C")+
  xlab("Longitude")+
  ylab("Latitiude")+
  labs(colour = "Temp (C)")

#####
##### Split into train/test
#####
#####
test_prop <- 0.6
test_size <- floor(nrow(argo2016)*test_prop)
id_test <- sample(1:nrow(argo2016), size = test_size)
id_train <- 1:nrow(argo2016)
id_train <- id_train[-id_test]
argo2016_test <- argo2016[id_test,]
argo2016_train <- argo2016[id_train,]

####
#### gpgp fit
####

y <- argo2016_train$temp100
X <- cbind(1, argo2016_train$lon, argo2016_train$lat, argo2016_train$lon^2, argo2016_train$lat^2)
locs = cbind(argo2016_train$lon, argo2016_train$lat)
ord <- GpGp::order_maxmin(locs = locs, lonlat = T)
NNarray <- GpGp::find_ordered_nn(locs = locs[ord,], m = 10, lonlat = T)
gpgp_fit <- GpGp::fit_model(y[ord], locs = locs[ord,], X = X[ord,], covfun_name = "matern_isotropic", NNarray = NNarray,
                      reorder = F, silent = F, m_seq = c(10))



#####
##### Now fit the SGRLD method
#####
n_epoch <- 200
n_batch <- 250
lr_no_gamma = 1e-4; lr_min_no_gamma = 5e-5
lr_rmsprop = 1e-6; lr_min_rmsprop = 1e-7
lr_msgld = 1e-6; lr_min_msgld = 1e-7
lr_adam = 1e-6; lr_min_adam = 1e-7
lr_fixed_G = 1e-6; lr_min_fixed_G = 1e-6
n_burn = 5000
thin = 1
covfun_name = "matern_isotropic"
covparams_prior_params <- cbind(c(.01, 100, 1, .1), c(.01, 2, 1, .1) )
covparms <- c(14, 50., .5, .03)

beta_c <- c(0., 0.)
# out_file <- file("~/NCSU/test_sgld_output/drifts_output/sgrld_argo_output.txt", open = "wt")
# sink(out_file ,type = "output", append = F)
# sink(out_file, type = "message", append = F)
#
# tic <- proc.time()
beta_c <- gpgp_fit$betahat; covparams0 <- gpgp_fit$covparms; initial_V <- diag(gpgp_fit$info)
# library(SGRLD)#library(testsgld)
# source("./../../test_sgld_output/sgmcmc_functions.R")
# source("./../sgmcmc_functions.R")
sgrld_fit <- spsgmcmc(y=y[ord], X = X[ord,], NNarray = NNarray, locs = locs[ord,], beta_0 = beta_c,
                      covparams0 = covparams0, covfun_name ="matern_isotropic", lr = lr_no_gamma,
                      lr_min = lr_min_no_gamma, n_epochs = n_epoch, n_batch = n_batch, n_burn = n_burn,
                      thin = thin, covparams_prior_params = covparams_prior_params, silent = F)
# toc <- proc.time()
# sink(type="message")
# sink(type="output")
# close(out_file)
# toc-tic

####
#### Plot chains
####
par(mfrow = c(2,2))
library(latex2exp)
names_sgrld_plots <- c("$\\sigma^2$", "$\\rho$", "$\\nu$", "$\\tau^2$")
for(j in 1:4) plot(sgrld_fit$theta_samples[,j], ylab = TeX(names_sgrld_plots[j]), type = "l")
names_sgrld_beta <- c("$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$")
for(j in 1:4) plot(sgrld_fit$beta_samples[,j+1], ylab = TeX(names_sgrld_beta[j]), type = "l")

sgrld_betahat <- colMeans(sgrld_fit$beta_samples); sgrld_thetahat <- colMeans(sgrld_fit$theta_samples)

####
#### Make predictions
####
predicted_sgrld_mean <- GpGp::predictions(locs_pred = locs_pred, X_pred = X_pred, y_obs = y[ord],
                               locs_obs = locs[ord,], X_obs = X[ord,], beta = gpgp_fit$betahat,
                               covparms = sgrld_thetahat, covfun_name = "matern_isotropic",
                               m =15, reorder = FALSE )
par(mfrow = c(1,1))
plot(y_pred, predicted_sgrld_mean, xlab = "True values", ylab = "Predicted values")
abline(a=0, b = 1, col = "red")
#Now predicted CI
library(parallel)
quantile_theta <- apply(sgrld_fit$theta_samples, MARGIN = 2, FUN = function(theta){
  return(quantile(theta, probs = seq(0.01, 0.99, 0.02)))})
predicted_sgrld <- apply(quantile_theta, MARGIN = 1, FUN = function(theta) {
  return(GpGp::cond_sim(locs_pred = locs_pred, X_pred = X_pred, y_obs = y[ord],
                                 locs_obs = locs[ord,], X_obs = X[ord,], beta = gpgp_fit$betahat,
                                 covparms = theta, covfun_name = "matern_isotropic",
                                 m =15, reorder = FALSE, nsims =1))
}
)
predicted_CI_sgrld <- apply(predicted_sgrld, 1, FUN = function(x){return(quantile(x, probs = c(0.025, 0.5, 0.975)))})
predicted_CI_sgrld <- t(predicted_CI_sgrld)

coverage_sgrld <- coverage_fun(predicted_CI_sgrld[, c(1,3)], y_pred)
cat(paste("Average coverage of SGRLD on the test set is: ", round(mean(coverage_sgrld),4), ".\n", sep = ""))
cat(paste("Mean squared error of SGRLD on the test set is: ", round( mean( (predicted_sgrld_mean - y_pred)^2 ), 4), ".\n", sep = ""))
cat(paste("Mean squared error of SGRLD median on the test set is: ",
          round( mean( (predicted_CI_sgrld[,2] - y_pred)^2 ), 4), ".\n", sep=""))
sgrld_ESS_beta <- coda::effectiveSize(x = sgrld_fit$beta_samples)/(sgrld_fit$elapsed_time/60)
sgrld_ESS_theta <- coda::effectiveSize(x = sgrld_fit$theta_samples)/(sgrld_fit$elapsed_time/60)
cat(paste("SGRLD ESS/min of $\\beta$ is:\n"));round(sgrld_ESS_beta, digits = 4)
cat(paste("SGRLD ESS/min of $\\theta$ is:\n"));round(sgrld_ESS_theta, digits = 4)



# #####
# ##### Now fit the SGRLD method
# #####
# n_epoch <- 400
# n_batch <- 250
# lr_no_gamma = 1e-4; lr_min_no_gamma = 1e-4
# lr_rmsprop = 1e-6; lr_min_rmsprop = 1e-7
# lr_msgld = 1e-6; lr_min_msgld = 1e-7
# lr_adam = 1e-6; lr_min_adam = 1e-7
# lr_fixed_G = 1e-6; lr_min_fixed_G = 1e-6
# n_burn = 5000
# thin = 1
# covfun_name = "matern_isotropic"
# covparams_prior_params <- cbind(c(.01, 100, 1, .1), c(.01, 2, 1, .1) )
# covparms <- c(14, 50., .5, .03)
#
# beta_c <- c(0., 0.)
# out_file <- file("~/NCSU/test_sgld_output/drifts_output/sgrld2_argo_output.txt", open = "wt")
# sink(out_file ,type = "output", append = F)
# sink(out_file, type = "message", append = F)
#
# tic <- proc.time()
# beta_c <- gpgp_fit$betahat; covparams0 <- gpgp_fit$covparms; initial_V <- diag(gpgp_fit$info)
# library(SGRLD)
# source("./../sgrld_sampling.R")
# sgrld2_fit <- tryCatch(
#   expr = {
#     sample_sgrld(y=y[ord], X = X[ord,], NNarray = NNarray, locs = locs[ord,], beta_0 = beta_c,
#                    covparams0 = covparams0, covfun_name ="matern_isotropic", lr = lr_no_gamma,
#                    lr_min = lr_min_no_gamma, lr_beta = 0.05, lr_beta_min = 0.01,
#                  n_epochs = n_epoch, n_batch = n_batch, n_burn = n_burn,
#                    thin = thin, covparams_prior_params = covparams_prior_params, silent = F)
#   },error=function(error_message) {
#     #message("This is my custom message.")
#     #message("And below is the error message from R:")
#     message(error_message)
#   }
# )
# toc <- proc.time()
# sink(type="message")
# sink(type="output")
# close(out_file)
#
#
# ####
# #### Plot chains
# ####
# par(mfrow = c(2,2))
# library(latex2exp)
# names_sgrld_plots <- c("$\\sigma^2$", "$\\rho$", "$\\nu$", "$\\tau^2$")
# for(j in 1:4) plot(sgrld2_fit$theta_samples[,j], ylab = TeX(names_sgrld_plots[j]), type = "l")
# names_sgrld_beta <- c("$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$")
# for(j in 1:4) plot(sgrld2_fit$beta_samples[,j+1], ylab = TeX(names_sgrld_beta[j]), type = "l")
#
# sgrld2_betahat <- colMeans(sgrld2_fit$beta_samples); sgrld2_thetahat <- colMeans(sgrld2_fit$theta_samples)
#
# ####
# #### Make predictions
# ####
# predicted_sgrld2_mean <- GpGp::predictions(locs_pred = locs_pred, X_pred = X_pred, y_obs = y[ord],
#                                           locs_obs = locs[ord,], X_obs = X[ord,], beta = gpgp_fit$betahat,
#                                           covparms = sgrld2_thetahat, covfun_name = "matern_isotropic",
#                                           m =15, reorder = FALSE )
# par(mfrow = c(1,1))
# plot(y_pred, predicted_sgrld2_mean, xlab = "True values", ylab = "Predicted values")
# abline(a=0, b = 1, col = "red")
# #Now predicted CI
# library(parallel)
# quantile_theta <- apply(sgrld2_fit$theta_samples, MARGIN = 2, FUN = function(theta){
#   return(quantile(theta, probs = seq(0.01, 0.99, 0.005)))})
# predicted_sgrld2 <- apply(quantile_theta, MARGIN = 1, FUN = function(theta) {
#   return(GpGp::cond_sim(locs_pred = locs_pred, X_pred = X_pred, y_obs = y[ord],
#                         locs_obs = locs[ord,], X_obs = X[ord,], beta = gpgp_fit$betahat,
#                         covparms = theta, covfun_name = "matern_isotropic",
#                         m =15, reorder = FALSE, nsims =1))
# }
# )
# predicted_CI_sgrld2 <- apply(predicted_sgrld2, 1, FUN = function(x){return(quantile(x, probs = c(0.025, 0.5, 0.975)))})
# predicted_CI_sgrld2 <- t(predicted_CI_sgrld2)
#
# coverage_sgrld2 <- coverage_fun(predicted_CI_sgrld2[, c(1,3)], y_pred)
# cat(paste("Average coverage of SGRLD 2 on the test set is: ", round(mean(coverage_sgrld2),4), ".\n", sep = ""))
# cat(paste("Mean squared error of SGRLD 2 on the test set is: ", round( mean( (predicted_sgrld2_mean - y_pred)^2 ), 4), ".\n", sep = ""))
# cat(paste("Mean squared error of SGRLD median on the test set is: ",
#           round( mean( (predicted_CI_sgrld2[,2] - y_pred)^2 ), 4), ".\n", sep=""))
# sgrld2_ESS_beta <- coda::effectiveSize(x = sgrld2_fit$beta_samples)/(sgrld2_fit$elapsed_time/60)
# sgrld2_ESS_theta <- coda::effectiveSize(x = sgrld2_fit$theta_samples)/(sgrld2_fit$elapsed_time/60)
# cat(paste("SGRLD ESS/min of $\\beta$ is:\n"));round(sgrld2_ESS_beta, digits = 4)
# cat(paste("SGRLD ESS/min of $\\theta$ is:\n"));round(sgrld2_ESS_theta, digits = 4)


#####
##### Conditionning 10
#####
NNarray_m10 <- GpGp::find_ordered_nn(locs = locs[ord,], m = 10, lonlat = T)
n_epoch <- 400
n_batch <- 250
lr_no_gamma = 1e-4; lr_min_no_gamma = 1e-4
lr_rmsprop = 1e-6; lr_min_rmsprop = 1e-7
lr_msgld = 1e-6; lr_min_msgld = 1e-7
lr_adam = 1e-6; lr_min_adam = 1e-7
lr_fixed_G = 1e-6; lr_min_fixed_G = 1e-6
n_burn = 5000
thin = 1
covfun_name = "matern_isotropic"
covparams_prior_params <- cbind(c(.01, 100, 1, .1), c(.01, 2, 1, .1) )
covparms <- c(14, 50., .5, .03)

beta_c <- c(0., 0.)
out_file <- file("~/NCSU/test_sgld_output/drifts_output/sgrld_m10_argo_output.txt", open = "wt")
sink(out_file ,type = "output", append = F)
sink(out_file, type = "message", append = F)

tic <- proc.time()
beta_c <- gpgp_fit$betahat; covparams0 <- gpgp_fit$covparms; initial_V <- diag(gpgp_fit$info)
library(SGRLD)
source("./../sgrld_sampling.R")
sgrld_m10_fit <- tryCatch(
  expr = {
    sample_sgrld(y=y[ord], X = X[ord,], NNarray = NNarray_m10, locs = locs[ord,], beta_0 = beta_c,
                 covparams0 = covparams0, covfun_name ="matern_isotropic", lr = lr_no_gamma,
                 lr_min = lr_min_no_gamma, lr_beta = 0.05, lr_beta_min = 0.01,
                 n_epochs = n_epoch, n_batch = n_batch, n_burn = n_burn,
                 thin = thin, covparams_prior_params = covparams_prior_params, silent = F)
  },error=function(error_message) {
    #message("This is my custom message.")
    #message("And below is the error message from R:")
    message(error_message)
  }
)
toc <- proc.time()
sink(type="message")
sink(type="output")
close(out_file)


####
#### Plot chains
####
par(mfrow = c(2,2))
library(latex2exp)
names_sgrld_plots <- c("$\\sigma^2$", "$\\rho$", "$\\nu$", "$\\tau^2$")
for(j in 1:4) plot(sgrld_m10_fit$theta_samples[,j], ylab = TeX(names_sgrld_plots[j]), type = "l", main = "Trace plots m=10")
names_sgrld_beta <- c("$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$")
for(j in 1:4) plot(sgrld_m10_fit$beta_samples[,j+1], ylab = TeX(names_sgrld_beta[j]), type = "l", main = "Trace plots m=10")

sgrld_m10_betahat <- colMeans(sgrld_m10_fit$beta_samples); sgrld_m10_thetahat <- colMeans(sgrld_m10_fit$theta_samples)

####
#### Make predictions
####
predicted_sgrld_m10_mean <- GpGp::predictions(locs_pred = locs_pred, X_pred = X_pred, y_obs = y[ord],
                                           locs_obs = locs[ord,], X_obs = X[ord,], beta = gpgp_fit$betahat,
                                           covparms = sgrld_m10_thetahat, covfun_name = "matern_isotropic",
                                           m =10, reorder = FALSE )
par(mfrow = c(1,1))
plot(y_pred, predicted_sgrld_m10_mean, xlab = "True values", ylab = "Predicted values", main="Prdictons m=10")
abline(a=0, b = 1, col = "red")
#Now predicted CI
library(parallel)
quantile_theta <- apply(sgrld_m10_fit$theta_samples, MARGIN = 2, FUN = function(theta){
  return(quantile(theta, probs = seq(0.01, 0.99, 0.005)))})
predicted_sgrld_m10 <- apply(quantile_theta, MARGIN = 1, FUN = function(theta) {
  return(GpGp::cond_sim(locs_pred = locs_pred, X_pred = X_pred, y_obs = y[ord],
                        locs_obs = locs[ord,], X_obs = X[ord,], beta = gpgp_fit$betahat,
                        covparms = theta, covfun_name = "matern_isotropic",
                        m =10, reorder = FALSE, nsims =1))
}
)
predicted_CI_sgrld_m10 <- apply(predicted_sgrld_m10, 1, FUN = function(x){return(quantile(x, probs = c(0.025, 0.5, 0.975)))})
predicted_CI_sgrld_m10 <- t(predicted_CI_sgrld_m10)

coverage_sgrld_m10 <- coverage_fun(predicted_CI_sgrld_m10[, c(1,3)], y_pred)
cat(paste("Average coverage of SGRLD 2 on the test set is: ", round(mean(coverage_sgrld_m10),4), ".\n", sep = ""))
cat(paste("Mean squared error of SGRLD 2 on the test set is: ", round( mean( (predicted_sgrld_m10_mean - y_pred)^2 ), 4), ".\n", sep = ""))
cat(paste("Mean squared error of SGRLD median on the test set is: ",
          round( mean( (predicted_CI_sgrld_m10[,2] - y_pred)^2 ), 4), ".\n", sep=""))
sgrld_m10_ESS_beta <- coda::effectiveSize(x = sgrld_m10_fit$beta_samples)/(sgrld_m10_fit$elapsed_time/60)
sgrld_m10_ESS_theta <- coda::effectiveSize(x = sgrld_m10_fit$theta_samples)/(sgrld_m10_fit$elapsed_time/60)
cat(paste("SGRLD ESS/min of $\\beta$ is:\n"));round(sgrld_m10_ESS_beta, digits = 4)
cat(paste("SGRLD ESS/min of $\\theta$ is:\n"));round(sgrld_m10_ESS_theta, digits = 4)
save.image(file = image_file_name)




####
#### Conditionning 30
####

NNarray_m30 <- GpGp::find_ordered_nn(locs = locs[ord,], m = 30, lonlat = T)
n_epoch <- 400
n_batch <- 250
lr_no_gamma = 1e-4; lr_min_no_gamma = 1e-4
lr_rmsprop = 1e-6; lr_min_rmsprop = 1e-7
lr_msgld = 1e-6; lr_min_msgld = 1e-7
lr_adam = 1e-6; lr_min_adam = 1e-7
lr_fixed_G = 1e-6; lr_min_fixed_G = 1e-6
n_burn = 5000
thin = 1
covfun_name = "matern_isotropic"
covparams_prior_params <- cbind(c(.01, 100, 1, .1), c(.01, 2, 1, .1) )
covparms <- c(14, 50., .5, .03)

beta_c <- c(0., 0.)
out_file <- file("~/NCSU/test_sgld_output/drifts_output/sgrld_m30_argo_output.txt", open = "wt")
sink(out_file ,type = "output", append = F)
sink(out_file, type = "message", append = F)

tic <- proc.time()
beta_c <- gpgp_fit$betahat; covparams0 <- gpgp_fit$covparms; initial_V <- diag(gpgp_fit$info)
library(SGRLD)
source("./../sgrld_sampling.R")
sgrld_m30_fit <- tryCatch(
  expr = {
    sample_sgrld(y=y[ord], X = X[ord,], NNarray = NNarray_m30, locs = locs[ord,], beta_0 = beta_c,
                 covparams0 = covparams0, covfun_name ="matern_isotropic", lr = lr_no_gamma,
                 lr_min = lr_min_no_gamma, lr_beta = 0.05, lr_beta_min = 0.01,
                 n_epochs = n_epoch, n_batch = n_batch, n_burn = n_burn,
                 thin = thin, covparams_prior_params = covparams_prior_params, silent = F)
  },error=function(error_message) {
    #message("This is my custom message.")
    #message("And below is the error message from R:")
    message(error_message)
  }
)
toc <- proc.time()
sink(type="message")
sink(type="output")
close(out_file)


####
#### Plot chains
####
par(mfrow = c(2,2))
library(latex2exp)
names_sgrld_plots <- c("$\\sigma^2$", "$\\rho$", "$\\nu$", "$\\tau^2$")
for(j in 1:4) plot(sgrld_m30_fit$theta_samples[,j], ylab = TeX(names_sgrld_plots[j]), type = "l", main = "Trace plots m=30")
names_sgrld_beta <- c("$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$")
for(j in 1:4) plot(sgrld_m30_fit$beta_samples[,j+1], ylab = TeX(names_sgrld_beta[j]), type = "l", main="Trace plots m=30")

sgrld_m30_betahat <- colMeans(sgrld_m30_fit$beta_samples); sgrld_m30_thetahat <- colMeans(sgrld_m30_fit$theta_samples)

####
#### Make predictions
####
predicted_sgrld_m30_mean <- GpGp::predictions(locs_pred = locs_pred, X_pred = X_pred, y_obs = y[ord],
                                              locs_obs = locs[ord,], X_obs = X[ord,], beta = gpgp_fit$betahat,
                                              covparms = sgrld_m30_thetahat, covfun_name = "matern_isotropic",
                                              m =30, reorder = FALSE )
par(mfrow = c(1,1))
plot(y_pred, predicted_sgrld_m30_mean, xlab = "True values", ylab = "Predicted values", main="Predictions m=30")
abline(a=0, b = 1, col = "red")
#Now predicted CI
library(parallel)
quantile_theta <- apply(sgrld_m30_fit$theta_samples, MARGIN = 2, FUN = function(theta){
  return(quantile(theta, probs = seq(0.01, 0.99, 0.005)))})
predicted_sgrld_m30 <- apply(quantile_theta, MARGIN = 1, FUN = function(theta) {
  return(GpGp::cond_sim(locs_pred = locs_pred, X_pred = X_pred, y_obs = y[ord],
                        locs_obs = locs[ord,], X_obs = X[ord,], beta = gpgp_fit$betahat,
                        covparms = theta, covfun_name = "matern_isotropic",
                        m = 30, reorder = FALSE, nsims =1))
}
)
predicted_CI_sgrld_m30 <- apply(predicted_sgrld_m30, 1, FUN = function(x){return(quantile(x, probs = c(0.025, 0.5, 0.975)))})
predicted_CI_sgrld_m30 <- t(predicted_CI_sgrld_m30)

coverage_sgrld_m30 <- coverage_fun(predicted_CI_sgrld_m30[, c(1,3)], y_pred)
cat(paste("Average coverage of SGRLD 2 on the test set is: ", round(mean(coverage_sgrld_m30),4), ".\n", sep = ""))
cat(paste("Mean squared error of SGRLD 2 on the test set is: ", round( mean( (predicted_sgrld_m30_mean - y_pred)^2 ), 4), ".\n", sep = ""))
cat(paste("Mean squared error of SGRLD median on the test set is: ",
          round( mean( (predicted_CI_sgrld_m30[,2] - y_pred)^2 ), 4), ".\n", sep=""))
sgrld_m30_ESS_beta <- coda::effectiveSize(x = sgrld_m30_fit$beta_samples)/(sgrld_m30_fit$elapsed_time/60)
sgrld_m30_ESS_theta <- coda::effectiveSize(x = sgrld_m30_fit$theta_samples)/(sgrld_m30_fit$elapsed_time/60)
cat(paste("SGRLD ESS/min of $\\beta$ is:\n"));round(sgrld_m30_ESS_beta, digits = 4)
cat(paste("SGRLD ESS/min of $\\theta$ is:\n"));round(sgrld_m30_ESS_theta, digits = 4)
save.image(file = image_file_name)



