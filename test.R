# install.packages("GpGp")

rm(list = ls())
library(ggplot2)
library(tidyverse)
library(viridis)
# library(GpGp)
library(spSGMCMC)
## Open and plot the data
argo2016 <- GpGp::argo2016
argo2016$lon <- ifelse(argo2016$lon>180,argo2016$lon-360,argo2016$lon)
world <- map_data("world")
argoplot <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "gray", fill = "white", linewidth = 0.01
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
#> Warning in geom_map(data = world, map = world, aes(long, lat, map_id = region),
#> : Ignoring unknown aesthetics: x and y
argoplot


##### Split into train/test
test_prop <- 0.6
test_size <- floor(nrow(argo2016)*test_prop)
id_test <- sample(1:nrow(argo2016), size = test_size)
id_train <- 1:nrow(argo2016)
id_train <- id_train[-id_test]
argo2016_test <- argo2016[id_test,]
argo2016_train <- argo2016[id_train,]

#### gpgp fit

y <- argo2016_train$temp100
X <- cbind(1, argo2016_train$lon, argo2016_train$lat, argo2016_train$lon^2, argo2016_train$lat^2)
locs = cbind(argo2016_train$lon, argo2016_train$lat)
ord <- GpGp::order_maxmin(locs = locs, lonlat = T)
NNarray <- GpGp::find_ordered_nn(locs = locs[ord,], m = 10, lonlat = T)


##### sgrld fit
n_epoch <- 200
n_batch <- 250
lr_sgrld = 1e-4; lr_min_sgrld = 5e-5
n_burn = 5000
thin = 1
covfun_name = "matern_isotropic"
covparams_prior_params <- cbind(c(.1, 100, 1, .1), c(.1, 2, 1, .1) )
# Note that the prior for range is Gamma(100,2), which gives us a posterior mean of 50. Apt for very smooth global data
# initial values from GpGp
aaa <- get_start_parms(y[ord],X[ord,],locs[ord,],'matern_isotropic')
beta_c <- aaa$betahat; covparams0 <- aaa$covparams
sgrld_fit <- fit_model_sgmcmc(y=y[ord], X = X[ord,], NNarray = NNarray, locs = locs[ord,], beta_0 = beta_c,
                              algorithm = 'SGRLD',
                        covparams0 = covparams0, covfun_name ="matern_isotropic", lr = lr_sgrld,
                        lr_min = lr_min_sgrld, n_epochs = n_epoch, n_batch = n_batch, n_burn = n_burn,
                        thin = thin, covparams_prior_params = covparams_prior_params, silent = F)

#### Total time taken
sgrld_fit$elapsed_time
apply(sgrld_fit$beta_samples,2,mean)
apply(sgrld_fit$theta_samples,2,mean)
