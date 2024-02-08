# This file is to compare first trial of predictive performance of spBART and BART
rm(list=ls())
library(dbarts)
library(mlbench)
library(purrr)
library(MOTRbart)
source("R/sim_functions.R")
source("R/main_function.R")
set.seed(42)

n_ <- 250
sd_ <- 1
train <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_) %>% as.data.frame()
test <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_) %>% as.data.frame()

# train <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame() %>% .[,c(1:5,11)]
# test <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame() %>% .[,c(1:5,11)]

# train <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
# test <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()

train <- mlbench.d1.break(n = n_,sd = sd_)  |> as.data.frame()
test <- mlbench.d1.break(n = n_,sd = sd_) |> as.data.frame()

# Getting the training elements
x_train <- train %>% dplyr::select(dplyr::starts_with("x"))
x_test <- test %>% dplyr::select(dplyr::starts_with("x"))
y_train <- train %>% dplyr::pull("y")

# Running the model
spBART <- rspBART(x_train = x_train,
        x_test = x_test,y_train = y_train,
        n_mcmc = 2000,node_min_size = 5,
        n_burn = 0,nIknots = 10,n_tree = 10,
        dif_order = 0,motrbart_bool = FALSE,
        update_tau_beta = TRUE,
        plot_preview = TRUE,scale_init = FALSE)

bartmod <- dbarts::bart(x.train = x_train,y.train = y_train,x.test = x_test)
softbartmod <- SoftBart::softbart(X = x_train,Y = y_train,X_test =  x_test)

motr_bart_mod <- motr_bart(x = x_train,y = y_train)
motrbart_pred <- predict_motr_bart(object = motr_bart_mod,newdata = x_test,type = "all")

# Summariseing the rmse
n_burn <- 500
rmse(x = colMeans(spBART$y_train_hat[(n_burn+1):spBART$mcmc$n_mcmc,]),y = y_train)
# plot(x = x_train$x,y = colMeans(spBART$y_train_hat))
# plot(x = x_test$x,y = colMeans(spBART$y_test_hat))

rmse(x = colMeans(spBART$y_test_hat[(n_burn+1):spBART$mcmc$n_mcmc,]),y = test$y)

# rmse(x = colMeans(motr_bart_mod$y_hat),y = train$y)
# rmse(x = motrbart_pred %>% colMeans(),y = test$y)

rmse(x = bartmod$yhat.train.mean,y = y_train)
rmse(x = bartmod$yhat.test.mean,y = test$y)
# plot(x = x_train$x,y = bartmod$yhat.train.mean)
# plot(x = x_test$x,y = bartmod$yhat.test.mean)

par(mfrow = c(3,1))
plot(spBART$all_tau, type = "l")
plot(spBART$all_tau_beta, type = "l")
plot(spBART$all_tau_gamma, type = "l")
par(mfrow = c(1,1))

par(mfrow = c(3,1))
plot(spBART$all_tau[(n_burn+1):spBART$mcmc$n_mcmc], type = "l")
plot(spBART$all_tau_beta[(n_burn+1):spBART$mcmc$n_mcmc], type = "l")
plot(spBART$all_tau_gamma[(n_burn+1):spBART$mcmc$n_mcmc], type = "l")
par(mfrow = c(1,1))

# Getting a counter for times that a variable is used in within a tree
p_counter <- numeric(ncol(x_train))

for(i in 501:spBART$mcmc$n_mcmc){

  for(t in 1:spBART$prior$n_tree){
      curr_tree <- spBART$mcmc$all_trees[[i]][[t]]
      terminals <- get_terminals(curr_tree)
      for(ell in 1:length(terminals)){
        p_counter[unique(curr_tree[[terminals[ell]]]$ancestors)] <- p_counter[unique(curr_tree[[terminals[ell]]]$ancestors)] + 1
      }
    }
}

round(p_counter/sum(p_counter),digits = 5)
rmse(x = colMeans(spBART$y_train_hat[(n_burn+1):spBART$mcmc$n_mcmc,]),y = y_train)
rmse(x = colMeans(spBART$y_test_hat[(n_burn+1):spBART$mcmc$n_mcmc,]),y = test$y)

rmse(x = bartmod$yhat.train.mean,y = y_train)
rmse(x = bartmod$yhat.test.mean,y = test$y)

rmse(x = softbartmod$y_hat_train_mean,y = y_train)
rmse(x = softbartmod$y_hat_test_mean,y = test$y)

rmse(x = colMeans(motr_bart_mod$y_hat),y = train$y)
rmse(x = motrbart_pred %>% colMeans(),y = test$y)
