library(mlbench)
library(tidyverse)
rm(list=ls())
devtools::load_all()
seed_ <- 42
set.seed(42)
n_ <- 250
sd_ <- 1
# sim_train <- mlbench.friedman1.nointeraction(n = n_,sd = sd_)  |> as.data.frame()
# sim_test <- mlbench.friedman1.nointeraction(n = n_,sd = sd_)  |> as.data.frame()

# sim_train <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_)  |> as.data.frame()
# sim_test <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_)  |> as.data.frame()

# sim_train <- mlbench.friedman1(n = n_,sd = sd_)  |> as.data.frame()
# sim_test <- mlbench.friedman1(n = n_,sd = sd_)  |> as.data.frame()

# sim_train <- smooth.main(n = n_,sd = sd_,p_ = 10) %>%  as.data.frame()
# sim_test <- smooth.main(n = n_,sd = sd_,p_ = 10)  %>% as.data.frame()

sim_train <- break.mlbench.friedman1(n = n_,sd = sd_)  |> as.data.frame()
sim_test <- break.mlbench.friedman1(n = n_,sd = sd_)  |> as.data.frame()
#

# sim_train <- mlbench.d1.break(n = n_,sd = 1)  |> as.data.frame()
# sim_test <- mlbench.d1.break(n = n_,sd = 1) |> as.data.frame()

# sim_train <- mlbench.d1(n = n_,sd = 1)  |> as.data.frame()
# sim_test <- mlbench.d1(n = n_,sd = 1) |> as.data.frame()


# Testing the model with only interactions
# sim_train <- mlbench.friedman1.interaction.only(n = n_,sd = sd_) %>% as.data.frame()
# sim_test <- mlbench.friedman1.interaction.only(n = n_,sd = sd_) %>% as.data.frame()


x_train <- sim_train |> dplyr::select(dplyr::starts_with("x"))
x_test <-  sim_test|> dplyr::select(dplyr::starts_with("x"))
y_train <- sim_train$y


#== Baltimore ==#
# baltimore<- baltimore %>% dplyr::rename(lon=X,lat=Y)
# baltimore<- baltimore[,-1,drop=FALSE] %>% dplyr::rename(y=PRICE)
# baltimore <- baltimore %>% select(lon,lat,SQFT,LOTSZ,y)
# x_test <- x_train <- baltimore[,1:4,drop = FALSE]
# y_train <- baltimore$y
# sim_train <- sim_test <- baltimore

# x_train <- x_train[,1:5]
# x_test <- x_test[,1:5]
n_tree <- 5
node_min_size = 25
n_mcmc = 5000
n_burn = 2500
alpha = 0.5
beta = 2
df = 3
sigquant = 0.9
kappa = 2
tau = 1
scale_bool = TRUE
stump = FALSE
numcut = 100L # Defining the grid of split rules
usequants = TRUE
linero_sampler <- FALSE
# Splines parameters
nIknots = 2
dif_order = 1
motrbart_bool <- FALSE
use_bs <- FALSE
plot_preview = FALSE
intercept <- FALSE
all_var <- FALSE
scale_init <- FALSE
update_tau_beta <- TRUE
main_effects_pred <- TRUE
# interaction_list <- interaction_list <- list(c(1,2))
interaction_list <- interaction_list <- NULL

# interaction_list <- NULL
store_tree_fit <- FALSE
interaction_term <- TRUE
cv_object_ <- kfold(data_ = sim_train,nfold_ = 10,seed_ = 42)
fold_ <- 1
cv_object_fold_ <- cv_object_[[fold_]]
varimportance_bool <- TRUE
use_D_bool <- TRUE
mle_prior <- FALSE
aux <- dbarts::bart(x.train = x_train,y.train = y_train,x.test = x_test)
library(SoftBart)
aux2 <- softbart(X = x_train,Y = y_train,X_test = x_test)
# plot(aux$sigma, type = 'l')
# plot(all_tau_norm^(-1/2), type = 'l'  , ylim = c(1,4), ylab = expression(sigma^2))
# lines(aux$sigma, type = 'l' , col = 'blue')
# softbart_mod <- SoftBart::softbart(X = x_train,Y = y_train,X_test = x_test)
# lines(softbart_mod$sigma, type = 'l' , col ="orange")
# =======================
# Doing some extra plots
# =======================
# par(mfrow = c(2,2))
# plot(x_test, cex = sim_test$y/5, pch = 20, main = paste("Test sample: 10*sin(pi*x1*x2)"))
# plot(x_test, cex = colMeans(all_y_hat_test_norm)/5, pch = 20,main = paste("Test \\hat sample: 10*sin(pi*x1*x2)"))
# rmse(sim_test$y, colMeans(all_y_hat_test_norm))

# #==========================
# # Plotting the main effects
# # =========================
# # Marginal plots
# selected_var_ <- 1
# plot(x_train[,selected_var_],colMeans(main_effects_train_list[[selected_var_]]),pch = 20, ylab = expression(f(x[1])),
# xlab = expression(x[1]), main = paste0("Main effect for: ",expression(x[1])))
#
# selected_var_ <- 2
# plot(x_train[,selected_var_],colMeans(main_effects_train_list[[selected_var_]]),pch = 20,
#      ylab = expression(f(x[2])), xlab = expression(x[2]),
#      main = paste0("Main effect for: ",expression(x[2])))
#
#
# # Comparing with BART and softBART
# softbart_mod <- SoftBart::softbart(X = x_train,Y = y_train,X_test = x_test)
# plot(x_test, cex = softbart_mod$y_hat_test_mean/10, pch = 20,main = paste("SoftBART \\hat sample: 10*sin(pi*x1*x2)"))
# rmse(sim_test$y, softbart_mod$y_hat_test_mean)
#
#
# bartmod <- dbarts::bart(x.train = x_train,y.train = y_train,x.test = x_test)
# rmse(sim_test$y,bartmod$yhat.test.mean)
#
#
# install.packages("BASS")
# bass_mod <- BASS::bass(xx = x_train,y = y_train,degree = 2)
# summary(bass_mod)
# predict_bass_mod <- predict(bass_mod,x_test)
# rmse(predict_bass_mod %>% colMeans(),sim_test$y)
