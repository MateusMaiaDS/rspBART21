# This file is to compare first trial of predictive performance of spBART and BART
rm(list=ls())
library(dbarts)
library(mlbench)
library(purrr)
library(MOTRbart)
library(doParallel)
source("R/sim_functions.R")
source("R/main_function.R")
set.seed(42)

n_ <- 100
sd_ <- 1
n_rep_ <- 10
nIknots_ <- 20
ntree_ <- 5
dif_order_ <- 2
use_bs_ <- FALSE
seed_ <- 42
y_scale_ <- TRUE
motr_bart_ <- FALSE
all_ <- FALSE
alpha_ <- 0.5
stump_ <- FALSE
scale_init_ <- FALSE
update_tau_beta_ <- TRUE
inter_ <- TRUE
mle_prior_ <- FALSE
n_mcmc_ <- 5000
n_burn_ <- 2000

# Selecting a simulated scenarion
# (1): "oned_break" one dimensionnal sin(2*x) with a break
# (2): "friedman_nointer_nonoise": four-dimensional friedmna setting with no interaction terms and no extra X noise variables
# (3): "interaction
# type_ <- c("friedman_break")
type_ <- c("friedman")
# type_ <- "smooth.main.formula"
# type_ <- "non.smooth.main.formula"
# type_ <- "non.and.smooth.main.formula"
# type_ <- 'mlbench.d1.break'
# ================
# Printing message
# ================


cv_ <- vector("list", n_rep_)

# Generating CV_ object
for( i in 1:n_rep_){


  # train <- mlbench.d1.break(n = n_,sd = sd_) %>% as.data.frame()
  # test <- mlbench.d1.break(n = n_,sd = sd_) %>% as.data.frame()

  if(type_ == "smooth.main.formula"){
    train_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = smooth.main.formula)
    test_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = smooth.main.formula)
    train <- data.frame(x = train_sim$x, y = train_sim$y)
    test <- data.frame(x = test_sim$x, y = test_sim$y)
  }

  if(type_ == "non.smooth.main.formula"){
    train_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = non.smooth.main.formula)
    test_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = non.smooth.main.formula)
    train <- data.frame(x = train_sim$x, y = train_sim$y)
    test <- data.frame(x = test_sim$x, y = test_sim$y)
  }

  if(type_ == "non.and.smooth.main.formula"){
    train_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = non.and.smooth.main.formula)
    test_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = non.and.smooth.main.formula)
    train <- data.frame(x = train_sim$x, y = train_sim$y)
    test <- data.frame(x = test_sim$x, y = test_sim$y)
  }

  if(type_ == "friedman_nointer_nonoise"){
    train <- mlbench.friedman1.nointeraction(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1.nointeraction(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == "friedman_break"){
    train <- break.mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
    test <- break.mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == "friedman_nointer_noise"){
    train <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_) %>% as.data.frame()
  }
  # train <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame() %>% .[,c(1:5,11)]
  # test <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame() %>% .[,c(1:5,11)]

  if(type_ == "friedman"){
    train <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == "friedman_interaction"){
    train <- mlbench.friedman1.interaction.only(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1.interaction.only(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == 'mlbench.d1.break'){
    train <- mlbench.d1.break(n = n_,sd = sd_)  |> as.data.frame()
    test <- mlbench.d1.break(n = n_,sd = sd_) |> as.data.frame()
  }
  cv_[[i]]$train <- train
  cv_[[i]]$test <- test
}



# Testing the simple n_tree
result <- foreach(i = 1:n_rep_, .packages = c("dbarts","SoftBart","MOTRbart","dplyr")) %dopar%{

  devtools::load_all()
  source("/users/research/mmarques/spline_bart_lab/rspBART21/R/sim_functions.R")
  source("/users/research/mmarques/spline_bart_lab/rspBART21/R/main_function.R")
  source("/users/research/mmarques/spline_bart_lab/rspBART21/R/cv_functions.R")
  # if(ntree_<50) {
  #   aux <- all_bart(cv_element = cv_[[i]],
  #                   nIknots_ = nIknots_,ntree_ = ntree_,seed_ = seed_,
  #                   use_bs_ = use_bs_,motr_bart_ = motr_bart_,rsp_bart_all_ = all_,
  #                   alpha_ = alpha_,stump = stump_)
  # } else {
    aux <- all_bart_lite_interaction(cv_element = cv_[[i]],
                         nIknots_ = nIknots_,ntree_ = ntree_,seed_ = seed_,
                         use_bs_ = use_bs_,alpha_ = alpha_,rsp_bart_all_ = all_,
                         j = i,motr_bart_ = motr_bart_, stump = stump_,dif_order_ = dif_order_,
                         scale_init = scale_init_,update_tau_beta_ = update_tau_beta_,
                         interaction_term_ = inter_,mle_prior_ = mle_prior_,y_scale_ = y_scale_,
                         n_mcmc_ = n_mcmc_,n_burn_ = n_burn_)
  # }

  aux
}


#
stopCluster(cl)



# Saving all
# if(all_){
# saveRDS(object = result,file = paste0("/localusers/research/mmarques/spline_bart_lab/preliminar_results/rspBART21/oned_n_",n_,
#                "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_bs_",use_bs_,"_motr_bart_",motr_bart_,"_allvar_",all_,".Rds"))
# } else {
#   saveRDS(object = result,file = paste0("/localusers/research/mmarques/spline_bart_lab/preliminar_results/rspBART21/oned_n_",n_,
#                                         "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_bs_",use_bs_,"_motr_bart_",motr_bart_,"_alpha_",alpha_,".Rds"))
# }



if(type_ == "friedman"){
      saveRDS(object = result,file = paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART21/friedman/v26_new_interaction_",inter_,"_oned_n_",n_,
                                            "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_allvar_",all_,
                                            "_alpha_",alpha_,"_uptaubeta_",update_tau_beta_,"_dif_",dif_order_,
                                            "_nmcmc_",n_mcmc_,"_nburn_",n_burn_,".Rds"))
}

if(type_ == "friedman_break"){
  saveRDS(object = result,file = paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART21/friedman_break/v26_new_interaction_",inter_,"_oned_n_",n_,
                                        "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_allvar_",all_,
                                        "_alpha_",alpha_,"_uptaubeta_",update_tau_beta_,"_dif_",dif_order_,
                                        "_nmcmc_",n_mcmc_,"_nburn_",n_burn_,".Rds"))
}


