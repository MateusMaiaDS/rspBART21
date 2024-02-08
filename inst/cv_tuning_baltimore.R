# This file is to compare first trial of predictive performance of spBART and BART
rm(list=ls())
library(dbarts)
library(mlbench)
library(purrr)
library(MOTRbart)
library(doParallel)
library(spData)
source("R/sim_functions.R")
source("R/main_function.R")
source("R/cv_functions.R")
set.seed(42)

n_ <- 250
sd_ <- 1
n_rep_ <- 10
nIknots_ <- 2
ntree_ <- 100
dif_order_ <- 1
use_bs_ <- FALSE
seed_ <- 42
motr_bart_ <- FALSE
all_ <- TRUE
alpha_ <- 0.5
stump_ <- FALSE
scale_init_ <- FALSE
update_tau_beta_ <- TRUE
inter_ <- TRUE
mle_prior_ <- TRUE


#== Baltimore ==#
baltimore<- baltimore %>% dplyr::rename(lon=X,lat=Y)
baltimore<- baltimore[,-1,drop=FALSE] %>% dplyr::rename(y=PRICE)
baltimore <- baltimore %>% select(lon,lat,SQFT,LOTSZ,y)


# ================
# Printing message
# ================

print(paste0("N: ",n_," SD: ", sd_, " nIknots: ", nIknots_,
             " Ntree: ",ntree_, " Seed: ",seed_, " Alpha:", alpha_,
             "Update \tau_\beta: ", update_tau_beta_))


cv_ <- vector("list", n_rep_)

cv_ <- kfold(data_ = baltimore,nfold_ = 10,seed_ = 42)


# Setting up the parallel simulation
number_cores <- n_rep_
cl <- parallel::makeCluster(number_cores)
doParallel::registerDoParallel(cl)


# Testing the simple n_tree
result <- foreach(i = 1:n_rep_, .packages = c("dbarts","SoftBart","MOTRbart","dplyr")) %dopar%{

  devtools::load_all()
  source("/users/research/mmarques/spline_bart_lab/rspBART13/R/sim_functions.R")
  source("/users/research/mmarques/spline_bart_lab/rspBART13/R/main_function.R")
  source("/users/research/mmarques/spline_bart_lab/rspBART13/R/cv_functions.R")
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
                                   interaction_term_ = inter_,mle_prior_ = mle_prior_)
  # }

  aux
}


#
stopCluster(cl)


