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

n_ <- 250
sd_ <- 5
n_rep_ <- 10
nIknots_ <- 5
ntree_ <- 10
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

# Selecting a simulated scenarion
# (1): "oned_break" one dimensionnal sin(2*x) with a break
# (2): "friedman_nointer_nonoise": four-dimensional friedmna setting with no interaction terms and no extra X noise variables
# (3): "interaction
type_ <- c("friedman_inter_noise")

# type_ <- c("friedman_inter_noise")

# ================
# Printing message
# ================

print(paste0("N: ",n_," SD: ", sd_, " nIknots: ", nIknots_,
             " Ntree: ",ntree_, " Seed: ",seed_, " Alpha:", alpha_,
             "Update \tau_\beta: ", update_tau_beta_))


cv_ <- vector("list", n_rep_)

# Generating CV_ object
for( i in 1:n_rep_){


  # train <- mlbench.d1.break(n = n_,sd = sd_) %>% as.data.frame()
  # test <- mlbench.d1.break(n = n_,sd = sd_) %>% as.data.frame()

  if(type_ == "friedman_nointer_nonoise"){
    train <- mlbench.friedman1.nointeraction(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1.nointeraction(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == "friedman_nointer_noise"){
    train <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_) %>% as.data.frame()
  }
  # train <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame() %>% .[,c(1:5,11)]
  # test <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame() %>% .[,c(1:5,11)]

  if(type_ == "friedman_inter_noise"){
    train <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == "friedman_interaction"){
    train <- mlbench.friedman1.interaction.only(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1.interaction.only(n = n_,sd = sd_) %>% as.data.frame()
  }
  # train <- mlbench.d1.break(n = n_,sd = sd_)  |> as.data.frame()
  # test <- mlbench.d1.break(n = n_,sd = sd_) |> as.data.frame()

  cv_[[i]]$train <- train
  cv_[[i]]$test <- test
}


# Setting up the parallel simulation
number_cores <- n_rep_
cl <- parallel::makeCluster(number_cores)
doParallel::registerDoParallel(cl)


# Testing the simple n_tree
result <- foreach(i = 1:n_rep_, .packages = c("dbarts","SoftBart","MOTRbart","dplyr")) %dopar%{

  devtools::load_all()
  source("/localusers/researchers/mmarques/spline_bart_lab/rspBART13/R/sim_functions.R")
  source("/localusers/researchers/mmarques/spline_bart_lab/rspBART13/R/main_function.R")
  source("/localusers/researchers/mmarques/spline_bart_lab/rspBART13/R/cv_functions.R")
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
                                   interaction_term_ = inter_)
  # }

  aux
}


#
stopCluster(cl)



# Saving all
# if(all_){
# saveRDS(object = result,file = paste0("/local/localusers/researchers/mmarques/spline_bart_lab/preliminar_results/rspBART13/oned_n_",n_,
#                "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_bs_",use_bs_,"_motr_bart_",motr_bart_,"_allvar_",all_,".Rds"))
# } else {
#   saveRDS(object = result,file = paste0("/local/localusers/researchers/mmarques/spline_bart_lab/preliminar_results/rspBART13/oned_n_",n_,
#                                         "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_bs_",use_bs_,"_motr_bart_",motr_bart_,"_alpha_",alpha_,".Rds"))
# }


if(type_ == "friedman_nointer_nonoise"){
  if(all_){
    if(!stump_){
      saveRDS(object = result,file = paste0("/localusers/researchers/mmarques/spline_bart_lab/preliminar_results/rspBART13/friedman_noint_nonoise/oned_n_",n_,
                                            "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_bs_",use_bs_,"_motr_bart_",motr_bart_,"_allvar_",all_,"_stump_",stump_,".Rds"))
    } else {
      saveRDS(object = result,file = paste0("/localusers/researchers/mmarques/spline_bart_lab/preliminar_results/rspBART13/friedman_noint_nonoise/oned_n_",n_,
                                            "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_bs_",use_bs_,
                                            "_motr_bart_",motr_bart_,"_allvar_",all_,"_stump_",stump_,
                                            "_sinit_",scale_init_,"_alpha_",alpha_,"_uptaubeta_",update_tau_beta_,"_dif_",dif_order_,".Rds"))
    }
  } else {
    saveRDS(object = result,file = paste0("/localusers/researchers/mmarques/spline_bart_lab/preliminar_results/rspBART13/friedman_noint_nonoise/oned_n_",n_,
                                          "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_bs_",use_bs_,"_motr_bart_",motr_bart_,"_alpha_",alpha_,"_stump_",stump_,"_dif_",dif_order_,".Rds"))
  }

}

if(type_ == "friedman_inter_noise"){
  saveRDS(object = result,file = paste0("/localusers/researchers/mmarques/spline_bart_lab/preliminar_results/rspBART13/friedman/v4_new_interaction_",inter_,"_oned_n_",n_,
                                        "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_bs_",use_bs_,
                                        "_motr_bart_",motr_bart_,"_allvar_",all_,"_stump_",stump_,
                                        "_sinit_",scale_init_,"_alpha_",alpha_,"_uptaubeta_",update_tau_beta_,"_dif_",dif_order_,".Rds"))
}

if(type_ == "friedman_nointer_noise"){
  saveRDS(object = result,file = paste0("/localusers/researchers/mmarques/spline_bart_lab/preliminar_results/rspBART13/friedman_noint_noise/v3_new_interaction_n_",n_,
                                        "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,"_bs_",use_bs_,
                                        "_motr_bart_",motr_bart_,"_allvar_",all_,"_stump_",stump_,
                                        "_sinit_",scale_init_,"_alpha_",alpha_,"_uptaubeta_",update_tau_beta_,"_dif_",dif_order_,".Rds"))
}


