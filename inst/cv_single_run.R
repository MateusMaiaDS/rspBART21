# This file is to compare first trial of predictive performance of spBART and BART
rm(list=ls())
library(dbarts)
library(mlbench)
library(purrr)
library(MOTRbart)
library(doParallel)
source("R/sim_functions.R")
source("R/main_function.R")
devtools::load_all()
set.seed(42)

# Simulation arguments
n_ <- 250
sd_ <- 1
n_rep_ <- 10


# Selecting a simulated scenarion
# (1): "oned_break" one dimensionnal sin(2*x) with a break
# (2): "friedman_nointer_nonoise": four-dimensional friedmna setting with no interaction terms and no extra X noise variables
# (3): "interaction
# type_ <- c("friedman_break")
# type_ <- c("friedman")
# type_ <- "smooth.main.formula"
type_ <- "non.smooth.main.formula"
# type_ <- "non.and.smooth.main.formula"

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
  # train <- mlbench.d1.break(n = n_,sd = sd_)  |> as.data.frame()
  # test <- mlbench.d1.break(n = n_,sd = sd_) |> as.data.frame()

  cv_[[i]]$train <- train
  cv_[[i]]$test <- test
}

# ======
# Setting the single model to be generated here.
# ======
selected_rep_ <- 1
selected_train <- cv_[[selected_rep_]]$train
selected_test <- cv_[[selected_rep_]]$test
sim_train <- selected_test
sim_test <- selected_test
# Defining the default parameters that are going to be used in the model, in this case I putting all
#of them as the same from function arguments to make it easier to save it.
x_train <- selected_train %>% dplyr::select(dplyr::starts_with("x"))
x_test <- selected_test %>% dplyr::select(dplyr::starts_with("x"))
y_train <- selected_train$y
y_test <- selected_test$y
n_tree <- 10
node_min_size <- 25
n_mcmc <- 5000
n_burn <- 2500
alpha <- 0.5
beta <- 2
df <- 1
sigquant <- 0.9
kappa <- 2
nIknots <- 20
dif_order <- 3
tau <- 1
scale_bool <- TRUE
stump <- FALSE
numcut <- 100
usequants <- FALSE
motrbart_bool <- FALSE
use_bs <- FALSE
plot_preview <- FALSE
all_var <- FALSE
scale_init <- FALSE
update_tau_beta <- TRUE
main_effects_pred <- TRUE
interaction_term <- TRUE
interaction_list <- NULL
store_tree_fit <- FALSE
mle_prior <- FALSE
linero_sampler <- FALSE
use_D_bool <- FALSE
varimportance_bool <- TRUE
seed_ <- 42
scale_basis_function <- FALSE
robust_prior <- FALSE
eta <- 1e-6
a_delta <- 0.5
d_delta <- 0.5
pen_basis <- TRUE

set.seed(seed_)

print(paste0("N: ",n_," SD: ", sd_, " nIknots: ", nIknots,
             " Ntree: ",n_tree, " Seed: ",seed_, " Alpha:", alpha,
             "Update \tau_\beta: ", update_tau_beta, " Dif.Order:", dif_order, "_type_", type_, "_nmcmc_", n_mcmc, "_nburn_",n_burn))

rsp_mod <- rspBART(x_train = x_train,
                   y_train = y_train,
                   x_test = x_test,
                   n_tree = n_tree,
                   node_min_size = node_min_size,
                   n_mcmc = n_mcmc,n_burn = n_burn,
                   alpha = alpha,beta = beta,
                   df = df,sigquant = sigquant,
                   kappa = kappa,nIknots = nIknots,
                   dif_order = dif_order,tau = tau,
                   scale_bool = scale_bool,stump = stump,
                   numcut = numcut,usequants = usequants,
                   motrbart_bool = motrbart_bool,use_bs = use_bs,
                   plot_preview = plot_preview,
                   all_var = all_var,scale_init = scale_init,
                   update_tau_beta = update_tau_beta,
                   main_effects_pred = main_effects_pred,
                   interaction_term = interaction_term,
                   interaction_list = interaction_list,
                   store_tree_fit = store_tree_fit,
                   mle_prior = mle_prior,
                   linero_sampler = linero_sampler,
                   use_D_bool = use_D_bool,
                   varimportance_bool = varimportance_bool,
                   scale_basis_function = scale_basis_function,
                   robust_prior = robust_prior,eta = eta,pen_basis = pen_basis)





saveRDS(object = rsp_mod,file = paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART20/",type_,"/single_run/v26_single_run_rep_",
                                       selected_rep_,"_n_",n_,
                                      "_sd_",sd_,"_nIknots_",nIknots,"_ntree_",n_tree,"_nodesize_",node_min_size,
                                      "_dif_",dif_order,"_scale_",scale_bool,"_sc_basis_",scale_basis_function,
                                      "_nmcmc_",n_mcmc,"_nburn_",n_burn,"_rb_prior_",robust_prior,"_bpen_",pen_basis,".Rds"))
