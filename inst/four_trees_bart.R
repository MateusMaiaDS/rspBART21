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
sd_ <- 1
n_rep_ <- 10
nIknots_ <- 10
ntree_ <- 2
use_bs_ <- FALSE
seed_ <- 42
motr_bart_ <- FALSE
all_ <- TRUE
alpha_ <- 0.95
stump_ <- TRUE


# Selecting a simulated scenarion
# (1): "oned_break" one dimensionnal sin(2*x) with a break
# (2): "friedman_noiter_nonoise": four-dimensional friedmna setting with no interaction terms and no extra X noise variables

type_ <- c("friedman_nointer_nonoise")

# ================
# Printing message
# ================

# Generating CV_ object
train <- mlbench.friedman1.nointeraction(n = n_,sd = sd_) %>% as.data.frame()
test <- mlbench.friedman1.nointeraction(n = n_,sd = sd_) %>% as.data.frame()

x_train <- train[,1:4]
x_test <- test[,1:4]
y_train <- test$y
# load_all()
four_tree <- rspBART(x_train = x_train,y_train = y_train,
                     x_test = x_test,n_tree = 4,all_var = FALSE,stump = FALSE)

source("R/cv_functions.R")
test <- tree_length_counter(result_ = four_tree,rep_ = 1,only_sp = TRUE)
test2 <- var_importance_counter(result_ = four_tree,rep_ = 1,only_sp_ = TRUE)

plot(four_tree$y_train_hat %>% colMeans(),y_train)
plot(four_tree$y_test_hat %>% colMeans(),test)

t <- 4
tree_one <- c()
df <- matrix(nrow = length(four_tree$mcmc$all_trees), ncol = 10) %>% as.data.frame()

# Iterating over all trees
for(i in 1000:length(four_tree$mcmc$all_trees)){
  
    t_nodes <- get_terminals(four_tree$mcmc$all_trees[[i]][[t]])
    for(ell in 1:length(t_nodes)){
      numbers <- gsub("[^0-9.]", "", t_nodes[ell])
      df[i, as.numeric(numbers)+1] <- paste(as.character(four_tree$mcmc$all_trees[[i]][[t]][[t_nodes[ell]]]$ancestors),collapse =  "-")
    }
    


}

cat(paste(" ------- Tree",t," --------:\n"))
for(i in 1:20){
node <- i
cat(paste("Node",i,":\n"))
print(df[,node][!is.na(df[,node])] %>% table() %>% sort(decreasing = TRUE))
}

