# Simulating the data
library(mlbench)
rm(list=ls())
seed_ <- 42
set.seed(42)
n_ <- 250
sd_ <- 1

# ======================================================================
# Auxiliar Functions
# ======================================================================

# Normalize BART function (Same way ONLY THE COVARIATE NOW)
normalize_covariates_bart <- function(y, a = NULL, b = NULL) {

  # Defining the a and b
  if( is.null(a) & is.null(b)){
    a <- min(y)
    b <- max(y)
  }
  # This will normalize y between -0.5 and 0.5
  y  <- (y - a)/(b - a)
  return(y)
}


# Normalize BART function (Same way ONLY THE COVARIATE NOW)
normalize_bart <- function(y, a = NULL, b = NULL) {

  # Defining the a and b
  if( is.null(a) & is.null(b)){
    a <- min(y)
    b <- max(y)
  }
  # This will normalize y between -0.5 and 0.5
  y  <- (y - a)/(b - a) - 0.5
  return(y)
}

# A function to create the penalty matrix P
P_gen <- function(D_train_, dif_order_,eta){

  P_train_ <- crossprod(diff(diag(NCOL(D_train_)),differences = dif_order_))

  if(dif_order_==1){
    if(nrow(P_train_)%%2==0){
      middle_ <- trunc(nrow(P_train_)/2)+1
    } else {
      middle_ <- trunc(nrow(P_train_)/2)
    }
    # middle_ <- 1
    P_train_[middle_,middle_] = P_train_[middle_,middle_] + eta
  } else if(dif_order_==2) {
    P_train_[1,1] = P_train_[1,1] + eta
    if(nrow(P_train_)%%2==0){
      middle_ <- trunc(nrow(P_train_)/2)+1
    } else {
      middle_ <- trunc(nrow(P_train_)/2)
    }
    P_train_[middle_,middle_] = P_train_[middle_,middle_] + eta
  } else if (dif_order_==3) {
    P_train_[1,1] = P_train_[1,1] + eta
    P_train_[nrow(P_train_),ncol(P_train_)] = P_train_[nrow(P_train_),ncol(P_train_)] + eta
    if(nrow(P_train_)%%2==0){
      middle_ <- trunc(nrow(P_train_)/2)+1
    } else {
      middle_ <- trunc(nrow(P_train_)/2)
    }
    P_train_[middle_,middle_] = P_train_[middle_,middle_] + eta

  } else {
    stop("Insert a lower order for the difference matrix")

  }

  return(P_train_)
}

# Creating the D (difference matrix)
D_gen <- function(p, n_dif){
  return(diff(diag(p),diff = n_dif))
}

# ==========================================
# End
# ==========================================

# Getting simulations
sim_train <- mlbench::mlbench.friedman1(n = n_,sd = sd_)  |> as.data.frame()
sim_test <- mlbench::mlbench.friedman1(n = n_,sd = sd_)  |> as.data.frame()

# Getting the x.scale()
x_train <- sim_train[,1:10,drop = FALSE]
x_test <- sim_test[,1:10,drop = FALSE]
y_train <- sim_train$y

# Getting the train and test set
x_train_scale <- as.matrix(x_train)
x_test_scale <- as.matrix(x_test)

# Scaling x
x_min <- apply(as.matrix(x_train_scale),2,min)
x_max <- apply(as.matrix(x_train_scale),2,max)

# Storing the original
x_train_original <- x_train
x_test_original <- x_test


# Normalising all the columns
for(i in 1:ncol(x_train)){
  x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
  x_test_scale[,i] <- normalize_covariates_bart(y = x_test_scale[,i],a = x_min[i], b = x_max[i])
}

# Creating the numcuts matrix of splitting rules
numcut <- 100
xcut_m <- matrix(NA,nrow = numcut,ncol = ncol(x_train_scale))
for(i in 1:ncol(x_train_scale)){

  if(nrow(x_train_scale)<numcut){
    xcut_m[,i] <- sort(x_train_scale[,i])
  } else {
    xcut_m[,i] <- seq(min(x_train_scale[,i]),
                      max(x_train_scale[,i]),
                      length.out = numcut+2)[-c(1,numcut+2)]
  }
}


p_var_basis <- 3 # Getting the Basis for X.3 for example
p_var_split <- 7 # Getting which gonna be the split rule for X3

# Generating parameters for the spline basis matrix
nIknots <- 2
ndx <- nIknots+1
ord_ <- 4
degree_ <- 3
x_min_sp <- apply(x_train_scale,2,min)
x_max_sp <- apply(x_train_scale,2,max)
dx <- (x_max_sp-x_min_sp)/ndx

# New_knots
new_knots <- matrix()
new_knots <- matrix(mapply(x_min_sp,x_max_sp,dx, FUN = function(MIN,MAX,DX){seq(from = MIN-(ord_-1)*DX, to = MAX+(ord_-1)*DX, by = DX)}), ncol = NCOL(x_train_scale)) # MIN and MAX are 0 and 1 respectively, because of the scale
colnames(new_knots) <- paste0('x.',1:NCOL(x_train_scale))


B_train <- as.matrix(splines::spline.des(x = x_train_scale[,p_var_basis, drop = FALSE],
                                        knots = new_knots[,p_var_basis],
                                        ord = ord_,
                                        derivs = 0*x_train_scale[,p_var_basis, drop = FALSE],outer.ok = TRUE)$design)
# Quick glance over the basis functions
plot(NA,ylim = range(B_train),xlim = range(x_train_scale), ylab = '',main = 'Basis Functions')
# Seeing the basis
for(i in 1:ncol(B_train)){
  points(x_train_scale[,p_var_basis],B_train[,i], col = i,pch=20)
}

# Getting left and right node
cutpoint <- 50 # choosing to split in the middle for example can be any value from 0 to 100
cutpoint_value <- xcut_m[cutpoint,p_var_split]
left_index <- which(x_train_scale[,p_var_split] <= cutpoint_value)
right_index <- which(x_train_scale[,p_var_split] > cutpoint_value)
curr_index <- 1:nrow(x_train_scale)
# Scaling y and setting other parameters
y_scale <- normalize_bart(y = y_train,a = min(y_train),b = max(y_train))
n_tree <- 1; k <- 2;
tau_beta <- 4*(k^2)*n_tree
tau <- (sd_^(-2))*diff(range(y_train))^2 # Getting the true tau value in the scaled version

dif_order_ <- 1

# Creating a function to calculate the loglike-lihood for a specific subset
logLikelihood_mvn <- function(y_scale,
                              B_train,
                              curr_index,
                              tau_beta,
                              tau,
                              dif_order){

  # Subsetting the basis matrix and residuals
  B_subset <- B_train[curr_index,,drop = FALSE]
  y_subset <- y_scale[curr_index]

  #
  mean_ <- rep(0,length(curr_index))
  if(dif_order==0){
    P <- diag(tau_beta, nrow =  NCOL(B_subset))
  } else {
    D_ <- D_gen(p = NCOL(B_train),n_dif = dif_order)
    P <- tau_beta*P_gen(D_train_ = D_,dif_order_ = dif_order,eta = 1)
  }
  cov_ <- diag(tau^(-1), nrow = length(curr_index)) + B_subset%*%solve(P,t(B_subset))


  return ( mvnfast::dmvn(X = y_subset,mu = mean_,sigma = cov_,log = TRUE))

}


# Calculating the loglikelihood associated with each move
g_node_ <- logLikelihood_mvn(y_scale = y_scale,B_train = B_train,curr_index = curr_index,
                  tau_beta = tau_beta,tau = tau,dif_order = dif_order_)
left_child_ <- logLikelihood_mvn(y_scale = y_scale,B_train = B_train,curr_index = left_index,
                  tau_beta = tau_beta,tau = tau,dif_order = dif_order_)
right_child_ <- logLikelihood_mvn(y_scale = y_scale,B_train = B_train,curr_index = right_index,
                  tau_beta = tau_beta,tau = tau,dif_order = dif_order_)

acceptance_loglike <- exp(-g_node_+left_child_+right_child_)
acceptance_loglike

