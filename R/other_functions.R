# Function to get the tensors
multiply_matrices_general <- function(A, B) {
  # Get the number of rows and columns for A and B
  nrow_A <- nrow(A)
  ncol_A <- ncol(A)
  nrow_B <- nrow(B)
  ncol_B <- ncol(B)
  # Initialize an empty matrix C with dimensions nrow_A x (ncol_A * ncol_B)
  C <- A[,1]*B
  # Check if both matrices have the same number of rows
  if (nrow_A != nrow_B) {
    stop("Matrices A and B must have the same number of rows")
  }
  # Loop to fill in the values of matrix C
  for (i in 2:ncol_A) {
    C <-  cbind(C,(A[, i]*B))
  }
  return(C)
}

# A function to create the penalty matrix P
P_gen_keefe <- function(D_train_, dif_order_,eta){

  if(dif_order_>0){
    P_train_ <- crossprod(diff(diag(NCOL(D_train_)),differences = dif_order_))
  } else {
    P_train_ <- diag(NCOL(D_train))
  }

  if(nrow(P_train_)%%2==0){
    middle_ <- trunc(nrow(P_train_)/2)+1
  } else {
    middle_ <- trunc(nrow(P_train_)/2)
  }
  # middle_ <- 1
  P_train_[middle_,middle_] = P_train_[middle_,middle_] + eta
  return(P_train_)

}

# A function to create the penalty matrix P
P_gen <- function(D_train_, dif_order_,eta){
  warning("Look at the actual values of where the eta is added.")
  if(dif_order_>0){
      P_train_ <- crossprod(diff(diag(NCOL(D_train_)),differences = dif_order_))
  } else {
      P_train_ <- diag(NCOL(D_train_))
  }
  if(dif_order_==1){
    if(nrow(P_train_)%%2==0){
      middle_ <- trunc(nrow(P_train_)/2)+1
    } else {
      middle_ <- trunc(nrow(P_train_)/2)
    }
    middle_ <- 1
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
    middle_ <- nrow(P_train)
    P_train_[middle_,middle_] = P_train_[middle_,middle_] + eta

  } else if(dif_order_==0){
    P_train_ <- diag(nrow = NCOL(D_train_))
  } else {
    stop("Insert a lower order for the difference matrix")

  }

  return(P_train_)
}

# Creating the D (difference matrix)
D_gen <- function(p, n_dif){
  if(n_dif>0){
    return(diff(diag(p),diff = n_dif))
  } else {
    return(diag(p))
  }
}

# In case where \mathbf{u} ~ MVN(Q^-1 %*% b, Q^-1)
keefe_mvn_sampler <-  function(b, Q) {
  p    <- NCOL(Q)
  Z    <- rnorm(p)
  if(p == 1) {
    U  <- sqrt(Q)
    drop((b/U + Z)/U)
  } else     {
    U  <- chol(Q)
    backsolve(U, backsolve(U, b, transpose=TRUE, k=p) + Z, k=p)
  }
}

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

# Getting back to the original scale
unnormalize_bart <- function(z, a, b) {
  # Just getting back to the regular BART
  y <- (b - a) * (z + 0.5) + a
  return(y)
}

unnormalize_bart_me <- function(z, a, b) {
  # Just getting back to the regular BART
  y <- (b - a) * (z)
  return(y)
}


# Half-cauchy log-density
log_hcauchy <- function(x, mu = 0, sigma) {
  if(x >= mu) {
    return(log(2) - log(pi*sigma) - log(1+(x-mu)/sigma^2))
  } else {
    return(-Inf)
  }
}

# Naive sigma_estimation
naive_sigma <- function(x,y){

  # Getting the valus from n and p
  n <- length(y)

  # Getting the value from p
  p <- ifelse(is.null(ncol(x)), 1, ncol(x))

  # Adjusting the df
  df <- data.frame(x,y)
  colnames(df)<- c(colnames(x),"y")

  # Naive lm_mod
  lm_mod <- stats::lm(formula = y ~ ., data =  df)

  # Getting sigma
  sigma <- summary(lm_mod)$sigma
  return(sigma)

}





# Function to create a vector of variables that being categorical will
#have the same code
recode_vars <- function(x_train, dummy_obj){

  vars <- numeric()
  j <- 0
  i <- 0
  c <- 1
  while(!is.na(colnames(x_train)[c])){
    if(colnames(x_train)[c] %in% dummy_obj$facVars){
      curr_levels <- dummy_obj$lvls[[colnames(x_train)[c]]]
      for(k in 1:length(curr_levels)){
        i = i+1
        vars[i] <- j
      }
    } else {

      i = i+1
      vars[i] <- j
    }
    j = j+1
    c = c+1
  }

  return(vars)
}

# Calculating the rmse
rmse <- function(x,y){
  return(sqrt(mean((y-x)^2)))
}


# Calculating CRPS from (https://arxiv.org/pdf/1709.04743.pdf)
crps <- function(y,means,sds){

  # scaling the observed y
  z <- (y-means)/sds

  crps_vector <- sds*(z*(2*stats::pnorm(q = z,mean = 0,sd = 1)-1) + 2*stats::dnorm(x = z,mean = 0,sd = 1) - 1/(sqrt(pi)) )

  return(list(CRPS = mean(crps_vector), crps = crps_vector))
}


# Calculating the rmse
rmse <- function(x,y){
  return(sqrt(mean((y-x)^2)))
}

# Calculating the rmse
mae <- function(x,y){
  return((mean(abs(y-x))))
}

pi_coverage <- function(y, y_hat_post, sd_post,only_post = FALSE, prob = 0.5,n_mcmc_replications = 1000){

  # Getting the number of posterior samples and columns, respect.
  np <- nrow(y_hat_post)
  nobs <- ncol(y_hat_post)

  full_post_draw <- list()

  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_mcmc_replications,
    style = 3, width = 50 )

  # Only post matrix
  if(only_post){
    post_draw <- y_hat_post
  } else {
    for(i in 1:n_mcmc_replications){
      utils::setTxtProgressBar(progress_bar, i)

      full_post_draw[[i]] <-(y_hat_post + replicate(sd_post,n = nobs)*matrix(stats::rnorm(n = np*nobs),
                                                                             nrow = np))
    }
  }

  if(!only_post){
    post_draw<- do.call(rbind,full_post_draw)
  }

  # CI boundaries
  low_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = prob/2)})
  up_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = 1-prob/2)})

  pi_cov <- sum((y<=up_ci) & (y>=low_ci))/length(y)

  return(pi_cov)
}

# A fucction to retrive the number which are the factor columns
base_dummyVars <- function(df) {
  num_cols <- sapply(df, is.numeric)
  factor_cols <- sapply(df, is.factor)

  return(list(continuousVars = names(df)[num_cols], facVars = names(df)[factor_cols]))
}
