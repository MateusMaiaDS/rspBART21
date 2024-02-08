rm(list=ls())
set.seed(42)
library(purrr)
# Generating a simple simulation scenario
beta1 <- 10
beta2 <- 2
n <- 1000
nmcmc <- 1000
p <- 3
x <- matrix(runif(p*n),ncol = p)
y <- beta1*x[,1] + beta2*x[,2] + rnorm(n = n)
# y <- beta1*x[,1] + rnorm(n = n)
# y <- y/10


# Creating a sampler

tau_beta_j_post <- matrix(0,nrow = nmcmc, ncol = NCOL(x))
tau <- 1
acceptance_rate <- 0


# Creating the basis splines
B1 <- bs(x[,1]) %>% as.matrix()
B2 <- bs(x[,2]) %>% as.matrix()
B3 <- bs(x[,3]) %>% as.matrix()
list_B <- list(B1 = B1, B2 = B2, B3 = B3)

betas <- vector("list", 3)
for(ii in 1:3){
  betas[[ii]] <- matrix(0,ncol = NCOL(B1))
}

seletected_vars <- 2
# a_tau_beta_j_prior <- rep(0.001,NCOL(x))
# d_tau_beta_j_prior <- rep(0.001,NCOL(x))
selected_vars_list <- list()
a_tau_beta_j_prior <- c(0.1,0.1,0.1)
d_tau_beta_j_prior <- c(0.1,0.1,0.1)

selected_vars <- 2
tau_beta_j <- rep(1,NCOL(x))

# Initialising Keefe Sampler;
for(mcmc_iter in 1:nmcmc){
  verb <- proposal <- sample(c("grow","change","prune"),size = 1)
  # verb <- "grow"

  sigma_likelihood <- matrix(0,nrow = n, ncol = n)
  sigma_likelihood <- sigma_likelihood + (tau^-1)*diag(nrow = n)

  for(i in 1:length(selected_vars)){
    sigma_likelihood <- sigma_likelihood + (tau_beta_j[selected_vars[i]]^-1)*tcrossprod(list_B[[selected_vars[i]]])
  }

  old <- mvnfast::dmvn(X = y,mu = rep(0,n),sigma =  sigma_likelihood ,log = TRUE)

  if(verb == "grow"){
    if(length(which(!(1:NCOL(x) %in% selected_vars)))!=0){
      new_selected_vars <- sort(unique(c(selected_vars,sample(x = which(!(1:NCOL(x) %in% selected_vars)),size = 1))))
    } else {
      new_selected_vars <- selected_vars
    }

  } else if (verb == "prune"){

    # Avoiding if there's only one variable
    if(length(selected_vars)>1){
      p_var <- sample(selected_vars,size = 1)
      new_selected_vars <- sort(unique(selected_vars[!(selected_vars %in% p_var)]))
    } else {
      new_selected_vars <- selected_vars
    }

  } else {
    if(length((1:NCOL(x))[!((1:NCOL(x)) %in% selected_vars)])!=0){
      c_var <- sample((1:NCOL(x))[!((1:NCOL(x)) %in% selected_vars)],size = 1)
      c_index <- sample(1:length(selected_vars),size = 1)
      new_selected_vars <- selected_vars
      new_selected_vars[c_index] <- c_var
      new_selected_vars <- sort(unique(new_selected_vars))
    } else {
      new_selected_vars <- selected_vars
    }
  }

  new_sigma_likelihood <- matrix(0,nrow = n, ncol = n)
  new_sigma_likelihood <- new_sigma_likelihood + (tau^-1)*diag(nrow = n)

  for(i in 1:length(new_selected_vars)){
    new_sigma_likelihood <- new_sigma_likelihood + (tau_beta_j[new_selected_vars[i]]^-1)*tcrossprod(list_B[[new_selected_vars[i]]])
  }

  new <- mvnfast::dmvn(X = y,mu = rep(0,n),sigma =  new_sigma_likelihood ,log = TRUE)
  new

  acceptance <- exp(new-old)

  if(runif(n = 1)<acceptance){
    selected_vars <- new_selected_vars
    acceptance_rate <- acceptance_rate + 1
  } else {
    # Do nothing
  }
  selected_vars_list[[mcmc_iter]] <- selected_vars


  # Updating betas according to the current model
  for(ii in 1:length(selected_vars)){
    # Other predictions
    residuals_ <- numeric(n)

    if(length(selected_vars)>1){
      other_pred <- selected_vars[-ii]
      for(jj in 1:length(other_pred)){
        residuals_  <- tcrossprod(list_B[[other_pred[jj]]],(betas[[other_pred[[jj]]]]))
      }
    } else {
      residuals_ <- rep(0,n)
    }
    betas[[selected_vars[ii]]] <- mvnfast::rmvn(n = 1,mu = solve((crossprod(list_B[[selected_vars[ii]]]) + diag(tau_beta_j[ii]*tau^-1, nrow = NCOL(list_B[[selected_vars[ii]]]))),
                                                          crossprod(list_B[[selected_vars[ii]]],(y - residuals_))), sigma = (tau^(-1))*solve(crossprod(list_B[[selected_vars[ii]]]) + diag(tau_beta_j[ii]*tau^-1, nrow = NCOL(list_B[[selected_vars[ii]]])) ) )

  }

  for(jj in 1:length(tau_beta_j)){
    a_shape <- 0.5 + a_tau_beta_j_prior[jj]
    d_rate <- 0.5*sum(betas[[jj]])^2 + d_tau_beta_j_prior[jj]
    tau_beta_j[jj] <- stats::rgamma(n = 1,shape = a_shape,rate = d_rate)
  }

  tau_beta_j_post[mcmc_iter,] <- tau_beta_j
  cat("Iteration number", mcmc_iter,"\n")
}

# Counting the number of acceptances
# par(mfrow=c(1,3))
# for(i in 1:3){
#   plot(betas[,i],type= "l", main = paste("Accepted",sum(betas[,i]!=0)))
# }
#
# for(i in 1:3){
#   plot(tau_beta_j_post[,i],type= "l", main = paste("Accepted",sum(betas[,i]!=0)))
# }

# Seeing frequency of selected vars
lapply(selected_vars_list,function(x)paste0(as.character(x),collapse = "-")) %>% unlist %>% table()
selected_vars_list_vec <- lapply(selected_vars_list,function(x)paste0(as.character(x),collapse = "-")) %>% unlist

# Visualising the frequency of all combination
selected_vars_list_vec <- ifelse(selected_vars_list_vec=="1",1,selected_vars_list_vec)
selected_vars_list_vec <- ifelse(selected_vars_list_vec=="2",2,selected_vars_list_vec)
selected_vars_list_vec <- ifelse(selected_vars_list_vec=="2",3,selected_vars_list_vec)
selected_vars_list_vec <- ifelse(selected_vars_list_vec=="1-2",4,selected_vars_list_vec)
selected_vars_list_vec <- ifelse(selected_vars_list_vec=="1-3",5,selected_vars_list_vec)
selected_vars_list_vec <- ifelse(selected_vars_list_vec=="2-3",6,selected_vars_list_vec)
selected_vars_list_vec <- ifelse(selected_vars_list_vec=="1-2-3",7,selected_vars_list_vec)
par(mfrow=c(1,1))
plot(1:1000,as.numeric(selected_vars_list_vec),axes = FALSE,ylab = "Parameters selected",pch=20)
axis(2,at = 1:7,labels = c("1","2","3","1-2","1-3","2-3","1-2-3"))
axis(1,at = seq(0,1000,by=250))
