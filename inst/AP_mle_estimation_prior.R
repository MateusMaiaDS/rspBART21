library(splines)
library(mvnfast)
library(mlbench)
rm(list=ls())
n <- 500
dat <- mlbench.friedman1(n, sd=1)
dat_df <- as.data.frame(dat)
b_var <- 4
B <- bs(dat$x[,b_var],intercept = TRUE)
y <- dat_df$y

# Have a first guess at an upper value for sigma_b
nll <- function(pars) {
  sigma <- pars[1]
  sigma_b <- pars[2]
  return(-mvnfast::dmvn(y, rep(0, n),
                        diag(sigma^2, n) + sigma_b^2 * B %*% t(B),
                        log = TRUE))
}
ans <- nlminb(c(1, 1), nll, lower = c(0, 0))
ans$par[2]

# Can we simulate beta?
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


# Can we fit the bigger model
interaction_list <- combn(1:NCOL(dat$x),2,simplify = TRUE)
B_all <- vector(mode = 'list', length = NCOL(dat$x) + NCOL(interaction_list))
for (j in 1:10) B_all[[j]] <- bs(dat$x[,j],intercept = FALSE)

# Adding interactions
for(jj in 1:(NCOL(interaction_list))){
  jj_ = jj + NCOL(dat$x)
  B_all[[jj_]] <- multiply_matrices_general(A = B_all[[interaction_list[1,jj]]],B = B_all[[interaction_list[2,jj]]])
}

# Number of all basis
p_all <- length(B_all)

nll2 <- function(pars) {
  sigma <- pars[1]
  sigma_b <- pars[2:(length(B_all)+1)]
  Big_var <- matrix(0, n, n)2+2

  for(j in 1:p_all) Big_var <- Big_var + sigma_b[j]^2 * B_all[[j]] %*% t(B_all[[j]])
  return(-mvnfast::dmvn(y, rep(0, n),
                        diag(sigma^2, n) + Big_var,
                        log = TRUE))
}

ans2 <- nlminb(rep(1, p_all+1), nll2, lower = rep(1e-8, p_all+1))
