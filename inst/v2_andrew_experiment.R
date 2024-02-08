library(splines)
library(mvnfast)
library(mlbench)
rm(list=ls())
set.seed(42)
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

# Can we fit the bigger model
B_all <- vector(mode = 'list', length = 10 + 2)
for (j in 1:10) B_all[[j]] <- bs(dat$x[,j],intercept = TRUE)

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

B_all[[11]] <- multiply_matrices_general(B_all[[1]],B_all[[2]])
B_all[[12]] <- multiply_matrices_general(B_all[[1]],B_all[[4]])


nll2 <- function(pars) {
  sigma <- pars[1]
  sigma_b <- pars[2:13]
  Big_var <- matrix(0, n, n)

  for(j in 1:12) Big_var <- Big_var + sigma_b[j] * B_all[[j]] %*% t(B_all[[j]])
  return(-mvnfast::dmvn(y, rep(0, n),
                        diag(sigma^2, n) + Big_var,
                        log = TRUE))
}
nll2(rep(1, 13))
ans2 <- nlminb(rep(1, 13), nll2, lower = c(0, 0))
