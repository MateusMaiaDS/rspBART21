rm(list=ls())
library(splines)
library(mvnfast)
library(mlbench)
set.seed(42)
n <- 500
dat <- mlbench.friedman1(n, sd=1)
dat_df <- as.data.frame(dat)
b1_var <- 4
b2_var <- 10
B1 <- bs(dat$x[,b1_var],intercept = TRUE)
B2 <- bs(dat$x[,b2_var],intercept = TRUE)
y_raw <- dat$y


# Creating a boolean to decide if scale or not
scale_y_bool <- FALSE
if(scale_y_bool){
  y <- (y_raw - min(y_raw))/(max(y_raw) - min(y_raw)) - 0.5
  sigma <- 1 / diff(range(y_raw))
  sigma_b <- (0.5/(2*sqrt(1))) # Sigma_b = 0.25/k*sqrt(T), where T is the number of trees
} else {
  sigma <- 1
  sigma_b <- (0.5/(2*sqrt(1)))/diff(range(y_raw))
  y <- y_raw
}

sigma_b <- 1

# Dealing with interactions

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

# Defining which interaction basis gonna be used.
b1_inter_ <- 1
b2_inter_ <- 2
B1_inter <- bs(dat$x[,b1_inter_],intercept = FALSE)
B2_inter <- bs(dat$x[,b2_inter_],intercept = FALSE)
B1_B2 <- multiply_matrices_general(A = as.matrix(B1_inter),B = as.matrix(B2_inter))

ll1 <- dmvn(y, rep(0, n), diag(sigma^2, n) + 0.1^2 * B1 %*% t(B1), log = TRUE) # Likelihood of only having the main effect
ll2 <- dmvn(y, rep(0, n), diag(sigma^2, n) + 1^2 * B1_B2 %*% t(B1_B2), log = TRUE) # Likelihood of having the interaction only
ll3 <- dmvn(y, rep(0, n), diag(sigma^2, n) + sigma_b^2 * B1 %*% t(B1) + sigma_b^2 * B1_B2 %*% t(B1_B2), log = TRUE) # Likelihood plus main effect + one interaction


print(ll1)
print(ll2)
print(ll3)

# Accepting the interaction only term
print(exp(ll2-ll1))
# Accepting interaction + main effect
print(exp(ll3-ll1))

