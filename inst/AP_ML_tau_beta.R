library(splines)
library(mvnfast)

set.seed(42)
n <- 500
dat <- mlbench.friedman1(n, sd=1)
dat_df <- as.data.frame(dat)
mle_tau_b <- numeric(NCOL(dat_df)-1)
for(jj in 1:length(mle_tau_b)){

  b_var <- jj
  B1 <- bs(dat$x[,b_var],intercept = FALSE)
  y <- dat_df$y

  # Have a first guess at an upper value for sigma_b
  nll <- function(pars) {
    sigma <- pars[1]
    sigma_b <- pars[2]
    return(-mvnfast::dmvn(y, rep(0, n),
                          diag(sigma^2, n) + sigma_b^2 * B1 %*% t(B1),
                          log = TRUE))
  }
  ans <- nlminb(c(1, 1), nll, lower = c(0, 0))
  mle_tau_b[jj] <- ans$par[2]

}
