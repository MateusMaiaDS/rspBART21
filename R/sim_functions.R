# Just running the default values so I can go through the function and debug all
#the things
# library(purrr)
# library(tidyverse)
#Auxliar object
# friedman_df <- std.mlbench.friedman1(n = 100,sd = 1) %>% as.data.frame()
# plot_nrow = 2
# break_friedman_df <- std.mlbench.friedman1(n = 250,sd = 1) %>% as.data.frame()


plot_spBART_me <- function(train_me,x_train, plot_nrow){


  # Getting the
  par(mfrow = c(4,3))
  fx1x2 <-colMeans(main_effects_train_list[[11]]) #+ main_effects_train_list[[1]] %>% colMeans()

  interaction_df <- cbind(x_train[,1:2],fx1x2) %>% dplyr::arrange(fx1x2) %>% mutate(color = gray.colors(nrow(x_train)))

  # Plotting interaction
  plot(interaction_df[,1:2],cex = interaction_df$fx1x2*20,col = interaction_df$color, pch=19,main = "f(x.1,x.2)")

  # Plotting x.3
  fx1 <- colMeans(main_effects_train_list[[1]])
  plot(x_train[,1],fx1,pch=20, xlab = "x.1", main = "f(x.1)")

  # Plotting x.3
  fx2 <- colMeans(main_effects_train_list[[2]])
  plot(x_train[,2],fx2,pch=20, xlab = "x.2", main = "f(x.2)")

  # Plotting x.3
  fx3 <- colMeans(main_effects_train_list[[3]])
  plot(x_train[,3],fx3,pch=20, xlab = "x.3", main = "f(x.3)")

  # Plotting x.4
  fx4 <- colMeans(main_effects_train_list[[4]])
  plot(x_train[,4],fx4,pch=20, xlab = "x.4", main = "f(x.4)")

  # Plotting x.5
  fx5 <- colMeans(main_effects_train_list[[5]])
  plot(x_train[,5],fx5,pch=20, xlab = "x.5", main = "f(x.5)")

  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  for(i in 6:10){
    plot(x_train[,i],colMeans(main_effects_train_list[[i]]),pch=20,
         xlab = paste0("x.",i,""), main = paste0("f(x.",i,")"),
         ylab = paste0("f(x.",i,")"))
  }


}

plot_break_friedman <- function(break_friedman_df, plot_nrow){


  # Getting the
  par(mfrow = c(plot_nrow,(NCOL(break_friedman_df)-1)/plot_nrow))
  fx1x2 <- 10* sin(pi * break_friedman_df[, 1] * break_friedman_df[, 2])
  fx1x2[break_friedman_df[,1]<0.5] <- fx1x2[break_friedman_df[,1]<0.5]  - 5
  fx1x2[break_friedman_df[,1]>0.5] <- fx1x2[break_friedman_df[,1]>0.5]  + 5

  interaction_df <- cbind(break_friedman_df[,1:2],fx1x2) %>% dplyr::arrange(fx1x2) %>% mutate(color = gray.colors(nrow(break_friedman_df)))

  # Plotting interaction
  plot(interaction_df[,1:2],cex = interaction_df$fx1x2*0.3,col = interaction_df$color, pch=19,main = "f(x.1,x.2)")

  # Plotting x.3
  fx3 <- numeric(nrow(break_friedman_df))
  fx3[break_friedman_df[,3]>0.5] <- fx3[break_friedman_df[,3]>0.5] + 20*(break_friedman_df$x.3[break_friedman_df$x.3>0.5]-0.5)^2 + 5
  fx3[break_friedman_df[,3]<0.5] <- fx3[break_friedman_df[,3]<0.5] - 20*(break_friedman_df$x.3[break_friedman_df$x.3<0.5]-0.5)^2 - 5

  plot(break_friedman_df$x.3,fx3,pch=20, xlab = "x.3", main = "f(x.3)")

  # Plotting x.4
  fx4 <- numeric(nrow(break_friedman_df))
  fx4[break_friedman_df[,4]>0.3] <- fx4[break_friedman_df[,4]>0.3] + 10*(break_friedman_df$x.4[break_friedman_df$x.4>0.3]) + 5
  fx4[break_friedman_df[,4]<0.3] <- fx4[break_friedman_df[,4]<0.3] - 15*(break_friedman_df$x.4[break_friedman_df$x.4<0.3]) - 5

  plot(break_friedman_df$x.4,fx4,pch=20, xlab = "x.4", main = "f(x.4)")

  # Plotting x.5
  fx5 <- 5*break_friedman_df$x.5
  plot(break_friedman_df$x.5,fx5,pch=20, xlab = "x.5", main = "f(x.5)")

  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  for(i in 6:10){
    plot(break_friedman_df[,i],rep(0,nrow(break_friedman_df)),pch=20,
         xlab = paste0("x.",i,""), main = paste0("f(x.",i,")"),
         ylab = paste0("f(x.",i,")"))
  }

}

# Plotting main effects  friedman1
plot_friedman1 <- function(friedman_df, plot_nrow = 2){

  # Getting the
  par(mfrow = c(plot_nrow,(NCOL(friedman_df)-1)/plot_nrow))
  fx1x2 <- 10* sin(pi * friedman_df[, 1] * friedman_df[, 2])
  interaction_df <- cbind(friedman_df[,1:2],fx1x2) %>% dplyr::arrange(fx1x2) %>% mutate(color = gray.colors(nrow(friedman_df)))

  # Plotting interaction
  plot(interaction_df[,1:2],cex = interaction_df$fx1x2*0.3,col = interaction_df$color, pch=19,main = "f(x.1,x.2)")

  # Plotting x.3
  fx3 <- 20*(friedman_df$x.3-0.5)^2
  plot(friedman_df$x.3,fx3,pch=20, xlab = "x.3", main = "f(x.3)")

  # Plotting x.4
  fx4 <- 10*friedman_df$x.4
  plot(friedman_df$x.4,fx4,pch=20, xlab = "x.4", main = "f(x.4)")

  # Plotting x.5
  fx5 <- 5*friedman_df$x.5
  plot(friedman_df$x.5,fx5,pch=20, xlab = "x.5", main = "f(x.5)")

  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  for(i in 6:10){
    plot(friedman_df[,i],rep(0,nrow(friedman_df)),pch=20,
         xlab = paste0("x.",i,""), main = paste0("f(x.",i,")"),
         ylab = paste0("f(x.",i,")"))
  }

}


break.mlbench.friedman1 <- function (n, sd = 1, scale_main = FALSE)
{
  x <- matrix(runif(10 * n), ncol = 10)
  y <- matrix(0, nrow = n, ncol = 4)

  # Adding a break
  y[,1] <- 10 * sin(pi * x[, 1] * x[, 2])
  # y[ x[,1] > 0.5] < -  y[ x[,1] > 0.5] + 5
  # y[ x[,1] < 0.5] < -  y[ x[,1] < 0.5] - 5

  # Adding a break for x.3
  y[x[,3] > 0.5,2] <-  + 20* (x[x[,3] > 0.5, 3] - 0.5)^2 + 5
  y[x[,3] < 0.5,2] <-  - 20* (x[x[,3] < 0.5, 3] - 0.5)^2 - 5

  y[x[,4] > 0.3,3] <-  + 10* (x[x[,4] > 0.3, 4]) + 5
  y[x[,4] < 0.3,3] <-  - 15* (x[x[,4] < 0.3, 4]) - 5

  y[,4] <- y[,4] + 5 * x[, 5]

  if(scale_main){
    y <- apply(y,2,normalize_main_effects)
  }

  if (sd > 0) {
    y_response <- rowSums(y) + rnorm(n, sd = sd)
  }
  list(x = x, y = y_response)
}

std.mlbench.friedman1 <- function (n, sd = 1)
{
  x <- matrix(runif(10 * n), ncol = 10)
  y <- 10 * sin(pi * x[, 1] * x[, 2])
  y <- y + 20 * (x[, 3] - 0.5)^2 + 10 * x[, 4] + 5 * x[, 5]
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}


mlbench.friedman1.nointeraction <- function (n, sd = 1)
{
  x <- matrix(runif(4 * n), ncol = 4)
  y <- 10 * sin(pi * x[, 1])
  y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.friedman1.interaction.only <- function (n, sd = 1)
{
  x <- matrix(runif(2* n), ncol = 2)
  y <- 10 * sin(pi * x[, 1]*x[,2])
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.friedman1.nointeraction.noise <- function (n, sd = 1)
{
  x <- matrix(runif(8 * n), ncol = 8)
  y <- 10 * sin(pi * x[, 1])
  y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.d1 <- function(n, sd = 1) {
  x <- matrix(runif(n,min = -pi,max = pi),ncol = 1)
  y <- sin(2*x)
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.d1.break <- function(n, sd = 1) {
  x <- matrix(runif(n,min = -pi,max = pi),ncol = 1)
  y <- sin(2*x)
  y[x<0] <- y[x<0] + 5
  y[x>=0] <- y[x>=0] - 5

  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.d1 <- function(n, sd = 1) {
  x <- matrix(runif(n,min = -pi,max = pi),ncol = 1)
  y <- sin(2*x)

  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}


# ====================================
# Generating new simulation functions
# ====================================
# n_ <- 250
# sd_ <- 1
# p_ <- 10

# Normalising the main effects to be between -1 and 1
normalize_main_effects <- function(vec_){
  a <- min(vec_)
  b <- max(vec_)
  if(a==b){
    norm_vec_ <- rep(0.0,length(vec_))
  } else {
    norm_vec_ <- (vec_-a)/(b-a) - 0.5
  }
  return(norm_vec_)
}

# Formula for the model with only main effects
smooth.main.formula <- function(x_){
  y.1 <- sin(2*pi*x_[,1])*10
  y.2 <- (2*(x_[,2])^3)*10
  y.3 <- (exp(-2*x_[,3]))*10
  y.4 <- (cos(x_[,4]))*10

  # Generating the main effect matrix
  main.effects <- cbind(y.1,y.2,y.3,y.4)
  colnames(main.effects) <- paste0("y.",1:4)

  # main.effects.norm <- apply(main.effects,2,normalize_main_effects)
  main.effects.norm <- main.effects

  return(main.effects.norm)
}

# Formula for the model with only main effects
non.smooth.main.formula <- function(x_){
  y.1 <- 5*x_[,1]
  y.2 <- rep(2,nrow(x_))
  y.2[x_[,2] < -0.3] <- 5
  y.2[(x_[,2] >= -0.3) & (x_[,2] < 0.5) ] <- 10.0

  y.3 <- rep(7,nrow(x_))
  y.3[x_[,3]<0] <- -3
  y.4 <- -5*x_[,4]

  # Generating the main effect matrix
  main.effects <- cbind(y.1,y.2,y.3,y.4)
  colnames(main.effects) <- paste0("y.",1:4)


  # main.effects.norm <- apply(main.effects,2,normalize_main_effects)
  main.effects.norm <- main.effects

  return(main.effects.norm)
}


# Formula for the model with only main effects
non.and.smooth.main.formula <- function(x_){
  y.1 <- sin(2*pi*x_[,1])
  y.2 <- rep(2,nrow(x_))
  y.2[x_[,2] < -0.3] <- 5
  y.2[(x_[,2] >= -0.3) & (x_[,2] < 0.5) ] <- 10.0

  y.3 <- rep(7,nrow(x_))
  y.3[x_[,3]<0] <- -3
  y.4 <- cos(x_[,4])

  # Generating the main effect matrix
  main.effects <- cbind(y.1,y.2,y.3,y.4)
  colnames(main.effects) <- paste0("y.",1:4)


  # main.effects.norm <- apply(main.effects,2,normalize_main_effects)
  main.effects.norm <- main.effects

  return(main.effects.norm)
}


sim.gen <- function(n_,
                    sd_,
                    p_,
                    formula_){

  # Generating the setting of covariates
  x_ <- matrix(runif(n = n_*p_,min = -1,max = 1),ncol = p_)

  # Generating the response
  f <- rowSums(formula_(x_ = x_))
  y <- f + rnorm(n = n_,sd = sd_)

  # Returning the data
  return(list(x = x_,
              y = y,
              f = f))
}


# ===============================
#     Setting plot functions
# ===============================

plot.smooth.main <- function(x_){

  # Recreating the main effects
  main.effects.norm_ <- smooth.main.formula(x_ = x_)

  # Setting the plot.grid.window
  par(mfrow = c(2,2))
  for(i in 1:4){
    plot(x_[,i],main.effects.norm_[,i], xlab = paste0("x.",i), ylab = paste0("f(x.",i,")"), pch = 20)
  }


}

plot.non.smooth.main <- function(x_){

  # Recreating the main effects
  main.effects.norm_ <- non.smooth.main.formula(x_ = x_)

  # Setting the plot.grid.window
  par(mfrow = c(2,2))
  for(i in 1:4){
    plot(x_[,i],main.effects.norm_[,i], xlab = paste0("x.",i), ylab = paste0("f(x.",i,")"), pch = 20)
  }


}

plot.non.and.smooth.main <- function(x_){

  # Recreating the main effects
  main.effects.norm_ <- non.and.smooth.main.formula(x_ = x_)

  # Setting the plot.grid.window
  par(mfrow = c(2,2))
  for(i in 1:4){
    plot(x_[,i],main.effects.norm_[,i], xlab = paste0("x.",i), ylab = paste0("f(x.",i,")"), pch = 20)
  }


}

# Testing generating
# set.seed(42)
sim_data <- sim.gen(n_ = 100,sd_ = 42,p_ = 10,formula_ = non.smooth.main.formula)
#
plot.non.smooth.main(x_ = sim_data$x)



