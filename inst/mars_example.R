library(mlbench)
rm(list=ls())
devtools::load_all()
set.seed(42)
n_ <- 250
sd_ <- 1
# sim_train <- mlbench.friedman1.nointeraction(n = n_,sd = sd_)  |> as.data.frame()
# sim_test <- mlbench.friedman1.nointeraction(n = n_,sd = sd_)  |> as.data.frame()
#
# sim_train <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_)  |> as.data.frame()
# sim_test <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_)  |> as.data.frame()

sim_train <- mlbench.friedman1(n = n_,sd = sd_)  |> as.data.frame()
sim_test <- mlbench.friedman1(n = n_,sd = sd_)  |> as.data.frame()

# sim_train <- mlbench.d1.break(n = n_,sd = sd_)  |> as.data.frame() %>% arrange(x)
# sim_test <- mlbench.d1.break(n = n_,sd = sd_) |> as.data.frame() %>% arrange(x)

# sim_train <- mlbench.d1(n = n_,sd = 1)  |> as.data.frame()
# sim_test <- mlbench.d1(n = n_,sd = 1) |> as.data.frame()

x_train <- sim_train |> dplyr::select(dplyr::starts_with("x"))
x_test <-  sim_test|> dplyr::select(dplyr::starts_with("x"))
y_train <- sim_train$y

mars1 <- earth(y_train~., data = cbind(x_train,y_train),degree = 2)
summary(mars1)
new_pred_mars <- predict(mars1,x_test)
rmse(new_pred_mars,sim_test$y)
# rmse(new_pred_mars,sim_test$y)
# summary(mars1)
# plot(sim_train,pch=20)
lines(x_train$x,mars1$fitted.values,col = "blue")
