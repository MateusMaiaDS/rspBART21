library(MOTRbart)
model <- motr_bart(x = x_train,y = y_train,sparse = FALSE,ancestors = TRUE)
# plot(x_train$x,y_train)
# points(x_train$x,model$y_hat |> colMeans(),pch = 20)
colMeans(model$vars_betas)
plot(model$sigma2^(-1),type = "l")

plot(model$y_hat %>% colMeans(),unnormalize_bart(z = colMeans(all_y_hat),a = min_y,b = max_y))
