# rm(list=ls())
library(ggplot2)
# rsp_mod <- readRDS("~/spline_bart_lab/preliminar_results/rspBART17/friedman/single_run/v20_single_run_rep_1_n_250_sd_1_nIknots_2_ntree_20_nodesize_15_dif_1_scale_TRUE_sc_basis_TRUE_nmcmc_5000_nburn_2500_rb_prior_FALSE.Rds")


# Main effect range difference
main_effects_train_list_norm <- rsp_mod$mcmc$main_effects_train
tree_main_effects <- rsp_mod$mcmc$tree_main_effects
x_train <- rsp_mod$data$x_train %>% as.matrix()
n_mcmc <- rsp_mod$mcmc$n_mcmc
n_burn <- rsp_mod$mcmc$n_burn

# n_mcmc <- 9200
# n_burn <- 5000

n_tree <- rsp_mod$prior$n_tree
n_burn_plot <- 4000
par(mfrow=c(1,1))
plot(rsp_mod$all_tau, type = 'l', main = expression(tau), ylab = expression(tau))

par(mfrow = c(2,floor(NCOL(x_train)/2)))

for(jj in 1:11){

  if(jj <= NCOL(x_train)){
    plot(x_train[,jj],colMeans(main_effects_train_list_norm[[jj]][n_burn_plot:n_mcmc,, drop = FALSE]),main = paste0('X',jj),
         ylab = paste0('G(X',jj,')'),pch=20,ylim = c(-10,10),xlab = paste0('x.',jj), col = alpha("black",1.0))
  }    else if(jj == 11 ) {
    # plot(x_train[,1], colMeans(main_effects_train_list_norm[[jj]][1:i,,drop = FALSE]), ylab = "G(x.1,x2)", xlab = "X.1")
    # plot(x_train[,2], colMeans(main_effects_train_list_norm[[jj]][1:i,,drop = FALSE]), ylab = "G(x.1,x2)", xlab = "X.1")
    # plot(x_train[,1], x_train[,2], cex = colMeans(main_effects_train_list_norm[[jj]][1000:i,,drop = FALSE]) + colMeans(main_effects_train_list_norm[[1]][1000:i,,drop = FALSE]) + colMeans(main_effects_train_list_norm[[2]][1000:i,,drop = FALSE]),
    #      ylab = "X.2", xlab = "X.1")
    plot(x_train[,1], x_train[,2], cex = abs(colMeans(main_effects_train_list_norm[[jj]][n_burn_plot:n_mcmc,,drop = FALSE])), ylab = "X.2", xlab = "X.1")
    # plot(x_train[,2], x_train[,1], cex = 10*sin(pi*x_train[,1]*x_train[,2]), col = "blue",  ylab = "X.2", xlab = "X.1")

  }


  for(tree_number in 1:n_tree){

        n_burn_plot_count <- 1
        # Create an auxiliar matrix for the main effect
        aux_main_effect_matrix <- matrix(0, nrow = (n_mcmc-n_burn_plot),ncol = nrow(x_train))

        for(mcmc_aux in (n_burn_plot+1):n_mcmc){
          aux_main_effect_matrix[n_burn_plot_count,] <- unnormalize_bart_me(tree_main_effects[[mcmc_aux]][,jj,tree_number],a = min(rsp_mod$data$y_train),b = max(rsp_mod$data$y_train))
          n_burn_plot_count <- n_burn_plot_count  + 1
        }

        if(jj <= NCOL(x_train)){
          if(all(colMeans(aux_main_effect_matrix)==0)){
            # print("NOTHING")
          } else {
            points(x_train[,jj],colMeans(aux_main_effect_matrix),main = paste0('X',jj),
                   col = ggplot2::alpha("darkred",0.5), pch = 20)
          }
        }
    }
    # range_tree_bais_j_predictions[tree_number, jj] <- diff(range(colMeans(aux_main_effect_matrix)))

}


par(mfrow=c(1,2))
burn_sample_ <- 4000
all_tau_beta <- rsp_mod$all_tau_beta
variable_importance_matrix <- rsp_mod$mcmc$variable_importance_matrix
plot(1:NCOL(variable_importance_matrix),variable_importance_matrix[burn_sample_:n_mcmc,,drop = FALSE] %>% colMeans(),
     ylab = "Prop. pred_var", xlab = "Predictor", main = c("Proportion Tree pred.vars"))
points((1:NCOL(variable_importance_matrix))[c(1:5,11)],variable_importance_matrix[burn_sample_:n_mcmc,c(1:5,11),drop = FALSE] %>% colMeans(),
       ylab = "Prop. pred_var", xlab = "Predictor/Basis", pch = 20)


plot(1:NCOL(variable_importance_matrix),all_tau_beta[burn_sample_:n_mcmc,,drop = FALSE] %>% colMeans(na.rm = TRUE),
     ylab = expression(bar(lambda[j])), xlab = "Predictor", main = c("Lambda_posterior_mean"))
points((1:NCOL(variable_importance_matrix))[c(1:5,11)],all_tau_beta[burn_sample_:n_mcmc,c(1:5,11),drop = FALSE] %>% colMeans(na.rm = TRUE),
       ylab = "mean_lambda_j", xlab = "Predictor/Basis", pch = 20)

par(mfrow = c(1,1))


# rsp_mod$all_tau_beta %>% apply(2,var) %>% plot
# rsp_mod$all_tau_beta[, c(1:5,11),drop = FALSE] %>% apply(2,var) %>% points(pch= 20)
#
rmse(x = rsp_mod$y_train_hat[4001:n_mcmc,,drop = FALSE] %>% colMeans(), rsp_mod$data$y_train)
rmse(x = rsp_mod$y_test_hat[4001:n_mcmc,,drop = FALSE] %>% colMeans(), y_test)

all_tau_beta[burn_sample_:n_mcmc,,drop = FALSE] %>% apply(2,sd) %>% plot
points((1:NCOL(variable_importance_matrix))[c(1:5,11)],apply(all_tau_beta[burn_sample_:n_mcmc,c(1:5,11),drop = FALSE],2,sd),
       ylab = "mean_tau_beta", xlab = "Predictor/Basis", pch = 20)

boxplot(all_tau_beta[4001:n_mcmc,,drop = FALSE], ylab = expression(lambda[j]), main = expression(lambda[j]))

