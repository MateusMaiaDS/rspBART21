all_bart <- function(cv_element,
                     nIknots_,
                     ntree_,
                     use_bs_,
                     seed_,
                     motr_bart_,
                     rsp_bart_all_,
                     alpha_ = 0.95,
                     stump_){

  # To replicate the results
  set.seed(seed_)
  train <- cv_element$train
  test <- cv_element$test

  # Getting the training elements
  x_train <- train %>% dplyr::select(dplyr::starts_with("x"))
  x_test <- test %>% dplyr::select(dplyr::starts_with("x"))
  y_train <- train %>% dplyr::pull("y")

  # Running the model
  spBART <- rspBART(x_train = x_train,
                    x_test = x_test,y_train = y_train,
                    n_mcmc = 2500,node_min_size = 25,alpha = alpha_,
                    n_burn = 0,nIknots = nIknots_,n_tree = ntree_,
                    use_bs = use_bs_,scale_bool = TRUE,plot_preview = FALSE,
                    motrbart_bool = motr_bart_)

  if(rsp_bart_all_){
    spBART_all <- rspBART(x_train = x_train,
                                    x_test = x_test,y_train = y_train,
                                    n_mcmc = 2500,node_min_size = 15, alpha = alpha_,
                                    n_burn = 500,nIknots = nIknots_,n_tree = ntree_,
                                    use_bs = use_bs_,scale_bool = TRUE,plot_preview = FALSE,
                                    motrbart_bool = motr_bart_,
                                    all_var = rsp_bart_all_)
  }

  bartmod <- dbarts::bart(x.train = x_train,y.train = y_train,x.test = x_test)

  # Special case when I would have univariate regression
  if(ncol(x_train)>=2){
    softbartmod <- SoftBart::softbart(X = x_train,Y = y_train,X_test =  x_test)
    motr_bart_mod <- motr_bart(x = x_train,y = y_train,ntrees = ntree_)
    # motrbart_pred <- predict_motr_bart(object = motr_bart_mod,newdata = x_test,type = "all")

    if(ntree_>1){
      motrbart_pred <- predict_motr_bart(object = motr_bart_mod,newdata = x_test,type = "all")
    } else {
      motrbart_pred <- NULL
    }

  } else {
    x_train_new <- cbind(x_train,x_train)
    x_test_new <- cbind(x_test,x_test)
    colnames(x_test_new) <- colnames(x_train_new) <- paste0("x.",1:2)

    softbartmod <- SoftBart::softbart(X = x_train_new,Y = y_train,X_test =  x_test_new)
    motr_bart_mod <- motr_bart(x = x_train_new,y = y_train,ntrees = ntree_)

    # motr_bart_mod$y_hat
    if(ntree_>1){
      motrbart_pred <- predict_motr_bart(object = motr_bart_mod,newdata = x_test_new,type = "all")
    } else {
      motrbart_pred <- NULL
    }
  }


  if(rsp_bart_all_){
    main_result_list <- list(spBART = spBART,
                             spBART_all = spBART_all,
                             bartmod = bartmod,
                             softbartmod = softbartmod,
                             motrbartmod = motr_bart_mod,
                             motrbart_pred  = motrbart_pred,
                             cv = cv_element)
  } else {
    main_result_list <- list(spBART = spBART,
                             bartmod = bartmod,
                             softbartmod = softbartmod,
                             motrbartmod = motr_bart_mod,
                             motrbart_pred  = motrbart_pred,
                             cv = cv_element)
  }

  return(main_result_list)

}


all_bart_lite <- function(cv_element,
                     nIknots_,
                     ntree_,
                     seed_,
                     use_bs_,
                     j,
                     motr_bart_,
                     alpha_,
                     rsp_bart_all_,
                     stump_,
                     scale_init_,
                     update_tau_beta_,
                     dif_order_,
                     interaction_term_){


  # Doing a warming for the case whichI don't have
  # if(ntree_<50){
  #   stop("Use the all_bart() function instead.")
  # }

  # To replicate the results
  set.seed(seed_)
  train <- cv_element$train
  test <- cv_element$test

  # Getting the training elements
  x_train <- train %>% dplyr::select(dplyr::starts_with("x"))
  x_test <- test %>% dplyr::select(dplyr::starts_with("x"))
  y_train <- train %>% dplyr::pull("y")
  y_test <- test %>% dplyr::pull("y")

  # Initialising df
  comparison_metrics <- data.frame(metric = NULL, value = NULL, model = NULL,fold = NULL)

  # # Running the model
  # spBART <- rspBART(x_train = x_train,
  #                   x_test = x_test,y_train = y_train,
  #                   n_mcmc = 2500,node_min_size = 15,alpha = alpha_,
  #                   n_burn = 0,nIknots = nIknots_,n_tree = ntree_,
  #                   use_bs = use_bs_,all_var = FALSE,stump = FALSE,
  #                   motrbart_bool = motr_bart_,
  #                   scale_init = scale_init_,update_tau_beta = update_tau_beta_)
  #
  #
  # n_burn_ <- 500
  # n_mcmc_ <- spBART$mcmc$n_mcmc
  #
  # # Calculating metrics for splinesBART
  # comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
  #                                                           value = rmse(x = colMeans(spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                        y = train$y),
  #                                                           model = "spBART",fold = j))
  #
  # comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
  #                                                           value = rmse(x = colMeans(spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                        y = test$y),
  #                                                           model = "spBART",fold = j))
  #
  # # Calculating the CRPS as well
  # comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
  #                                                           value = crps(y = train$y ,
  #                                                                        means = colMeans(spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                        sds = rep(mean(spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(train$y)))$CRPS,
  #                                                           model = "spBART",fold = j))
  #
  # comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
  #                                                           value = crps(y = test$y ,
  #                                                                        means = colMeans(spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                        sds = rep(mean(spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(test$y)))$CRPS,
  #                                                           model = "spBART",fold = j))
  #
  # # Removing the model
  # rm(spBART)

  # if(rsp_bart_all_){
    # Running the model

    if(interaction_term_){
    spBART <- rspBART(x_train = x_train,
                      x_test = x_test,y_train = y_train,
                      n_mcmc = 2500,node_min_size = 15,alpha = alpha_,
                      n_burn = 0,nIknots = nIknots_,n_tree = ntree_,
                      use_bs = use_bs_,all_var = rsp_bart_all_,
                      stump = stump_,dif_order = dif_order_,
                      motrbart_bool = motr_bart_,
                      scale_init = scale_init_,
                      interaction_term = interaction_term_,
                      interaction_list = list(c(1,2)),
                      update_tau_beta = update_tau_beta_)

    } else {
      spBART <- rspBART(x_train = x_train,
                        x_test = x_test,y_train = y_train,
                        n_mcmc = 2500,node_min_size = 15,alpha = alpha_,
                        n_burn = 0,nIknots = nIknots_,n_tree = ntree_,
                        use_bs = use_bs_,all_var = rsp_bart_all_,
                        stump = stump_,dif_order = dif_order_,
                        motrbart_bool = motr_bart_,
                        scale_init = scale_init_,
                        interaction_term = FALSE,
                        interaction_list = NULL,
                        update_tau_beta = update_tau_beta_)
    }
    n_burn_ <- 500
    n_mcmc_ <- spBART$mcmc$n_mcmc

    # Calculating metrics for splinesBART
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(x = colMeans(spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           y = train$y),
                                                              model = "spBART_all",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(x = colMeans(spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           y = test$y),
                                                              model = "spBART_all",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = train$y ,
                                                                           means = colMeans(spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           sds = rep(mean(spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(train$y)))$CRPS,
                                                              model = "spBART_all",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = test$y ,
                                                                           means = colMeans(spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           sds = rep(mean(spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(test$y)))$CRPS,
                                                              model = "spBART_all",fold = j))

    # Removing the model
    rm(spBART)
  # }


  #   if(rsp_bart_all_){
  #     # Running the model
  #     spBART <- rspBART(x_train = x_train,
  #                       x_test = x_test,y_train = y_train,
  #                       n_mcmc = 2500,node_min_size = 15,alpha = alpha_,
  #                       n_burn = 0,nIknots = nIknots_,n_tree = ntree_,
  #                       use_bs = use_bs_,all_var = rsp_bart_all_,dif_order = 1,
  #                       motrbart_bool = motr_bart_,stump = stump_,
  #                       scale_init = scale_init_,update_tau_beta = update_tau_beta_)
  #
  #
  #     n_burn_ <- 500
  #     n_mcmc_ <- spBART$mcmc$n_mcmc
  #
  #     # Calculating metrics for splinesBART
  #     comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
  #                                                               value = rmse(x = colMeans(spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                            y = train$y),
  #                                                               model = "spBART_stump",fold = j))
  #
  #     comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
  #                                                               value = rmse(x = colMeans(spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                            y = test$y),
  #                                                               model = "spBART_stump",fold = j))
  #
  #     # Calculating the CRPS as well
  #     comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
  #                                                               value = crps(y = train$y ,
  #                                                                            means = colMeans(spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                            sds = rep(mean(spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(train$y)))$CRPS,
  #                                                               model = "spBART_stump",fold = j))
  #
  #     comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
  #                                                               value = crps(y = test$y ,
  #                                                                            means = colMeans(spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                            sds = rep(mean(spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(test$y)))$CRPS,
  #                                                               model = "spBART_stump",fold = j))
  #
  #     # Removing the model
  #     rm(spBART)
  #   }
  # Initializing the modelling for the BART model.
  bartmod <- dbarts::bart(x.train = x_train,y.train = y_train,x.test = x_test)

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = bartmod$yhat.train.mean,
                                                                         y = train$y),
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = bartmod$yhat.test.mean,
                                                                         y = test$y),
                                                            model = "BART",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = bartmod$yhat.train.mean,
                                                                         sds = rep(mean(bartmod$sigma), length(train$y) ))$CRPS,
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = bartmod$yhat.test.mean,
                                                                         sds = rep(mean(bartmod$sigma), length(test$y) ))$CRPS,
                                                            model = "BART",fold = j))
  rm(bartmod)


  # Since SOFTBART and MOTR-BART dont do it properly for 1-d
  if(ncol(x_train)>1){
      # Doing for SoftBART
      softbartmod <- SoftBart::softbart(X = x_train,Y = y_train,X_test =  x_test)

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(x = softbartmod$y_hat_train_mean,
                                                                             y = train$y),
                                                                model = "softBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                                value = rmse(x = softbartmod$y_hat_test_mean,
                                                                             y = test$y),
                                                                model = "softBART",fold = j))

      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = train$y ,
                                                                             means = softbartmod$y_hat_train_mean,
                                                                             sds = rep(mean(softbartmod$sigma), length(train$y) ))$CRPS,
                                                                model = "softBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                                value = crps(y = test$y ,
                                                                             means = softbartmod$y_hat_test_mean,
                                                                             sds = rep(mean(softbartmod$sigma), length(test$y) ))$CRPS,
                                                                model = "softBART",fold = j))

      rm(softbartmod)

      # Doing for MOTR-BART
      motrbartmod <- motr_bart(x = x_train,y = y_train,ancestors = TRUE,ntrees = ntree_)
      if(ntree_!=1){
      motrbart_pred <- predict_motr_bart(object = motrbartmod,newdata = x_test,type = "all")

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(x = colMeans(motrbartmod$y_hat),
                                                                             y = train$y),
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                                value = rmse(x = colMeans(motrbart_pred),
                                                                             y = test$y),
                                                                model = "motrBART",fold = j))

      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = train$y ,
                                                                             means = colMeans(motrbartmod$y_hat),
                                                                             sds = rep(mean(sqrt(motrbartmod$sigma2)), length(train$y) ))$CRPS,
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                                value = crps(y = test$y ,
                                                                             means = colMeans(motrbart_pred),
                                                                             sds = rep(mean(sqrt(motrbartmod$sigma2)), length(test$y) ))$CRPS,
                                                                model = "motrBART",fold = j))
      }
      rm(motrbartmod)

  } else {

    # Doing for SoftBART
    x_train_new <- cbind(x_train,x_train)
    x_test_new <- cbind(x_test,x_test)
    colnames(x_train_new) <- colnames(x_test_new) <- c("x.1","x.2")

    softbartmod <- SoftBart::softbart(X = x_train_new,Y = y_train,X_test =  x_test_new)

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(x = softbartmod$y_hat_train_mean,
                                                                           y = train$y),
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(x = softbartmod$y_hat_test_mean,
                                                                           y = test$y),
                                                              model = "softBART",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = train$y ,
                                                                           means = softbartmod$y_hat_train_mean,
                                                                           sds = rep(mean(softbartmod$sigma), length(train$y) ))$CRPS,
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = test$y ,
                                                                           means = softbartmod$y_hat_test_mean,
                                                                           sds = rep(mean(softbartmod$sigma), length(test$y) ))$CRPS,
                                                              model = "softBART",fold = j))

    rm(softbartmod)

    # Doing for MOTR-BART
    motrbartmod <- motr_bart(x = x_train_new,y = y_train,ancestors = TRUE,ntrees = ntree_)

    if(ntree_!=1){
      motrbart_pred <- predict_motr_bart(object = motrbartmod,newdata = x_test_new,type = "all")

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(x = colMeans(motrbartmod$y_hat),
                                                                             y = train$y),
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                                value = rmse(x = colMeans(motrbart_pred),
                                                                             y = test$y),
                                                                model = "motrBART",fold = j))

      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = train$y ,
                                                                             means = colMeans(motrbartmod$y_hat),
                                                                             sds = rep(mean(sqrt(motrbartmod$sigma2)), length(train$y) ))$CRPS,
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                                value = crps(y = test$y ,
                                                                             means = colMeans(motrbart_pred),
                                                                             sds = rep(mean(sqrt(motrbartmod$sigma2)), length(test$y) ))$CRPS,
                                                                model = "motrBART",fold = j))
    }
    rm(motrbartmod)

  }

  return(comparison_metrics)

}

# All BART-lite interaction
all_bart_lite_interaction <- function(cv_element,
                          nIknots_,
                          ntree_,
                          seed_,
                          use_bs_,
                          j,
                          motr_bart_,
                          alpha_,
                          rsp_bart_all_,
                          stump_,
                          scale_init_,
                          update_tau_beta_,
                          dif_order_,
                          interaction_term_,
                          mle_prior_,
                          y_scale_,
                          n_mcmc_,
                          n_burn_){


  if(!interaction_term_){
    warning("Do not run this if you dont want interactions")
  }


  # Doing a warming for the case whichI don't have
  # if(ntree_<50){
  #   stop("Use the all_bart() function instead.")
  # }

  # To replicate the results
  set.seed(seed_)
  train <- cv_element$train
  test <- cv_element$test

  # Getting the training elements
  x_train <- train %>% dplyr::select(dplyr::starts_with("x"))
  x_test <- test %>% dplyr::select(dplyr::starts_with("x"))
  y_train <- train %>% dplyr::pull("y")
  y_test <- test %>% dplyr::pull("y")

  # Initialising df
  comparison_metrics <- data.frame(metric = NULL, value = NULL, model = NULL,fold = NULL)


  spBART_interaction <- rspBART(x_train = x_train,
                    x_test = x_test,y_train = y_train,
                    n_mcmc = n_mcmc_,node_min_size = 25,alpha = alpha_,
                    n_burn = 0,nIknots = nIknots_,n_tree = ntree_,
                    use_bs = use_bs_,all_var = rsp_bart_all_,
                    stump = FALSE,dif_order = dif_order_,scale_bool = y_scale_,
                    motrbart_bool = motr_bart_,
                    scale_init = scale_init_,
                    interaction_term = interaction_term_,main_effects_pred = FALSE,
                    update_tau_beta = update_tau_beta_,mle_prior = mle_prior_,linero_sampler = FALSE,
                    use_D_bool = FALSE,scale_basis_function = FALSE,store_tree_fit = FALSE,varimportance_bool = TRUE,
                    robust_prior = FALSE)


  n_burn_ <- n_burn_
  n_mcmc_ <- spBART_interaction$mcmc$n_mcmc

  # Calculating metrics for splinesBART
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = colMeans(spBART_interaction$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         y = train$y),
                                                            model = "spBART_inter",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = colMeans(spBART_interaction$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         y = test$y),
                                                            model = "spBART_inter",fold = j))

  # Calculating metrics for splinesBART
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                            value = mae(x = colMeans(spBART_interaction$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         y = train$y),
                                                            model = "spBART_inter",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                            value = mae(x = colMeans(spBART_interaction$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         y = test$y),
                                                            model = "spBART_inter",fold = j))


  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = colMeans(spBART_interaction$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         sds = rep(mean(spBART_interaction$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(train$y)))$CRPS,
                                                            model = "spBART_inter",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = colMeans(spBART_interaction$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                         sds = rep(mean(spBART_interaction$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(test$y)))$CRPS,
                                                            model = "spBART_inter",fold = j))


  rm(spBART_interaction)

  # spBART <- rspBART(x_train = x_train,
  #                   x_test = x_test,y_train = y_train,
  #                   n_mcmc = 2500,node_min_size = 15,alpha = alpha_,
  #                   n_burn = 0,nIknots = nIknots_,n_tree = ntree_,
  #                   use_bs = use_bs_,all_var = rsp_bart_all_,
  #                   stump = FALSE,dif_order = dif_order_,
  #                   motrbart_bool = motr_bart_,
  #                   scale_init = scale_init_,
  #                   interaction_term = FALSE,
  #                   interaction_list = NULL,
  #                   update_tau_beta = update_tau_beta_)
  #
  #
  # n_burn_ <- 500
  # n_mcmc_ <- spBART$mcmc$n_mcmc
  #
  # # Calculating metrics for splinesBART
  # comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
  #                                                           value = rmse(x = colMeans(spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                        y = train$y),
  #                                                           model = "spBART_all",fold = j))
  #
  # comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
  #                                                           value = rmse(x = colMeans(spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                        y = test$y),
  #                                                           model = "spBART_all",fold = j))
  #
  # # Calculating the CRPS as well
  # comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
  #                                                           value = crps(y = train$y ,
  #                                                                        means = colMeans(spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                        sds = rep(mean(spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(train$y)))$CRPS,
  #                                                           model = "spBART_all",fold = j))
  #
  # comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
  #                                                           value = crps(y = test$y ,
  #                                                                        means = colMeans(spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
  #                                                                        sds = rep(mean(spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(test$y)))$CRPS,
  #                                                           model = "spBART_all",fold = j))

  #  Removing the model
  # rm(spBART)

  # Initializing the modelling for the BART model.
  bartmod <- dbarts::bart(x.train = x_train,y.train = y_train,x.test = x_test)

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = bartmod$yhat.train.mean,
                                                                         y = train$y),
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = bartmod$yhat.test.mean,
                                                                         y = test$y),
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                            value = mae(x = bartmod$yhat.train.mean,
                                                                         y = train$y),
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                            value = mae(x = bartmod$yhat.test.mean,
                                                                         y = test$y),
                                                            model = "BART",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = bartmod$yhat.train.mean,
                                                                         sds = rep(mean(bartmod$sigma), length(train$y) ))$CRPS,
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = bartmod$yhat.test.mean,
                                                                         sds = rep(mean(bartmod$sigma), length(test$y) ))$CRPS,
                                                            model = "BART",fold = j))
  rm(bartmod)


  # Since SOFTBART and MOTR-BART dont do it properly for 1-d
  if(ncol(x_train)>1){
    # Doing for SoftBART
    softbartmod <- SoftBart::softbart(X = x_train,Y = y_train,X_test =  x_test)

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(x = softbartmod$y_hat_train_mean,
                                                                           y = train$y),
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(x = softbartmod$y_hat_test_mean,
                                                                           y = test$y),
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                              value = mae(x = softbartmod$y_hat_train_mean,
                                                                           y = train$y),
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                              value = mae(x = softbartmod$y_hat_test_mean,
                                                                           y = test$y),
                                                              model = "softBART",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = train$y ,
                                                                           means = softbartmod$y_hat_train_mean,
                                                                           sds = rep(mean(softbartmod$sigma), length(train$y) ))$CRPS,
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = test$y ,
                                                                           means = softbartmod$y_hat_test_mean,
                                                                           sds = rep(mean(softbartmod$sigma), length(test$y) ))$CRPS,
                                                              model = "softBART",fold = j))

    rm(softbartmod)

    # Doing for MOTR-BART
    motrbartmod <- motr_bart(x = x_train,y = y_train,ancestors = TRUE,ntrees = ntree_)
    if(ntree_!=1){
      motrbart_pred <- predict_motr_bart(object = motrbartmod,newdata = x_test,type = "all")

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(x = colMeans(motrbartmod$y_hat),
                                                                             y = train$y),
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                                value = rmse(x = colMeans(motrbart_pred),
                                                                             y = test$y),
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                                value = mae(x = colMeans(motrbartmod$y_hat),
                                                                             y = train$y),
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                                value = mae(x = colMeans(motrbart_pred),
                                                                             y = test$y),
                                                                model = "motrBART",fold = j))

      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = train$y ,
                                                                             means = colMeans(motrbartmod$y_hat),
                                                                             sds = rep(mean(sqrt(motrbartmod$sigma2)), length(train$y) ))$CRPS,
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                                value = crps(y = test$y ,
                                                                             means = colMeans(motrbart_pred),
                                                                             sds = rep(mean(sqrt(motrbartmod$sigma2)), length(test$y) ))$CRPS,
                                                                model = "motrBART",fold = j))
    }
    rm(motrbartmod)

  } else {

    # Doing for SoftBART
    x_train_new <- cbind(x_train,x_train)
    x_test_new <- cbind(x_test,x_test)
    colnames(x_train_new) <- colnames(x_test_new) <- c("x.1","x.2")

    softbartmod <- SoftBart::softbart(X = x_train_new,Y = y_train,X_test =  x_test_new)

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(x = softbartmod$y_hat_train_mean,
                                                                           y = train$y),
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(x = softbartmod$y_hat_test_mean,
                                                                           y = test$y),
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                              value = mae(x = softbartmod$y_hat_train_mean,
                                                                           y = train$y),
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                              value = mae(x = softbartmod$y_hat_test_mean,
                                                                           y = test$y),
                                                              model = "softBART",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = train$y ,
                                                                           means = softbartmod$y_hat_train_mean,
                                                                           sds = rep(mean(softbartmod$sigma), length(train$y) ))$CRPS,
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = test$y ,
                                                                           means = softbartmod$y_hat_test_mean,
                                                                           sds = rep(mean(softbartmod$sigma), length(test$y) ))$CRPS,
                                                              model = "softBART",fold = j))

    rm(softbartmod)

    # Doing for MOTR-BART
    motrbartmod <- motr_bart(x = x_train_new,y = y_train,ancestors = TRUE,ntrees = ntree_)

    if(ntree_!=1){
      motrbart_pred <- predict_motr_bart(object = motrbartmod,newdata = x_test_new,type = "all")

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(x = colMeans(motrbartmod$y_hat),
                                                                             y = train$y),
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                                value = rmse(x = colMeans(motrbart_pred),
                                                                             y = test$y),
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_train",
                                                                value = mae(x = colMeans(motrbartmod$y_hat),
                                                                             y = train$y),
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "mae_test",
                                                                value = mae(x = colMeans(motrbart_pred),
                                                                             y = test$y),
                                                                model = "motrBART",fold = j))

      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = train$y ,
                                                                             means = colMeans(motrbartmod$y_hat),
                                                                             sds = rep(mean(sqrt(motrbartmod$sigma2)), length(train$y) ))$CRPS,
                                                                model = "motrBART",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                                value = crps(y = test$y ,
                                                                             means = colMeans(motrbart_pred),
                                                                             sds = rep(mean(sqrt(motrbartmod$sigma2)), length(test$y) ))$CRPS,
                                                                model = "motrBART",fold = j))
    }
    rm(motrbartmod)

  }

  return(comparison_metrics)

}



# Summarising all the metrics and results
wrapping_comparison <- function(result_){

  # Initialising df
  comparison_metrics <- data.frame(metric = NULL, value = NULL, model = NULL,fold = NULL)

  for(j in 1:length(result_)){


    n_burn_ <- 500
    n_mcmc_ <- result_[[j]]$spBART$mcmc$n_mcmc

    if(!is.null(result_[[j]]$spBART_all)){
      # Calculating metrics for splinesBART
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(x = colMeans(result_[[j]]$spBART_all$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                             y = result_[[j]]$cv$train$y),
                                                                model = "spBART_all",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                                value = rmse(x = colMeans(result_[[j]]$spBART_all$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                             y = result_[[j]]$cv$test$y),
                                                                model = "spBART_all",fold = j))

      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = result_[[j]]$cv$train$y ,
                                                                             means = colMeans(result_[[j]]$spBART_all$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                             sds = rep(mean(result_[[j]]$spBART_all$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(result_[[j]]$cv$train$y)))$CRPS,
                                                                model = "spBART_all",fold = j))

      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                                value = crps(y = result_[[j]]$cv$test$y ,
                                                                             means = colMeans(result_[[j]]$spBART_all$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                             sds = rep(mean(result_[[j]]$spBART_all$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(result_[[j]]$cv$test$y)))$CRPS,
                                                                model = "spBART_all",fold = j))
    }

    # Calculating metrics for splinesBART
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(x = colMeans(result_[[j]]$spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           y = result_[[j]]$cv$train$y),
                                                              model = "spBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(x = colMeans(result_[[j]]$spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           y = result_[[j]]$cv$test$y),
                                                              model = "spBART",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = result_[[j]]$cv$train$y ,
                                                                           means = colMeans(result_[[j]]$spBART$y_train_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           sds = rep(mean(result_[[j]]$spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(result_[[j]]$cv$train$y)))$CRPS,
                                                              model = "spBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = result_[[j]]$cv$test$y ,
                                                                           means = colMeans(result_[[j]]$spBART$y_test_hat[(n_burn_+1):n_mcmc_,,drop = FALSE]),
                                                                           sds = rep(mean(result_[[j]]$spBART$all_tau[(n_burn_+1):n_mcmc_])^(-1/2), length(result_[[j]]$cv$test$y)))$CRPS,
                                                              model = "spBART",fold = j))

    # ============================
    # Calculating metrics for BART
    # ============================

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(x = result_[[j]]$bartmod$yhat.train.mean,
                                                                           y = result_[[j]]$cv$train$y),
                                                              model = "BART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(x = result_[[j]]$bartmod$yhat.test.mean,
                                                                           y = result_[[j]]$cv$test$y),
                                                              model = "BART",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = result_[[j]]$cv$train$y ,
                                                                           means = result_[[j]]$bartmod$yhat.train.mean,
                                                                           sds = rep(mean(result_[[j]]$bartmod$sigma), length(result_[[j]]$cv$train$y) ))$CRPS,
                                                              model = "BART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = result_[[j]]$cv$test$y ,
                                                                           means = result_[[j]]$bartmod$yhat.test.mean,
                                                                           sds = rep(mean(result_[[j]]$bartmod$sigma), length(result_[[j]]$cv$test$y) ))$CRPS,
                                                              model = "BART",fold = j))


    # ============================
    # Calculating metrics for softBART
    # ============================

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                              value = rmse(x = result_[[j]]$softbartmod$y_hat_train_mean,
                                                                           y = result_[[j]]$cv$train$y),
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                              value = rmse(x = result_[[j]]$softbartmod$y_hat_test_mean,
                                                                           y = result_[[j]]$cv$test$y),
                                                              model = "softBART",fold = j))

    # Calculating the CRPS as well
    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                              value = crps(y = result_[[j]]$cv$train$y ,
                                                                           means = result_[[j]]$softbartmod$y_hat_train_mean,
                                                                           sds = rep(mean(result_[[j]]$softbartmod$sigma), length(result_[[j]]$cv$train$y) ))$CRPS,
                                                              model = "softBART",fold = j))

    comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                              value = crps(y = result_[[j]]$cv$test$y ,
                                                                           means = result_[[j]]$softbart$y_hat_test_mean,
                                                                           sds = rep(mean(result_[[j]]$softbartmod$sigma), length(result_[[j]]$cv$test$y) ))$CRPS,
                                                              model = "softBART",fold = j))
    # ============================
    # Calculating metrics for MOTRBART
    # ============================


    if(result_[[j]]$spBART$prior$n_tree>1){
        comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                  value = rmse(x = colMeans(result_[[j]]$motrbartmod$y_hat),
                                                                               y = result_[[j]]$cv$train$y),
                                                                  model = "motrBART",fold = j))

        comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                                  value = rmse(x = colMeans(result_[[j]]$motrbart_pred),
                                                                               y = result_[[j]]$cv$test$y),
                                                                  model = "motrBART",fold = j))

        # Calculating the CRPS as well
        comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                  value = crps(y = result_[[j]]$cv$train$y ,
                                                                               means = colMeans(result_[[j]]$motrbartmod$y_hat),
                                                                               sds = rep(mean(sqrt(result_[[j]]$motrbartmod$sigma2)), length(result_[[j]]$cv$train$y) ))$CRPS,
                                                                  model = "motrBART",fold = j))

        comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                                  value = crps(y = result_[[j]]$cv$test$y ,
                                                                               means = colMeans(result_[[j]]$motrbart_pred),
                                                                               sds = rep(mean(sqrt(result_[[j]]$motrbartmod$sigma2)), length(result_[[j]]$cv$test$y) ))$CRPS,
                                                                  model = "motrBART",fold = j))
    } else {
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                                value = rmse(x = colMeans(result_[[j]]$motrbartmod$y_hat),
                                                                             y = result_[[j]]$cv$train$y),
                                                                model = "motrBART",fold = j))


      # Calculating the CRPS as well
      comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                                value = crps(y = result_[[j]]$cv$train$y ,
                                                                             means = colMeans(result_[[j]]$motrbartmod$y_hat),
                                                                             sds = rep(mean(sqrt(result_[[j]]$motrbartmod$sigma2)), length(result_[[j]]$cv$train$y) ))$CRPS,
                                                                model = "motrBART",fold = j))

  }

  }

  return(comparison_metrics)

}


# Getting a model to evaluate variable importance
var_importance_counter <- function(result_,rep_, only_sp_ = FALSE){

  if(!only_sp_){
    # Getting a counter for times that a variable is used in within a tree
    p_counter <- numeric(ncol(result_[[rep_]]$cv$train)-1)
    spBART <- result_[[rep_]]$spBART
  } else {
    spBART <- result_
    p_counter <- numeric(NCOL(spBART$data$x_train))
  }
  for(i in 501:spBART$mcmc$n_mcmc){

    for(t in 1:spBART$prior$n_tree){
      curr_tree <- spBART$mcmc$all_trees[[i]][[t]]
      terminals <- get_terminals(curr_tree)
      for(ell in 1:length(terminals)){
        p_counter[unique(curr_tree[[terminals[ell]]]$ancestors)] <- p_counter[unique(curr_tree[[terminals[ell]]]$ancestors)] + 1
      }
    }
  }

  return(round(p_counter/sum(p_counter),digits = 5))

}

# Getting a model to evaluate variable importance
tree_length_counter <- function(result_,rep_, only_sp = FALSE){

  # Getting a counter for times that a variable is used in within a tree
  if(!only_sp){
    spBART <- result_[[rep_]]$spBART
  } else {
    spBART <- result_
  }
  matrix_tree <- matrix(0,ncol = spBART$prior$n_tree, nrow = 2000)
  curr <- 0
  for(i in 501:spBART$mcmc$n_mcmc){
    curr <- curr + 1
    for(t in 1:spBART$prior$n_tree){
        curr_tree <- spBART$mcmc$all_trees[[i]][[t]]
        terminals <- get_terminals(curr_tree)
        matrix_tree[curr,t] <- length(terminals)
      }
  }

  return(matrix_tree)

}

# =======
# k_fold
# =======
kfold <- function(data_,
                  nfold_ = 10,
                  seed_){

  # Always set a seed
  if(is.null(seed_)){
    stop("Insert a valid seed.")
  }

  # Setting a seed
  set.seed(seed_)

  # Getting the cv
  cv_object_ <- vector("list",length = nfold_)

  # Getting the indentifier
  fold_index <- dismo::kfold(x = data_,k = nfold_)

  for(i in 1:nfold_){

    # Creating a list to store all information
    data_train <- data_[fold_index!=i,]
    data_test <- data_[fold_index==i,]
    x_train <- data_train[,which(colnames(data_train)!="y"), drop = FALSE]
    x_test <- data_test[,which(colnames(data_train)!="y"), drop = FALSE]
    y_train <- data_train[["y"]]
    y_test <- data_test[["y"]]

    # List witha all elements
    obj_ <- list(data_train = data_train,
                 data_test = data_test,
                 x_train = x_train,
                 x_test = x_test,
                 y_train = y_train,
                 y_test = y_test,
                 seed = seed_)

    cv_object_[[i]] <- obj_

  }

  # Return the list with all cv_objects
  return(cv_object_)

}

# spBART competitors
competitors_comparison_ <- function(cv_object_fold_,
                                    fold_,
                                    seed_,
                                    interaction_term_,
                                    interaction_list_,
                                    return_models = TRUE){

  # Setting a seed
  set.seed(cv_object_fold_$seed)

  # Initialising df
  comparison_metrics <- data.frame(metric = NULL,
                                   value = NULL,
                                   model = NULL,
                                   fold = NULL)


  # Running all the models
  x_train <- cv_object_fold_$x_train
  x_test <- cv_object_fold_$x_test
  y_train <- cv_object_fold_$y_train
  y_test <- cv_object_fold_$y_test
  train <- cv_object_fold_$data_train
  test <- cv_object_fold_$data_test
  j <- fold_

  # Loading packages
  library(dbarts)
  library(SoftBart)
  library(earth)
  library(mgcv)
  library(MOTRbart)
  library(BASS)

  # Running all models
  bartmod <- dbarts::bart(x.train = x_train,y.train = y_train,x.test = x_test)
  softbartmod <- SoftBart::softbart(X = x_train,Y = y_train,X_test = x_test)
  motrbartmod <- MOTRbart::motr_bart(x = x_train,y = y_train)
  motrbart_pred <- predict_motr_bart(motrbartmod,x_test,type = "mean")
  bassmod <- BASS::bass(xx = x_train,y = y_train)
  bass_pred <- predict(bassmod,x_test)

  if(interaction_term_){

    # MARS
    marsmod <- earth(y~., data = train,degree = length(interaction_list_))
    marsmod_pred <- predict(marsmod,x_test)

    # Fitting mgcv
    main_effects_formula <- paste0(paste0("s(",colnames(x_train),", bs = 'cr')"),collapse = "+")
    interactions <- character(0)
    for(i in 1:length(interaction_list_)){
      interactions[i] <- paste0("te(",paste0(paste0("",colnames(x_train)[interaction_list_[[i]]],""),collapse = ","),", bs = 'cr')")

    }

    interactions <- paste0(interactions,collapse = " + ")

    # Getting the final formula
    mgcv_formula <- formula(paste0("y~",paste0(c(main_effects_formula,interactions),collapse = "+")))

    gammod <- gam(mgcv_formula,data = train)
    gammod_pred <- predict(gammod,x_test)

  } else {

    marsmod <- earth(y~., data = train,degree = 1)
    marsmod_pred <- predict(marsmod,x_test)

    main_effects_formula <- paste0(paste0("s(",colnames(x_train),", bs = 'cr')"),collapse = "+")
    mgcv_formula <- formula(paste0("y~",paste0(c(main_effects_formula),collapse = "+")))

    gammod <- gam(mgcv_formula,data = train)
    gammod_pred <- predict(gammod,x_test)

  }

  # ==================================
  # Calculating all metrics from now
  # ==================================

  # ============================
  # Calculating metrics for BART
  # ============================

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = bartmod$yhat.train.mean,
                                                                         y = train$y),
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = bartmod$yhat.test.mean,
                                                                         y = test$y),
                                                            model = "BART",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = bartmod$yhat.train.mean,
                                                                         sds = rep(mean(bartmod$sigma), length(train$y) ))$CRPS,
                                                            model = "BART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = bartmod$yhat.test.mean,
                                                                         sds = rep(mean(bartmod$sigma), length(test$y) ))$CRPS,
                                                            model = "BART",fold = j))


  # ============================
  # Calculating metrics for softBART
  # ============================

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = softbartmod$y_hat_train_mean,
                                                                         y = train$y),
                                                            model = "softBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = softbartmod$y_hat_test_mean,
                                                                         y = test$y),
                                                            model = "softBART",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = softbartmod$y_hat_train_mean,
                                                                         sds = rep(mean(softbartmod$sigma), length(train$y) ))$CRPS,
                                                            model = "softBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = softbartmod$y_hat_test_mean,
                                                                         sds = rep(mean(softbartmod$sigma), length(test$y) ))$CRPS,
                                                            model = "softBART",fold = j))
  # ============================
  # Calculating metrics for MOTRBART
  # ============================
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = colMeans(motrbartmod$y_hat),
                                                                         y = train$y),
                                                            model = "motrBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = colMeans(motrbart_pred),
                                                                         y = test$y),
                                                            model = "motrBART",fold = j))

    # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = colMeans(motrbartmod$y_hat),
                                                                         sds = rep(mean(sqrt(motrbartmod$sigma2)), length(train$y) ))$CRPS,
                                                            model = "motrBART",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = colMeans(motrbart_pred),
                                                                         sds = rep(mean(sqrt(motrbartmod$sigma2)), length(test$y) ))$CRPS,
                                                            model = "motrBART",fold = j))

  # ============================
  # Calculating metrics for MARS
  # ============================
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = c(marsmod$fitted.values),
                                                                         y = train$y),
                                                            model = "MARS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = c(marsmod_pred),
                                                                         y = test$y),
                                                            model = "MARS",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = c(marsmod$fitted.values),
                                                                         sds = rep(sigma(marsmod), length(train$y) ))$CRPS,
                                                            model = "MARS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = c(marsmod_pred),
                                                                         sds = rep(sigma(marsmod), length(test$y) ))$CRPS,
                                                            model = "MARS",fold = j))

  # ============================
  # Calculating metrics for GAM
  # ============================
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = c(gammod$fitted.values),
                                                                         y = train$y),
                                                            model = "GAM",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = c(gammod_pred),
                                                                         y = test$y),
                                                            model = "GAM",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = c(gammod$fitted.values),
                                                                         sds = rep(sigma(marsmod), length(train$y) ))$CRPS,
                                                            model = "GAM",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = c(gammod_pred),
                                                                         sds = rep(sigma(marsmod), length(test$y) ))$CRPS,
                                                            model = "GAM",fold = j))


  # =======================================
  #     Doing a comparison with BASS
  # =======================================

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_train",
                                                            value = rmse(x = c(bassmod$yhat.mean),
                                                                         y = train$y),
                                                            model = "BASS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "rmse_test",
                                                            value = rmse(x = colMeans(bass_pred),
                                                                         y = test$y),
                                                            model = "BASS",fold = j))

  # Calculating the CRPS as well
  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_train",
                                                            value = crps(y = train$y ,
                                                                         means = c(bassmod$yhat.mean),
                                                                         sds = rep(mean(bassmod$s2), length(train$y) ))$CRPS,
                                                            model = "BASS",fold = j))

  comparison_metrics <- rbind(comparison_metrics,data.frame(metric = "crps_test",
                                                            value = crps(y = test$y ,
                                                                         means = colMeans(bass_pred),
                                                                         sds = rep(mean(bassmod$s2), length(test$y) ))$CRPS,
                                                            model = "BASS",fold = j))

  # Returning a list with all elements
  if(return_models){
    return(list(comparison_metrics = comparison_metrics,
                bartmod = bartmod,
                softbartmod = softbartmod,
                motrbartmod = motrbartmod,
                motrbart_pred = motrbart_pred,
                marsmod = marsmod,
                marsmod_pred = marsmod_pred,
                gammod = gammod,
                gammod_pred = gammod_pred,
                bassmod = bassmod,
                bass_pred = bass_pred))
  } else {
    return(list(comparison_metrics = comparison_metrics))
  }

}
