# source("R/debugging_rspBART.R")
# rm(list=ls())
# source("R/other_functions.R")
# source("R/sim_functions.R")
# source("R/tree_functions.R")
# source("inst/debugging_rspBART.R")
# devtools::load_all()
set.seed(42)

# Creating the main function from the rspBART
rspBART <- function(x_train,
                    y_train,
                    x_test,
                    n_tree = 10,
                    node_min_size = 15,
                    n_mcmc = 2000,
                    n_burn = 500,
                    alpha = 0.95,
                    beta = 2,
                    df = 3,
                    sigquant = 0.9,
                    kappa = 2,
                    # Splines parameters
                    nIknots = 3,
                    dif_order = 1,
                    tau = 100,
                    scale_bool = TRUE,
                    stump = FALSE,
                    numcut = 100L, # Defining the grid of split rules
                    usequants = FALSE,
                    motrbart_bool = FALSE,
                    use_bs = FALSE,
                    plot_preview = FALSE,
                    all_var = FALSE,
                    scale_init = FALSE,
                    update_tau_beta = FALSE,
                    main_effects_pred = FALSE,
                    interaction_term =  TRUE,
                    interaction_list = NULL,
                    store_tree_fit = FALSE,
                    mle_prior = FALSE,
                    linero_sampler = FALSE,
                    use_D_bool = TRUE,
                    varimportance_bool = FALSE,
                    scale_basis_function = TRUE,
                    robust_prior = FALSE,
                    eta = 1e-6,
                    a_delta = 0.5,
                    d_delta = 0.5,
                    pen_basis = FALSE,
                    centered_basis = FALSE
) {


  # # Another constrain
  # if(interaction_term & is.null(interaction_list)){
  #   stop("Define the interaction list.")
  # }

  # Vector to intialise the stumps
  stump_init_indicator <- rep(1:NCOL(x_train), each = 1)
  if(n_tree > length(stump_init_indicator)){
    stump_init_indicator <- sample(1:NCOL(x_train),size = n_tree, replace = TRUE)
    stump_init_indicator[1:NCOL(x_train)] <- 1:NCOL(x_train)
  }

  # if(length(stump_init_indicator)!=n_tree){
  #   n_tree <- length(stump_init_indicator)
  #   warning("initialise with the correct number of trees")
  # }


  if(robust_prior){
    stop("Do not use the robust prior approach for the penalised splines")
  }
  # Making interaction FALSE for univariate
  if(NCOL(x_train)==1){
    interaction_term <- FALSE
  }

  if(isFALSE(interaction_term) & !is.null(interaction_list)){
    stop("Is a model without the interaction but a interaction list is defined")
  }

  # Verifying if x_train and x_test are matrices
  if(!is.data.frame(x_train) || !is.data.frame(x_test)){
    stop("Insert valid data.frame for both data and xnew.")
  }


  # Getting the valid
  dummy_x <- base_dummyVars(x_train)

  # Create a list
  if(length(dummy_x$facVars)!=0){
    for(i in 1:length(dummy_x$facVars)){
      # See if the levels of the test and train matches
      if(!all(levels(x_train[[dummy_x$facVars[i]]])==levels(x_test[[dummy_x$facVars[i]]]))){
        levels(x_test[[dummy_x$facVars[[i]]]]) <- levels(x_train[[dummy_x$facVars[[i]]]])
      }
      df_aux <- data.frame( x = x_train[,dummy_x$facVars[i]],y)
      formula_aux <- stats::aggregate(y~x,df_aux,mean)
      formula_aux$y <- rank(formula_aux$y)
      x_train[[dummy_x$facVars[i]]] <- as.numeric(factor(x_train[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

      # Doing the same for the test set
      x_test[[dummy_x$facVars[i]]] <- as.numeric(factor(x_test[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

    }
  }


  # Getting the train and test set
  x_train_scale <- as.matrix(x_train)
  x_test_scale <- as.matrix(x_test)

  # Scaling x
  x_min <- apply(as.matrix(x_train_scale),2,min)
  x_max <- apply(as.matrix(x_train_scale),2,max)

  # Storing the original
  x_train_original <- x_train
  x_test_original <- x_test



  # Normalising all the columns
  for(i in 1:ncol(x_train)){
    x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
    x_test_scale[,i] <- normalize_covariates_bart(y = x_test_scale[,i],a = x_min[i], b = x_max[i])
  }



  # Creating the numcuts matrix of splitting rules
  xcut_m <- matrix(NA,nrow = numcut,ncol = ncol(x_train_scale))
  for(i in 1:ncol(x_train_scale)){

    if(nrow(x_train_scale)<numcut){
      xcut_m[,i] <- sort(x_train_scale[,i])
    } else {
      xcut_m[,i] <- seq(min(x_train_scale[,i]),
                        max(x_train_scale[,i]),
                        length.out = numcut+2)[-c(1,numcut+2)]
    }
  }

  # =========================================================================================================
  # Getting the Splines Basis functions
  # =========================================================================================================
  # Creating a list of basis functions
  B_train_obj <- B_test_obj <- vector("list",length = length(dummy_x$continuousVars) + NCOL(interaction_list))

  # Setting new parameters for the spline
  ndx <- nIknots+1
  ord_ <- 4
  degree_ <- 3
  x_min_sp <- apply(x_train_scale,2,min)
  x_max_sp <- apply(x_train_scale,2,max)
  dx <- (x_max_sp-x_min_sp)/ndx

  # Creating the knots based on the DALSM package
  if(isFALSE(centered_basis)){
      new_knots <- matrix()
      new_knots <- matrix(mapply(x_min_sp,x_max_sp,dx, FUN = function(MIN,MAX,DX){seq(from = MIN-(1*ord_-1)*DX, to = MAX+(1*ord_-1)*DX, by = DX)}), ncol = length(dummy_x$continuousVars)) # MIN and MAX are 0 and 1 respectively, because of the scale
      colnames(new_knots) <- dummy_x$continuousVars

      colnames(new_knots) <- paste0("x.",1:length(dummy_x$continuousVars))
  } else {

    # Creating the knots based on the DALSM:: package
    new_knots <- list()
    for(i in 1:length(dummy_x$continuousVars)){
      if(dif_order==0){
        dif_order = 1
      }
      new_knots[[i]] <- DALSM::qknots(x = x_train_scale[,i],equid.knots = TRUE,
                                      xmin = min(min(x_test_scale[,i]),min(x_train_scale[,i])),
                                      xmax = max(max(x_test_scale[,i]),max(x_train_scale[,i])),
                                      pen.order = dif_order,K = nIknots)
      # Simplification in case dif_order = 0
      if(dif_order==0){
        new_knots[[i]]$Pd <- diag(nrow = nrow(new_knots[[i]]$Pd))
        new_knots[[i]]$pen.order = 0
        dif_order = 0
      }
    }

    # Setting the names of the basis
    names(new_knots) <- paste0("x.",1:length(dummy_x$continuousVars))

  }


  # New_knots
  if(interaction_term){

    if(is.null(interaction_list)){
      interaction_list <- combn(1:NCOL(x_train),2,simplify = TRUE)
      interaction_size <- NCOL(interaction_list)
    } else {
      interaction_list <- interaction_list
      if(is.matrix(interaction_list)){
        interaction_size <- NCOL(interaction_list)
      } else {
        interaction_list_matrix <- matrix(NA, nrow = 2, ncol = length(interaction_list))
        for(iter_inter_list in 1:length(interaction_list)){
          interaction_list_matrix[,iter_inter_list] <- interaction_list[[iter_inter_list]]
        }
        interaction_size <- NCOL(interaction_list_matrix)
        interaction_list <- interaction_list_matrix
        warning("Double check the interaction list")
      }
    }
  } else {
    stop("The default is to use interaction.")
  }


  # Creating the natural B-spline for each predictor
  for(i in 1:length(dummy_x$continuousVars)){

    if(isFALSE(centered_basis)){
      train_B <- splines::spline.des(x = x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],
                          knots = new_knots[,dummy_x$continuousVars[i]],
                          ord = ord_,
                          derivs = 0*x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = TRUE)$design

      # Choosing or not scaling the basis
      train.B.ref <- train_B

      intercept_train <- 1
      train_basis_intercept <- cbind(intercept_train,train.B.ref)


      B_train_obj[[i]] <- train_basis_intercept

      # Doing the same for the test
      test_B <- splines::spline.des(x = x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                    knots = new_knots[,dummy_x$continuousVars[i]],
                                    ord = ord_,
                                    derivs = 0*x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = TRUE)$design

      test_basis_intercept <- cbind(intercept_train, test_B)
      # CHANGE THIS LATER
      B_test_obj[[i]] <- train_basis_intercept

    } else { ## Accounting for the case wiht centered-basis


          # Modify the basis only with respect to the main effects at the moment
          centered_basis_aux <- DALSM::centeredBasis.gen(x = x_train_scale[,dummy_x$continuousVars[i]],
                                                         knots = new_knots[[i]]$knots,
                                                         pen.order = dif_order)

          train.B.ref <- centered_basis_aux$B
          intercept_train <- 1
          train_basis_intercept <- cbind(intercept_train,train.B.ref)


          B_train_obj[[i]] <- train_basis_intercept


          # Doing the same for the test samples
          centered_basis_aux_test <-DALSM::centeredBasis.gen(x = x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                                             knots = new_knots[[i]]$knots,pen.order = dif_order)
          test.B.ref <- centered_basis_aux_test$B
          intercept_test <- 1
          test_basis_intercept <- cbind(intercept_test,test.B.ref)

          B_test_obj[[i]] <- test_basis_intercept

    }

  }




  # Visualzing the main effects
  if(plot_preview){
    par(mfrow = c(2,5))
    # Visualizing some basis functions
    for(visu_basis in 1:10){
      plot(x  = 0, xlim = c(0,1), ylim = c(-1,1),
           type= 'n', ylab= paste0("B.",visu_basis), main = paste0("B(",visu_basis,")"))
      for(basis_col in 1:NCOL(B_train_obj[[visu_basis]])){
        points(x_train_scale[,visu_basis], B_train_obj[[visu_basis]][,basis_col],
               col = basis_col, pch = 20, xlab = paste0("x.",visu_basis))
      }
    }
  }

  # Interaction matrix list to be used in the penalised
  interaction_matrix_list <- vector("list",length = interaction_size)

  # Recreating the interaction list in the same format as the combn() function
  if(is.list(interaction_list)){
    interaction_list_aux <- matrix(NA, nrow = 2, ncol = interaction_size)
    for(col_iter in 1:NCOL(interaction_list_aux)){
      interaction_list_aux[,col_iter] <- interaction_list[[col_iter]]
    }
    interaction_list <- interaction_list_aux
  }


  # Adding the interaction basis
  if(interaction_term){
    jj_ = length(dummy_x$continuousVars)
    for (jj in 1:NCOL(interaction_list)) {
      jj_ = jj_ +1
      B_train_obj[[jj_]] <- multiply_matrices_general(A = B_train_obj[[interaction_list[1,jj]]][,-1],B = B_train_obj[[interaction_list[2,jj]]][,-1])
      B_test_obj[[jj_]] <- multiply_matrices_general(A = B_test_obj[[interaction_list[1,jj]]][,-1],B = B_test_obj[[interaction_list[2,jj]]][,-1])
    }
  }



  # Renaming the columns of inthe basis subindex
  main_effects_names <- paste0(1:length(dummy_x$continuousVars))

  # Getting a string of all pairwise possibilities
  concatenate_columns <- function(matrix) {
    result <- apply(matrix, 2, function(col) paste(col, collapse = ""))
    return(result)
  }

  # Renaming only the interactions
  if(interaction_term){
    interaction_names <- concatenate_columns(interaction_list)
  }




  # Scaling the y
  min_y <- min(y_train)
  max_y <- max(y_train)

  # Getting the min and max for each column
  min_x <- apply(x_train_scale,2,min)
  max_x <- apply(x_train_scale, 2, max)

  # Scaling "y"
  if(scale_bool){
    y_scale <- normalize_bart(y = y_train,a = min_y,b = max_y)

    # New update
    # tau_mu <- 4*n_tree*(kappa^2)
    tau_mu <- (n_tree/(NCOL(x_train)))
  } else {
    y_scale <- y_train

    # New parameter update
    # tau_mu <- (4*n_tree*(kappa^2))/((max_y-min_y)^2)
    tau_mu <- (n_tree/(NCOL(x_train)))


  }


  # Getting the naive sigma value
  nsigma <- naive_sigma(x = x_train_scale,y = y_scale)

  # Calculating tau hyperparam
  a_tau <- 3/2
  df_tau <- 3
  # Calculating lambda
  qchi <- stats::qchisq(p = 1-sigquant,df = df_tau,lower.tail = 1,ncp = 0)
  lambda <- (nsigma*nsigma*qchi)/df_tau
  d_tau <- (lambda*df_tau)/2


  # Simplying a_tau and d_tau priors
  # a_tau <- d_tau <- 1e-6

  # Call the bart function
  tau_init <- nsigma^(-2)

  mu_init <- mean(y_scale)

  # Creating the vector that stores all trees
  all_tree_post <- vector("list",length = round(n_mcmc-n_burn))


  # =====================================================================
  # ========= From here I gonna initialise the BART function itself =====
  # =====================================================================
  n_post <- (n_mcmc-n_burn)
  all_trees <- vector("list", n_mcmc)
  all_betas <- vector("list",n_mcmc)
  if(interaction_term){
    tau_beta <- rep(tau_mu,length(dummy_x$continuousVars)+NCOL(interaction_list))
    tau_beta <- rep(1,length(dummy_x$continuousVars)+NCOL(interaction_list))

    # tau_beta <- rep(tau_mu, NCOL(x_train_scale))
  } else {
    tau_beta <- rep(tau_mu,length(dummy_x$continuousVars))
  }

  # all_df <- cbind(x_train,y_train)
  # mod <- mgcv::gam(y_train ~ bs(x.10), data = all_df)
  # mod$coefficients

  # In this first scenario we are going to work with a single value of \tau
  if(interaction_term){
    all_tau_beta <- matrix(NA, nrow = (n_mcmc), ncol = NCOL(x_train_scale)+NCOL(interaction_list))
  } else {
    all_tau_beta <- matrix(NA, nrow = (n_mcmc), ncol = NCOL(x_train_scale))
  }
  # all_delta <- numeric(n_mcmc)
  all_tau <- numeric(n_mcmc)
  all_delta <- matrix(NA, nrow = n_mcmc, ncol = NCOL(x_train_scale)+NCOL(interaction_list))

  all_y_hat <- matrix(NA,nrow = n_mcmc,ncol = nrow(x_train_scale))
  all_y_hat_test <- matrix(NA, nrow = n_mcmc, ncol = nrow(x_test_scale))
  all_trees_fit <- vector("list",n_mcmc)
  all_trees <- vector("list",n_mcmc)
  forest <- vector("list",n_tree)

  # Partial component pieces
  # (OPTIONAL)
  if(store_tree_fit){
    partial_train_fits <- vector("list", n_tree)
  }

  proposal_outcomes <- setNames(data.frame(matrix(nrow = 0, ncol = 6)),
                                c("tree_number" , "proposal", "status","mcmc_iter", "new_tree_loglike", "old_tree_loglike"))
  all_train_indexes <- data.frame(matrix(data = NA,nrow = nrow(xcut_m),ncol = ncol(xcut_m)))

  # Creating a list to store all the main effects for the sum of trees
  if(main_effects_pred){

    # Calculating the main effects
    if(interaction_term){

      # Interacting
      main_effects_train_list_norm <- main_effects_test_list_norm <- main_effects_train_list <- main_effects_test_list <- vector("list", length = length(dummy_x$continuousVars)+NCOL(interaction_list))

      for(list_size in 1:length(main_effects_train_list)){
        main_effects_train_list[[list_size]] <- main_effects_train_list_norm[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_train))
        main_effects_test_list[[list_size]] <- main_effects_test_list_norm[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_test))
      }

      # Renaming the main list
      if(interaction_term){
        internames_ <- apply(interaction_list,2,function(vars_){paste0("x",paste0(vars_,collapse = ":"))})
        names(main_effects_train_list_norm) <- names(main_effects_test_list_norm) <- names(main_effects_train_list) <- names(main_effects_test_list) <- c(dummy_x$continuousVars, internames_)
      } else {
        names(main_effects_train_list_norm) <- names(main_effects_test_list_norm) <- names(main_effects_train_list) <- names(main_effects_test_list) <-dummy_x$continuousVars
      }


    } else { # In case there are no interactions
      main_effects_train_list_norm <- main_effects_test_list_norm <- main_effects_train_list <- main_effects_test_list <- vector("list", length = length(dummy_x$continuousVars))

      for(list_size in 1:length(dummy_x$continuousVars)){
        main_effects_train_list[[list_size]] <- main_effects_train_list_norm[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_train))
        main_effects_test_list[[list_size]] <- main_effects_train_list_norm[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_test))
      }

      names(main_effects_train_list_norm) <- names(main_effects_test_list_norm) <-  names(main_effects_train_list) <- names(main_effects_test_list) <- dummy_x$continuousVars
    }
  } else {
    main_effects_train_list_norm <- main_effects_test_list_norm <- main_effects_train_list <- main_effects_test_list <- NULL
  }

  # Gonna create a list of lists to store all the indexes for all split rules and cutpoints
  all_var_splits <- vector("list",ncol(x_train_scale))
  names(all_var_splits) <- colnames(x_train_scale)



  # Iterating over all possible x.columns
  for(i in 1:length(all_var_splits)){

    # Creating the dummy for a list of index to store all numeric split values
    all_cut_points <- vector("list", nrow(xcut_m))


    for(j in 1:length(all_cut_points)){

      # Getting the node indexes object
      left_train_list <- vector("list",length = 1L)
      names(left_train_list) <- "left_train"
      right_train_list <- vector("list",length = 1L)
      names(right_train_list) <- "right_train"
      left_test_list <- vector("list",length = 1L)
      names(left_test_list) <- "left_test"
      right_test_list <- vector("list",length = 1L)
      names(right_test_list) <- "right_test"

      node_index <- append(left_train_list, right_train_list) |>
        append(left_test_list) |> append(right_test_list)

      all_cut_points[[j]]$left_train <-  which(x_train_scale[,i] < xcut_m[j,i])
      all_cut_points[[j]]$right_train <-  which(x_train_scale[,i] >= xcut_m[j,i])
      all_cut_points[[j]]$left_test <-  which(x_test_scale[,i] < xcut_m[j,i])
      all_cut_points[[j]]$right_test <-  which(x_test_scale[,i] >= xcut_m[j,i])

    }

    all_var_splits[[i]] <- all_cut_points

  }

  # Creating the penalty matrix
  if(dif_order==0){
    P <- diag(nrow = NCOL(B_train_obj[[1]])-1)
    P_train_main <- diag(nrow = NCOL(B_train_obj[[1]]))
  } else {
    if(centered_basis){
      D <- centered_basis_aux$Dd
      P <- crossprod(centered_basis_aux$Dd)
      if(dif_order==1) {
        P[1,1] <- P[1,1] + eta
      } else if(dif_order==2){
        P[1,1] <- P[1,1] + eta
        P[2,2] <- P[2,2] + eta
      }
    } else {
      D <- (diff(diag(NCOL(B_train_obj[[1]])-1),differences = dif_order))
      P <- crossprod(D)
    }
    auxP <- as.matrix(Matrix::bdiag(1,P))
    P_train_main <- as.matrix(Matrix::nearPD(auxP)$mat)
  }

  # Generating the interaction matrix
  if(dif_order==0){
    P_train_interaction <- kronecker(P,P)
  } else {
    P_train_interaction <- as.matrix(Matrix::nearPD(kronecker(P,P))$mat)
  }

  # Tau_beta_priors
  # Getting hyperparameters for \tau_beta_j
  tau_beta_df <- 3
  a_tau_beta_j <- tau_beta_df/2

  # Setting the new value considering the paper
  d_tau_beta_j <- 0.01

  par(mfrow=c(1,1))
  # Visualising the prior
  if(plot_preview){
    rgamma(n = 1000,shape = a_tau_beta_j,rate = d_tau_beta_j) |> density() |> plot(main = "density prior \tau_beta_j")
    mean(rgamma(n = 1000,shape = a_tau_beta_j,rate = d_tau_beta_j))
  }

  # Getting the inverse of the P matrix
  P_inv <- chol2inv(chol(P_train_main))
  # P_inv_two <- chol2inv(chol(P_train_main_two))
  P_inv_interaction <- chol2inv(chol(P_train_interaction))


  # For the case of penalised basis I transform back the P into a diagonal matrix
  if(pen_basis){
    if(isFALSE(centered_basis)){
      stop("Don't use this approach at the moment")
    }

    # Getting the D matrix
    D <- centered_basis_aux$Dd
    Diff_term <- crossprod(D,solve(tcrossprod(D)))
    for(main_iter in 1:length(dummy_x$continuousVars)){
        B.train_aux <- B_train_obj[[main_iter]][,-1]%*%Diff_term
        B_train_obj[[main_iter]] <- cbind(1,B.train_aux)
        B.test_aux <- B_test_obj[[main_iter]][,-1]%*%Diff_term
        B_test_obj[[main_iter]] <- cbind(1,B.test_aux)
    }

    P_train_main <- P_inv <- diag(nrow = NCOL(B_train_obj[[1]]))
    P_train_interaction <- P_inv_interaction <- diag(nrow = (NCOL(B_train_obj[[1]])-1)^2)

    # Adding the interaction basis
    if(interaction_term){
      jj_ = length(dummy_x$continuousVars)
      for (jj in 1:NCOL(interaction_list)) {
        jj_ = jj_ +1
        B_train_obj[[jj_]] <- multiply_matrices_general(A = B_train_obj[[interaction_list[1,jj]]][,-1],B = B_train_obj[[interaction_list[2,jj]]][,-1])
        B_test_obj[[jj_]] <- multiply_matrices_general(A = B_test_obj[[interaction_list[1,jj]]][,-1],B = B_test_obj[[interaction_list[2,jj]]][,-1])
      }
    }

  }


  # ======== Defining index for each coefficient ============== #

  # Checking all possible interactions
  basis_size <-  NCOL(B_train_obj[[1]])     # Change this value to the desired size of each sublist

  # Setting the indexes for the betas
  D_seq <- 1:(basis_size*length(dummy_x$continuousVars))  # Replace this with the columns of D

  # Creating a vector
  basis_inter_size <- NCOL(P_train_interaction)
  basis_subindex_main <- split(D_seq, rep(1:((length(D_seq) %/% basis_size)), each = basis_size, length.out = length(D_seq)))


  D_int_seq <- (length(D_seq)+1):(basis_size*length(dummy_x$continuousVars)+(basis_inter_size)*interaction_size)
  interactions_subindex <- split(D_int_seq, rep(1:((length(D_int_seq) %/% basis_inter_size)), each = basis_inter_size, length.out = length(D_int_seq)))

  # Getting the final basis subindex
  if(interaction_term){
    basis_subindex <- append(basis_subindex_main,interactions_subindex)
  } else {
    basis_subindex <- basis_subindex_main
  }

  # Renaming only the interactions
  if(interaction_term){
    names(basis_subindex)[(length(dummy_x$continuousVars)+1):length(basis_subindex)] <- concatenate_columns(interaction_list)
    names(B_train_obj) <- names(B_test_obj) <- c(main_effects_names,interaction_names)
  } else {
    names(B_train_obj) <- names(B_test_obj) <- main_effects_names
  }

  # Geneating initial values for delta
  robust_delta <- rep(3, length(basis_subindex))

  # P_int_inv <-
  #most of the functions
  data <- list(x_train = x_train_scale,
               x_test = x_test_scale,
               y_train = y_scale,
               xcut_m = xcut_m,
               B_train = B_train_obj,
               B_test = B_test_obj,
               dif_order = dif_order,
               alpha = alpha,
               beta = beta,
               nu = df,
               all_var_splits = all_var_splits,
               n_tree = n_tree,
               tau_mu = tau_mu,
               tau = nsigma^(-2),
               a_tau = a_tau,
               d_tau = d_tau,
               tau_beta = tau_beta,
               basis_subindex = basis_subindex,
               P = P_train_main,
               P_interaction  = P_train_interaction,
               P_inv = P_inv,
               P_inv_interaction = P_inv_interaction,
               node_min_size = node_min_size,
               all_var = all_var,
               stump = stump,
               a_tau_beta_j = a_tau_beta_j,
               d_tau_beta_j = d_tau_beta_j,
               interaction_term = interaction_term,
               interaction_list = interaction_list,
               dummy_x = dummy_x,
               linero_sampler = linero_sampler,
               robust_prior = robust_prior,
               robust_delta = robust_delta,
               a_delta = a_delta,
               d_delta = d_delta)

  #   So to simply interepret the element all_var_splits each element correspond
  #to each variable. Afterwards each element corresponds to a cutpoint; Finally,
  #inside that level we would have the index for the the left and right nodes;



  # Create a list of arrays to store the contribution of each tree to the main effects
  if(main_effects_pred){
      tree_main_effects <- vector("list",n_mcmc)
      for( iii_ in 1:n_mcmc){
        tree_main_effects[[iii_]] <- array(0,dim = c(nrow(data$x_train),length(data$basis_subindex),data$n_tree))
      }
  }




  # Visualzing the main effects
  if(plot_preview){
    par(mfrow = c(2,5))
    # Visualizing some basis functions
    for(visu_basis in 1:10){
      plot(x  = 0, xlim = c(0,1), ylim = c(-1,1),
           type= 'n', ylab= paste0("B.",visu_basis), main = paste0("B(",visu_basis,")"))
      for(basis_col in 1:NCOL(B_train_obj[[visu_basis]])){
        points(x_train_scale[,visu_basis], B_train_obj[[visu_basis]][,basis_col],
               col = basis_col, pch = 20, xlab = paste0("x.",visu_basis))
      }
    }
  }

  # Initialing for storing post samples
  post <- 0


  # Initialising all the stumps
  for(k in 1:data$n_tree){
      if(data$all_var){
        forest[[k]] <- init_stump(data = data,initial_var = NULL)
      } else {
        # if(k <= length(dummy_x$continuousVars)){
          forest[[k]] <- init_stump(data = data,initial_var = stump_init_indicator[k])
        # } else {
        #   forest[[k]] <- init_stump(data = data, initial_var = sample(1:length(dummy_x$continuousVars),size = 1))
        # }
      }
  }

  # This warning was to avoid to initiliase an unbalanced number of trees regarding some variables.
  # if(length(forest %>% lapply(function(x)x$node0$master_var) %>% unlist() %>% table() %>% unique())>1){
  #   stop("Error on the stumps initalisations")
  # }

  # and tree predictions
  trees_fit <- matrix(0,nrow = n_tree,ncol = nrow(x_train_scale))
  trees_fit_test <- matrix(0,nrow = n_tree, ncol  = nrow(x_test_scale))

  # For cases where the tree is greater than one;
  if(scale_init){
    if(n_tree>1){
      # Initial prediction
      for(i in 1:n_tree){
        trees_fit[i,] <- y_scale/n_tree
      }
    }
  }


  # Setting a matrix to store the frequency that a variable appear within terminal nodes
  variable_importance_matrix <- matrix(0,nrow = n_mcmc, ncol = length(data$basis_subindex))


  # Initialsing the loop
  for(i in 1:n_mcmc){

    # Initialising the partial train tree fits
    if(store_tree_fit){
      partial_train_fits[[i]] <- vector("list",data$n_tree)
      names(partial_train_fits[[i]]) <- paste0("tree",1:data$n_tree)
    }

    # Initialising orogress bar
    progress <- i / n_mcmc * 100

    x1_pred <- numeric(nrow(x_train))

    for(t in 1:data$n_tree){


      # Calculating the partial residuals
      if(n_tree>1){
        partial_residuals <- y_scale-colSums(trees_fit[-t,,drop = FALSE])
      } else {
        partial_residuals <- y_scale
      }


      # Sample a verb
      if(isFALSE(all_var)){
        # verb <- sample(c("grow","add_variable","prune", "remove_variable", "change"), prob = c(0.15,0.15,0.15,0.15,0.4),size = 1)
        verb <- sample(c("grow","add_variable","prune", "remove_variable", "change","change_variable"), prob = c(0.15,0.15,0.15,0.15,0.2,0.2),size = 1)
      } else {
        verb <- sample(c("grow", "prune", "change"), prob = c(0.3,0.3,0.4))
      }
      # verb <- sample(c("grow","prune"),size = 1)

      # Verifuy this line here, if does really make change
      # if(i<500){
      #   verb <- "change"
      # }

      # Checking the trees variables

      # Sampling a verb
      if(verb == "grow"){
          forest[[t]] <- grow(tree = forest[[t]],
                              curr_part_res = partial_residuals,
                              data = data)

      } else if (verb == "add_variable"){

          forest[[t]] <- add_variable(tree = forest[[t]],
                                      curr_part_res = partial_residuals,
                                      data = data)
      } else if (verb == "prune"){
          forest[[t]] <- prune(tree = forest[[t]],
                               curr_part_res = partial_residuals,
                               data = data)

      } else if(verb == "remove_variable"){
          forest[[t]] <- remove_variable(tree = forest[[t]],
                                         curr_part_res = partial_residuals,
                                         data = data)
      } else if (verb == "change"){
          forest[[t]] <- change(tree = forest[[t]],
                                curr_part_res = partial_residuals,
                                data = data)

      } else if(verb == "change_variable") {
         forest[[t]] <- change_stump(tree = forest[[t]],
                                         curr_part_res = partial_residuals,
                                         data = data)
      }
      # Updating the betas
      update_betas_aux <- updateBetas(tree = forest[[t]],
                                      curr_part_res = partial_residuals,
                                      data = data,
                                      trees_fit = trees_fit,tree_number = t)

      forest[[t]] <- update_betas_aux$tree

      # Getting the predictions
      # tree_predictions <- getPredictions(tree = forest[[t]],
      #                                    data = data)
      if(main_effects_pred){
          tree_main_effects[[i]][,,t] <- update_betas_aux$y_hat_train
      }

      trees_fit[t,] <- update_betas_aux$tree_hat_train
      trees_fit_test[t,] <- update_betas_aux$tree_hat_test


      # Visualzing the update for the basis function
      if(plot_preview){
        par(mfrow =  c(1,2))
        plot(NULL,type = 'n', xlim = c(0,1), ylim = c(-1,1), ylab = "",
             main = "Individual Basis")
        selected_basis <- 1
        final_sum <- numeric(nrow(data$x_train))
        for(col_basis in 1:NCOL(data$B_train[[selected_basis]])){
          current_beta <- forest[[t]]$node0$betas_vec[data$basis_subindex[[selected_basis]]]
          # current_beta <- rep(1, NCOL(data$B_train[[selected_basis]])) # If I just want to see the basis
          basis_pred <- (data$B_train[[selected_basis]][,col_basis, drop = FALSE]%*%matrix(current_beta[col_basis]))
          final_sum <- final_sum + basis_pred
          points(x_train[,selected_basis],
                 basis_pred,
                 col = ggplot2::alpha(col_basis,0.5), pch = 20)
        }
        plot(x_train[,selected_basis], ylim = c(-1,1), final_sum , pch = 4, main = "Sum all basis.")
      }

      if(store_tree_fit){
        partial_train_fits[[t]] <- update_betas_aux$y_hat_train
      }

      # Adding up the contribution for each tree with respect to the covariate (i)
      if(main_effects_pred){
        for(ii in 1:length(main_effects_train_list)){
          main_effects_train_list[[ii]][i,] <- main_effects_train_list[[ii]][i,] + update_betas_aux$y_hat_train[,ii]
          # main_effects_test_list[[ii]][i,] <- main_effects_test_list[[ii]][i,] + update_betas_aux$y_hat_test[,ii]
        }

      }

    }


    # Counting the number of times that a vriable was selected
    if(varimportance_bool){
        variable_importance_matrix[i, ] <- varimportance(forest = forest,data = data)
    }

    # Getting final predcition
    y_hat <- colSums(trees_fit)
    y_hat_test <- colSums(trees_fit_test)


    # Seeing the results for the unidimensional cases.
    if(NCOL(data$x_train)==1){
      if(plot_preview){
        plot(x_train_scale,y_scale)
        for(plot_i in 1:n_tree){
          points(x_train_scale,trees_fit[plot_i,],pch=20,col = ggplot2::alpha(plot_i,0.2))
        }
        points(x_train_scale,y_hat,col = "blue",pch=20)
      }
    }


    # Updating delta
    data$robust_delta <- update_delta(data = data)

    # Updating all other parameters
    if(update_tau_beta){
        data$tau_beta <- update_tau_betas_j(forest = forest,data = data)
    }


    # Getting tau
    data$tau <- update_tau(y_train_hat = y_hat,
                           data = data, forest = forest)


    # Storing all predictions and (TREES depending on the boolean)
    if(store_tree_fit){
      all_trees[[i]] <- forest
      all_trees_fit[[i]] <- partial_train_fits
    }

    all_tau[[i]] <- data$tau
    all_delta[i,] <- data$robust_delta
    all_y_hat[i,] <- y_hat
    all_y_hat_test[i,] <- y_hat_test
    all_tau_beta[i,] <- data$tau_beta
    # all_delta[i] <- data$delta


    # Print progress bar
    cat("\rProgress: [", paste(rep("=", floor(progress / 5)), collapse = ""),
        paste(rep(" ", floor((100 - progress) / 5)), collapse = ""),
        "] ", sprintf("%.1f%%", progress))

    # Flush the output
    flush.console()

    # Simulate some work
    Sys.sleep(0.1)

  }



  # Normalising elements
  all_tau_norm <- numeric(n_mcmc)

  if(store_tree_fit){
    all_trees_fit_norm <- vector("list",n_mcmc)
  }

  all_y_hat_norm <- matrix(NA,nrow = nrow(all_y_hat),ncol = ncol(all_y_hat))
  all_y_hat_test_norm <- matrix(NA,nrow = nrow(all_y_hat_test),ncol = ncol(all_y_hat_test))

  # Returning to the original scale
  if(scale_bool){

    all_tau_norm <- all_tau/((max_y-min_y)^2)

    for(post_iter in 1:n_mcmc){

      if(store_tree_fit){
        for(tree_number in 1:n_tree){
          all_trees_fit_norm[[post_iter]][[tree_number]] <- unnormalize_bart(z = all_trees_fit[[post_iter]][[tree_number]],a = min_y,b = max_y)
        }
      }
      all_y_hat_norm[post_iter,] <- unnormalize_bart(z = all_y_hat[post_iter,],a = min_y,b = max_y)
      all_y_hat_test_norm[post_iter, ] <- unnormalize_bart(z = all_y_hat_test[post_iter,],a = min_y,b = max_y)

      if(main_effects_pred){
        for(ii in 1:length(main_effects_train_list)){
          main_effects_train_list_norm[[ii]][post_iter,] <- unnormalize_bart_me(z = main_effects_train_list[[ii]][post_iter,],a = min_y,b = max_y)
          # main_effects_test_list_norm[[ii]][post_iter,] <- unnormalize_bart_me(z = main_effects_test_list[[ii]][post_iter,],a = min_y,b = max_y)
        }
      }
    }
  } else {
    all_tau_norm <- all_tau
    all_trees_fit_norm <-all_trees_fit

    all_y_hat_norm <- all_y_hat
    all_y_hat_test_norm <- all_y_hat_test


    for(post_iter in 1:n_mcmc){

      if(main_effects_pred){
        for(ii in 1:length(main_effects_train_list)){
          main_effects_train_list_norm[[ii]][post_iter,] <- main_effects_train_list[[ii]][post_iter,]
        }
      }
    }
  }

  # Visualizing some plots
  if(plot_preview){


    # Looking the effect of the basis calculations and its coefficients
    selected_node_one <- forest[[3]]$node1
    selected_node_two <- forest[[3]]$node2
    betas_indexes <- data$basis_subindex[[selected_node_one$master_var]]

    B_tree_one <-  data$B_train[["3"]][selected_node_one$train_index,]

    par(mfrow=c(1,1))
    plot(NA, type = "n", ylim =c(-20,20), xlim = c(0,1))
    betas_ <- selected_node_one$betas_vec[betas_indexes]
    # betas_ <- rep(1,length(betas_indexes))
    for(basis_num in 1:length(betas_)){
      points(x_train_scale[selected_node_one$train_index,3],
             unnormalize_bart_me(B_tree_one[,basis_num]*betas_[basis_num],a = min_y,b = max_y),
             col = "blue", pch = 20)
    }

    B_tree_two <-  data$B_train[["3"]][selected_node_two$train_index,]
    betas_ <- selected_node_two$betas_vec[betas_indexes]
    # betas_ <- rep(1,length(betas_indexes))

    for(basis_num in 1:length(betas_)){
      points(x_train_scale[selected_node_two$train_index,3],
             unnormalize_bart_me(B_tree_two[,basis_num]*betas_[basis_num],a = min_y,b = max_y),
             col = "red", pch = 20)
    }

    points(x_train_scale[,3],
           unnormalize_bart_me(main_effects_train_list[[3]][i,],a = min_y,b = max_y),
           col = "black", pch = 20)

    if(interaction_term){

      # Main effect range difference
      range_basis_j_predictions <- numeric(length = length(main_effects_train_list_norm))
      range_tree_bais_j_predictions <- matrix(NA, ncol = n_tree, nrow = length(main_effects_train_list_norm))
      n_burn_plot <- 100

      # points(x_train[,1],y_train)
      par(mfrow = c(2,floor(NCOL(data$x_train)/2)))
      for(jj in 1:12){

        if(jj <= NCOL(data$x_train)){
            plot(x_train[,jj],colMeans(main_effects_train_list_norm[[jj]][n_burn_plot:i,, drop = FALSE]),main = paste0('X',jj),
                 ylab = paste0('G(X',jj,')'), ylim = c(-10,10) ,pch=20,xlab = paste0('x.',jj), col = ggplot2::alpha("black",1.0))
        }    else if(jj == 11 ) {
            # plot(x_train[,1], colMeans(main_effects_train_list_norm[[jj]][1:i,,drop = FALSE]), ylab = "G(x.1,x2)", xlab = "X.1")
            # plot(x_train[,2], colMeans(main_effects_train_list_norm[[jj]][1:i,,drop = FALSE]), ylab = "G(x.1,x2)", xlab = "X.1")
            # plot(x_train[,1], x_train[,2], cex = colMeans(main_effects_train_list_norm[[jj]][1000:i,,drop = FALSE]) + colMeans(main_effects_train_list_norm[[1]][1000:i,,drop = FALSE]) + colMeans(main_effects_train_list_norm[[2]][1000:i,,drop = FALSE]),
            #      ylab = "X.2", xlab = "X.1")
            plot(x_train[,1], x_train[,2], cex = colMeans(main_effects_train_list_norm[[jj]][1000:i,,drop = FALSE]), ylab = "X.2", xlab = "X.1")
            # plot(x_train[,2], x_train[,1], cex = 10*sin(pi*x_train[,1]*x_train[,2]), col = "blue",  ylab = "X.2", xlab = "X.1")

        }


        range_basis_j_predictions[jj] <- abs(diff(range(colMeans(main_effects_train_list_norm[[jj]][1:i,, drop = FALSE]))))

        for(tree_number in 1:data$n_tree){

            n_burn_plot_count <- 1
            # Create an auxiliar matrix for the main effect
            aux_main_effect_matrix <- matrix(0, nrow = (n_mcmc-n_burn_plot),ncol = nrow(data$x_train))


            if(main_effects_pred){
                for(mcmc_aux in (n_burn_plot+1):n_mcmc){
                  aux_main_effect_matrix[n_burn_plot_count,] <- unnormalize_bart_me(tree_main_effects[[mcmc_aux]][,jj,tree_number],a = min_y,b = max_y)
                  n_burn_plot_count <- n_burn_plot_count  + 1

                }
            }

            if(jj <= NCOL(data$x_train)){
              if(all(colMeans(aux_main_effect_matrix)==0)){
              } else {
                  points(x_train[,jj],colMeans(aux_main_effect_matrix),main = paste0('X',jj),
                       col = ggplot2::alpha("darkred",0.1), pch = 20)
              }
            }
            # range_tree_bais_j_predictions[tree_number, jj] <- diff(range(colMeans(aux_main_effect_matrix)))

        }

        # points(x_train[,jj],colMeans(main_effects_train_list[[jj]][1:i,, drop = FALSE]),main = paste0('X',jj),
        #      ylab = paste0('G(X',jj,')'),ylim = c(-0.5,0.5),pch=20,xlab = paste0('x.',jj))
        # if( jj ==3 ){
        #   points(x_train[,jj],20*(x_train[,jj]-0.5)^2,pch = 20, col = "blue")
        # } else if ( jj == 4){
        #   points(x_train[,jj],10*x_train[,jj],pch = 20, col = "blue")
        # } else if ( jj == 5){
        #   points(x_train[,jj],5*x_train[,jj],pch = 20, col = "blue")
        # }

      }
      par(mfrow = c(1,1))
    }

  }

  # Quick check over the RMSE
  if(plot_preview){
    # par(mrow=c(1,3))
    rmse(all_y_hat_test_norm[1:i,] %>% colMeans(na.rm = TRUE), y = sim_test$y)
    rmse(sim_test$y, aux$yhat.test.mean)
    rmse(sim_test$y, aux2$y_hat_test_mean)

    # Getting some plots
    par(mfrow=c(1,3))
    plot(all_y_hat_test_norm[2000:i,] %>% colMeans(na.rm = TRUE), y = sim_test$y,
         xlab = "spBART_test", ylab = "y_test", main = "spBART")
    plot(sim_test$y, aux$yhat.test.mean,
         xlab = "BART_test", ylab = "y_test", main = "BART")
    plot(sim_test$y, aux2$y_hat_test_mean,
         xlab = "softBART_test", ylab = "y_test", main = "softBART")
  }


  # Plotting vaiable importance
  if(plot_preview){
    par(mfrow=c(1,2))
    burn_sample_ <- 100
    plot(1:NCOL(variable_importance_matrix),variable_importance_matrix[burn_sample_:i,,drop = FALSE] %>% colMeans(),
         ylab = "Prop. pred_var", xlab = "Predictor", main = c("Proportion Tree pred.vars"))
    points((1:NCOL(variable_importance_matrix))[c(1:5,11)],variable_importance_matrix[burn_sample_:i,c(1:5,11),drop = FALSE] %>% colMeans(),
         ylab = "Prop. pred_var", xlab = "Predictor/Basis", pch = 20)


    plot(1:NCOL(variable_importance_matrix),all_tau_beta[burn_sample_:i,,drop = FALSE] %>% colMeans(na.rm = TRUE),
         ylab = expression(bar(tau[beta])), xlab = "Predictor", main = c("Tau_beta_posterior_mean"))
    points((1:NCOL(variable_importance_matrix))[c(1:5,11)],all_tau_beta[burn_sample_:i,c(1:5,11),drop = FALSE] %>% colMeans(na.rm = TRUE),
           ylab = "mean_tau_beta", xlab = "Predictor/Basis", pch = 20)

  }
  # importance_var_ <- variable_importance_matrix[501:i,,drop = FALSE] %>% colMeans()

  # plot(1:55,log((importance_var_*tau_beta_post_mean_)/sum(importance_var_*tau_beta_post_mean_)))
  # currr_ <- numeric(n_mcmc)
  # get_which <- 0
  # for(test_mcmc_ in 1:n_mcmc){
  #     for(test_ in 1:data$n_tree){
  #       if(any(tree_main_effects[[test_mcmc_]][,11,test_]>0)){
  #         currr_[test_mcmc_] <- currr_[test_mcmc_] + 1
  #         get_which <- test_
  #       }
  #     }
  # }
  # plot(currr_)
  # # Some extra analysis
  # par(mfrow=c(2,2))
  # plot(all_tau_norm, type = 'l', ylab = expression(tau), xlab = "MCMC_iter",main = expression(tau))
  # abline(v = 1318, lty = 'dashed', col = 'blue')
  #
  # plot(all_tau_beta[,11], type = 'l', ylab = expression(tau[11]), xlab = "MCMC_iter",main = expression(tau[11]))
  # abline(v = 1318, lty = 'dashed', col = 'blue')
  #
  # plot(all_tau_beta[,1], type = 'l', ylab = expression(tau[1]), xlab = "MCMC_iter",main = expression(tau[1]))
  # abline(v = 1318, lty = 'dashed', col = 'blue')
  #
  # plot(all_tau_beta[,2], type = 'l', ylab = expression(tau[2]), xlab = "MCMC_iter",main = expression(tau[2]))
  # abline(v = 1318, lty = 'dashed', col = 'blue')


  # plot(colMeans(all_y_hat_norm),y_train)
  # sqrt(crossprod((colMeans(all_y_hat_norm)-y_train))/n_)

  # plot(10*sin(pi*x_train_scale[,1]*x_train_scale[,2]),main_effects_train_list[[11]][1:i,] %>% colMeans(), ylab = "y")

  # ====== Few analyses from the results ======
  #           (Uncomment to run those)
  # ===========================================
  #

  # plot(all_tau,type = "l")
  # plot(all_tau,type = "l")
  #
  # y1_hat <- matrix(0,nrow = n_post,ncol = nrow(data$x_train))
  #
  # for(i in 1:nrow(y1_hat)){
  #   for(t in 1:n_tree){
  #     y1_hat[i,] <- all_trees_fit[[i]][[t]][,1]
  #   }
  # }
  #
  # plot(x_train_scale[,1],colMeans(y1_hat[1801:3000,,drop = FALSE]))
  # plot(colMeans(y1_hat[1801:3000,,drop = FALSE]),y_train)
  # # plot(x_train[,1],10 * sin(pi * x_train[, 1])) # For x1
  # # plot(x_train[,2],20 * (x_train[, 2] - 0.5)^2) # For x1
  #
  # plot(all_tau_beta, type = "l")
  # plot(all_delta, type = "l")
  # plot(x_train_scale,y_scale)
  # points(x_train_scale,colMeans(all_y_hat[1501:2000,]),pch= 20, col = "blue")
  # # ============================================
  #
  # curr <- 0
  # tree_lengths <- numeric()
  #
  # for(i in 1:length(all_trees)){
  #
  #     curr <- curr + 1
  #
  #     for(j in 1:n_tree){
  #         tree_lengths[curr] <- length(all_trees[[i]][[j]])
  #     }
  #
  # }
  # rmse(x = sim_test$y,y = all_y_hat_test_norm %>% colMeans())

  # Deleting var_importance_matrix
  if(isFALSE(varimportance_bool)){
    variable_importance_matrix <- NULL
  }

  # Setting tree_main_effects as null if necessary
  if(isFALSE(main_effects_pred)){
    tree_main_effects <- NULL
  }


  # Return the list with all objects and parameters
  return(list(y_train_hat = all_y_hat_norm,
              y_test_hat = all_y_hat_test_norm,
              all_tau = all_tau_norm,
              all_tau_beta = all_tau_beta,
              all_delta = all_delta,
              prior = list(n_tree = n_tree,
                           alpha = alpha,
                           beta = beta,
                           tau_mu = tau_mu,
                           a_tau = a_tau,
                           d_tau = d_tau,
                           a_tau_beta = a_tau_beta_j,
                           d_tau_beta = d_tau_beta_j,
                           a_delta = a_delta,
                           d_delta = d_delta,
                           nu = df,
                           centered_basis = centered_basis,
                           eta = eta),
              mcmc = list(n_mcmc = n_mcmc,
                          n_burn = n_burn,
                          all_trees = all_trees,
                          main_effects_train = main_effects_train_list_norm,
                          main_effects_test = main_effects_test_list_norm,
                          tree_main_effects = tree_main_effects,
                          variable_importance_matrix = variable_importance_matrix),
              data = list(x_train = x_train,
                          y_train = y_train,
                          B_train = B_train_obj,
                          x_test = x_test,
                          B_test = B_test_obj,
                          basis_subindex = basis_subindex)))

}




