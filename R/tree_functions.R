# A function to get which variables are used within each terminal node
varimportance <- function(forest, data){

  # Getting the auxilinar vector
  var_counter <- var_counter_aux <- numeric(length = length(data$basis_subindex))

  for(ii_ in 1:data$n_tree){
    tree_ <- forest[[ii_]]
    terminal_names <- get_terminals(tree = tree_)
    n_terminal <- length(terminal_names)
    var_counter_aux <- numeric(length = length(data$basis_subindex))

    for(jj_ in 1:length(terminal_names)){
      var_counter_aux[tree_[[terminal_names[jj_]]]$pred_vars] <- var_counter_aux[tree_[[terminal_names[jj_]]]$pred_vars] + 1
    }
    var_counter <- var_counter + var_counter_aux/n_terminal
  }

  return(var_counter/data$n_tree)
}
# Creating a stump for a tree
init_stump <- function(data,
                  initial_var){

  # Creating the base node
  node <- list()

      # Create a boolean to evaluate if all variables are gonna to be used or not
      if(data$all_var){
          node[["node0"]] <- list(
            # Creating the node number
            node_number = 0,
            master_var = initial_var,
            pred_vars = 1:length(data$basis_subindex),
            inter = NA,
            isRoot = TRUE,
            # Creating a vector with the tranining index
            train_index = 1:nrow(data$x_train),
            test_index = 1:nrow(data$x_test),
            depth_node = 0,
            node_var = NA,
            node_cutpoint_index = NA,
            left = NA,
            right = NA,
            parent_node = NA,
            ancestors = NA,
            terminal = TRUE,
            betas_vec = rep(0, length(unlist(data$basis_subindex)))
        )
      } else {
        node[["node0"]] <- list(
          # Creating the node number
          node_number = 0,
          master_var = initial_var,
          pred_vars = initial_var,
          inter = NA,
          isRoot = TRUE,
          # Creating a vector with the tranining index
          train_index = 1:nrow(data$x_train),
          test_index = 1:nrow(data$x_test),
          depth_node = 0,
          node_var = NA,
          node_cutpoint_index = NA,
          left = NA,
          right = NA,
          parent_node = NA,
          ancestors = NA,
          terminal = TRUE,
          betas_vec = rep(0, length(unlist(data$basis_subindex)))
        )
    }


  # Returning the node
  return(node)

}

# Get all the terminal nodes
get_terminals <- function(tree){

  # Return the name of the termianl nodes
  return(names(tree)[unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)])
}

# Get nog terminal nodes
get_nogs <- function(tree){

  # Return the name of the termianl nodes
  non_terminal <- names(tree)[!unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)]

  # In case there are non nonterminal nondes
  if(length(non_terminal)==0){
    return(non_terminal)
  }

  bool_nog <- vector("logical",length = length(non_terminal))
  for(i in 1:length(bool_nog)){
    # Checking if both children are terminal
    if( tree[[tree[[non_terminal[i]]]$left]]$terminal & tree[[tree[[non_terminal[i]]]$right]]$terminal) {
      bool_nog[i] <- TRUE
    }
  }

  return(  non_terminal[bool_nog])
}

# Getting the maximum node index number
get_max_node <- function(tree){

  # Return the name of the termianl nodes
  return(max(unlist(lapply(tree, function(x){x$node_number}),use.names =  TRUE)))
}




# A function to calculate the loglikelihood
nodeLogLike <- function(curr_part_res,
                        j_,
                        index_node,
                        data){

  # Subsetting the residuals
  curr_part_res_leaf <- curr_part_res[index_node]

  result <- numeric()
  # for(hh in 1:55){
  # j_ <- hh
  # Getting the number of observationsin the terminal node
  n_leaf <- length(index_node)
  d_basis <- length(j_)
  ones <- matrix(1,nrow = n_leaf)


  if(length(j_)==0){
    stop(" Node Log-likelihood: No variables")
  }

  # Using the Andrew's approach I would have
  mean_aux <- rep(0,length(curr_part_res_leaf))
  cov_aux <- matrix(0,nrow = length(index_node),ncol = length(index_node))
  # diag_tau_beta_inv <- diag(x = 1/unique(data$tau_beta), nrow = )

  for(jj in 1:length(j_)){
      # Adding the quantities with respect to the interaction
      if(j_[jj] <= length(data$dummy_x$continuousVars)){
        cov_aux <- cov_aux + (data$tau_beta[j_[jj]]^(-1))*data$B_train[[j_[jj]]][index_node,,drop = FALSE]%*%tcrossprod(data$P_inv,(data$B_train[[j_[jj]]][index_node,,drop = FALSE]))
      } else {
        cov_aux <- cov_aux + (data$tau_beta[j_[jj]]^(-1))*data$B_train[[j_[jj]]][index_node,,drop = FALSE]%*%tcrossprod(data$P_inv_interaction,(data$B_train[[j_[jj]]][index_node,,drop = FALSE]))
      }
  }

  # Adding the main diagonal
  cov_aux <- (diag(nrow = n_leaf) + cov_aux)/data$tau

  # Defining Keefe suggestion
  pivot_chol <- function(x) {
    x <- chol(x, pivot=TRUE)
    r <- attr(x, "rank")
    p <- order(attr(x, "pivot"))
    x[-seq_len(r),-seq_len(r)] <- 0
    x <- x[,p]
    x
  }


  chol_cov_aux <- tryCatch(chol(cov_aux), error=function(e) pivot_chol(cov_aux))

  result <- mvnfast::dmvn(X = curr_part_res_leaf,mu = mean_aux,
                          sigma = chol_cov_aux, log = TRUE, isChol=TRUE)
  # result <- mvnfast::dmvn(X = curr_part_res_leaf,mu = mean_aux,
  #                         sigma = (data$tau^(-1))*cov_aux ,log = TRUE)


  # plot(1:55,result)

  return(c(result))

}

# # ========================================
# # ATTENTION HERE: THIS IS THE STABLE VERSIOMN OF THE LEGIT FUNTION FOR NODE
# # LOGLIKEBIGD , THE OTHER FUNCTION HAS THE SAME NAME BUT USE THE SUM OF COVARIANCE
# #FUNCTIONS TO CALCULATE THE MAIN VARIANCE
# # ======================================
# nodeLogLikeBigD <- function(curr_part_res,
#                         j_,
#                         index_node,
#                         data){
#
# # Subsetting the residuals
# curr_part_res_leaf <- curr_part_res[index_node]
#
# # Getting the number of observationsin the terminal node
# n_leaf <- length(index_node)
# d_basis <- length(j_)
# ones <- matrix(1,nrow = n_leaf)
# D_subset_index <- unlist(data$basis_subindex[j_])
# D_leaf <- data$D_train[index_node,D_subset_index, drop = FALSE]
#
# if(ncol(D_leaf)==0){
#   stop(" Node Log-likelihood: No variables")
# }
#
# # Using the Andrew's approach I would have
# mean_aux <- rep(0,length(curr_part_res_leaf))
# cov_aux <- matrix(0,nrow = length(index_node),ncol = length(index_node))
#
# # Setting the dif_order==0
# if(data$dif_order==0){
#   stop( "Work with penalised splines")
# } else {
#   P_aux <- matrix(0, nrow = nrow(data$P), ncol = ncol(data$P))
#   for(jj in 1:length(j_)){
#     P_aux[data$basis_subindex[[j_[jj]]],data$basis_subindex[[j_[jj]]]] <- (data$tau_beta[jj])*data$P[data$basis_subindex[[j_[jj]]],data$basis_subindex[[j_[jj]]]]
#   }
#   P_aux <- P_aux[D_subset_index,D_subset_index]
#
#   cov_aux <-  diag(x = (data$tau^(-1)),nrow = n_leaf) + D_leaf%*%solve(P_aux)%*%t(D_leaf)
#
# }
#
#
# result <- mvnfast::dmvn(X = curr_part_res_leaf,mu = mean_aux,
#                         sigma = cov_aux ,log = TRUE)
#
#
#   return(c(result))
#
# }

# Grow a tree
grow <- function(tree,
                 curr_part_res,
                 data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0

  # acceptance_grid <- numeric(100)
  # for(kk in 1:100){
  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2

    # Sample a split var
    # ===== Uncomment this line below after ========
    p_var <- sample(1:NCOL(data$x_train),size = 1)
    # ==============================================
    # p_var <- 7

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[g_node$train_index,p_var])

    # Case of invalid range
    if(length(valid_range_grow)==0){
      return(tree)
    }

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)
    # sample_cutpoint <- valid_cutpoint[kk]

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% g_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% g_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% g_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% g_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(g_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(g_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }


    # === Uncomment those lines after
    if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }

  # Calculating loglikelihood for the grown node, the left and the right node
  # Recover the g_node index

  node_index_var <- g_node$pred_vars

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = node_index_var,
                           index_node = g_node$train_index,
                           data = data)


  left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               j_ = node_index_var,
                               index_node = left_index,
                               data = data)

  right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                j_ = node_index_var,
                                index_node = right_index,
                                data = data)

  # Calculating the prior
  prior_loglike <- log(data$alpha*(1+g_node$depth_node)^(-data$beta)) + # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+g_node$depth_node+1)^(-data$beta)) - # plus the prior of the two following nodes being terminal
    log(1-data$alpha*(1+g_node$depth_node)^(-data$beta)) # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_nog_nodes+1))-log(0.3/n_t_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+left_loglike+right_loglike+prior_loglike+log_trasition_prob)
  # acceptance_grid[kk] <- acceptance



  # par(mfrow=c(1,2))
  # plot(data$xcut_m[,2],(acceptance_grid), main = "Acceptance to split on X2", xlab = "X2", ylab = "Prob. Acceptance")

  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    if(any(is.na(g_node$ancestors))){
      new_ancestors <- p_var
    } else {
      new_ancestors <- c(g_node$ancestors,p_var)
    }

    left_node <- list(node_number = max_index+1,
                      master_var = g_node$master_var,
                      pred_vars = g_node$pred_vars,
                      inter = g_node$inter,
                      isRoot = FALSE,
                      train_index = left_index,
                      test_index = left_test_index,
                      depth_node = g_node$depth_node+1,
                      node_var = p_var,
                      node_cutpoint_index = sample_cutpoint,
                      left = NA,
                      right = NA,
                      parent_node = g_node_name,
                      ancestors = new_ancestors,
                      terminal = TRUE,
                      betas_vec = g_node$betas_vec)

    right_node <- list(node_number = max_index+2,
                       master_var = g_node$master_var,
                       pred_vars = g_node$pred_vars,
                       inter = g_node$inter,
                       isRoot = FALSE,
                       train_index = right_index,
                       test_index = right_test_index,
                       depth_node = g_node$depth_node+1,
                       node_var = p_var,
                       node_cutpoint_index = sample_cutpoint,
                       left = NA,
                       right = NA,
                       parent_node = g_node_name,
                       ancestors = new_ancestors,
                       terminal = TRUE,
                       betas_vec = g_node$betas_vec)

    # Modifying the current node
    tree[[g_node_name]]$left = paste0("node",max_index+1)
    tree[[g_node_name]]$right = paste0("node",max_index+2)
    tree[[g_node_name]]$terminal = FALSE

    tree[[paste0("node",max_index+1)]] <- left_node
    tree[[paste0("node",max_index+2)]] <- right_node


  } else {

    # Do nothing

  }

  # Return the new tree
  return(tree)
}


# Add variable
add_variable <- function(tree,
                 curr_part_res,
                 data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  g_node_name <- sample(terminal_nodes,size = 1)
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0

  # Getting all possible candidates regarding the 'master variable' of that tree

  available_interactions <- which(apply(data$interaction_list, 2 , function(x) {g_node$master_var %in% x}))

  names_available_interactions <- apply(data$interaction_list[,available_interactions, drop = FALSE], 2, function(x){paste0(x,collapse = "")})
  index_available_interactions <- c(g_node$master_var,which(names(data$basis_subindex) %in% names_available_interactions))

  actual_index_available_interactions <-  index_available_interactions[!(index_available_interactions %in% g_node$pred_vars)]# Removing the ones that are already in the pred_vars vector
  # Merging the index of the avaible candidates and the candidations that can be add removing the ones that are in the prediction var

  # Selecting the variable to be added into the prediction vars
  if(length(actual_index_available_interactions)>1){
      p_var <- sample(actual_index_available_interactions,size = 1)
  } else {
      p_var <- actual_index_available_interactions
  }


  if(length(actual_index_available_interactions) == 0) {
    warning("No available interaction or variables")
    return(tree)
  }

  new_node_index_var <- unique(sort(c(g_node$pred_vars,p_var)))

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                               j_ = g_node$pred_vars,
                               index_node = g_node$train_index,
                               data = data)


  new_g_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                       j_ = new_node_index_var,
                                       index_node = g_node$train_index,
                                       data = data)

  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+new_g_loglike)


  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # In this case I will only modify the pred_vars component in the terminal nodes
    tree[[g_node_name]]$pred_vars <- sort(unique(new_node_index_var))

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)
}


# Pruning a tree
prune <- function(tree,
                  curr_part_res,
                  data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)

  # Just in case to avoid errors
  if(n_nog_nodes==0){
    return(tree)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(nog_nodes,size = 1)
  p_node <- tree[[p_node_name]]

  # Getting the indexes from the left and right children from the pruned node
  children_left_index <- tree[[p_node$left]]$train_index
  children_right_index <- tree[[p_node$right]]$train_index
  children_left_ancestors <- tree[[p_node$left]]$ancestors
  children_right_ancestors <- tree[[p_node$right]]$ancestors

  # Calculating loglikelihood for the grown node, the left and the right node

  node_index_var <- p_node$pred_vars

  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           index_node = p_node$train_index,
                           j_ = node_index_var,
                           data = data)


  p_left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                 index_node =  children_left_index,
                                 j_ = node_index_var,
                                 data = data)

  p_right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                  index_node = children_right_index,
                                  j_ = node_index_var,
                                  data = data)

  # Calculating the prior
  prior_loglike <- log(1-data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the new terminal node
    log(data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+p_node$depth_node+1)^(-data$beta))  # plus the prior of the two following nodes being terminal
  # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_t_nodes))-log(0.3/n_nog_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(p_loglike-p_left_loglike-p_right_loglike+prior_loglike+log_trasition_prob)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[p_node$left]] <- NULL
    tree[[p_node$right]] <- NULL

    # Modifying back the pruned node
    tree[[p_node_name]]$left <- NA
    tree[[p_node_name]]$right <- NA
    tree[[p_node_name]]$terminal <- TRUE

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

# Remove variable
remove_variable <- function(tree,
                            curr_part_res,
                            data) {

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  p_node_name <- sample(terminal_nodes,size = 1)
  p_node <- tree[[p_node_name]]


  # Pruning a node
  if(length(p_node$pred_vars)==1){
    return(tree)
  }

  if(identical(p_node$pred_vars,p_node$master_var)){
    return(tree)
  }

  valid_terminal_node <- TRUE
  valid_count <- 0

  # Removing one variable from the candidates

  # Change this so you can change any interaction
  # Comment this line is not used to anything
  interactions_available <- p_node$pred_vars[!(p_node$pred_vars %in% p_node$master_var)]

  # Attention on this line!!
  # remove_var <- sample(interactions_available,size = 1)
  remove_var <- sample(p_node$pred_vars,size = 1)


  new_node_index_var <- unique(sort(p_node$pred_vars[!(p_node$pred_vars %in% remove_var)]))

  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                               j_ = p_node$pred_vars,
                               index_node = p_node$train_index,
                               data = data)


  new_p_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                    j_ = new_node_index_var,
                                    index_node = p_node$train_index,
                                    data = data)

  # Calculating the acceptance probability
  acceptance <- exp(-p_loglike+new_p_loglike)


  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # In this case I will only modify the pred_vars component in the terminal nodes
    tree[[p_node_name]]$pred_vars <- sort(unique(new_node_index_var))

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}


# Change a tree
change <- function(tree,
                   curr_part_res,
                   data){

  # # Changing the stump
  if(length(tree)==1 & isFALSE(data$all_var)){
    change_stump_obj <- change_stump(tree = tree,
                                     curr_part_res = curr_part_res,
                                     data = data)
    return(change_stump_obj)
  }


  # For the seocnd case
  if(length(tree)==1){
    return(tree)
  }

  # Sampling a terminal node
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  c_node_name <- sample(nog_nodes,size = 1)
  c_node <- tree[[c_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0


  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2
    # Sample a split var
    p_var <- sample(1:ncol(data$x_train),size = 1)

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[c_node$train_index,p_var])

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% c_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% c_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% c_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% c_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(c_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(c_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }

    # Avoiding having terminal nodes with just one observation
    if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }


  # Getting the node_index var
  node_index_var <- c_node$pred_vars


  # Calculating loglikelihood for the new changed nodes and the old ones
  c_loglike_left <- nodeLogLike(curr_part_res = curr_part_res,
                                    index_node = tree[[c_node$left]]$train_index,
                                    j_ = node_index_var,
                                    data = data)


  c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                      index_node = tree[[c_node$right]]$train_index,
                                      j_ =  node_index_var,
                                      data = data)

  # Calculating a new ancestors left and right
  old_p_var <- tree[[c_node$left]]$node_var

  # Storing new left and right ancestors
  new_left_ancestors <- tree[[c_node$left]]$ancestors
  new_left_ancestors[length(new_left_ancestors)] <- p_var

  new_right_ancestors <- tree[[c_node$right]]$ancestors
  new_right_ancestors[length(new_right_ancestors)] <- p_var


  new_c_loglike_left <-  nodeLogLike(curr_part_res = curr_part_res,
                                         index_node = left_index,
                                         j = node_index_var,
                                         data = data)

  new_c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                          index_node = right_index,
                                          j =  node_index_var,
                                          data = data)


  # Calculating the acceptance probability
  acceptance <- exp(new_c_loglike_left+new_c_loglike_right-c_loglike_left-c_loglike_right)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1,min = 0,max = 1)<acceptance){

    # Updating the left and the right node
    # === Left =====
    tree[[c_node$left]]$node_var <- p_var
    tree[[c_node$left]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$left]]$train_index <- left_index
    tree[[c_node$left]]$test_index <- left_test_index
    tree[[c_node$left]]$ancestors <- new_left_ancestors

    #==== Right ====
    tree[[c_node$right]]$node_var <- p_var
    tree[[c_node$right]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$right]]$train_index <- right_index
    tree[[c_node$right]]$test_index <- right_test_index
    tree[[c_node$right]]$ancestors <- new_right_ancestors

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

# Change a tree ( this gonna allow to change the 'master variable' within the stump)
change_stump <- function(tree,
                         curr_part_res,
                         data){


  # while(TRUE){ # A while to check how to function ins working --- interesting that get stuck on pred_var == 2

  if(length(tree)!=1){
    # stop("This function does not work on something beyond stumps.")
    # Sampling a terminal node
    nog_nodes <- get_nogs(tree)
    n_nog_nodes <- length(nog_nodes)
    c_node_name <- sample(nog_nodes,size = 1)
    c_node <- tree[[c_node_name]]
  } else {
    c_node_name <- "node0"
    c_node <- tree$node0
  }

  # Sampling a terminal node

  # Selecting if change the master_variable or the interaction as well
  if(stats::runif(1) < 0.5){
    change_type <- "master"
  } else {
    change_type <- "interaction"
  }

  # Recalling the vector of the interactions
  all_interactions_index_ <- (length(data$dummy_x$continuousVars)+1):length(data$basis_subindex)


  # # --- Not valid maipulation
  # Checking if there's any interaction available
  # if(isFALSE(any(c_node$pred_vars %in% all_interactions_index_))){ # I still need to allow to change one master for a interaction
  #   change_type <- "master"
  # }

  # Selecting available interactions and which interaction is gonna be replaced
  available_interactions <- which(apply(data$interaction_list, 2 , function(x) {c_node$master_var %in% x}))
  names_available_interactions <- apply(data$interaction_list[,available_interactions, drop = FALSE], 2, function(x){paste0(x,collapse = "")})

  # Should remove the MASTER-VAR from that line???
  index_available_interactions <- c(c_node$master_var,which(names(data$basis_subindex) %in% names_available_interactions))
  index_available_interactions_ONLY <- which(names(data$basis_subindex) %in% names_available_interactions) # This ONLY contains interactions


  if(change_type == "interaction"){

    # Seeing which interactions can be changed
    replaced_interaction_candidates <- c_node$pred_vars[c_node$pred_vars %in% index_available_interactions_ONLY]

    if(length(replaced_interaction_candidates)==0){# If there is no interactions to be replaced simply change the master_variable
      # IN THIS CASE I SWAP THE MAIN EFFECT BY AN INTERACTION
      available_changes <- index_available_interactions_ONLY
      change_master_pred <- TRUE
    } else {
      if(length(replaced_interaction_candidates)>1){
        replace_interaction_selected <- sample(replaced_interaction_candidates, size = 1)
      } else {
        replace_interaction_selected <- replaced_interaction_candidates
      }
      change_master_pred <- FALSE
    }

    actual_index_available_interactions <-  index_available_interactions_ONLY[!(index_available_interactions_ONLY %in% c_node$pred_vars)]# Removing the ones that are already in the pred_vars vector


    # Testing if there is any main effect in the actual_index_available_interactions, here we want to avoid that
    # ====== ACTUALLY YOU BY COMMENT THE IF() BELOW WE ALLOW TO RETURN TO THE MAIN_PRED() in the preditors ======
    # if(c_node$master_var %in% actual_index_available_interactions){ # Maybe you can and that line isn't necessary
    #   stop("The master variable should not be included into the available interactions vector.")
    # }

    available_changes <- actual_index_available_interactions

    # All possible interactions are already inclued => Try to do the CHANGE_type
    if(length(available_changes)==0){
      change_type <- 'master'
      available_changes <- (1:NCOL(data$x_train))[!((1:NCOL(data$x_train)) %in% c_node$master_var)]
      if( !(c_node$master_var %in% c_node$pred_vars) ){
        available_changes <- sort(unique(c(available_changes,c_node$master_var)))
      }
      # return(tree)
    }
  } else if (change_type == "master"){
    available_changes <- (1:NCOL(data$x_train))[!((1:NCOL(data$x_train)) %in% c_node$master_var)]
    if( !(c_node$master_var %in% c_node$pred_vars) ){
      available_changes <- sort(unique(c(available_changes,c_node$master_var)))
    }
  }

  # This only replaces for a new MAIN effect. If I want to replace for an interaction?
  if(length(available_changes) > 1) {
    p_var <- sample(available_changes,size = 1)
  } else {
    p_var <- available_changes
  }


  # Setting the new pred_vars
  if(change_type == "master"){
    new_node_index <- p_var

    # Extra error message
    if(isFALSE(new_node_index %in% 1:(NCOL(data$x_train))) ){
      stop("Selecting an interaction to be master-change type")
    }

  } else if(isTRUE(change_master_pred) & change_type == "interaction"){
    new_node_index <- p_var

    # Extra error message
    if(isFALSE(new_node_index %in% all_interactions_index_) ){
      stop("Not a main effect to be the new variable")
    }

  } else if( isFALSE(change_master_pred) & change_type == "interaction") {
    new_node_index <- c_node$pred_vars
    new_node_index[c_node$pred_vars %in% replace_interaction_selected] <- p_var
    new_node_index <- sort(unique(new_node_index))
  }

  # cat("New Node indexes: ",new_node_index,"\n")


  # Calculating loglikelihood for the new changed nodes and the old ones
  c_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                               index_node = c_node$train_index,
                               j_ = c_node$pred_vars,
                               data = data)


  c_new_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                    index_node = c_node$train_index,
                                    j_ =  new_node_index,
                                    data = data)

  # Calculating the acceptance probability
  acceptance <- exp(-c_loglike + c_new_loglike)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1,min = 0,max = 1)<acceptance){

    #  Updating the variables accordingly
    if(change_type == 'master'){
      tree[[c_node_name]]$master_var <- new_node_index
      tree[[c_node_name]]$pred_vars <- new_node_index

      if(!identical(new_node_index,c_node$pred_vars)){
          tree[[c_node_name]]$betas_vec[data$basis_subindex[[p_var]]] < 0.0  # Setting the "old value" of beta as zero.
      }

    } else if (change_type == 'interaction'){
      tree[[c_node_name]]$pred_vars <- new_node_index

      if(!identical(new_node_index,c_node$pred_vars)){
        tree[[c_node_name]]$betas_vec[data$basis_subindex[[p_var]]] < 0.0  # Setting the "old value" of beta as zero.
      }

    } else {
      error("No valid change_type into the change_stump() function")

    }

    # cat("Prediction vars selected:",(tree[["node0"]]$pred_vars),"\n")

  } else {
    # Do nothing
  }
  # }

  # Getting main-effect if there's any // Getting the interaction if there's any
  main_var_pred <- tree[[c_node_name]]$pred_vars[tree[[c_node_name]]$pred_vars %in% (1:(length(data$dummy_x$continuousVars)))]
  inter_var_pred <- tree[[c_node_name]]$pred_vars[tree[[c_node_name]]$pred_vars %in% index_available_interactions_ONLY]

  if(length(main_var_pred)==0 & length(inter_var_pred)==0){
    stop("There's a main effect and a interaction not related to the master variable")
  }

  if(!identical(main_var_pred,tree[[c_node_name]]$master_var) & length(inter_var_pred)==0){
    stop("Either the new main effect do not correspond to the the master variable! Error type 2")
  }


  # Return the new tree
  return(tree)

}

# Change a interaction in a tree, this gonna allow to change the interaction of a
#terminal node isn't exclusive to terminal nodes
change_variable <- function(tree,
                         curr_part_res,
                         data){


  if(length(tree)==1){
    new_tree <- change_stump(tree = tree,
                             curr_part_res = curr_part_res,
                             data = data)
    return(new_tree)
  }
  # Sampling a terminal node
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  c_node_name <- sample(nog_nodes,size = 1)
  c_node <- tree[[c_node_name]]


  # Recalling the vector of the interactions
  all_interactions_index_ <- (length(data$dummy_x$continuousVars)+1):length(data$basis_subindex)

  # Selecting available interactions and which interaction is gonna be replaced
  available_interactions <- which(apply(data$interaction_list, 2 , function(x) {c_node$master_var %in% x}))
  names_available_interactions <- apply(data$interaction_list[,available_interactions, drop = FALSE], 2, function(x){paste0(x,collapse = "")})
  index_available_interactions <- c(c_node$master_var,which(names(data$basis_subindex) %in% names_available_interactions))
  index_available_interactions_ONLY <- which(names(data$basis_subindex) %in% names_available_interactions) # This ONLY contains interactions

  # Seeing which interactions can be changed in the c_node$pred_vars (here I'm allowing to change the main effect as well)
  # ATTENTION HERE ONLY INTERACTIONS SHOULD BE ALLOWED AT THE MOMENT LATER YOU CAN ALLOW
  # replaced_interaction_candidates <- c_node$pred_vars[c_node$pred_vars %in% index_available_interactions_ONLY]
  replaced_interaction_candidates <- c_node$pred_vars[c_node$pred_vars %in% index_available_interactions]

  if(length(replaced_interaction_candidates)==0){# If there is no interactions/variable to be replaced simply change the master_variable
    return(tree)
    warning("All possible variables were included in the terminal node, there is no possible change.")
  } else {
    if(length(replaced_interaction_candidates>1)){ # This is just to avoid the sample bug when there's only one variable available
      replace_interaction_selected <- sample(replaced_interaction_candidates, size = 1)
    } else {
      replace_interaction_selected <- replaced_interaction_candidates
    }
  }

  actual_index_available_interactions <-  index_available_interactions[!(index_available_interactions %in% c_node$pred_vars)]# Removing the ones that are already in the pred_vars vector

  # Seeing which changes are available
  available_changes <- actual_index_available_interactions

  # All possible interactions are already inclued => Try to do the CHANGE_type
  if(length(available_changes)==0){
    warning("All possible variables were included in the terminal node, there is no possible change.")
    return(tree)
  }


  # This only replaces for a new MAIN effect. If I want to replace for an interaction?
  if(length(available_changes) > 1) {
    p_var <- sample(available_changes,size = 1)
  } else {
    p_var <- available_changes
  }


  # Setting the new pred_vars
  new_node_index <- c_node$pred_vars
  new_node_index[c_node$pred_vars %in% replace_interaction_selected] <- p_var
  new_node_index <- sort(unique(new_node_index))

  # cat("New Node indexes: ",new_node_index,"\n")


  # Calculating loglikelihood for the new changed nodes and the old ones
  c_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                               index_node = c_node$train_index,
                               j_ = c_node$pred_vars,
                               data = data)


  c_new_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                    index_node = c_node$train_index,
                                    j_ =  new_node_index,
                                    data = data)

  # Calculating the acceptance probability
  acceptance <- exp(-c_loglike + c_new_loglike)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1,min = 0,max = 1)<acceptance){

    #  Updating the variables accordingly
      tree[[c_node_name]]$pred_vars <- sort(unique(new_node_index))

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}



# ============
# Update Betas
# ============
updateBetas <- function(tree,
                        curr_part_res,
                        data,
                        trees_fit,
                        tree_number){


  # Getting the terminals
  t_nodes_names <- get_terminals(tree)


  # Getting the current prediction for that tree
  y_hat_train <- matrix(0,nrow = nrow(data$x_train),ncol = length(data$basis_subindex))
  y_hat_test <- matrix(0,nrow = nrow(data$x_test),ncol = length(data$basis_subindex))

  for(i in 1:length(t_nodes_names)){


    # Select the current terminal node
    cu_t <- tree[[t_nodes_names[i]]]

    # The lines above are summarised here
    node_index_var <- cu_t$pred_vars

    # Selecting the actually parameters subsetting
    basis_dim <- NCOL(data$P)
    basis_dim_interaction <- NCOL(data$P_interaction)
    n_leaf <- length(cu_t$train_index)
    diag_leaf <- diag(nrow = n_leaf)
    diag_basis <- diag(nrow = basis_dim)

    #

    #  Calculating the quantities need to the posterior of \beta
    # == Starting to iterate over those coefficients ==========#
    for(jj in 1:length(node_index_var)){


      leaf_basis_subindex <- unlist(data$basis_subindex[node_index_var[jj]]) # Recall to the unique() here too
      other_leaf_basis_subindex <- node_index_var[-jj]

      # RES_LEAF also need to updated here from the new_curr_part_res
      pred_minus_jj <- numeric(n_leaf)
      if(length(other_leaf_basis_subindex)!=0){
        for(iter_minus_jj in 1:length(other_leaf_basis_subindex)){
          current_betas <- matrix(tree[[t_nodes_names[i]]]$betas_vec[data$basis_subindex[[other_leaf_basis_subindex[iter_minus_jj]]]],nrow= 1)
          pred_minus_jj <- pred_minus_jj + tcrossprod(data$B_train[[other_leaf_basis_subindex[iter_minus_jj]]][cu_t$train_index,,drop=FALSE],current_betas)
        }

      }

      # old_betas <- matrix(tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex],nrow = 1) # This can be problematic

      res_leaf <- matrix(curr_part_res[cu_t$train_index], ncol=1) - pred_minus_jj

      # Adding a new step where I remove the mean of the res_leaf
      # res_leaf_mean <- mean(res_leaf)
      # res_leaf_norm <- res_leaf-res_leaf_mean

      # Getting the index for the vector of betas
      b_ <- crossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],res_leaf)


      if(node_index_var[jj]<=length(data$dummy_x$continuousVars)){
        U_ <- data$P*data$tau_beta[node_index_var[jj]]
      } else {
        U_ <- data$P_interaction*data$tau_beta[node_index_var[jj]]
      }

      Q_ <- (crossprod(data$B_train[[node_index_var[jj]]]) + U_)
      Q_inv_ <- chol2inv(chol(Q_))

      # Storing the old betas
      # See that I also creating a vector with the new betas
      new_betas <- mvnfast::rmvn(n = 1,mu = Q_inv_%*%b_,sigma = (data$tau^(-1))*Q_inv_)
      tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex] <- new_betas
      new_betas <- matrix(new_betas,nrow = 1)
      # Updating the residuals
      new_partial_pred <- tcrossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],new_betas) #+ res_leaf_mean
      # Need to update the trees fit!
      # trees_fit[tree_number,cu_t$train_index] <- trees_fit[tree_number,cu_t$train_index] - tcrossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],old_betas) + new_partial_pred

      y_hat_train[cu_t$train_index,node_index_var[jj]] <- new_partial_pred
      y_hat_test[cu_t$test_index,node_index_var[jj]] <- tcrossprod(data$B_test[[node_index_var[jj]]][cu_t$test_index,,drop=FALSE],new_betas) #+ res_leaf_mean
    }

  }


  # Getting tree predictions
  tree_hat_train <- rowSums(y_hat_train)
  tree_hat_test <- rowSums(y_hat_test)

  # Returning the tree
  return(list(tree = tree,
              y_hat_train = y_hat_train,
              tree_hat_train = tree_hat_train,
              tree_hat_test = tree_hat_test,
              y_hat_test = y_hat_test))

}


# Generating betas using try
sample_betas_try <- function(b_,Q_inv_, data){

  new_betas_ <- try(mvnfast::rmvn(n = 1,mu = data$tau*Q_inv_%*%b_,sigma = Q_inv_),silent = TRUE)

  if(inherits(new_betas_,'try-error')){
    new_betas_add <- mvnfast::rmvn(n = 1,mu = data$tau*Q_inv_%*%b_,sigma = Q_inv_ + diag(1e-6,nrow = length(b_)))
    return(new_betas_add)
  } else {
    return(new_betas_)
  }
}

updateBetasBigD <- function(tree,
                        curr_part_res,
                        data){


  # Getting the terminals
  t_nodes_names <- get_terminals(tree)

  # Getting the current prediction for that tree
  tree_hat_train <- numeric(nrow(data$x_train))
  tree_hat_test <- numeric(nrow(data$x_test))
  y_hat_train <- matrix(0, nrow = nrow(data$x_train), ncol = length(data$basis_subindex))

  for(i in 1:length(t_nodes_names)){

    # Generating the new betas_vec
    betas_vec_node = rep(0, length(unlist(data$basis_subindex)))

    # Select the current terminal node
    cu_t <- tree[[t_nodes_names[i]]]

    # The lines above are summarised here
    node_index_var <- cu_t$pred_vars
    res_leaf <- matrix(curr_part_res[cu_t$train_index], ncol=1)


    # Selecting the actually parameters subsetting
    leaf_basis_subindex <- unlist(data$basis_subindex[unique(node_index_var)])
    basis_dim <- NCOL(data$P)
    basis_dim_interaction <- NCOL(data$P_interaction)
    n_leaf <- length(cu_t$train_index)
    diag_leaf <- diag(nrow = n_leaf)
    diag_basis <- diag(nrow = basis_dim)

    # Need to define the subset for the D
    D_leaf <- data$D_train[cu_t$train_index, leaf_basis_subindex, drop = FALSE]

    # == Starting to iterate over those coefficients ==========#

    # Getting the index for the vector of betas
    b_ <- crossprod(D_leaf,res_leaf)
    U_ <- matrix(0, nrow = nrow(data$P), ncol = ncol(data$P))
    for(k in 1:length(unique(node_index_var))){
      aux_P_indexes <- unlist(data$basis_subindex[node_index_var[k]])
      U_[aux_P_indexes,aux_P_indexes] <- data$P[aux_P_indexes,aux_P_indexes]*(data$tau_beta[node_index_var[k]])
    }

    # Using another parametrization
    U_ <- U_[leaf_basis_subindex,leaf_basis_subindex, drop = FALSE]
    Q_ <- data$tau*crossprod(D_leaf) + U_
    Q_inv_ <- MASS::ginv(X = Q_)


    # Storing the old betas
    # See that I also creating a vector with the new betas

    new_betas <- sample_betas_try(b_ = b_,Q_inv_ = Q_inv_,data = data)
    betas_vec_node[leaf_basis_subindex] <- new_betas
    tree[[t_nodes_names[i]]]$betas_vec <- betas_vec_node
    new_betas <- matrix(new_betas,nrow = 1)
    # Updating the residuals
    new_partial_pred <- tcrossprod(D_leaf,new_betas)
    new_partial_pred_test <- tcrossprod(data$D_test[cu_t$test_index,leaf_basis_subindex, drop = FALSE], new_betas)
    # Need to update the trees fit!
    tree_hat_train[cu_t$train_index] <- new_partial_pred
    tree_hat_test[cu_t$test_index] <- new_partial_pred_test



    # To get the main effects
    for(jj in 1:length(cu_t$pred_vars)){
        y_hat_train[cu_t$train_index,cu_t$pred_vars[jj]] <- tcrossprod(data$D_train[cu_t$train_index,data$basis_subindex[[cu_t$pred_vars[jj]]], drop = FALSE], matrix(tree[[t_nodes_names[i]]]$betas_vec[data$basis_subindex[[tree[[t_nodes_names[i]]]$pred_vars[jj]]]],nrow = 1))
    }
  }


  # Doing some few experiments (check after some iterations how this behave)
  # main_effect_aux_ <-
  # plot(data$x_train[,main_effect_aux_], 0*normalize_bart(10*sin(pi*data$x_train[,1]),a = min_y, b = max_y),ylim = c(-0.5,0.1))
  # for(basis_n in 1:6){
  #   points(data$x_train[,main_effect_aux_,drop = FALSE], D_train[,data$basis_subindex[[main_effect_aux_]][basis_n],drop = FALSE]%*%new_betas[,data$basis_subindex[[main_effect_aux_]][basis_n]], col = basis_n)
  # }
  # points(data$x_train[,main_effect_aux_], y_hat_train[,main_effect_aux_], pch=  20)



  # Returning the tree
  return(list(tree = tree,
              y_hat_train = y_hat_train,
              tree_hat_train = tree_hat_train,
              tree_hat_test = tree_hat_test))

}





# =================
# Update \tau_betas_gammas from the prior
# =================
update_delta <- function(data){


  # Updating the robust delta
  for(i in 1:length(data$robust_delta)){
    # Pay attention here that (data$a_tau_beta_j is the same as 0.5*nu)
    data$robust_delta[i] <- rgamma(n = 1,
                                   shape = data$a_delta + 0.5*data$nu,
                                   rate = data$d_delta + 0.5*data$nu*data$tau_beta[i])
  }
  # Returning hte rhe robust delta
  return(data$robust_delta)
}

# =================
# Update \tau_betas
# =================
update_tau_betas_j <- function(forest,
                               data){


  # This is the same regardless the approach
  a_tau_beta <- data$a_tau_beta_j
  d_tau_beta <- data$d_tau_beta_j

  tau_b_shape <- 0.0
  tau_b_rate <- 0.0


  if(data$interaction_term){
    tau_b_shape <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
    tau_b_rate <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
    tau_beta_vec_aux <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
  } else{
    tau_b_shape <- numeric(NCOL(data$x_train))
    tau_b_rate <- numeric(NCOL(data$x_train))
    tau_beta_vec_aux_proposal <- tau_beta_vec_aux <- numeric(NCOL(data$x_train))
  }

  # Iterating over all trees
  for(i in 1:length(forest)){

    # Getting terminal nodes
    t_nodes_names <- get_terminals(forest[[i]])
    n_t_nodes <- length(t_nodes_names)

    # Iterating over the terminal nodes
    for(j in 1:length(t_nodes_names)){

      cu_t <- forest[[i]][[t_nodes_names[j]]]


      # All the information from var_ now is summarised inside the element from ht enode pred_vars
      var_ <- cu_t$pred_vars


      # Getting ht leaf basis
      for(kk in 1:length(var_)){


        leaf_basis_subindex <- unlist(data$basis_subindex[var_[kk]]) # Recall to the unique() function here
        p_ <- length(leaf_basis_subindex)
        betas_mat_ <- matrix(cu_t$betas_vec[leaf_basis_subindex],nrow = p_)
        tau_b_shape[var_[kk]] <- tau_b_shape[var_[kk]] + p_
        if(var_[[kk]] <= length(data$dummy_x$continuousVars)){
          tau_b_rate[var_[kk]] <- tau_b_rate[var_[kk]] + c(crossprod(betas_mat_,crossprod(data$P,betas_mat_)))
        } else {
          tau_b_rate[var_[kk]] <- tau_b_rate[var_[kk]] + c(crossprod(betas_mat_,crossprod(data$P_interaction,betas_mat_)))
        }
      }

    }


  }

  if(data$interaction_term){

      # Getting the sample from tau_beta
      for(j in 1:(NCOL(data$x_train)+NCOL(data$interaction_list)) ){

        # if( j== 8){
        #   stop("")
        # }
        if(length(d_tau_beta)>1){
          tau_beta_vec_aux_proposal <-rgamma(n = 1,
                                             shape = 0.5*tau_b_shape[j] + 0.5*data$nu,
                                             rate = 0.5*data$tau*tau_b_rate[j] + 0.5*data$nu*data$robust_delta[j])
        } else {
          tau_beta_vec_aux_proposal <- rgamma(n = 1,
                                              shape = 0.5*tau_b_shape[j] + 0.5*data$nu,
                                              rate = 0.5*data$tau*tau_b_rate[j] + 0.5*data$nu*data$robust_delta[j])

            # Just checking any error with tau_beta sampler
            if(tau_beta_vec_aux_proposal > 100){
              tau_beta_vec_aux_proposal <- 100
              warning("Warning: modified value for tau_beta to avoid numerical issues")
            }
        }

        tau_beta_vec_aux[j] <- tau_beta_vec_aux_proposal
      }


    } else { # This else is if there is no interaction

      for(j in 1:(NCOL(data$x_train)) ){


        if(length(d_tau_beta)>1){
          tau_beta_vec_aux_proposal <-rgamma(n = 1,
                                             shape = 0.5*tau_b_shape[j] + 0.5*data$nu,
                                             rate = 0.5*data$tau*tau_b_rate[j] + 0.5*data$nu*data$robust_delta[j])
        } else {
          tau_beta_vec_aux_proposal <- rgamma(n = 1,
                                              shape = 0.5*tau_b_shape[j] + 0.5*data$nu,
                                              rate = 0.5*data$tau*tau_b_rate[j] + 0.5*data$nu*data$robust_delta[j])

          # Just checking any error with tau_beta sampler
          if(tau_beta_vec_aux_proposal > 100){
            tau_beta_vec_aux_proposal <- 100
            warning("Warning: modified value for tau_beta to avoid numerical issues")
          }
        }

        tau_beta_vec_aux[j] <- tau_beta_vec_aux_proposal
      }
    }


  return(tau_beta_vec_aux)

}

# ===================
# Updating the \delta
# ===================

# A function to get predictions
getPredictions <- function(tree,
                           data){

  # Creating the vector to hold the values of the prediction
  if(data$interaction_term){
    y_hat <- matrix(0, nrow = nrow(data$x_train), ncol = NCOL(data$x_train)+NCOL(data$interaction_list))
    y_hat_test <- matrix(0,nrow(data$x_test), ncol = NCOL(data$x_test)+NCOL(data$interaction_list))
  } else {
    y_hat <- matrix(0, nrow = nrow(data$x_train), ncol = ncol(data$x_train))
    y_hat_test <- matrix(0,nrow(data$x_test), ncol = ncol(data$x_test))
  }

  # Getting terminal nodes
  t_nodes <- get_terminals(tree = tree)
  n_t_nodes <- length(t_nodes)

  for(i in 1:n_t_nodes){


    # Getting the current terminal node
    cu_t <- tree[[t_nodes[[i]]]]
    leaf_train_index <- cu_t$train_index
    leaf_test_index <- cu_t$test_index

    # Getting the variables used in the model
    node_index_var <- cu_t$pred_vars

    leaf_ancestors <- node_index_var # here isnt really the ancestors, but the variables that are being used

    leaf_basis_subindex <- data$basis_subindex[leaf_ancestors]

    # This test doesn't make sense anymore
    # # Test unit
    # if(length(leaf_ancestors)!=length(leaf_basis_subindex)){
    #   stop("Error on the getPredictions function")
    # }

    # Only add the marginal effects if the variables are within that terminal node
    if(length(leaf_basis_subindex)!=0){
      for(k in 1:length(leaf_basis_subindex)){

        y_hat[leaf_train_index,leaf_ancestors[k]] <- y_hat[leaf_train_index,leaf_ancestors[k]] + data$D_train[leaf_train_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]
        y_hat_test[leaf_test_index,leaf_ancestors[k]] <- y_hat_test[leaf_test_index,leaf_ancestors[k]] + data$D_test[leaf_test_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]

      }
    }

  }

  # Returning both training and test set predictions
  return(list(y_train_hat = y_hat,
              y_hat_test = y_hat_test))

}

# Updating tau
update_tau <- function(y_train_hat,
                       data,
                       forest){

  # Sampling a tau value
  n_ <- nrow(data$x_train)
  k_vec <- numeric(length = length(data$basis_subindex))
  j_vec <- numeric(length = length(data$basis_subindex))

  # Doing the parameter sum
  for(tree in 1:length(forest)){
    terminal_nodes_names <- get_terminals(forest[[tree]])
    for(curr_node in terminal_nodes_names) {
      pred_vars <- forest[[tree]][[curr_node]]$pred_vars

      for(k in pred_vars){
        k_vec[k] <- k_vec[k] + length(data$basis_subindex[[k]])
        # Getting the current betas
        current_betas <- forest[[tree]][[curr_node]]$betas_vec[data$basis_subindex[[k]]]
        if(k <= length(data$dummy_x$continuousVars)){
          j_vec[k] <- j_vec[k] + data$tau_beta[k]*current_betas%*%crossprod(data$P,current_betas)
        } else {
          j_vec[k] <- j_vec[k] + data$tau_beta[k]*current_betas%*%crossprod(data$P_interaction,current_betas)
        }
      }

    }
  }

  tau_sample <- stats::rgamma(n = 1,shape = 0.5*n_ + 0.5*sum(k_vec) + data$a_tau,
                              rate = 0.5*crossprod((data$y_train-y_train_hat))+ 0.5*sum(j_vec) + data$d_tau)

  return(tau_sample)

}


