library(dbarts)
library(dplyr)
library(purrr)
library(mlbench)

 # -------------------------------------------------------------------------
myData <- mlbench.friedman1(n = 250,sd = 1) %>% as.data.frame()

 # Create simple dbarts model:
set.seed(99)

bartFit <- model <- bart(myData[,1:10],
                               myData[,11],
                               ntree = 20,
                               keeptrees = TRUE,
                               nskip = 500,
                               ndpost = 1500)

# -------------------------------------------------------------------------
# Initial Setup
# -------------
# Get the total number of trees and iterations from the model
treesTotal <- model$call$ntree
iteration  <- model$call$ndpost
#
# # Extract trees based on the total number of trees and iterations
trees <- model$fit$getTrees(treeNums = 1:treesTotal, sampleNums = 1:iteration)
# trees$isLeaf <- ifelse(trees$var == -1,TRUE,FALSE)
# trees <- trees %>% group_by(sample,tree) %>% mutate(node = dplyr::row_number()) %>% mutate(ancestors = "")
# trees
#
#
#
# for(i in 1:nrow(trees)){
#
#   curr_tree <- trees$tree[i]
#   curr_sample <- trees$sample[i]
#   curr_tree <- trees[trees$sample==curr_sample & trees$tree==curr_tree,]
#
#   ancestors <- character(0)
#
#   for(j in 1:nrow(curr_tree)){
#
#     if(curr_tree$isLeaf[j]){
#
#       # Checking if it's a left node or a right node
#       if(curr_tree$isLeaf[j-1]){
#         ancestors <- c(ancestors,curr_tree$var[j-2])
#       }
#
#       if(!curr_tree$isLeaf[j-1]){
#         ancestors <- c(ancestors,curr_tree$var[j-1])
#       }
#     }
#
#     trees$ancestors[i] <- paste0(ancestors,collapse = "-")
#   }
#
#
#   }
# }


mapOverNodes <- function(tree, f, ...) {
  mapOverNodesRecurse <- function(tree, depth, f, ...) {
    node <- list(
      value = tree$value[1],
      n = tree$n[1],
      depth = depth
    )
    if (tree$var[1] == -1) {
      node$n_nodes <- 1
      node$f.x <- f(node, ...)
      return(node)
    }
    node$var <- tree$var[1]
    node$f.x <- f(node, ...)
    headOfLeftBranch <- tree[-1,]
    left <- mapOverNodesRecurse(headOfLeftBranch, depth + 1, f, ...)
    n_nodes.left <- left$n_nodes
    left$n_nodes <- NULL
    node$left <- left
    headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
    right <- mapOverNodesRecurse(headOfRightBranch, depth + 1, f, ...)
    n_nodes.right <- right$n_nodes
    right$n_nodes <- NULL
    node$right <- right
    node$n_nodes <- 1 + n_nodes.left + n_nodes.right
    return(node)
  }
  result <- mapOverNodesRecurse(tree, 1, f, ...)
  result$n_nodes <- NULL
  return(result)
}

observeInteractions <- function(node, ...) {
  if (is.null(node$var)) return(NULL)
  interactionData <- list(...)$interactionData
  # Make the current node visibile inside the environment.
  interactionData$node <- node
  with(interactionData, {
    if (node$depth <= currentDepth) {
      # If true, we have backtracked to go down the right branch, so we
      # remove the variables from the left branch.
      currentVariables <- currentVariables[seq_len(node$depth - 1)]
    }
    if (length(interactionData$currentVariables) > 0) {
      # This is a brute-force way of updating the following indices,
      # relying on the column-major storage order that R uses:
      # hasInteraction[currentVariables,,drop = FALSE][,node$var]
      updateIndices <- currentVariables +
        (node$var - 1) * nrow(hasInteraction)
      hasInteraction[updateIndices] <- TRUE
    }
    currentVariables <- c(currentVariables, node$var)
    currentDepth <- node$depth
  })
  rm("node", envir = interactionData)
  # Since the function is used for its side effects, there isn't a return
  # value.
  return(NULL)
}


numVariables <- ncol(bartFit$fit$data@x)
variableNames <- colnames(bartFit$fit$data@x)



interaction_table <- matrix(0,nrow = NCOL(x_train),ncol = NCOL(x_train))

for(jj in 1:max(trees$sample)){

  for(ii in 1:max(trees$tree)){


      treeOfInterest <- subset(trees,  sample == jj & tree == ii)

      # Define this as an environment as they are mutable
      interactionData <- list2env(list(
        currentDepth = 0,
        currentVariables = integer(),
        hasInteraction = matrix(
          data = FALSE,
          ncol = numVariables, nrow = numVariables,
          dimnames = list(ancestor = variableNames, descendant = variableNames)
        )
      ))

      invisible(mapOverNodes(
        tree = treeOfInterest,
        observeInteractions,
        interactionData = interactionData
      ))

      interaction_table <- interaction_table + interactionData$hasInteraction

  }

  print(jj)

}

# Creating the vector with names
# norm_interaction_table <- interaction_table/(max(trees$sample)*max(trees$tree))
norm_interaction_table <- interaction_table/(sum(norm_interaction_table))

lower_tri <- norm_interaction_table[lower.tri(norm_interaction_table)]
upper_tri <- norm_interaction_table[upper.tri(norm_interaction_table)]

lower_triangular_names <- combn(rownames(norm_interaction_table), 2, FUN = function(x) paste(x, collapse = "-"))

all_proportion <- lower_tri+upper_tri
names(all_proportion) <- lower_triangular_names
sort(all_proportion,decreasing = TRUE)
plot(all_proportion, xaxt = "n", xlab = "",pch=20, main = "BART interaction proportion of ancestors pairs in Friedman")
axis(1,at = 1:45,labels = names(all_proportion),las = 2)
