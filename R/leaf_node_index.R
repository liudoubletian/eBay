#' @export
###an indicator table that represents the leaf nodes of each internal node##
leaf_node_index=function (tree){
  p <- length(tree$tip.label)
  total.d <- p+tree$Nnode
  otu_mat <- matrix(0, total.d,total.d)
  colnames(otu_mat) <- as.character(1:total.d)
  rownames(otu_mat) <- as.character(1:total.d)

  for (i in 1:total.d) {
    node <- i
    while(node %in% tree$edge[, 2]) {
      index <- which(tree$edge[, 2] == node)
      parent <- tree$edge[index, 1]
      otu_mat[i, which(colnames(otu_mat) == as.character(parent))] <- 1
      node <- parent
    }
  }
  return(otu_mat)
}


