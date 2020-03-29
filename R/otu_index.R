#' @export
###an indicator table that represents the parent nodes of each leaf node##
otu_index=function (tree){
  p <- length(tree$tip.label)
  otu_mat <- matrix(0, p,tree$Nnode)
  colnames(otu_mat) <- as.character((p + 1):(p + tree$Nnode))
  rownames(otu_mat) <- as.character(1:p)

  for (i in 1:p) {
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
