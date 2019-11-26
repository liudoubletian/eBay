#' @title simulate data from DTM
#'
#' @description
#' @examples
#' @param n the number of taxa
#' @export
simulate_tree = function (n) {
  foo <- function(n, pos) {
    n1 <- sample.int(n - 1, 1, FALSE, NULL)
    n2 <- n - n1
    po2 <- pos + 2 * n1 - 1
    edge[c(pos, po2), 1] <<- nod
    nod <<- nod + 1L
    if (n1 > 2) {
      edge[pos, 2] <<- nod
      foo(n1, pos + 1)
    }
    else if (n1 == 2) {
      edge[c(pos + 1, pos + 2), 1] <<- edge[pos, 2] <<- nod
      nod <<- nod + 1L
    }
    if (n2 > 2) {
      edge[po2, 2] <<- nod
      foo(n2, po2 + 1)
    }
    else if (n2 == 2) {
      edge[c(po2 + 1, po2 + 2), 1] <<- edge[po2, 2] <<- nod
      nod <<- nod + 1L
    }
  }
  nbr <- 2 * n - 3 + 1
  edge <- matrix(NA, nbr, 2)
  n <- as.integer(n)
  nod <- n + 1L
  foo(n, 1)
  i <- which(is.na(edge[, 2]))
  edge[i, 2] <- 1:n
  phy <- list(edge = edge)
  phy$tip.label <- as.character(1: n)
  phy$node.label <- as.character((n+1):(2*n-1))
  phy$edge.length <- rep(1, length.out = nbr)
  phy$Nnode <- n - 2L + 1
  class(phy) <- "phylo"
  attr(phy, "order") <- "cladewise"
  phy
}
