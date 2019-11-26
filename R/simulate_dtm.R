#' @title simulate data from DTM
#'
#' @description
#' @examples
#' @param p the number of taxa
#' @param seed a random seed
#' @param N sample size of each group
#' @param tree a phylogenetic tree
#' @param control_pi a probability vector for each taxon which were sum to 1 in control group
#' @param case_pi a probability vector for each taxon which were sum to 1 in case group
#' @param control_theta the over-dispersion parameter
#' @param case_theta the over-dispersion parameter
#' @return dtm_table the otu table
#' @export
simulation_dtm <- function(p,tree, seed = 1, N, tree,control_pi, case_pi,control_theta,case_theta){
  library(dirmult)
  dtm_table <- matrix(NA, 2*N, p + tree$Nnode - 1)
  colnames(dtm_table) <- as.character(tree$edge[, 2])
  for(j in (p + 1) : (p + tree$Nnode)){
    control_pi[which(tree$edge[, 1] == j)] -> pi0
    case_pi[which(tree$edge[, 1] == j)] -> pi1

    for(i in 1:N){
      if(j == (p+1)){
        set.seed(i*j*seed)
        parent1 <- round(runif(1, 5000,50000))
        set.seed(i*j*(seed+1))
        parent2 <- round(runif(1, 5000,50000))
      }else{
        parent1 <- dtm_table[i, which(tree$edge[, 2] == j)]
        parent2 <- dtm_table[i + N, which(tree$edge[, 2] == j)]
      }
      set.seed(i*j*seed)

      dtm_table[i, which(tree$edge[, 1] == j)] <- simPop(J=1, K=2, n=parent1, pi = pi0, theta = control_theta[j-p])$data
      set.seed(i*j*(seed+1))

      dtm_table[i+N, which(tree$edge[, 1] == j)] <- simPop(J=1, K=2, n=parent2, pi = pi1, theta =case_theta[j-p])$data
    }
  }

  return(dtm_table)
}
