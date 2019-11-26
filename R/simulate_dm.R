#' @title simulate data from DM
#'
#' @description
#' @examples
#' @param p the number of taxa
#' @param seed a random seed
#' @param N sample size of each group
#' @param control_pi a probability vector for each taxon which were sum to 1 in control group
#' @param case_pi a probability vector for each taxon which were sum to 1 in case group
#' @param control_theta the over-dispersion parameter
#' @param case_theta the over-dispersion parameter
#' @return dtm_table the otu table
#' @export
simulation_dm=function(p,seed, N,control_pi, case_pi,control_theta,case_theta){
  library(dirmult)
  ntree_table = matrix(NA, 2*N, p)
  for(i in 1:N){
  set.seed(i*seed)
  parent1 = round(runif(1, 5000, 50000))
  parent2 = round(runif(1, 5000, 50000))
  set.seed(i*seed)
  ntree_table[i, ] <- simPop(J=1, K=p, n=parent1, pi = control_pi, theta = control_theta)$data
  set.seed(i*(seed+1))
  ntree_table[i+N,] <- simPop(J=1, K=p, n=parent2, pi = case_pi, theta = case_theta)$data
  }
colnames(ntree_table)=as.character(1:p)
return(ntree_table)
}

