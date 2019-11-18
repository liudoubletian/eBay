DTM_simulation_sam=function(p,seed = 1, N, tree,  control_pi, case_pi,control_theta,case_theta,m){
  dtm_alltab = matrix(NA, 2*N, p + tree$Nnode - 1)
  colnames(dtm_alltab) <- as.character(tree$edge[, 2])
  for(j in (p + 1) : (p + tree$Nnode)){
    control_pi[which(tree$edge[, 1] == j)] -> pi0
    case_pi[which(tree$edge[, 1] == j)] -> pi1
    
    
    #theta0 = ifelse(sum(taxa_index[, colnames(taxa_index) == as.character(j)]) > p/5, 1e-8, theta)
    for(i in 1:N){
      if(j == (p+1)){
        set.seed(i*j*m*seed)
        parent1 = round(runif(1, 5000,50000))
        set.seed(i*j*m*(seed+1))
        parent2 =round(runif(1, 5000,50000))
      }
      else{
        parent1 = dtm_alltab[i, which(tree$edge[, 2] == j)]
        parent2 = dtm_alltab[i + N, which(tree$edge[, 2] == j)]
      }
      set.seed(i*j*m*seed)
      
      dtm_alltab[i, which(tree$edge[, 1] == j)] <- simPop(J=1, K=2, n=parent1, pi = pi0, theta = control_theta[j-p])$data
      set.seed(i*j*m*(seed+1))
      
      dtm_alltab[i+N, which(tree$edge[, 1] == j)] <- simPop(J=1, K=2, n=parent2, pi = pi1, theta =case_theta[j-p])$data
    }
  }
  return(dtm_alltab)
}  
Tree_Simulate = function (n) {
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
inter_func=function(g=g,tree=tree,inter_node=inter_node,tree_table=tree_table,case_samp,con_samp,samps,leaf=leaf){
  library(MGLM)
  sam_glm=sam_glm_wil=c()
  sam_glm_ubay=sam_glm_wil_ubay=c()
  
  
  lit_tree=tree$edge[which(tree$edge[,1]==inter_node[g]),2]
  lit_tree_table=tree_table[,match(lit_tree,colnames(tree_table))]
  colnames(lit_tree_table)=as.character(lit_tree)
  
  
  
  fit_glm_con=MGLMfit(data.frame(lit_tree_table),dist = "DM")
  alpha_con_glm=fit_glm_con@estimate
  
 
  
  pi_exp_glm=matrix(rep(0,samps*length(lit_tree)),ncol=length(lit_tree))
  pi_exp_glm_ubay=matrix(rep(0,samps*length(lit_tree)),ncol=length(lit_tree))
  
  
  alpha_matrix_glm = rbind(matrix(rep(alpha_con_glm, samps),byrow = TRUE,ncol = length(lit_tree)))
  
  
  for (n in 1:samps) {
    pi_exp_glm[n,]= unlist((lit_tree_table[n,] + alpha_matrix_glm[n,]) / (sum(lit_tree_table[n, ]) + sum(alpha_matrix_glm[n,])))
    pi_exp_glm_ubay[n,]= unlist((lit_tree_table[n,] + 0.5) / (sum(lit_tree_table[n, ])+sum(rep(1/2,length(lit_tree)))))
  }
  
  
  cij_expi_glm=apply(t(pi_exp_glm), 2, function(x){log2(x) - mean(log2(x))})
  cij_expi_glm_ubay=apply(t(pi_exp_glm_ubay), 2, function(x){log2(x) - mean(log2(x))})
  
  for( t in 1:length(lit_tree)){
    
    sam_glm[t]= t.test(cij_expi_glm[t,1:case_samp], cij_expi_glm[t,(case_samp+1):samps])$p.value
    sam_glm_wil[t]= wilcox.test(cij_expi_glm[t,1:case_samp], cij_expi_glm[t,(case_samp+1):samps])$p.value
    
    sam_glm_ubay[t]= t.test(cij_expi_glm_ubay[t,1:case_samp], cij_expi_glm_ubay[t,(case_samp+1):samps])$p.value
    sam_glm_wil_ubay[t]= wilcox.test(cij_expi_glm_ubay[t,1:case_samp], cij_expi_glm_ubay[t,(case_samp+1):samps])$p.value
  }
  return(rbind(sam_glm,sam_glm_wil,sam_glm_ubay,sam_glm_wil_ubay,pi_exp_glm,pi_exp_glm_ubay))
}
pw_fdr = function(p.adj, dif = dif, cutf = 0.05){
  reject = which(p.adj <= cutf)
  TP = length(intersect(reject, dif))
  FP = length(reject) - TP
  power = TP/length(dif)
  fdr = ifelse(length(reject) == 0, 0, FP/length(reject))
  return(c(power, fdr))
}


Taxa.index=function (p, phy_tree) 
{
  taxa_index <- matrix(0, p, phy_tree$Nnode)
  colnames(taxa_index) <- as.character((p + 1):(p + phy_tree$Nnode))
  for (i in 1:p) {
    temp <- i
    for (j in 1:1000) {
      index <- which(phy_tree$edge[, 2] == temp)
      if (length(index) != 0) {
        parent <- phy_tree$edge[index, 1]
        taxa_index[i, which(colnames(taxa_index) == as.character(parent))] <- 1
        temp <- parent
      }
      else {
        break
      }
    }
  }
  return(taxa_index)
}