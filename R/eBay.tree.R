#' @title An method for phylogeny-aware differential abundance detection
#' the eBay incorporates the tree into the normalization process and implement
#' differential abuance testing.
#' @param otu.data an OTU table with n rows (samples) and m columns (taxa)
#' @param tree a phylogenetic tree among m taxa
#' @param group  a n-vector of group indicators
#' @param test.method  t-test or Wilcoxon rank sum test
#' @param cutf  level of significance
#' @param adj.m the adjustment methods for p-values
#' @return final.p the adjusted p values
#' @return dif.otus  the set of differentially abundant OTUs
#' @examples
#' ####generate data####
#' p <- 40
#' set.seed(1)
#' tree <- simulate_tree(p)
#' set.seed(1)
#' control_pi = case_pi = c()
#' for(j in (p+1):(p+tree$Nnode)){
#' set.seed(j)
#' random_pi <- runif(1,0.2,0.4)
#' control_pi[which(tree$edge[,1]==j)] <- c(random_pi, 1-random_pi)
#' case_pi[which(tree$edge[,1]==j)] <- c(random_pi, 1-random_pi)}
#' control_theta = case_theta = rep(0.1, tree$Nnode)
#' group <- rep(c(0,1),each =20)
#' tree_table <- simulation_dtm(p=40,tree, seed=1, N=20,control_pi, case_pi,
#' control_theta,case_theta)
#' #####differential abundance testing based on a phylogenetic tree###
#' ebay_tree.res <- eBay_tree(otu.data=tree_table, tree=tree, group=group,
#' test.method="t", cutf=0.05,adj.m="BH")
#' @export
#' @importFrom MGLM  MGLMfit
#' @import foreach
#' @import stats
#' @import doParallel
#' @export
eBay_tree =function(otu.data,tree,group,test.method,cutf,adj.m){
  otu.data <- otu.data
  sample.s<- nrow(otu.data)
  otu.n <- ncol(otu.data)
  case.s <- length(which(group == 0))
  con.s <- length(which(group == 1))
  case <- which(group == 0)
  con <- which(group == 1)


  taxa.p <- ncol(otu.data)
  colnames(otu.data) <- as.character(1:taxa.p)


  pru_tree <- tree
  taxa_index <- otu_index(pru_tree)
  taxa_table <- as.matrix(otu.data) %*% as.matrix(taxa_index)
  total_table <- cbind(taxa_table, otu.data)
  tree_table <- total_table[,-1]

  inter_node <- unique(pru_tree$edge[,1]) ##internal nodes of the tree
  leaf=c()
  for (i in 1:length(inter_node)){
    leaf_node=pru_tree$edge[which(pru_tree$edge[,1]==inter_node[i]),2]
    leaf=append(leaf,leaf_node)
  } ### the child node of each internal node
  p_table <- foreach(g=1:length(inter_node),.combine='cbind') %dopar% inter_func(g=g,tree=pru_tree,tree_table=tree_table, group,test.method=test.method)
  ###obtain the p value of each internal node
  p.vector <-  p.adjust(p_table[1,],"BH")

  names(p.vector) <- as.character(leaf)


  detec_result=detec_otus(p.val=p.vector,tree_table=tree_table,tree=pru_tree,group=group,test.method=test.method,cutf=cutf,adj.m=adj.m)


  return(detec_result)
}







