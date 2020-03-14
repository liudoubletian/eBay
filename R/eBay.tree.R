#' @title An method for phylogeny-aware differential abundance detection
#' the eBay incorporates the tree into the normalization process and implement
#' differential abuance testing.
#' @param otu.data the data set where each row was a sample and each column was a taxa,
#' @param tree a phylogenetic tree for otus
#' @param test.method a test including "t" and "wilcoxon"
#' @param cutf the significance level
#' @return dif.node the detected differential nodes
#' @return p.value the p value for each taxa
#' @examples
#' p <- 40
#' set.seed(1)
#' tree <- simulate_tree(p)
#' set.seed(1)
#' control_pi = case_pi = c()
#' for(j in (p+1):(p+tree$Nnode)){
#' set.seed(j)
#' random_pi <- runif(1,0.2,0.4)
#' control_pi[which(tree$edge[,1]==j)] <- c(random_pi, 1-random_pi)
#' case_pi[which(tree$edge[,1]==j)] <- c(random_pi, 1-random_pi)
#' }
#' control_theta = case_theta = rep(0.1, tree$Nnode)
#' group <- rep(c(0,1),each =20)
#' tree_table <- simulation_dtm(p=40,tree, seed=1, N=20,control_pi, case_pi,control_theta,case_theta)
#' ebay_tree.res <- eBay_tree(otu.data=tree_table, tree=tree, group=group, test.method="t", cutf=0.05)
#' @export
#' @importFrom MGLM  MGLMfit
#' @import foreach
#' @import stats
#' @import doParallel
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
#' @export
#### estimate the parameter alpha on each internal node and conduct normalization####
inter_func=function(g,tree,inter_node,tree_table,case,con,test.method){
  sample.s <- nrow(tree_table)
  sam_glm=sam_glm_wil=c()
  sam_glm_ubay=sam_glm_wil_ubay=c()


  lit_tree <- tree$edge[which(tree$edge[,1]==inter_node[g]),2]
  lit_tree_table <- tree_table[,match(lit_tree,colnames(tree_table))]
  colnames(lit_tree_table) <- as.character(lit_tree)
  taxa.p <- ncol(lit_tree_table)

  fit_glm <- MGLMfit(data.frame(lit_tree_table),dist = "DM")
  tree_para <- fit_glm@estimate

  tree_para_mat <- rbind(matrix(rep(tree_para, sample.s),byrow = TRUE,ncol = taxa.p))

  exp_norm <- matrix(NA,ncol=taxa.p,nrow=sample.s)

  for (n in 1:sample.s) {
    exp_norm[n, ] = unlist((lit_tree_table[n, ] + tree_para_mat[n,]) / (sum(lit_tree_table[n, ]) +sum(tree_para_mat[n,])))
  }

  exp_clr <- apply(exp_norm, 1, function(x){log2(x) - mean(log2(x))})



  if (test.method == "t"){
    exp_test <-  apply(exp_clr, 1, function(input){ t.test(input[case], input[con])$p.value})
    final.p <- exp_test
  }
  if (test.method == "wilcoxon"){
    exp_test <- apply(exp_clr, 1, function(input){ wilcox.test(input[case], input[con])$p.value})
    final.p <- exp_test
  }

  return(rbind(final.p,exp_norm))
}
#' @export
eBay_tree =function(otu.data,tree,group,test.method,cutf,adj.m){
  otu.data <- otu.data
  sample.s<- nrow(otu.data)
  otu.n <- ncol(otu.data)
  case.s <- length(which(group == 0))
  con.s <- length(which(group == 1))
  case <- which(group == 0)
  con <- which(group == 1)


  tree_table<- otu.data### final otu table
  taxa.p <- ncol(tree_table)
  colnames(tree_table) <- as.character(1:taxa.p)


  pru_tree <- tree
  taxa_index <- otu_index(pru_tree)
  taxa_table <- as.matrix(tree_table) %*% as.matrix(taxa_index)
  total_table <- cbind(taxa_table, tree_table)
  tree_table <- total_table[,-1]

  inter_node <- unique(pru_tree$edge[,1]) ##internal nodes of the tree
  leaf=c()
  for (i in 1:length(inter_node)){
    leaf_node=pru_tree$edge[which(pru_tree$edge[,1]==inter_node[i]),2]
    leaf=append(leaf,leaf_node)
  } ### the child node of each internal node
  p_table <- foreach(g=1:length(inter_node),.combine='cbind') %dopar% inter_func(g=g,tree=pru_tree,inter_node=inter_node,tree_table=tree_table,case=case,con=con,test.method=test.method)
  ###obtain the p value of each internal node
  p.vector <- p_table[1,]

  names(p.vector) <- as.character(leaf)


  detec_result=detec_otus(p.val=p.vector,tree_table=tree_table,tree=pru_tree,group=group,test.method=test.method,cutf=cutf,adj.m=adj.m)


  return(detec_result)
}







