#' @title An method for phylogeny-aware differential abundance detection
#'
#' @description
#' @param otu.data the data set where each row was a sample and each column was a taxa,
#' @param tree a phylogenetic tree for otus
#' @param test.method a test including "t" and "wilcoxon"
#' @param cutf the significance level
#' @return dif.node the detected differential nodes
#' @return p.value the p value for each taxa
#' @export
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

inter_func <- function(g,tree,inter_node,tree_table,case,con,test.method){
  library(MGLM)
  sample.s <- nrow(tree_table)

  sam_glm = sam_glm_wil = c()
  sam_glm_ubay = sam_glm_wil_ubay = c()


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
eBay_tree=function(otu.data,tree,group,test.method,cutf){

  otu.data <- otu.data
  sample.s<- nrow(otu.data)
  otu.n <- ncol(otu.data)
  case.s <- length(which(group == 0))
  con.s <- length(which(group == 1))
  case <- which(group == 0)
  con <- which(group == 1)


  phy_data <- phyloseq(otu_table(as.matrix(t(otu.data)),taxa_are_rows=TRUE), phy_tree(tree))##transfrom data to phyloseq structure
  minobs<- 0
  prevalence <- apply(as(otu_table(phy_data), "matrix"), 1, function(x, minobs) {
    return(sum(x > minobs))
  }, minobs)/sample.s

  kep_pre1 <- apply(as(otu_table(phy_data)[,case], "matrix"), 1, function(x) {
    return(sum(x > 0))
  })/case.s

  kep_pre2 <- apply(as(otu_table(phy_data)[,con], "matrix"), 1, function(x) {
    return(sum(x > 0))
  })/con.s

  keepOTUs <-  prevalence > 0.2 & taxa_sums(t(phy_data)) > (0.5 *sample.s) & kep_pre1>0 & kep_pre2>0
  ###filtering the otus


  pru_data <- prune_taxa(keepOTUs, phy_data) ###prune the data set

  tree_table<- t(otu_table(pru_data))### final otu table
  colnames(tree_table) <- as.character(1:taxa.p)
  taxa.p <- ncol(tree_table)

  pru_tree <- phy_tree(pru_data)
  taxa_index <- otu_index(pru_tree)

  taxa_table <- tree_table %*% taxa_index
  total_table <- cbind(taxa_table, tree_table)
  tree_table <- total_table[,-1]

  inter_node <- unique(pru_tree$edge[,1]) ##internal nodes of the tree
  leaf=c()
  for (i in 1:length(inter_node)){
    leaf_node=pru_tree$edge[which(pru_tree$edge[,1]==inter_node[i]),2]
    leaf=append(leaf,leaf_node)
  } ### the child node of each internal node



  p_table <- foreach(g=1:length(inter_node),.combine='cbind') %dopar% inter_func(g=g,tree=pru_tree,inter_node=inter_node,tree_table=tree_table,case=case,con=con,test.method)
###obtain the p value of each internal node
  p.vector <- p_table[1,]

  names(p.vector) <- as.character(leaf)


  detec_result=detec_otus(p.val=p.vector,tree_table=tree_table,tree=pru_tree,group,test_method,cutf)


  return(detec_result)
}








