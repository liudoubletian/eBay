#' @export
#### estimate the parameter alpha on each internal node and conduct normalization####
inter_func=function(g,tree,tree_table,group,test.method){
  sample.s <- nrow(tree_table)
  case <- which(group == 0)
  con <- which(group == 1)
  inter_node <- unique(tree$edge[,1])
  sam_glm=sam_glm_wil=c()
  sam_glm_ubay=sam_glm_wil_ubay=c()


  lit_tree <- tree$edge[which(tree$edge[,1]==inter_node[g]),2]
  lit_tree_table <- tree_table[,match(lit_tree,colnames(tree_table))]
  colnames(lit_tree_table) <- as.character(lit_tree)
  taxa.p <- ncol(lit_tree_table)

  resp <- as(lit_tree_table,"matrix")
    coe<-matrix(0,1,ncol(resp))
    B_e <- try(MGLMreg(resp~1, dist="DM")@coefficients, silent=TRUE)
    if(inherits(B_e,"try-error")){
    B_e <- MGLMfit(lit_tree_table, dist="DM")@estimate
    tree_para_mat <- rbind(matrix(rep(B_e, sample.s),byrow = TRUE,ncol = ncol(lit_tree_table)))
    }
    else{
      coe=B_e
      gr <- matrix(rep(1,nrow(lit_tree_table)))
      alpha_e=exp(gr%*%coe)
      tree_para_mat <- alpha_e
    }

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
    exp_test <- apply(exp_norm, 2, function(input){ wilcox.test(input[case], input[con])$p.value})
    final.p <- exp_test
  }

  return(rbind(final.p,exp_norm))
}
