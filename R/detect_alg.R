#' An algrithm for phylogeny-aware differential abundance detection
#'
#' @param p.val the p value for each node
#' @param tree_table a data set
#' @param tree a phylogenetic tree among m taxa
#' @param group  a n-vector of group indicators
#' @param test.method t-test or Wilcoxon rank sum test
#' @param cutf level of significance
#' @param adj.m the adjustment methods for p-values
#' @return final.p the adjusted p values
#' @return dif.otus  the set of differentially abundant OTUs
#' @examples
#' p <- 40
#' set.seed(1)
#' tree <- simulate_tree(p)
#' leaf_ind <- leaf_node_index(tree)
#' @export
#' @importFrom MGLM  MGLMfit
#' @import stats
#' @export
##################the algrithm for detecting differential abundant taxa########
detec_otus= function(p.val,tree_table,tree,group,test.method,cutf,adj.m){
  taxa.p <- length(tree$tip.label)
  sample.s <- nrow(tree_table)
  case <- which(group == 0)
  con <- which(group == 1)


  reject <- names(which(p.val<= cutf))
  child_index <- leaf_node_index(tree)

  if (length(reject)==0) {
    final.p <- p.adjust(p.val,adj.m)
    dif.otus <- NULL
  }

  if (length(reject)>=1) {
    edge_dif=lapply(1:taxa.p,function(x) {
      edge <- names(which(child_index[x,]==1))
      return(intersect(edge,reject))})
    site <- which(lengths(edge_dif)>=2)  ###a potential degenerate node
    nest_node <- sapply(edge_dif[site], function(v) return(min(v)))
    nest_node <- unique(nest_node)

    parent_nest <- tree$edge[match(nest_node,tree$edge[,2]),1]


    leaf_node <- lapply(parent_nest,function(x) {
      sub_node <- which(child_index[,parent_nest]==1)
      return(sub_node)})

    leaf_node <- unique(unlist(leaf_node))
    lone_node <- setdiff(reject,leaf_node)
    if(length(lone_node)==0){ p.val <- p.val}
    else{
      for (i in 1:length(lone_node)){
        x <- lone_node[i]
        if(x%in% 1: taxa.p){
          p.val <- p.val
        }else{
          leaf_sub <- which(child_index[1: taxa.p,x]==1)
          p.val[match(leaf_sub,names(p.val))]=p.val[match(x,names(p.val))]
        }
      }
    }
    if(length(parent_nest)==0){
      final.p  <- p.adjust(p.val[match(1:taxa.p,names(p.val))],adj.m)
      dif.otus <-  names(which(final.p<= cutf))
    }
    else{

      all_qian=lapply(parent_nest,function(x) {
        if(x%in%1:taxa.p){
          all_qian_sub=x
        }else{
          all_qian_sub <- which(child_index[1:taxa.p,x]==1)
        }
        return(all_qian_sub)})

      all_qian <- unique(unlist(all_qian))
      exp_test <- c()
      for(w in 1:length(parent_nest)){

        parent_nest_node <- which(child_index[1:taxa.p,parent_nest[w]]==1)
        tree_table_sub <- tree_table[,match(parent_nest_node,colnames(tree_table))]


        fit_glm <- MGLMfit(data.frame(tree_table_sub),dist = "DM")
        tree_para <- fit_glm@estimate



        tree_para_mat <- rbind(matrix(rep(tree_para, sample.s),byrow = TRUE,ncol = length(parent_nest_node)))
        exp_norm <- matrix(NA,ncol= length(parent_nest_node),nrow=sample.s)

        for (n in 1:sample.s) {
          exp_norm[n, ] = unlist((tree_table_sub[n, ] + tree_para_mat[n,]) / (sum(tree_table_sub[n, ]) +sum(tree_para_mat[n,])))
        }

        exp_clr <- apply(exp_norm, 1, function(x){log2(x) - mean(log2(x))})


        if (test.method == "t"){
          exp_test_sub <-  apply(exp_clr, 1, function(input){ t.test(input[case], input[con])$p.value})
        }
        if (test.method == "wilcoxon"){
          exp_test_sub <- apply(exp_norm, 2, function(input){ wilcox.test(input[case], input[con])$p.value})
        }

        names(exp_test_sub) <- as.character(parent_nest_node)
        exp_test <- append(exp_test,exp_test_sub)
      }

      left_node <- setdiff(1:taxa.p,all_qian)
      left_pval <- p.val[match(left_node,names(p.val))]
      exp_test_final <- exp_test[match(all_qian,names(exp_test))]
      final.p <- p.adjust(c(exp_test_final,left_pval),adj.m)
      dif.otus <- names(which(final.p < cutf))
    }
  }
  return(list(final.p=final.p,dif.otus=dif.otus))
}
