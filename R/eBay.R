#' An method for differential abundance detection
#'
#'
#' @param otu.data  an OTU table with n rows (samples) and m columns (taxa)
#' @param group a n-vector of group indicators
#' @param cutf level of significance
#' @param test.method t-test or Wilcoxon rank sum test
#' @param adj.m the adjustment methods for p-values
#' @return final.p    the adjusted p values
#' @return dif.otus   the set of differentially abundant OTUs
#' @examples
#' ####generate data####
#' library(eBay)
#' set.seed(1)
#' rand_pi <- runif(20)
#' control_pi = case_pi = rand_pi/sum(rand_pi)
#' control_theta = case_theta = 0.1
#' group <- rep(c(0,1),each =20)
#' ntree_table <- simulation_dm(p=20,seed=1, N=20,control_pi, case_pi,
#' control_theta,case_theta)
#' #####differential abundance testing###
#' ebay.res <- eBay(otu.data=ntree_table, group=group, test.method="t",
#' cutf=0.05,adj.m="BH")
#' @export
#' @importFrom MGLM  MGLMfit
#' @importFrom MGLM  MGLMreg
#' @import stats
eBay=function(otu.data,group,cutf,test.methods,adj.m){
  library(MGLM)
  sample.s<- nrow(otu.data)
  otu.n <- ncol(otu.data)
  case.s <- length(which(group == 0))
  con.s <- length(which(group == 1))
  case <- which(group == 0)
  con <- which(group == 1)

  ntree_table<- otu.data
  taxa.p <- ncol(ntree_table)
  rownames(ntree_table) <- as.character(1:sample.s)

  resp <- as(ntree_table,"matrix")
  coe<-matrix(0,1,ncol(resp))
  B_e <- try(MGLMreg(resp~1, dist="DM")@coefficients, silent=TRUE)
  coe <- B_e
  gr <- matrix(rep(1,nrow(ntree_table)))
  alpha_e=exp(gr%*%coe)
  ntree_para_mat <- alpha_e

  exp_norm <- matrix(NA,ncol=taxa.p,nrow=sample.s)

  for (n in 1:sample.s) {
    exp_norm[n, ] = unlist((ntree_table[n, ] + ntree_para_mat[n,]) / (sum(ntree_table[n, ]) +sum(ntree_para_mat[n,])))
  }

  colnames(exp_norm)= colnames(ntree_table)

  exp_clr <- apply(exp_norm, 1, function(x){log2(x) - mean(log2(x))})
 if(test.methods=="t"){
  exp_test_t <- apply(exp_clr, 1, function(input) {
    t.test(input[case], input[con])$p.value
  })
  final.p.t <- p.adjust(exp_test_t,adj.m)
  dif.otus.t <- which( final.p.t< cutf)
  return(list(final.p=final.p.t, dif.otus=dif.otus.t))

 }
  if(test.methods=="wilcoxon"){
  exp_test_wil <- apply(exp_norm, 2, function(input) {
    wilcox.test(input[case], input[con])$p.value
  })
  final.p.wil <- p.adjust(exp_test_wil,adj.m)
  dif.otus.wil <-which(final.p.wil < cutf)
  return(list(final.p=final.p.wil,dif.otus=dif.otus.wil))
}
}
