#' An method for differential abundance detection
#'
#' @param cutf the significance level
#' @param otu.data the data set where each row was a sample and each column was a taxa,
#' @param test.method a test including "t" and "wilcoxon"
#' @return dif.otus   the detected differential otus
#' @return final.p    the p value for each otu
#' @examples
#' library(eBay)
#' set.seed(1)
#' rand_pi <- runif(20)
#' control_pi = case_pi = rand_pi/sum(rand_pi)
#' control_theta = case_theta = 0.1
#' group <- rep(c(0,1),each =20)
#' ntree_table <- simulation_dm(p=20,seed=1, N=20,control_pi, case_pi,
#' control_theta,case_theta)
#' ebay.res <- eBay(otu.data=ntree_table, group=group, test.method="t", cutf=0.05)
#' @export
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

  fit_glm <- MGLMfit(ntree_table, dist = "DM")  ####estimate the parameter alpha
  ntree_para <- fit_glm@estimate

  ntree_para_mat <- matrix(rep(ntree_para, sample.s),byrow = TRUE,ncol = taxa.p)

  exp_norm <- matrix(NA,ncol=taxa.p,nrow=sample.s)
  exp_ubay <- matrix(NA,ncol=taxa.p,nrow=sample.s)

  for (n in 1:sample.s) {
    exp_norm[n, ] = unlist((ntree_table[n, ] + ntree_para_mat[n,]) / (sum(ntree_table[n, ]) +sum(ntree_para_mat[n,])))
    exp_ubay[n, ] = unlist((ntree_table[n, ] +0.5) / (sum(ntree_table[n, ]) +0.5*taxa.p))

  }

  colnames(exp_norm)=colnames(exp_ubay) = colnames(ntree_table)

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
  exp_test_wil <- apply(exp_clr, 1, function(input) {
    wilcox.test(input[case], input[con])$p.value
  })
  final.p.wil <- p.adjust(exp_test_wil,adj.m)
  dif.otus.wil <-which(final.p.wil < cutf)
  return(list(final.p=final.p.wil,dif.otus=dif.otus.wil))
}
}
