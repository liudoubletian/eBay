#' An method for differential abundance detection
#'
#' @param cutf the significance level
#' @param otu.data the data set where each row was a sample and each column was a taxa,
#' @param test.method a test including "t" and "wilcoxon"
#' @return dif.otus the detected differential otus
#' @return final.p the p value for each otu
#' @examples
#' set.seed(1)  
#' rand_pi <- runif(20)   
#' control_pi = case_pi = rand_pi/sum(rand_pi)   
#' control_theta = case_theta = 0.1  
#' group <- rep(c(0,1),each =20)  
#' ntree_table <- simulation_dm(p=20,seed=1, N=20,control_pi, case_pi,control_theta,case_theta) 
#' ebay.res <- eBay(otu.data=ntree_table, group=group, test.method="t", cutf=0.05)  
#' @export
eBay=function(otu.data,group,cutf){

  sample.s<- nrow(otu.data)
  otu.n <- ncol(otu.data)
  case.s <- length(which(group == 0))
  con.s <- length(which(group == 1))
  case <- which(group == 0)
  con <- which(group == 1)


  ntree_table<- otu.data
  
  taxa.p <- ncol(ntree_table)

  rownames(ntree_table) <- as.character(1:sample.s)

  fit_glm <- MGLMfit(ntree_table, dist = "DM") ######using MGLMfit to estimate the parameter alpha

  ntree_para <- fit_glm@estimate

  ntree_para_mat <- rbind(matrix(rep(ntree_para, sample.s),byrow = TRUE,ncol = taxa.p))

  exp_norm <- matrix(NA,ncol=taxa.p,nrow=sample.s)

  for (n in 1:sample.s) {
    exp_norm[n, ] = unlist((ntree_table[n, ] + ntree_para_mat[n,]) / (sum(ntree_table[n, ]) +sum(ntree_para_mat[n,])))
  }

  colnames(exp_norm) <- colnames(ntree_table)

  exp_clr <- apply(exp_norm, 1, function(x){log2(x) - mean(log2(x))}) ###CLR transformation

  exp_test_t <- apply(exp_clr, 1, function(input) {
      t.test(input[case], input[con])$p.value
    })
  final.p.t <- p.adjust(exp_test_t,"BH")
  dif.otus.t <- names(which(final.p.t < cutf))
 
  exp_test_wil <- apply(exp_clr, 1, function(input) {
      wilcox.test(input[case], input[con])$p.value
    })
  final.p.wil <- p.adjust(exp_test_wil,"BH")
  dif.otus.wil <- names(which(final.p.wil < cutf))

  return(list(final.p.t=final.p.t, dif.otus.t=dif.otus.t,
              final.p.wil=final.p.wil,dif.otus.wil=dif.otus.wil))

}



