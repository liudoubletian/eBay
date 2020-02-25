#' An method for differential abundance detection
#'
#' @param cutf the significance level
#' @param otu.data the data set where each row was a sample and each column was a taxa,
#' @param test.method a test including "t" and "wilcoxon"
#' @return dif.otus the detected differential otus
#' @return final.p the p value for each otu
#' @examples
#' group <- rep(c(0,1),each=46)
#' ebay.res <- eBay(otu.table,group,test.method="t",cutf=0.05)
#' @export
eBay=function(otu.data,group,test.method,cutf){

  sample.s<- nrow(otu.data)
  otu.n <- ncol(otu.data)
  case.s <- length(which(group == 0))
  con.s <- length(which(group == 1))
  case <- which(group == 0)
  con <- which(group == 1)


  phy_data <- phyloseq(otu_table(as.matrix(t(otu.data)),taxa_are_rows=TRUE))##transfrom data to phyloseq structure
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

  keepOTUs <-  prevalence > 0.2& taxa_sums(t(phy_data)) > (0.5 *sample.s) #& kep_pre1>0 & kep_pre2>0
  ###filtering the otus


  pru_data <- prune_taxa(keepOTUs, phy_data) ###prune the data set

  ntree_table<- t(otu_table(pru_data))### final otu table

  taxa.p <- ncol(ntree_table)

  rownames(ntree_table) <- as.character(1:sample.s)

  fit_glm <- MGLMfit(ntree_table, dist = "DM")

  ntree_para <- fit_glm@estimate

  ntree_para_mat <- rbind(matrix(rep(ntree_para, sample.s),byrow = TRUE,ncol = taxa.p))

  exp_norm <- matrix(NA,ncol=taxa.p,nrow=sample.s)

  for (n in 1:sample.s) {
    exp_norm[n, ] = unlist((ntree_table[n, ] + ntree_para_mat[n,]) / (sum(ntree_table[n, ]) +sum(ntree_para_mat[n,])))
  }

  colnames(exp_norm) <- colnames(ntree_table)

  exp_clr <- apply(exp_norm, 1, function(x){log2(x) - mean(log2(x))})

  if (test.method == "t"){
    exp_test <-  apply(exp_clr, 1, function(input){ t.test(input[case], input[con])$p.value})
    final.p <- p.adjust(exp_test,"BH")
    dif.otus <-  names(which(final.p < cutf))
  }
  if (test.method == "wilcoxon"){
    exp_test <- apply(exp_clr, 1, function(input){ wilcox.test(input[case], input[con])$p.value})
    final.p <- p.adjust(exp_test,"BH")
    dif.otus <-  names(which(final.p < cutf))
  }


  return(list(final.p=final.p,dif.otus=dif.otus))
}




