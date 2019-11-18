

library(biomformat)
library(ape)
library(phyloseq)
library(dirmult)
library(MGLM)
library(foreach)
library(doParallel)

source("/export/home/yangsheng/double/emp_norm/alg_total.R")
source("/export/home/yangsheng/double/emp_norm/inter_func_total.R")
real_control_pi=read.table("/export/home/yangsheng/double/emp_norm/real_control_pi.txt")
the_eff=read.table("/export/home/yangsheng/double/emp_norm/the_eff.txt")
real_control_pi_nt=read.table("/export/home/yangsheng/double/emp_norm/real_control_pi_nt.txt")
the_eff_nt=read.table("/export/home/yangsheng/double/emp_norm/the_eff_nt.txt")



real_control_pi_nt=unlist(real_control_pi_nt)
the_eff_nt=unlist(the_eff_nt)



new_tree_power_samll_effect=matrix(rep(0,8*12),ncol=12)
para_count_t=para_rare_t=para_rel_t=para_count_wil=para_rare_wil=para_rel_wil=matrix(ncol = 2,nrow=100)
para_t_alpha=para_wil_alpha=para_t_alpha_nt=para_wil_alpha_nt=para_t_alpha_ubay_nt=para_wil_alpha_ubay_nt=matrix(ncol = 2,nrow=100)
para_t_alpha_ubay=para_wil_alpha_ubay=matrix(ncol = 2,nrow=100)
power=fdr=matrix(rep(0,100*3),ncol=3)

N.sam=50
samps=100
p=40
lotime=100
set.seed(2016)
pre_tree=Tree_Simulate(p)

pi_exp_glm_nt=pi_0.5_exp_nt=matrix(ncol = p,nrow=100)


for (h in 1: 6){ 
  p=40
  delta=0.005
  cc=c(0.1,1.5,2,2.5,3,3.5)
  
  set.seed(2019)
  sam_pi=sample(real_control_pi_nt,replace = FALSE,p)
  case_pi=control_pi=sam_pi/sum(sam_pi)
  names(case_pi)=names(control_pi)=1:p
  
  eta=cc[h]
  control_pi[1]=control_pi[1]-eta*delta;control_pi[7]=control_pi[7]+eta*delta;
  control_pi[22]=control_pi[22]-eta*delta;control_pi[11]=control_pi[11]+eta*delta
  control_pi[23]=control_pi[23]-eta*delta;control_pi[40]=control_pi[40]+eta*delta
  
  
  case_pi[1]=case_pi[1]+eta*delta;case_pi[7]=case_pi[7]-eta*delta;
  case_pi[22]=case_pi[22]+eta*delta;case_pi[11]=case_pi[11]-eta*delta
  case_pi[23]=case_pi[23]+eta*delta;case_pi[40]=case_pi[40]-eta*delta
 
  for ( m in  1:lotime){
    p=40
    seed=1
    ntree_table = matrix(NA, samps, p)
    for(i in 1:N.sam){
      set.seed(i*m*seed)
      parent1 = round(runif(1, 5000, 50000))
      parent2 = round(runif(1, 5000, 50000))
      set.seed(i*m*seed)
      ntree_table[i, ] <- simPop(J=1, K=p, n=parent1, pi = control_pi, theta = 0.15)$data
      set.seed(i*m*(seed+1))
      ntree_table[i+N.sam,] <- simPop(J=1, K=p, n=parent2, pi = case_pi, theta = 0.15)$data
    }
    colnames(ntree_table)=as.character(1:p)
    
    group <- c(rep(0,N.sam),rep(1,N.sam))
    
    sim_data = phyloseq(otu_table(as.matrix(t(ntree_table)),taxa_are_rows=TRUE), phy_tree(pre_tree))
    minobs=0
    prevalence = apply(as(t(otu_table(sim_data)), "matrix"), 2, function(x, minobs) {
      return(sum(x > minobs))
    }, minobs)/(2*N.sam)
    
    kep_pre1 = apply(as(t(otu_table(sim_data))[1:N.sam,], "matrix"), 2, function(x) {
      return(sum(x > 0))
    })/(N.sam)
    
    kep_pre2 = apply(as(t(otu_table(sim_data))[(N.sam+1):(2*N.sam),], "matrix"), 2, function(x) {
      return(sum(x > 0))
    })/(N.sam)
    
    keepOTUs = prevalence > 0.2 & taxa_sums(t(sim_data)) > (0.5 *2*N.sam)&kep_pre1>0&kep_pre2>0
    
    sim_data_pru=prune_taxa(keepOTUs, sim_data)
    tree=phy_tree(sim_data_pru)
    p=length(tree$tip.label)
    ntree_table=t(otu_table(sim_data_pru))
    
    
    
    taxa_index <- Taxa.index(p, tree)
    rownames(taxa_index)=1:p
    colnames(ntree_table)=as.character(1:p)
    rownames(ntree_table)=as.character(1:samps)
    
    taxa.tab <- ntree_table %*% taxa_index
    taxa.alltab <- cbind(taxa.tab, ntree_table)
    tree_table=taxa.alltab[,-1]
    
    inter_node=unique(tree$edge[,1])
    leaf=c()
    for (i in 1:length(inter_node)){
      leaf_node=tree$edge[which(tree$edge[,1]==inter_node[i]),2]
      leaf=append(leaf,leaf_node)
    }
    
    
    func_dtm=foreach(g=1:length(inter_node),.combine='cbind') %dopar%  inter_func(g=g,tree=tree,inter_node=inter_node,tree_table=tree_table,case_samp=N.sam,con_samp=N.sam,samps=samps,leaf=leaf)
    
    exp_te_glm=func_dtm[1,]
    exp_te_wil_glm=func_dtm[2,]
    exp_te_glm_ubay=func_dtm[3,]
    exp_te_wil_glm_ubay=func_dtm[4,]
    
    
    names(exp_te_glm)=names(exp_te_wil_glm)=names(exp_te_glm_ubay)=names(exp_te_wil_glm_ubay)=as.character(leaf)
    
    dif=colnames(ntree_table)[tree$tip.label%in%c(1,7,11,22,23,40)]
    
    para1=fdr_power_repair(exp_te_glm,cutf=0.05,dif=dif,leaf=leaf,func_dtm=func_dtm,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=p,case_samp=N.sam,con_samp=N.sam,test_method="t.test",normal_method="glm")
    para_t_alpha[m,]=para1
    
    
    
    para2=fdr_power_repair((exp_te_wil_glm),cutf=0.05,dif=dif,leaf=leaf,func_dtm=func_dtm,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=p,case_samp=N.sam,con_samp=N.sam,test_method="wilcox",normal_method="glm")
    para_wil_alpha[m,]=para2
    
    
    para3=fdr_power_repair(exp_te_glm_ubay,cutf=0.05,dif=dif,leaf=leaf,func_dtm=func_dtm,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=p,case_samp=N.sam,con_samp=N.sam,test_method="t.test",normal_method="prior")
    para_t_alpha_ubay[m,]=para3
    
    
    
    para4=fdr_power_repair((exp_te_wil_glm_ubay),cutf=0.05,dif=dif,leaf=leaf,func_dtm=func_dtm,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=p,case_samp=N.sam,con_samp=N.sam,test_method="wilcox",normal_method="prior")
    para_wil_alpha_ubay[m,]=para4
    
  


fit_control_glm <- MGLMfit(ntree_table, dist = "DM")
ntree_alpha_con_glm=fit_control_glm@estimate

alpha_ntree_glm=rbind(matrix(rep(ntree_alpha_con_glm, 2*N.sam),byrow = TRUE,ncol = p))

pi_exp_glm_nt=matrix(NA,ncol=p,nrow=(2*N.sam))
pi_exp_ubay_nt=matrix(NA,ncol=p,nrow=(2*N.sam))
for (n in 1:(2*N.sam)) {
  pi_exp_glm_nt[n, ] = unlist((ntree_table[n, ] + alpha_ntree_glm[n,]) / (sum(ntree_table[n, ]) +sum(alpha_ntree_glm[n,])))
  pi_exp_ubay_nt[n,]= unlist((ntree_table[n,] + 0.5) / (sum(ntree_table[n, ])+sum(rep(1/2,p))))
  
}
colnames(pi_exp_glm_nt)=colnames(pi_exp_ubay_nt)=colnames(ntree_table)

cij_expi_glm_nt=apply(t(pi_exp_glm_nt), 2, function(x){log2(x) - mean(log2(x))})
cij_expi_ubay_nt=apply(t(pi_exp_ubay_nt), 2, function(x){log2(x) - mean(log2(x))})


exp_te_glm_nt=exp_te_wil_glm_nt=exp_te_glm_ubay_nt=exp_te_wil_glm_ubay_nt=c()
count_t=count_rare_t=count_rel_t=count_wil=count_rare_wil=count_rel_wil=c()

for( t in 1:p){
  
  
  exp_te_glm_nt[t] = t.test(cij_expi_glm_nt[t,1:N.sam], cij_expi_glm_nt[t,(N.sam+1):(2*N.sam)])$p.value
  exp_te_wil_glm_nt[t] = wilcox.test(cij_expi_glm_nt[t,1:N.sam], cij_expi_glm_nt[t,(N.sam+1):(2*N.sam)])$p.value
  
  exp_te_glm_ubay_nt[t] = t.test(cij_expi_ubay_nt[t,1:N.sam], cij_expi_ubay_nt[t,(N.sam+1):(2*N.sam)])$p.value
  exp_te_wil_glm_ubay_nt[t] = wilcox.test(cij_expi_ubay_nt[t,1:N.sam], cij_expi_ubay_nt[t,(N.sam+1):(2*N.sam)])$p.value
  
 
}


para_t_alpha_nt[m,]=pw_fdr(p.adjust(exp_te_glm_nt),dif=dif,cutf = 0.05)


para_wil_alpha_nt[m,]=pw_fdr(p.adjust(exp_te_wil_glm_nt),dif=dif,cutf = 0.05)


para_t_alpha_ubay_nt[m,]=pw_fdr(p.adjust(exp_te_glm_ubay_nt),dif=dif,cutf = 0.05)


para_wil_alpha_ubay_nt[m,]=pw_fdr(p.adjust(exp_te_wil_glm_ubay_nt),dif=dif,cutf = 0.05)


print(m)





  }
  
  
  pwf_glm_t=apply(para_t_alpha,2,mean)
  pwf_glm_wil=apply(para_wil_alpha,2,mean)
  
  pwf_ubay_t=apply(para_t_alpha_ubay,2,mean)
  pwf_ubay_wil=apply(para_wil_alpha_ubay,2,mean)
  
  
  
  pwf_nt_t=apply(para_t_alpha_nt,2,mean)
  pwf_nt_wil=apply(para_wil_alpha_nt,2,mean)
  
  pwf_nt_ubay_t=apply(para_t_alpha_ubay_nt,2,mean)
  pwf_nt_ubay_wil=apply(para_wil_alpha_ubay_nt,2,mean)
  
  
  
  te_pwf=matrix(c(pwf_glm_t,pwf_ubay_t,pwf_nt_t,pwf_nt_ubay_t),byrow=TRUE,ncol=2)
  wil_pwf=matrix(c(pwf_glm_wil,pwf_ubay_wil,pwf_nt_wil,pwf_nt_ubay_wil),byrow=TRUE,ncol=2)
  final_n_tr_pwf=rbind(te_pwf,wil_pwf)
  new_tree_power_samll_effect[,(2*h-1):(2*h)]=final_n_tr_pwf
}

write.table(new_tree_power_samll_effect,"/export/home/yangsheng/double/emp_norm/algtotal_ntree_effect.txt")
