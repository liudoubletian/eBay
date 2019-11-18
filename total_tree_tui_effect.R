

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



new_tree_power_samll_effect=matrix(rep(0,8*12),ncol=12)
para_count_t=para_rare_t=para_rel_t=para_count_wil=para_rare_wil=para_rel_wil=matrix(ncol = 2,nrow=100)
para_t_alpha=para_wil_alpha=para_t_alpha_nt=para_wil_alpha_nt=para_t_alpha_ubay_nt=para_wil_alpha_ubay_nt=matrix(ncol = 2,nrow=100)
para_t_alpha_ubay=para_wil_alpha_ubay=matrix(ncol = 2,nrow=100)
power=fdr=matrix(rep(0,100*3),ncol=3)




N.sam=50
lotime=100
p=50
set.seed(2018)
pre_tree=Tree_Simulate(p)

set.seed(1)
real_pi=sample(unlist(real_control_pi),replace = FALSE,pre_tree$Nnode)
case_pi=control_pi=c()

for ( j in (p + 1) : (p + pre_tree$Nnode)){
  
  random_pi=real_pi[j-p]
  control_pi[which(pre_tree$edge[, 1] == j)]= c(random_pi,1-random_pi)
  case_pi[which(pre_tree$edge[, 1] == j)]= c(random_pi,1-random_pi)
}


for (h in 1: 6){
  p=50
  cc=c(0.1,2,4,6,7,8)
  delta=0.05
  eta=cc[h]
  
  control_theta=case_theta=rep( mean(unlist(the_eff)),pre_tree$Nnode)
  case_pi=control_pi=c()
  
  for ( j in (p + 1) : (p + pre_tree$Nnode)){
    
    random_pi=real_pi[j-p]
    control_pi[which(pre_tree$edge[, 1] == j)]= c(random_pi,1-random_pi)
    case_pi[which(pre_tree$edge[, 1] == j)]= c(random_pi,1-random_pi)
  }
  control_pi[9]=control_pi[8];case_pi[9]=case_pi[8]
  control_pi[12]=control_pi[5];case_pi[12]=case_pi[5]
  
  control_pi[5]=control_pi[5]+eta*delta;control_pi[8]=control_pi[8]-eta*delta;
  control_pi[9]=control_pi[9]+eta*delta;control_pi[12]=control_pi[12]-eta*delta;
  
  case_pi[5]=case_pi[5]-eta*delta;case_pi[8]=case_pi[8]+eta*delta
  case_pi[9]=case_pi[9]-eta*delta;case_pi[12]=case_pi[12]+eta*delta
  
  
  
  for ( m in  1:lotime){
    p=50
    pre_tree_table=DTM_simulation_sam(p=p,seed=1,N=N.sam,tree=pre_tree, control_pi=control_pi, 
                                      case_pi=case_pi,control_theta=control_theta,case_theta=case_theta,m=m)
    
    colnames(pre_tree_table)=as.character(pre_tree$edge[, 2])
    rownames(pre_tree_table)=as.character(1:(2*N.sam))
    group <- c(rep(0,N.sam),rep(1,N.sam))
    
    pre_tree_table_aft=pre_tree_table[,colnames(pre_tree_table)%in%(1:p)]
    
    sim_data = phyloseq(otu_table(as.matrix(t(pre_tree_table_aft)),taxa_are_rows=TRUE), phy_tree(pre_tree))
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
    
    
    
    
    gen_tree=phy_tree(sim_data_pru)
    gen_p=length(gen_tree$tip.label)
    taxa_index <- Taxa.index(gen_p, gen_tree)
    rownames(taxa_index)=1:gen_p
    
    gen_otu=t(otu_table(sim_data_pru))
    sum(gen_otu==0)/(dim(gen_otu)[1]*dim(gen_otu)[2])
    colnames(gen_otu)=as.character(1:gen_p)
    taxa.tab <- gen_otu %*% taxa_index
    taxa.alltab <- cbind(taxa.tab, gen_otu)
    tree_table=taxa.alltab[,-1]
    tree=gen_tree
    # tree=real_tree
    # gen_p=length(tree$tip.label)
    # taxa_index <- Taxa.index(gen_p, tree)
    # rownames(taxa_index)=1:gen_p
    # tree_table=pre_tree_table
    
    
    inter_node=unique(tree$edge[,1])
    leaf=c()
    for (i in 1:length(inter_node)){
      leaf_node=tree$edge[which(tree$edge[,1]==inter_node[i]),2]
      leaf=append(leaf,leaf_node)
    }
    
    pri_node_new=alpha_node_new=alpha_glm_node_new=c()
    exp_te_wil_glm=exp_te_glm=c()
    exp_te_wil_glm_ubay=exp_te_glm_ubay=c()
    samps=2*N.sam
    
    
    
    
    
    func_dtm=foreach(g=1:length(inter_node),.combine='cbind') %dopar%  inter_func(g=g,tree=tree,inter_node=inter_node,tree_table=tree_table,case_samp=N.sam,con_samp=N.sam,samps=samps,leaf=leaf)
    
    exp_te_glm=func_dtm[1,]
    exp_te_wil_glm=func_dtm[2,]
    exp_te_glm_ubay=func_dtm[3,]
    exp_te_wil_glm_ubay=func_dtm[4,]
    
    
    p=gen_p
    names(exp_te_glm)=names(exp_te_wil_glm)=names(exp_te_glm_ubay)=names(exp_te_wil_glm_ubay)=as.character(leaf)
    dif=colnames(gen_otu)[tree$tip.label%in%c(2,3,6,7,8)]
 
    
    para1=fdr_power_repair(exp_te_glm,cutf=0.05,dif=dif,leaf=leaf,func_dtm=func_dtm,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=gen_p,case_samp=N.sam,con_samp=N.sam,test_method="t.test",normal_method="glm")
    para_t_alpha[m,]=para1
    
    
    
    para2=fdr_power_repair((exp_te_wil_glm),cutf=0.05,dif=dif,leaf=leaf,func_dtm=func_dtm,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=gen_p,case_samp=N.sam,con_samp=N.sam,test_method="wilcox",normal_method="glm")
    para_wil_alpha[m,]=para2
    
    
    para3=fdr_power_repair(exp_te_glm_ubay,cutf=0.05,dif=dif,leaf=leaf,func_dtm=func_dtm,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=gen_p,case_samp=N.sam,con_samp=N.sam,test_method="t.test",normal_method="prior")
    para_t_alpha_ubay[m,]=para3
    
    
    
    para4=fdr_power_repair((exp_te_wil_glm_ubay),cutf=0.05,dif=dif,leaf=leaf,func_dtm=func_dtm,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=gen_p,case_samp=N.sam,con_samp=N.sam,test_method="wilcox",normal_method="prior")
    para_wil_alpha_ubay[m,]=para4
    
    order_tree=match(1:gen_p,colnames(tree_table))
    
    tree_table_leaf=tree_table[,order_tree]
    
    
    
    ntree_table=tree_table_leaf
    fit_control_glm <- MGLMfit(ntree_table, dist = "DM")
    
    ntree_alpha_con_glm=fit_control_glm@estimate
    
    alpha_ntree_glm=rbind(matrix(rep(ntree_alpha_con_glm, samps),byrow = TRUE,ncol = gen_p))
    
    pi_exp_glm_nt=matrix(NA,ncol=gen_p,nrow=(2*N.sam))
    pi_exp_ubay_nt=matrix(NA,ncol=gen_p,nrow=(2*N.sam))
    for (n in 1:(2*N.sam)) {
      pi_exp_glm_nt[n, ] = (ntree_table[n, ] + alpha_ntree_glm[n,]) / (sum(ntree_table[n, ]) +sum(alpha_ntree_glm[n,]))
      pi_exp_ubay_nt[n,]= unlist((ntree_table[n,] + 0.5) / (sum(ntree_table[n, ])+sum(rep(1/2,gen_p))))
      
    }
    colnames(pi_exp_glm_nt)=colnames(pi_exp_ubay_nt)=colnames(ntree_table)
    
    cij_expi_glm_nt=apply(t(pi_exp_glm_nt), 2, function(x){log2(x) - mean(log2(x))})
    cij_expi_ubay_nt=apply(t(pi_exp_ubay_nt), 2, function(x){log2(x) - mean(log2(x))})
    
    
    exp_te_glm_nt=exp_te_wil_glm_nt=exp_te_glm_ubay_nt=exp_te_wil_glm_ubay_nt=c()
    count_t=count_rare_t=count_rel_t=count_wil=count_rare_wil=count_rel_wil=c()
    
    for( t in 1:gen_p){
      
      
      exp_te_glm_nt[t] = t.test(cij_expi_glm_nt[t,1:N.sam], cij_expi_glm_nt[t,(N.sam+1):(2*N.sam)])$p.value
      exp_te_wil_glm_nt[t] = wilcox.test(cij_expi_glm_nt[t,1:N.sam], cij_expi_glm_nt[t,(N.sam+1):(2*N.sam)])$p.value
      
      exp_te_glm_ubay_nt[t] = t.test(cij_expi_ubay_nt[t,1:N.sam], cij_expi_ubay_nt[t,(N.sam+1):(2*N.sam)])$p.value
      exp_te_wil_glm_ubay_nt[t] = wilcox.test(cij_expi_ubay_nt[t,1:N.sam], cij_expi_ubay_nt[t,(N.sam+1):(2*N.sam)])$p.value
      
      
      
    }
    # 
    # 
    
    
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

write.table(new_tree_power_samll_effect,"/export/home/yangsheng/double/emp_norm/algtotal_tree_tui_effect.txt")


