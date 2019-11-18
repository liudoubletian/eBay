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


N.sam=50
new_tree_power_samll_effect=matrix(rep(0,8*12),ncol=12)
para_count_t=para_rare_t=para_rel_t=para_count_wil=para_rare_wil=para_rel_wil=matrix(ncol = 2,nrow=100)
para_t_alpha=para_wil_alpha=para_t_alpha_nt=para_wil_alpha_nt=para_t_alpha_ubay_nt=para_wil_alpha_ubay_nt=matrix(ncol = 2,nrow=100)
para_t_alpha_ubay=para_wil_alpha_ubay=matrix(ncol = 2,nrow=100)
power=fdr=matrix(rep(0,100*3),ncol=3)



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

delta=0.05
eta=4

control_pi[5]=control_pi[5]+eta*delta;control_pi[8]=control_pi[8]-eta*delta;

case_pi[5]=case_pi[5]-eta*delta;case_pi[8]=case_pi[8]+eta*delta


for (h in 1: 6){ 
  
  cc=c(0.05,0.1,0.15,0.2,0.25,0.3)
  
  control_theta=rep( cc[h],pre_tree$Nnode)#runif( tree$Nnode,0,0.1)
  case_theta=rep(cc[h],pre_tree$Nnode)#r
  
  
  
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
    
    p=gen_p
    
    func_dtm=foreach(g=1:length(inter_node),.combine='cbind') %dopar% inter_func(g=g,tree=tree,inter_node=inter_node,tree_table=tree_table,leaf=leaf)
    
    pi_exp_glm=(func_dtm[5:104,])
    pi_exp_glm_ubay=(func_dtm[105:204,])
    
    colnames(pi_exp_glm)=colnames(pi_exp_glm_ubay)=as.character(leaf)
    
    chil_node_fun=function(j, tree,pi_exp,real_qian) {
      pa = tree$edge[which(tree$edge[, 2] == j), 1]
      while (pa[length(pa)] > as.numeric(real_qian)) {
        pa = append(pa, tree$edge[which(tree$edge[, 2] == pa[length(pa)]), 1], length(pa))
      }
      pa = c(pa[-which(pa==as.numeric(real_qian))], j)
      pi_co=c()
      for( n in 1: nrow(pi_exp)){
        pro=1
        for( k in 1: length(pa)){
          pro=pi_exp[n,which(colnames(pi_exp)==pa[k])]*pro
        }
        pi_co[n]=pro
      }
      return(pi_co)
    }
    dif=colnames(gen_otu)[tree$tip.label%in%c(2:8)]
    
    pa_pi_glm=t(foreach(x=1:p,.combine='rbind') %dopar% chil_node_fun(j = x,tree = tree,pi_exp = pi_exp_glm,real_qian=p+1))
    pa_pi_ubay=t(foreach(x=1:p,.combine='rbind') %dopar% chil_node_fun(j = x,tree = tree,pi_exp = pi_exp_glm_ubay,real_qian=p+1))
    
    colnames(pa_pi_glm)=colnames(pa_pi_ubay)=as.character(1:p)
    
    cij_expi_glm=apply(t(pa_pi_glm), 2, function(x){log2(x) - mean(log2(x))})
    cij_expi_glm_ubay=apply(t(pa_pi_ubay), 2, function(x){log2(x) - mean(log2(x))})
    
    
    sam_glm <- apply(cij_expi_glm, 1, function(t.input){as.numeric(t.test(x=t.input[1:N.sam],y=t.input[(N.sam+1):(2*N.sam)])$p.value)})
    sam_glm_ubay <- apply(cij_expi_glm_ubay, 1, function(t.input){as.numeric(t.test(x=t.input[1:N.sam],y=t.input[(N.sam+1):(2*N.sam)])$p.value)})
    sam_glm_wil <- apply(cij_expi_glm, 1, function(t.input){as.numeric(wilcox.test(x=t.input[1:N.sam],y=t.input[(N.sam+1):(2*N.sam)])$p.value)})
    sam_glm_wil_ubay <- apply(cij_expi_glm_ubay, 1, function(t.input){as.numeric(wilcox.test(x=t.input[1:N.sam],y=t.input[(N.sam+1):(2*N.sam)])$p.value)})
    
    
    para_t_alpha_nt[m,]=pw_fdr(p.adjust(sam_glm),dif=dif,cutf = 0.05)
    para_wil_alpha_nt[m,]=pw_fdr(p.adjust(sam_glm_wil),dif=dif,cutf = 0.05)
    para_t_alpha_ubay_nt[m,]=pw_fdr(p.adjust(sam_glm_ubay),dif=dif,cutf = 0.05)
    para_wil_alpha_ubay_nt[m,]=pw_fdr(p.adjust(sam_glm_wil_ubay),dif=dif,cutf = 0.05)
    print(m)
    
  }
  
  pwf_nt_t=apply(para_t_alpha_nt,2,mean)
  pwf_nt_wil=apply(para_wil_alpha_nt,2,mean)
  
  pwf_nt_ubay_t=apply(para_t_alpha_ubay_nt,2,mean)
  pwf_nt_ubay_wil=apply(para_wil_alpha_ubay_nt,2,mean)
  
  
  
  te_pwf=matrix(c(pwf_nt_t,pwf_nt_ubay_t),byrow=TRUE,ncol=2)
  wil_pwf=matrix(c(pwf_nt_wil,pwf_nt_ubay_wil),byrow=TRUE,ncol=2)
  final_n_tr_pwf=rbind(te_pwf,wil_pwf)
  new_tree_power_samll_effect[1:4,(2*h-1):(2*h)]=final_n_tr_pwf
}

write.table(new_tree_power_samll_effect[1:4,],"/export/home/yangsheng/double/emp_norm/algtotal_theta_global.txt")




