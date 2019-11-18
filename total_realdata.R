library(biomformat)
library(ape)
library(phyloseq)
library(metagenomeSeq)
library(DESeq2)
library(MGLM)
library(foreach)
library(parallel)
library(plyr)
library(dirmult)
library(VennDiagram)

pers=read_biom("/Users/tiantian/Documents/data/gordon lab/Subramanian_et_al_otu_table.biom")

#imp_tree=read.tree(file = "/Users/tiantian/Documents/data/gordon lab/rep_phylo.tre")

imp_txt=read.table(file = "/Users/tiantian/Documents/data/gordon lab/imp.txt",sep="\t",header=FALSE)
imp_id=imp_txt[,2]
meta_sam=read.table("/Users/tiantian/Documents/data/gordon lab/sam.txt",sep="\t",header=FALSE)
meta_sam=as.matrix(meta_sam)
meta_healthy=read.table("/Users/tiantian/Documents/data/gordon lab/healthy_twins_samp.txt",sep="\t",header=FALSE)

imp_tree=read.tree(file = "/Users/tiantian/Documents/data/gordon lab/UPGMA.tre")

Taxa.index=function (p, phy_tree) 
{
  taxa_index <- matrix(0, p, phy_tree$Nnode)
  colnames(taxa_index) <- as.character((p + 1):(p + phy_tree$Nnode))
  for (i in 1:p) {
    temp <- i
    for (j in 1:1000) {
      index <- which(phy_tree$edge[, 2] == temp)
      if (length(index) != 0) {
        parent <- phy_tree$edge[index, 1]
        taxa_index[i, which(colnames(taxa_index) == as.character(parent))] <- 1
        temp <- parent
      }
      else {
        break
      }
    }
  }
  return(taxa_index)
}

inter_func=function(g=g,tree=tree,inter_node=inter_node,tree_table=tree_table,case_samp,con_samp,samps,leaf=leaf){
  library(MGLM)
  sam_glm=sam_glm_wil=c()
  sam_glm_ubay=sam_glm_wil_ubay=c()
  
  
  lit_tree=tree$edge[which(tree$edge[,1]==inter_node[g]),2]
  lit_tree_table=tree_table[,match(lit_tree,colnames(tree_table))]
  colnames(lit_tree_table)=as.character(lit_tree)
  
  
  
  fit_glm_con=MGLMfit(data.frame(lit_tree_table[1:case_samp,]),dist = "DM")
  alpha_con_glm=fit_glm_con@estimate
  
  fit_glm_case=MGLMfit(data.frame(lit_tree_table[(case_samp+1):(case_samp+con_samp),]),dist = "DM")
  alpha_case_glm=fit_glm_case@estimate
  
  pi_exp_glm=matrix(rep(0,samps*length(lit_tree)),ncol=length(lit_tree))
  pi_exp_glm_ubay=matrix(rep(0,samps*length(lit_tree)),ncol=length(lit_tree))
  
  
  alpha_matrix_glm = rbind(matrix(rep(alpha_con_glm, case_samp),byrow = TRUE,ncol = length(lit_tree)),
                           matrix(rep(alpha_case_glm, con_samp),byrow = TRUE, ncol = length(lit_tree)))
  
  
  for (n in 1:samps) {
    pi_exp_glm[n,]= unlist((lit_tree_table[n,] + alpha_matrix_glm[n,]) / (sum(lit_tree_table[n, ]) + sum(alpha_matrix_glm[n,])))
    pi_exp_glm_ubay[n,]= unlist((lit_tree_table[n,] + 0.5) / (sum(lit_tree_table[n, ])+sum(rep(1/2,length(lit_tree)))))
  }
  
  
  cij_expi_glm=apply(t(pi_exp_glm), 2, function(x){log2(x) - mean(log2(x))})
  cij_expi_glm_ubay=apply(t(pi_exp_glm_ubay), 2, function(x){log2(x) - mean(log2(x))})
  
  for( t in 1:length(lit_tree)){
    
    sam_glm[t]= t.test(cij_expi_glm[t,1:case_samp], cij_expi_glm[t,(case_samp+1):samps])$p.value
    sam_glm_wil[t]= wilcox.test(cij_expi_glm[t,1:case_samp], cij_expi_glm[t,(case_samp+1):samps])$p.value
    
    sam_glm_ubay[t]= t.test(cij_expi_glm_ubay[t,1:case_samp], cij_expi_glm_ubay[t,(case_samp+1):samps])$p.value
    sam_glm_wil_ubay[t]= wilcox.test(cij_expi_glm_ubay[t,1:case_samp], cij_expi_glm_ubay[t,(case_samp+1):samps])$p.value
  }
  return(rbind(sam_glm,sam_glm_wil,sam_glm_ubay,sam_glm_wil_ubay,pi_exp_glm,pi_exp_glm_ubay))
}
meta_healthy=as.matrix(meta_healthy)

sam_samp=meta_sam[grep("Acute Phase",meta_sam[, 4]),]
sam_samp1=sam_samp#[which(as.numeric(sam_samp[,8])< -3),]
sam_mon=as.numeric(sam_samp1[,5])
sam_samp=sam_samp1[which(sam_mon>12&sam_mon<18),]
sam_id=sam_samp[,2]
meta_healthy=meta_healthy[which(as.numeric(meta_healthy[,7])>-2),]
meta_healthy=meta_healthy[1:92,]#meta_healthy[1:274,]#
healthy_id=meta_healthy[which(as.numeric(meta_healthy[,6])>12&as.numeric(meta_healthy[,6])<18),4]
sam_persis=c(sam_id,healthy_id)##70,251
samp_to=colnames(biom_data(pers))
otu_to=rownames(biom_data(pers))
mam_data=biom_data(pers)[,match(sam_persis,samp_to)]

mam_da=mam_data[match(imp_tree$tip.label,rownames(mam_data)),]
#mam_fam=observation_metadata(pers)[match(imp_tree$tip.label,names(observation_metadata(pers)))]
mam_da=t(as(mam_da,"matrix"))
healthy_id=rownames(mam_da)[-(1:length(sam_id))]
samps=length(c(sam_id,healthy_id))


sim_data = phyloseq(otu_table(as.matrix(t(mam_da)),taxa_are_rows=TRUE), phy_tree(imp_tree))

prevalence = apply(as(t(otu_table(sim_data)), "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(samps)



kep_pre1 = apply(as(t(otu_table(sim_data))[1:length(sam_id),], "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(length(sam_id))

kep_pre2 = apply(as(t(otu_table(sim_data))[(length(sam_id)+1):samps,], "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(length(healthy_id))

keepOTUs = prevalence > 0.2 & taxa_sums(t(sim_data)) > (0.5 *samps)&kep_pre1>0&kep_pre2>0

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

grouplabel <- c(rep(0,length(sam_id)),rep(1,length(healthy_id)))



order_tree=match(1:gen_p,colnames(tree_table))

tree_table_leaf=tree_table[,order_tree]

group= c(rep(0,length(sam_id)),rep(1,length(healthy_id)))
otu <- otu_table(as.matrix(t(tree_table_leaf)), taxa_are_rows=TRUE)
meta_otu=data.frame(group=c(rep(0,length(sam_id)),rep(1,length(healthy_id))))
sample_data<- sample_data(meta_otu)
rownames(sample_data)=colnames(otu)
otu_nontree <- phyloseq(otu, sample_data)
rare_otu=rarefy_even_depth(otu_nontree, sample.size = min(sample_sums(otu_nontree))*0.9,rngseed = 2018)
otu_rare <- otu_table((rare_otu), taxa_are_rows=TRUE)
tree_table_otu=as.matrix(as.data.frame(t(otu_rare)))

#


##########metagenomeseq_raw
rownames(meta_otu)=rownames(tree_table_otu)
phenotypeData = AnnotatedDataFrame(meta_otu)
obj_raw = newMRexperiment(t(tree_table_leaf),phenoData=phenotypeData)
p_obj_raw = cumNormStatFast(obj_raw)
obj_raw = cumNorm(obj_raw, p = p_obj_raw)
normFactor_raw = normFactors(obj_raw)
normFactor_raw = log2(normFactor_raw/median(normFactor_raw) + 1)
mod_raw = model.matrix(~group,data = environment(obj_raw))
settings = zigControl(maxit = 100, verbose = TRUE)
res_raw = fitZig(obj = obj_raw, mod = mod_raw, useCSSoffset = FALSE,
                 control = settings)

meta_raw_res=rownames(MRcoefs(res_raw))[which(MRcoefs(res_raw)$adjPvalues<0.05)]

#############DESeq2_raw
deseq_group=factor(c(rep(0,length(sam_id)),rep(1,length(healthy_id))))
dds_raw <- DESeqDataSetFromMatrix(t(tree_table_leaf+1), DataFrame(deseq_group), ~ deseq_group)
dds_raw <- estimateSizeFactors(dds_raw)
dds_raw <- estimateDispersions(dds_raw)
de_result_raw=nbinomWaldTest(dds_raw)
re_des_raw=results(de_result_raw)
deseq_raw_res=which(re_des_raw$padj<0.05)

source("/Users/tiantian/Documents/data/gordon lab/ancom.R_1.1-3/ancom.R/R/ancom_functions.r")

ancom.test=ANCOM(OTUdat=data.frame((tree_table_leaf),Group=c(rep(0,length(sam_id)),rep(1,length(healthy_id)))),
                 sig=0.05, multcorr=3, tau=0.02, theta=0.1, repeated=FALSE)$detected

ancom.test=gsub("X","",ancom.test)



ntree_table=tree_table_leaf
fit_control_glm <- MGLMfit(ntree_table[1:length(sam_id),], dist = "DM")
fit_case_glm <- MGLMfit(ntree_table[(length(sam_id)+1):(samps),], dist = "DM")
ntree_alpha_con_glm=fit_control_glm@estimate
ntree_alpha_ca_glm=fit_case_glm@estimate

alpha_ntree_glm=rbind(matrix(rep(ntree_alpha_con_glm, length(sam_id)),byrow = TRUE,ncol = gen_p),
                      matrix(rep(ntree_alpha_ca_glm, length(healthy_id)),byrow = TRUE, ncol = gen_p))

pi_exp_glm_nt=matrix(NA,ncol=gen_p,nrow=samps)
pi_exp_ubay_nt=matrix(NA,ncol=gen_p,nrow=samps)
for (n in 1:samps) {
  pi_exp_glm_nt[n, ] = (ntree_table[n, ] + alpha_ntree_glm[n,]) / (sum(ntree_table[n, ]) +sum(alpha_ntree_glm[n,]))
  pi_exp_ubay_nt[n,]= unlist((ntree_table[n,] + 0.5) / (sum(ntree_table[n, ])+sum(rep(1/2,gen_p))))
  
}
colnames(pi_exp_glm_nt)=colnames(pi_exp_ubay_nt)=colnames(ntree_table)

cij_expi_glm_nt=apply(t(pi_exp_glm_nt), 2, function(x){log2(x) - mean(log2(x))})
cij_expi_ubay_nt=apply(t(pi_exp_ubay_nt), 2, function(x){log2(x) - mean(log2(x))})


exp_te_glm_nt=exp_te_wil_glm_nt=exp_te_glm_ubay_nt=exp_te_wil_glm_ubay_nt=c()
count_t=count_rare_t=count_rel_t=count_wil=count_rare_wil=count_rel_wil=c()

for( t in 1:gen_p){
  
  
  exp_te_glm_nt[t] = t.test(cij_expi_glm_nt[t,1:length(sam_id)], cij_expi_glm_nt[t,(length(sam_id)+1):(samps)])$p.value
  exp_te_wil_glm_nt[t] = wilcox.test(cij_expi_glm_nt[t,1:length(sam_id)], cij_expi_glm_nt[t,(length(sam_id)+1):(samps)])$p.value
  
  exp_te_glm_ubay_nt[t] = t.test(cij_expi_ubay_nt[t,1:length(sam_id)], cij_expi_ubay_nt[t,(length(sam_id)+1):(samps)])$p.value
  exp_te_wil_glm_ubay_nt[t] = wilcox.test(cij_expi_ubay_nt[t,1:length(sam_id)], cij_expi_ubay_nt[t,(length(sam_id)+1):(samps)])$p.value
  
  
  count=tree_table_leaf[,t]
  count_t[t]=t.test(count[1:length(sam_id)],count[(length(sam_id)+1):(samps)])$p.value
  
  count_rare=tree_table_otu[,t]
  count_rare_t[t]=t.test(count_rare[1:length(sam_id)],count_rare[(length(sam_id)+1):(samps)])$p.value
  
  count_rel=tree_table_leaf[,t]/rowSums(tree_table_leaf)
  count_rel_t[t]=t.test(count_rel[1:length(sam_id)],count_rel[(length(sam_id)+1):(samps)])$p.value
  
  
  
  
  count_wil[t]=wilcox.test(count[1:length(sam_id)],count[(length(sam_id)+1):(samps)])$p.value
  
  count_rare_wil[t]=wilcox.test(count_rare[1:length(sam_id)],count_rare[(length(sam_id)+1):(samps)])$p.value
  count_rel_wil[t]=wilcox.test(count_rel[1:length(sam_id)],count_rel[(length(sam_id)+1):(samps)])$p.value
  
}


names(exp_te_glm_nt)=names(exp_te_wil_glm_nt)=names(exp_te_glm_ubay_nt)=names(exp_te_wil_glm_ubay_nt)=as.character(1:gen_p)

names(count_t)=names(count_rare_t)=names(count_rel_t)=names(count_wil)=names(count_rare_wil)=names(count_rel_wil)=as.character(1:gen_p)




inter_node=unique(tree$edge[,1])
leaf=c()
for (i in 1:length(inter_node)){
  leaf_node=tree$edge[which(tree$edge[,1]==inter_node[i]),2]
  leaf=append(leaf,leaf_node)
}

pri_node_new=alpha_node_new=alpha_glm_node_new=c()
exp_te_wil_glm=exp_te_glm=c()
exp_te_wil_glm_ubay=exp_te_glm_ubay=c()



func_dtm=foreach(g=1:length(inter_node),.combine='cbind') %dopar% inter_func(g=g,tree=tree,inter_node=inter_node,tree_table=tree_table,leaf=leaf,case_samp=length(sam_id),con_samp=length(healthy_id),samps=samps)

exp_te_glm=func_dtm[1,]
exp_te_wil_glm=func_dtm[2,]
exp_te_glm_ubay=func_dtm[3,]
exp_te_wil_glm_ubay=func_dtm[4,]


p=gen_p
names(exp_te_glm)=names(exp_te_wil_glm)=names(exp_te_glm_ubay)=names(exp_te_wil_glm_ubay)=as.character(leaf)
tree$tip.label=1:p

find_dif=function(p.val,cutf,leaf=leaf,taxa_index,tree_table,tree,p,case_samp,con_samp,test_method,normal_method){
  library(doParallel)
  MGLMfit=MGLMfit
  degene=function(j, tree,root.node) {
    pa = tree$edge[which(tree$edge[, 2] == j), 1]
    while (pa[length(pa)] > as.numeric(root.node)) {
      pa = append(pa, tree$edge[which(tree$edge[, 2] == pa[length(pa)]), 1], length(pa))
    }
    pa = c(j,pa)
    return(pa)
  }
  
  child_fun=function (p, phy_tree){
    taxa_index <- matrix(0, length(phy_tree$tip.label)+phy_tree$Nnode,length(phy_tree$tip.label)+phy_tree$Nnode)
    colnames(taxa_index) <- as.character(1:(p + phy_tree$Nnode))
    for (i in 1:(p + phy_tree$Nnode)) {
      temp <- i
      while(temp %in% phy_tree$edge[, 2]) {
        index <- which(phy_tree$edge[, 2] == temp)
        if (length(index) != 0) {
          parent <- phy_tree$edge[index, 1]
          taxa_index[i, which(colnames(taxa_index) == as.character(parent))] <- 1
          temp <- parent
        }
      }
    }
    return(taxa_index)
  }
  
  judg_fun=function(node,tree,tree_table){
    
    pa_node=tree$edge[match(node,tree$edge[,2]),1]
    final_node=c()
    for(i in 1:length(pa_node)){
      
      pa_count_sub=tree_table[,colnames(tree_table)%in%pa_node[i]]
      nam_pa=pa_node[i]
      name_set=c()
      if(sum(pa_count_sub==0)/100<=0.1){
        nam_pa=node[i]
      }
      else{
        while(sum(tree_table[,colnames(tree_table)%in%nam_pa]==0)/100>0.1){
          
          nam_pa=tree$edge[tree$edge[,2]%in%nam_pa,1]
          
          
        }
      }
      final_node=append(final_node,nam_pa)
    }
    return(final_node)
  }
  
  
  est_fun=function(node,tree,tree_table){
    
    
    final_node=c()
    for(i in 1:length(node)){
      pa_count_sub=tree_table[,colnames(tree_table)%in%node[i]]
      nam_pa=node[i]
      if(sum(pa_count_sub==0)/100<=0.1){
        nam_pa=node[i]
      }
      else{
        while(sum(tree_table[,colnames(tree_table)%in%nam_pa]==0)/100>0.1){
          nam_pa=tree$edge[tree$edge[,2]%in%nam_pa,1]
        }
      }
      final_node=append(final_node,nam_pa)
    }
    return(final_node)
  }
  
  
  reject =names(which(p.val<= cutf))
  child_index=child_fun(p, tree)
  
  pri_node_in=pri_node_out=c()
  node_in=node_out=c()
  
  if (length(reject)==0) {
    TP = length(intersect(reject, dif))
    FP = length(reject) - TP
    power_test = TP / length(dif)
    fdr_test = ifelse(length(reject) == 0, 0, FP / length(reject))
  }
  
  if (length(reject)>=1) {
    
    
    edge_dif=lapply(1:p,function(x) {
      edge=degene(j=x,tree=tree,root.node=p+1)
      return(intersect(edge,reject))}
    )
    site1=which(lengths(edge_dif)>=2)
    real_qian_do1=sapply(edge_dif[site1], function(v) return(tail(v,1))) 
    real_qian_do1=unique(real_qian_do1)
    
    
    # 
    real_qian1=c()
    if(length(real_qian_do1)==0){real_qian1=NULL}
    else{
      for(w in 1:length(real_qian_do1)){
        real_node=tree$edge[which(tree$edge[,2]==real_qian_do1[w]),1]
        real_qian_sub=est_fun(real_node,tree=tree,tree_table=tree_table)
        real_qian1=append(real_qian1,real_qian_sub)
      }
    }  
    
    qian_sub1=lapply(tree$edge[match(real_qian_do1,tree$edge[,2]),1],function(x) {
      all_qian_sub=which(child_index[,x]==1)
      return(all_qian_sub)
    }
    )
    qian_sub1=unique(unlist(qian_sub1))
    real_qian_do2=setdiff(reject,qian_sub1)
    
    real_qian2=c()
    if(length(real_qian_do2)==0){real_qian2=real_qian2}
    else{
      for(l in 1:length(real_qian_do2)){
        real_node=judg_fun(real_qian_do2[l],tree=tree,tree_table=tree_table)
        real_qian2=append(real_qian2,real_node)
      }
    }
    
    
    node_intec=intersect(real_qian_do2,unique(real_qian2))
    
    real_qian=c(unique(real_qian1),setdiff(real_qian2,node_intec))
    
    if(length(node_intec)==0){p.val=p.val}
    else{
      for (i in 1:length(node_intec)){
        x=node_intec[i]
        if(x%in%1:p){p.val=p.val}
        else{
          all_qian_sub=which(taxa_index[,which(colnames(taxa_index)==x)]==1)
          p.val[match(all_qian_sub,names(p.val))]=p.val[match(x,names(p.val))]
        }
      }
    }
    
    if(length(real_qian)==0){
      reje =names(which(p.adjust(p.val[match(1:p,names(p.val))])<= cutf))
      TP = length(intersect(reje, dif))
      FP = length(reje) - TP
      power_test = TP/length(dif)
      fdr_test = ifelse(length(reje) == 0, 0, FP/length(reje))
    }
    
    else{
      
      
      all_qian=lapply(c(tree$edge[match(real_qian_do1,tree$edge[,2]),1],setdiff(real_qian_do2,node_intec)),function(x) {
        if(x%in%1:p){all_qian_sub=x}
        else{
          all_qian_sub=which(taxa_index[,which(colnames(taxa_index)==x)]==1)
        }
        return(all_qian_sub)}
      )
      
      all_qian=unique(unlist(all_qian))
      
      
      
      
      
      exp_te_wil_glm_re=exp_te_glm_re=exp_te_wil_glm_re_ubay=exp_te_glm_re_ubay=c()
      reject_test_alpha=reject_test_glm=reject_test_alpha_wil=reject_test_glm_wil=c()
      
      
      for(w in 1:length(real_qian)){
        
        real_qian_node=which(taxa_index[,which(colnames(taxa_index)==real_qian[w])]==1)
        
        
        
        
        
        
        exp_te_wil_glm_re_sub=exp_te_glm_re_sub=exp_te_wil_glm_re_sub_ubay=exp_te_glm_re_sub_ubay=c()
        pi_exp_glm=matrix(rep(0,(case_samp+con_samp)*length(real_qian_node)),ncol=length(real_qian_node))
        pi_exp_glm_ubay=matrix(rep(0,(case_samp+con_samp)*length(real_qian_node)),ncol=length(real_qian_node))
        
        
        
        fit_glm_con=MGLMfit(data.frame(tree_table[1:case_samp,match(real_qian_node,colnames(tree_table))]),dist = "DM")
        alpha_con_glm=fit_glm_con@estimate
        
        fit_glm_case=MGLMfit(data.frame(tree_table[(case_samp+1):((case_samp+con_samp)),match(real_qian_node,colnames(tree_table))]),dist = "DM")
        alpha_case_glm=fit_glm_case@estimate
        
        alpha_matrix_glm = rbind(matrix(rep(alpha_con_glm, case_samp),byrow = TRUE,ncol = length(real_qian_node)),
                                 matrix(rep(alpha_case_glm, con_samp),byrow = TRUE, ncol = length(real_qian_node)))
        for (n in 1:(case_samp+con_samp)) {
          pi_exp_glm[n,]= unlist((tree_table[n,match(real_qian_node,colnames(tree_table))] + alpha_matrix_glm[n,]) / (sum(tree_table[n, match(real_qian_node,colnames(tree_table))]) + sum(alpha_matrix_glm[n,])))
          pi_exp_glm_ubay[n,]= unlist((tree_table[n,match(real_qian_node,colnames(tree_table))] +0.5) / (sum(tree_table[n, match(real_qian_node,colnames(tree_table))]) + +sum(rep(1/2,length(real_qian_node)))))
          
        }
        
        cij_expi_glm_re=apply(t(pi_exp_glm), 2, function(x){log2(x) - mean(log2(x))})
        cij_expi_glm_ubay=apply(t(pi_exp_glm_ubay), 2, function(x){log2(x) - mean(log2(x))})
        
        for( t in 1:length(real_qian_node)){
          exp_te_glm_re_sub[t] = t.test(cij_expi_glm_re[t,1:case_samp], cij_expi_glm_re[t,(case_samp+1):(case_samp+con_samp)])$p.value
          exp_te_wil_glm_re_sub[t] = wilcox.test(cij_expi_glm_re[t,1:case_samp], cij_expi_glm_re[t,(case_samp+1):(case_samp+con_samp)])$p.value
          
          
          
          exp_te_glm_re_sub_ubay[t] = t.test(cij_expi_glm_ubay[t,1:case_samp], cij_expi_glm_ubay[t,(case_samp+1):(case_samp+con_samp)])$p.value
          exp_te_wil_glm_re_sub_ubay[t] = wilcox.test(cij_expi_glm_ubay[t,1:case_samp], cij_expi_glm_ubay[t,(case_samp+1):(case_samp+con_samp)])$p.value
          
        } 
        names(exp_te_glm_re_sub)=names(exp_te_wil_glm_re_sub)=names(exp_te_glm_re_sub_ubay)=names(exp_te_wil_glm_re_sub_ubay)=as.character(real_qian_node)
        exp_te_glm_re=append(exp_te_glm_re,exp_te_glm_re_sub)
        exp_te_wil_glm_re=append(exp_te_wil_glm_re,exp_te_wil_glm_re_sub)
        exp_te_glm_re_ubay=append(exp_te_glm_re_ubay,exp_te_glm_re_sub_ubay)
        exp_te_wil_glm_re_ubay=append(exp_te_wil_glm_re_ubay,exp_te_wil_glm_re_sub_ubay)
        
      }
      
      
      reject_c_pri_t=reject_l_pri_t=reject_c_alpha_t=reject_l_alpha_t=c()
      reject_c_pri_wil=reject_l_pri_wil=reject_c_alpha_wil=reject_l_alpha_wil=c()
      
      left_node=setdiff(1:p,all_qian)
      left_pval=p.val[match(left_node,names(p.val))]
      
      exp_te_wil_glm_re=exp_te_wil_glm_re[match(all_qian,names(exp_te_wil_glm_re))]
      
      exp_te_glm_re=exp_te_glm_re[match(all_qian,names(exp_te_glm_re))]
      
      exp_te_wil_glm_re_ubay=exp_te_wil_glm_re_ubay[match(all_qian,names(exp_te_wil_glm_re_ubay))]
      
      exp_te_glm_re_ubay=exp_te_glm_re_ubay[match(all_qian,names(exp_te_glm_re_ubay))]
      
      
      if(test_method=="t.test"&normal_method=="glm"){
        padj_alpha_t=names(which(p.adjust(c(exp_te_glm_re,left_pval))<=cutf))
        final.p=p.adjust(c(exp_te_glm_re,left_pval))
      }
      
      
      if(test_method=="wilcox"&normal_method=="glm"){
        padj_alpha_t=names(which(p.adjust(c(exp_te_wil_glm_re,left_pval))<=cutf))
        final.p=p.adjust(c(exp_te_wil_glm_re,left_pval))
        
      }
      
      
      if(test_method=="t.test"&normal_method=="prior"){
        padj_alpha_t=names(which(p.adjust(c(exp_te_glm_re_ubay,left_pval))<=cutf))
        
        final.p=p.adjust(c(exp_te_glm_re_ubay,left_pval))
        
      }
      
      if(test_method=="wilcox"&normal_method=="prior"){
        padj_alpha_t=names(which(p.adjust(c(exp_te_wil_glm_re_ubay,left_pval))<=cutf))
        final.p=p.adjust(c(exp_te_wil_glm_re_ubay,left_pval))
        
        
      }
    }
    
    
    
  }
  return(c(final.p))
}




rank_glm_t=find_dif(exp_te_glm,cutf=0.05,leaf=leaf,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=gen_p,case_samp=length(sam_id),con_samp=length(healthy_id),test_method="t.test",normal_method="glm")



rank_glm_wil=find_dif((exp_te_wil_glm),cutf=0.05,leaf=leaf,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=gen_p,case_samp=length(sam_id),con_samp=length(healthy_id),test_method="wilcox",normal_method="glm")


rank_pri_t=find_dif(exp_te_glm_ubay,cutf=0.05,leaf=leaf,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=gen_p,case_samp=length(sam_id),con_samp=length(healthy_id),test_method="t.test",normal_method="prior")



rank_pri_wil=find_dif((exp_te_wil_glm_ubay),cutf=0.05,leaf=leaf,taxa_index=taxa_index,tree_table=tree_table,tree=tree,p=gen_p,case_samp=length(sam_id),con_samp=length(healthy_id),test_method="wilcox",normal_method="prior")













diff.otu1=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(count_t))<0.05)))]
mat1_mon612=match(diff.otu1,imp_id)
imp_txt[mat1_mon612,3]

diff.otu2=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(count_rel_t))<0.05)))]
mat2_mon612=(match(diff.otu2,imp_id))
imp_txt[mat2_mon612,3]

diff.otu3=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(count_rare_t))<0.05)))]
mat3_mon612=(match(diff.otu3,imp_id))
imp_txt[mat3_mon612,3]

diff.otu4=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(exp_te_glm_ubay_nt))<0.05)))]
mat4_mon612=(match(diff.otu4,imp_id))
imp_txt[mat4_mon612,3]

diff.otu5=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(exp_te_glm_nt))<0.05)))]
mat5_mon612=(match(diff.otu5,imp_id))
imp_txt[mat5_mon612,3]

diff.otu6=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(count_wil))<0.05)))]
mat6_mon612=(match(diff.otu6,imp_id))
imp_txt[mat6_mon612,3]

diff.otu7=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(count_rel_wil))<0.05)))]
mat7_mon612=(match(diff.otu7,imp_id))
imp_txt[mat7_mon612,3]

diff.otu8=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(count_rare_wil))<0.05)))]
mat8_mon612=(match(diff.otu8,imp_id))
imp_txt[mat8_mon612,3]

diff.otu9=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(exp_te_wil_glm_ubay_nt))<0.05)))]
mat9_mon612=(match(diff.otu9,imp_id))
imp_txt[mat9_mon612,3]

diff.otu10=gen_tree$tip.label[as.numeric(names(which(sort(p.adjust(exp_te_wil_glm_nt))<0.05)))]
mat10_mon612=(match(diff.otu10,imp_id))
imp_txt[mat10_mon612,3]




pas=MRcoefs(res_raw)[,4]
names(pas)=rownames(MRcoefs(res_raw))
diff.otu11=gen_tree$tip.label[as.numeric(names(which(sort(pas)<0.05)))]
mat11_mon612=(match(diff.otu11,imp_id))
imp_txt[mat11_mon612,3]




names(re_des_raw$padj)=1:p
diff.otu12=gen_tree$tip.label[as.numeric(names(which(sort(re_des_raw$padj)<0.05)))]
mat12_mon612=(match(diff.otu12,imp_id))
imp_txt[mat12_mon612,3]



anc_ra=ANCOM(OTUdat=data.frame(tree_table_leaf,Group=c(rep(0,length(sam_id)),rep(1,length(healthy_id)))),
             sig=0.05, multcorr=3, tau=0.02, theta=0.1, repeated=FALSE)$W
diff.otu13=gen_tree$tip.label[as.numeric(ancom.test)[order(anc_ra[as.numeric(ancom.test)],decreasing = TRUE)]]
mat13_mon612=(match(diff.otu13,imp_id))
imp_txt[mat13_mon612,3]

diff.otu21=gen_tree$tip.label[as.numeric(names(which(sort((rank_pri_t))<0.05)))]
mat21_mon612=(match(diff.otu21,imp_id))
imp_txt[mat12_mon612,3]

diff.otu22=gen_tree$tip.label[as.numeric(names(which(sort((rank_glm_t))<0.05)))]
mat22_mon612=(match(diff.otu22,imp_id))
imp_txt[mat22_mon612,3]


diff.otu23=gen_tree$tip.label[as.numeric(names(which(sort((rank_pri_wil))<0.05)))]
mat23_mon612=(match(diff.otu23,imp_id))
imp_txt[mat23_mon612,3]


diff.otu24=gen_tree$tip.label[as.numeric(names(which(sort((rank_glm_wil))<0.05)))]
mat24_mon612=(match(diff.otu24,imp_id))
imp_txt[mat24_mon612,3]


real_ran=10
pg_mean1=data.frame(number=c(sum((mat1_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat2_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat3_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat4_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat5_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat13_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat12_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat11_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat21_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat22_mon612)[1:real_ran] %in% (1:real_ran))),
                    group=c("none-t","tss-t","rarefying-t","uBay-t","eBay-t","ANCOM","DESeq2","metagenomeSeq","uBay-tree-t","eBay-tree-t"))





real_ran=15

pg_mean2=data.frame(number=c(sum((mat1_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat2_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat3_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat4_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat5_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat13_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat12_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat11_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat21_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat22_mon612)[1:real_ran] %in% (1:real_ran))),
                    group=c("none-t","tss-t","rarefying-t","uBay-t","eBay-t","ANCOM","DESeq2","metagenomeSeq","uBay-tree-t","eBay-tree-t"))






real_ran=20

pg_mean3=data.frame(number=c(sum((mat1_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat2_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat3_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat4_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat5_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat13_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat12_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat11_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat21_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat22_mon612)[1:real_ran] %in% (1:real_ran))),
                    group=c("none-t","tss-t","rarefying-t","uBay-t","eBay-t","ANCOM","DESeq2","metagenomeSeq","uBay-tree-t","eBay-tree-t"))







real_ran=25

pg_mean4=data.frame(number=c(sum((mat1_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat2_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat3_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat4_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat5_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat13_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat12_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat11_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat21_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat22_mon612)[1:real_ran] %in% (1:real_ran))),
                    group=c("none-t","tss-t","rarefying-t","uBay-t","eBay-t","ANCOM","DESeq2","metagenomeSeq","uBay-tree-t","eBay-tree-t"))







pg_mean=data.frame(rbind(pg_mean1,pg_mean2,pg_mean3,pg_mean4),rank=rep(c("top 10","top 15","top 20","top 25"),each=10))
pg_mean=pg_mean[-which(pg_mean[,2]=="uBay-tree-t"),]
pg_mean<-within(pg_mean,{
  group<-factor(group,levels=c("none-t","tss-t","rarefying-t","DESeq2","ANCOM","metagenomeSeq","uBay-t","eBay-t","eBay-tree-t"))
})
pg_mean[,3]=as.factor(pg_mean[,3])


pic1=ggplot(pg_mean, aes(x=group, y=number,fill=rank))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+
  theme_bw()+
  scale_fill_manual(values=brewer.pal( 9, 'PuBu')[c(4,6,7,8)])+
  labs(x = "Methods",y = "Number of matches")+                        
  guides(fill = guide_legend(reverse = F))+                 
  theme(plot.title = element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),   
        legend.title = element_blank(),                    
        legend.text = element_text(size = 13, face = "bold"),        
        legend.position = 'right',               
        legend.key.size=unit(0.8,'cm')) +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(plot.margin=unit(c(1,0,1,1),'lines'))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17,face="bold"))

#    上 右 下 左
pic1




real_ran=10
pg_mean1=data.frame(number=c(sum((mat6_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat7_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat8_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat9_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat10_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat13_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat12_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat11_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat23_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat24_mon612)[1:real_ran] %in% (1:real_ran))),
                    group=c("none-Wilcoxon","tss-Wilcoxon","rarefying-Wilcoxon","uBay-Wilcoxon","eBay-Wilcoxon","ANCOM","DESeq2","metagenomeSeq","uBay-tree-Wilcoxon","eBay-tree-Wilcoxon"))




real_ran=15

pg_mean2=data.frame(number=c(sum((mat6_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat7_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat8_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat9_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat10_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat13_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat12_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat11_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat23_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat24_mon612)[1:real_ran] %in% (1:real_ran))),
                    group=c("none-Wilcoxon","tss-Wilcoxon","rarefying-Wilcoxon","uBay-Wilcoxon","eBay-Wilcoxon","ANCOM","DESeq2","metagenomeSeq","uBay-tree-Wilcoxon","eBay-tree-Wilcoxon"))





real_ran=20

pg_mean3=data.frame(number=c(sum((mat6_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat7_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat8_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat9_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat10_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat13_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat12_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat11_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat23_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat24_mon612)[1:real_ran] %in% (1:real_ran))),
                    group=c("none-Wilcoxon","tss-Wilcoxon","rarefying-Wilcoxon","uBay-Wilcoxon","eBay-Wilcoxon","ANCOM","DESeq2","metagenomeSeq","uBay-tree-Wilcoxon","eBay-tree-Wilcoxon"))






real_ran=25

pg_mean4=data.frame(number=c(sum((mat6_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat7_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat8_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat9_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat10_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat13_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat12_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat11_mon612)[1:real_ran] %in% (1:real_ran)),
                             sum((mat23_mon612)[1:real_ran] %in% (1:real_ran)),sum((mat24_mon612)[1:real_ran] %in% (1:real_ran))),
                    group=c("none-Wilcoxon","tss-Wilcoxon","rarefying-Wilcoxon","uBay-Wilcoxon","eBay-Wilcoxon","ANCOM","DESeq2","metagenomeSeq","uBay-tree-Wilcoxon","eBay-tree-Wilcoxon"))



pg_mean=data.frame(rbind(pg_mean1,pg_mean2,pg_mean3,pg_mean4),rank=rep(c("top 10","top 15","top 20","top 25"),each=10))
pg_mean=pg_mean[-which(pg_mean[,2]=="uBay-tree-Wilcoxon"),]
pg_mean<-within(pg_mean,{
  group<-factor(group,levels=c("none-Wilcoxon","tss-Wilcoxon","rarefying-Wilcoxon","DESeq2","ANCOM","metagenomeSeq","uBay-Wilcoxon","eBay-Wilcoxon","eBay-tree-Wilcoxon"))
})
pg_mean[,3]=as.factor(pg_mean[,3])


pic2=ggplot(pg_mean, aes(x=group, y=number,fill=rank))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+
  theme_bw()+
  scale_fill_manual(values=brewer.pal( 9, 'PuBu')[c(4,6,7,8)])+
  labs(x = "Methods",y = "Number of matches")+                        
  guides(fill = guide_legend(reverse = F))+                  
  theme(plot.title = element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),  
        legend.title = element_blank(),                    
        legend.text = element_text(size = 12.5, face = "bold"),        
        legend.position = 'right',               
        legend.key.size=unit(0.8,'cm')) +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(plot.margin=unit(c(1,0,1,1),'lines'))+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=17,face="bold"))


pic2



plot(rnorm,lty=0,bty='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(-1,2),ylim=c(-1,1))

T=venn.diagram(x= list("none-t" = as.character(mat1_mon612),"tss-t" = as.character(mat2_mon612),"rarefying-t" = as.character(mat3_mon612),"uBay-t" = as.character(mat4_mon612),
                       "eBay-t" = as.character(mat5_mon612)),
               filename = NULL,imagetype="png", col="transparent",
               fill=c("cornflowerblue","green","yellow","darkorchid1","darkorange1"),alpha = 0.50,
               cex = 1.6,
               fontfamily = "serif",
               fontface = "bold",
               cat.cex = 1.6,
               cat.default.pos = "text",
               cat.dist = c(0.05, 0.06, 0.06,0.06, 0.05),
               cat.pos = 0.2,
               cat.fontfamily = "serif",
               rotation.degree =360,
               margin = 0)
grid.draw(T)


plot(rnorm,lty=0,bty='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(-1,2),ylim=c(-1,1))

T=venn.diagram(x= list(ANCOM = as.character(mat13_mon612),DESeq2 = as.character(mat12_mon612),metagenomeSeq = as.character(mat11_mon612),
                       "eBay-t" = as.character(mat5_mon612),"eBay-tree-t" = as.character(mat22_mon612)),
               filename = NULL,imagetype="png", col="transparent",
               fill=c("cornflowerblue","green","yellow","darkorchid1","darkorange1"),alpha = 0.50,
               cex = 1.6,
               fontfamily = "serif",
               fontface = "bold",
               cat.cex = 1.6,
               cat.default.pos = "text",
               cat.dist = c(0.05, 0.06, 0.06,0.06, 0.05),
               cat.pos = 0.2,
               cat.fontfamily = "serif",
               rotation.degree =360,
               margin = 0)
grid.draw(T)





























plot(rnorm,lty=0,bty='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(-1,2),ylim=c(-1,1))

T=venn.diagram(x= list("none-Wilcoxon" = as.character(mat6_mon612),"tss-Wilcoxon" = as.character(mat7_mon612),"rarefying-Wilcoxon" = as.character(mat8_mon612),
                       "uBay-Wilcoxon" = as.character(mat9_mon612),"eBay-Wilcoxon" = as.character(mat10_mon612)),
               filename = NULL,imagetype="png", col="transparent",
               fill=c("cornflowerblue","green","yellow","darkorchid1","darkorange1"),alpha = 0.50,
               cex = 1.6,
               fontfamily = "serif",
               fontface = "bold",
               cat.cex = 1.6,
               cat.default.pos = "text",
               cat.dist = c(0.05, 0.06, 0.06,0.06, 0.05),
               cat.pos = 0.2,
               cat.fontfamily = "serif",
               rotation.degree =360,
               margin = 0)
grid.draw(T)


plot(rnorm,lty=0,bty='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(-1,2),ylim=c(-1,1))

T=venn.diagram(x= list(ANCOM = as.character(mat13_mon612),DESeq2 = as.character(mat12_mon612),metagenomeSeq = as.character(mat11_mon612),
                       "eBay-Wilcoxon" = as.character(mat10_mon612),"eBay-tree-Wilcoxon" = as.character(mat24_mon612)),
               filename = NULL,imagetype="png", col="transparent",
               fill=c("cornflowerblue","green","yellow","darkorchid1","darkorange1"),alpha = 0.50,
               cex = 1.6,
               fontfamily = "serif",
               fontface = "bold",
               cat.cex = 1.6,
               cat.default.pos = "text",
               cat.dist = c(0.05, 0.06, 0.06,0.06, 0.05),
               cat.pos = 0.2,
               cat.fontfamily = "serif",
               rotation.degree =360,
               margin = 0)
grid.draw(T)








