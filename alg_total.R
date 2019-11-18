fdr_power_repair=function(p.val,cutf,dif,leaf=leaf,func_dtm=func_dtm,taxa_index,tree_table,tree,p,case_samp,con_samp,test_method,normal_method){
  library(doParallel)
  dirmult=dirmult
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
        real_qian_sub=real_node
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
        real_node=real_qian_do2[l]
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
        pi_exp_glm=matrix(rep(0,samps*length(real_qian_node)),ncol=length(real_qian_node))
        pi_exp_glm_ubay=matrix(rep(0,samps*length(real_qian_node)),ncol=length(real_qian_node))
        
        
        
        fit_glm_con=MGLMfit(data.frame(tree_table[,match(real_qian_node,colnames(tree_table))]),dist = "DM")
        alpha_con_glm=fit_glm_con@estimate
        
       
        
        alpha_matrix_glm = rbind(matrix(rep(alpha_con_glm, case_samp+con_samp),byrow = TRUE,ncol = length(real_qian_node)))
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
        TP1 = length(intersect(padj_alpha_t, dif))
        FP1 = length(padj_alpha_t) - TP1
        power_test = TP1/length(dif)
        fdr_test= ifelse(length(padj_alpha_t) == 0, 0, FP1/length(padj_alpha_t))
      }
      
      
      if(test_method=="wilcox"&normal_method=="glm"){
        padj_alpha_wil=names(which(p.adjust(c(exp_te_wil_glm_re,left_pval))<=cutf))
        
        TP2 = length(intersect(padj_alpha_wil, dif))
        FP2 = length(padj_alpha_wil) - TP2
        power_test= TP2/length(dif)
        fdr_test = ifelse(length(padj_alpha_wil) == 0, 0, FP2/length(padj_alpha_wil))
      }
      
      
      if(test_method=="t.test"&normal_method=="prior"){
        padj_alpha_t_ubay=names(which(p.adjust(c(exp_te_glm_re_ubay,left_pval))<=cutf))
        
        
        TP3 = length(intersect(padj_alpha_t_ubay, dif))
        FP3 = length(padj_alpha_t_ubay) - TP3
        power_test = TP3/length(dif)
        fdr_test= ifelse(length(padj_alpha_t_ubay) == 0, 0, FP3/length(padj_alpha_t_ubay))
      }
      
      if(test_method=="wilcox"&normal_method=="prior"){
        padj_alpha_wil_ubay=names(which(p.adjust(c(exp_te_wil_glm_re_ubay,left_pval))<=cutf))
        
        TP4 = length(intersect(padj_alpha_wil_ubay, dif))
        FP4 = length(padj_alpha_wil_ubay) - TP4
        power_test= TP4/length(dif)
        fdr_test = ifelse(length(padj_alpha_wil_ubay) == 0, 0, FP4/length(padj_alpha_wil_ubay))
      }
    }
    
    
    
  }
  return(c(power_test,fdr_test))
}
