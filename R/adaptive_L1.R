adaptive_L1 = function(allgene_res,df_true = NULL){
  
  allgene_res <- remove_na_mrk(allgene_res)$keep
  
  phat_all <- unlist(lapply(allgene_res, getphat),recursive = FALSE)
  names(phat_all) <-gsub("TRUE.","TRUE_",names(phat_all))
  names(phat_all) <-gsub(".","_",names(phat_all), fixed = TRUE)
  names(phat_all) <-gsub("p_value","p.value",names(phat_all))
  if(!is.null(df_true)){
    phat_all <- lapply(phat_all,ord_name,df_true = df_true)
  }else{
    phat_all <- lapply(phat_all,ord_name,df_true = phat_all[[1]])
  }
  
  phat_allgene <- lapply(phat_all,sum_to_one)
  
  mean_phat <-  Reduce(`+`, phat_allgene) / length(phat_allgene)
  
  l1_e_all = list()
  k_all = c()
  gre = c()
  
  
  l1_0 = 1/l1_prep(pre = mean_phat,phatlist = phat_allgene)
  pre_w = l1_0/sum(l1_0)
  
  
  PL1 = 0
  for (j in 1:length(phat_allgene)) {
    PL1 = phat_allgene[[j]] * as.vector(pre_w[j]) + PL1
  }
  
  err_f = l1_prep(PL1,phat_allgene)
  
  i = 1
  e=Inf
  while (e >1e-6) {
    
    new_w=(1/err_f)/sum((1/err_f))
    i = i+1
    new_p = 0
    for (j in 1:length(phat_allgene)) {
      new_p = phat_allgene[[j]] * as.vector(new_w[j]) + new_p
    }
    new_p = new_p/rowSums(new_p)
    #e = mean(new_w-pre_w)
    e = mean(abs(c (as.matrix(new_w-pre_w))))
    
    pre_w = new_w
    err_f = l1_prep(pre=new_p,phat_allgene)
    
  }
  
  return(list(ensemble_p = new_p, tmp_weight = new_w,i=i,mean_phat = mean_phat))
}

l1_prep = function(pre,phatlist){
  tmp  = sapply(phatlist, function(x){
    frobenius.norm(pre-x)
  })
  return(tmp)
}