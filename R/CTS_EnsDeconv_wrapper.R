CTS_EnsDeconv_wrapper = function(allgene_res,df_true = NULL,maxit= 1000){
  
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
  mean_phat <- sum_to_one(mean_phat)
  CTS_EnsDeconv_p_LS <- CTS_EnsDeconv_LS(What = phat_allgene,avar = F,maxit = maxit)
  
  return(list(ensemble_p = CTS_EnsDeconv_p_LS$W, mean_phat = mean_phat,ind = CTS_EnsDeconv_p_LS$ind))
}
