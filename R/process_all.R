#This function is used for process output
#
# param  allfold_res List of output from gen_all_res or gen_all_res_list of two fold
#
# Each element of list is result for one scenario, which contains element \code{a}, \code{p}, and \code{ensemble}.
#
# \itemize{
# \item{'a'}{  Contains p_hat (list, deconvolution results of each fold), p_truth (list, true cell type proportion), n (list, number of markers of each fold), time, markers (list, markers gene of each fold), all_time (computation time)}
# \item{'p'}{  Parameters}
# \item{'ensemble'}{  Testing data for each fold}
# }
#
# @param  allgene_res List of output from gen_all_res or gen_all_res_list of one fold
#
# Each element of list is result for one scenario, which contains element \code{a}, \code{p}, and \code{ensemble}.
#
# \itemize{
# \item{'a'}{  Contains p_hat (list of deconvolution results), p_truth (list, true cell type proportion), n (list, number of markers of each fold), time, markers (list, markers gene of each fold), all_time (computation time)}
# \item{'p'}{  Parameters}
# \item{'ensemble'}{  Testing data for}
# }
#
# @param df_true Dataframe, True cell type expression, samples*cell type
# @param trueMet Character, name of bulk data
# @param feature_selection Logical, conducting feature selection or not
# @param lambda (Optional)  tuning parameter for feature selection
# @param get_metric Logical, calculating metric or not
#param rm_dmeths Logical, remove several deconvolution methods or not when ensemble
#


process_all <- function(allfold_res,allgene_res,df_true = NULL,trueMet = NULL,feature_selection = TRUE,lambda = NULL,incl_lasso = FALSE,get_metric = TRUE,rm_dmeths = TRUE,r2 = FALSE,r2_ind = NULL,scale = "log",transformation = "CPM",optimization = FALSE){


  # change previous scale name
  allfold_res <- lapply(allfold_res, change_scale_name)
  allgene_res <- lapply(allgene_res, change_scale_name)
  # remove incomplete markers scenarios.
  allfold_res <- remove_na_mrk(allfold_res, cv = T)$keep
  allgene_res <- remove_na_mrk(allgene_res)$keep

  phat_all <- unlist(lapply(allgene_res, getphat),recursive = FALSE)
  names(phat_all) <-gsub("TRUE.","TRUE_",names(phat_all))
  names(phat_all) <-gsub(".","_",names(phat_all), fixed = TRUE)
  names(phat_all) <-gsub("p_value","p.value",names(phat_all))
  # phat_all = lapply(phat_all,function(x){
  #   colnames(x) = gsub("log(sig_matrix)","",colnames(x),fixed = T)
  #   colnames(x) = gsub("sig_matrix","",colnames(x))
  #   return(x)
  # })
  phat_all <- lapply(phat_all,ord_name,df_true = df_true)
  phat_all <- lapply(phat_all,sum_to_one)


  if(is.null(df_true)){
    df_true <- allgene_res[[1]][["a"]][["p_truth"]]
  }


  output <- process_basis(allfold_res,allgene_res,phat_all = phat_all,feature_selection = feature_selection,lambda = lambda,incl_lasso = incl_lasso,rm_dmeths = rm_dmeths,scale = scale,transformation = transformation,r2 =r2, r2_ind = r2_ind,optimization = optimization)


  mean_phat <-  Reduce(`+`, phat_all) / length(phat_all)

  # append the ensemble result to original onefold result
  allgene_res[[1]] <- append_ensemblewt(allgene_res[[1]],phat = output[[1]])
  allgene_res[[1]] <- append_ensemblewt(allgene_res[[1]],phat =list(Mean = mean_phat))

  if(is.null(trueMet)){
    trueMet <- "none"
  }

  if(get_metric){
    cl <- makeCluster(50)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    re_org_out_all <- foreach(i = 1:length(allgene_res)) %dopar% {
      ress <- get_metric_wrap(allgene_res[[i]], true = df_true,trueMet = trueMet)
      ress <- ExportRes_all_new(ress)
      return(ress)
    }
    stopCluster(cl)
    re_org_out_all <- bind_rows(re_org_out_all)
    # re_org_res_new <- suppressMessages(lapply(allgene_fin,get_metric_wrap, true = df_true,trueMet = trueMet))
    # re_org_out <-suppressMessages(lapply(re_org_res_new,ExportRes_all_new))
    # re_org_out_all <- bind_rows(re_org_out)



    re_org_out_all$NormalizationTC <- paste0(re_org_out_all$TNormalization,"_",re_org_out_all$CNormalization)
    re_org_out_all$ind <- paste0(re_org_out_all$NormalizationTC,"_",re_org_out_all$Scale)

    neworder <- c("2fold_ensemble_lasso","2fold_ensemble_isse","2fold_ensemble_nnls","2fold_ensemble_isae","2fold_ensemble_rmsd","2fold_ensemble_solnp","2fold_ensemble_mean","1fold_ensemble_lasso","1fold_ensemble_isse","1fold_ensemble_nnls","1fold_ensemble_isae","1fold_ensemble_rmsd","1fold_ensemble_solnp","1fold_ensemble_mean","Mean","BayesPrism_GEP","BayesPrism_scRNA","dtangle_scRNA","dtangle2","dtangle_GEP","dtangle","CIBERSORT","MuSiC","BisqueRNA","LS Fit","EPIC","TOASTP","logRegression","linearRegression","deconf", "ssFrobenius", "ssKL", "DSA","Q Prog","ICeDT","DeconRNASeq")
    re_org_out_all$Method <- factor(re_org_out_all$Method, levels=neworder)

    re_org_out_all <- re_org_out_all %>%
      mutate(miss_cat = ifelse(miss == "NA","No miss",ifelse(str_count(miss, ',')==0,"one", ifelse(str_count(miss, ',') ==1, "two", ifelse(is.na(miss),"No miss","three")) )),miss_cat = as.character(miss_cat))


    res_clean <- re_org_out_all %>% ungroup() %>%select(Scale,Method,RMSE_Total,MAE_Total,Spearman_median,Marker.Method,data_name,NormalizationTC,n_markers,ind) %>% unique()


    return(list(res_all = re_org_out_all,res_clean = res_clean,weight = output$weight))
  }else{
    return(allgene_fin)
  }


}

######### process_basis ########
process_basis = function(allfold_res,allgene_res,phat_all,feature_selection,lambda ,incl_lasso,rm_dmeths = TRUE,scale = scale,transformation = transformation,r2 =r2, r2_ind = r2_ind,optimization = optimization){

  res_all = list()

  # extract p_hat and ensemble test_bulk
  for (i in 1:2) {
    res_all[[i]]=lapply(allfold_res,ex_fold,i)
  }

  res_all_gene = lapply(allgene_res, ex_fold,1)


  # remove dmeths : TOASTP, dtangle2, LS Fit
  if(rm_dmeths){
    res_all = lapply(res_all,rm_dmeths)
    res_all_gene = rm_dmeths(res_all_gene)
  }


  res_all_up <- lapply(res_all, get_all_phatyhat,scale = scale,transformation = transformation,r2 = r2, r2_ind)
  res_all_gene_up <- get_all_phatyhat(res_all_gene,scale = scale,transformation = transformation,r2 = r2, r2_ind)
  # res_allgene_up = get_all_phat(res_allgene)

  # get ensemble weights
  res_all_weight = lapply(res_all_up, get_weight_cv,feature_selection = feature_selection,lambda = lambda,incl_lasso = incl_lasso,optimization = optimization)
  tmp = nrow(res_all_weight[[1]][[1]])+1
  res_all_weight[[1]][[1]][tmp,] = 1/ncol(res_all_weight[[1]][[1]])
  rownames(res_all_weight[[1]][[1]])[tmp] = "mean_sub"
  res_all_weight[[2]][[1]][tmp,] = 1/ncol(res_all_weight[[2]][[1]])
  rownames(res_all_weight[[2]][[1]])[tmp] = "mean_sub"
  res_all_weight_allgene = get_weight_cv(res_all_gene_up,feature_selection = feature_selection,lambda = lambda,incl_lasso = incl_lasso,optimization = optimization)
  res_all_weight_allgene[[1]][tmp,] = 1/ncol(res_all_weight_allgene[[1]])
  rownames(res_all_weight_allgene[[1]])[tmp] = "mean_sub"
  # calculate final ensemble phat
  phat = get_phat_fin(res_all_weight,phat_all,fold = 2)
  names(phat) = paste0("2fold","_",names(phat))
  phat_allgene = get_phat_fin(res_all_weight_allgene,phat_all,fold = 1)
  names(phat_allgene) = paste0("1fold","_",names(phat_allgene))
  phat = append(phat,phat_allgene)
  return(list(phat = phat, weight = res_all_weight))
}



######## append_ensemblewt ########

append_ensemblewt = function(obj,phat){
  obj[["a"]][["p_hat"]][[1]] = append(obj[["a"]][["p_hat"]][[1]], phat)
  return(obj)
}



####### rm_dmeths #########

rm_dmeths = function(x){
  # remove dmeths
  for (i in 1:length(x)) {
    x[[i]][["p_hat"]]$TOASTP = NULL
    x[[i]][["p_hat"]]$'LS Fit' = NULL
    x[[i]][["p_hat"]]$DSA = NULL
    x[[i]][["p_hat"]]$deconf = NULL
    x[[i]][["p_hat"]]$'Q Prog' = NULL
    x[[i]][["p_hat"]]$ssFrobenius = NULL
    x[[i]][["p_hat"]]$ssKL = NULL
    x[[i]][["p_hat"]]$linearRegression = NULL
    x[[i]][["p_hat"]]$logRegression = NULL
    x[[i]][["p_hat"]]$dtangle = NULL
  }
  return(x)
}

####### ex_fold #########
ex_fold = function(x,i){
  # function to extract p_hat and ensemble
  res = list()
  res$p_hat = x[["a"]][["p_hat"]][[i]]
  res$p = x[["p"]]
  #res$ensemble = x[["ensemble"]][[i]]
  return(res)
}

######## remove_na_mrk ########

remove_na_mrk = function(res_list, cv =FALSE){
  if(cv){
    check = lapply(res_list , function(obj) append(obj[["a"]][["markers"]][[1]],obj[["a"]][["markers"]][[2]]))
  }else{
    check = lapply(res_list , function(obj) obj[["a"]][["markers"]][[1]])
  }

  keep =res_list[which(sapply(check , function(x) sum(is.na(unlist(x, recursive = F)))) == 0)]

  re_names = names(res_list)[which(sapply(check , function(x) sum(is.na(unlist(x, recursive = F)))) != 0)]
  return(list(keep= keep,re_names = re_names))
}

####### remove_na_phat #########

remove_na_phat = function(res_list, cv =FALSE){
  if(cv){
    check = lapply(res_list , function(obj) append(obj[["a"]][["p_hat"]][[1]],obj[["a"]][["p_hat"]][[2]]))
  }else{
    check = lapply(res_list , function(obj) obj[["a"]][["p_hat"]][[1]])
  }

  keep =res_list[which(sapply(check , function(x) sum(is.na(unlist(x, recursive = F)))) == 0)]

  re_names = names(res_list)[which(sapply(check , function(x) sum(is.na(unlist(x, recursive = F)))) != 0)]
  return(list(keep= keep,re_names = re_names))
}

#### getphat ########

getphat = function(res){
  res = res[["a"]][["p_hat"]][[1]]
  return(res)
}

#### ExportRes_all_new ########

ExportRes_all_new = function(ori_res){
  params <- expand.grid(ori_res[[5]], stringsAsFactors = FALSE)
  params = bind_cols(lapply(params,unlist))
  res = ori_res[[2]]
  params_new = params
  for (i in 2:nrow(res)) {
    params_new[i,] = params
  }

  res_mis = ori_res[[3]]
  mis_dif = data.frame(Method = names(res_mis))
  d2 = NULL
  b = sapply(res_mis, paste, collapse=",")
  for (i in 1:length(res_mis)) {
    dd = b[i]
    d2 = c(d2,dd)
  }
  mis_dif$miss = d2


  output = cbind(params_new,res) %>% full_join(mis_dif)
  output$ind = ori_res[[4]]
  if("TScaling" %in% colnames(output)){
    output = output %>% group_by(TScaling,CScaling,Method) %>%
      fill(Spearman_min_celltype, Spearman_max_celltype,.direction ="downup") %>% unique()
  }else{
    output = output %>% group_by(TNormalization,CNormalization,Method) %>%
      fill(Spearman_min_celltype, Spearman_max_celltype,.direction ="downup") %>% unique()
  }

  return(output)
}

change_scale_name = function(obj){
  colnames(obj[["p"]])[1] = "TNormalization"
  colnames(obj[["p"]])[2] = "CNormalization"
  return(obj)

}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
