# This function is used for generate metric
# param  res_list : list, each element has a, p, and ensemble
# param true : true fraction
# param trueMet : Bulk data name
# param cmap : named cell type
# param cal_mean : logical, calculate mean

get_metric_wrap <- function(res_list,true,trueMet,cal_mean = FALSE){

  res <- getOutput_ensemble_cv(res_list,true = true,trueMet = trueMet,cal_mean = FALSE)
  param <- paste0(res_list[["p"]][["data_name"]],"_",unlist(res_list[["p"]][["TNormalization"]]),"_",unlist(res_list[["p"]][["CNormalization"]]),"_",unlist(res_list[["p"]][["Scale"]]),"_",unlist(res_list[["p"]][["Marker.Method"]]))
  res[['param']] <- param
  res[['param_list']] <- res_list[["p"]]

  return(res)
}

getOutput_ensemble_cv <- function(res,true,trueMet,cal_mean = FALSE){
  p <- res[["p"]]
  res <- res$a$p_hat
  if(length(res) == 1){
    res <- res[[1]]
  }
  all_names <- names(res)
  res <- lapply(res, ord_name, df_true = true)
  res <- lapply(res, sum_to_one)
  K  <-length(res)

  df <- list()
  for (i in 1:length(res)) {
    u <- intersect(rownames(true),rownames(res[[i]]))
    d <- data.frame(Subject = u,Method = names(res)[i],res[[i]][u,]/rowSums(res[[i]][u,]))
    df[[i]] <- melt(d,id.vars = c('Subject','Method'))
  }
  d <- data.frame(Subject = u,Method =trueMet ,true[u,])
  true <- melt(d,id.vars = c('Subject','Method'))
  names(true) <- c("Subject","Method","CellType","p_true")
  true <- true %>% select(Subject, CellType, p_true)
  df <- bind_rows(df) %>% filter()
  names(df) <- c("Subject","Method","CellType","p_hat")
  df <- df%>% full_join(true)%>% filter(!is.na(p_hat))

  df <- df %>% filter(!is.na(p_true))

  if(cal_mean){

    # df_a <- df %>% filter(!Method%in%c("2fold_ensemble_lasso","2fold_ensemble_isse","2fold_ensemble_nnls","2fold_ensemble_isae","2fold_ensemble_rmsd","2fold_ensemble_mean","1fold_ensemble_lasso","1fold_ensemble_isse","1fold_ensemble_nnls","1fold_ensemble_isae","1fold_ensemble_rmsd","1fold_ensemble_mean")) %>%
      df_a <- df %>% filter(! str_detect(Method,"fold")) %>%filter(! str_detect(Method,"allgene")) %>%filter(! str_detect(Method,"new_alg")) %>%
      group_by(Subject,CellType) %>%
      mutate(p_mean = mean(p_hat),
             Method = "Mean") %>%
      select(-p_hat) %>%
      select(p_hat = p_mean, everything())
    df_a <- df_a[,colnames(df)]

    df <- rbind(df,df_a)
  }



  df1 <- suppressWarnings(df %>% group_by(CellType,Method) %>%
    mutate(RMSE= sqrt(mean((p_hat-p_true)^2))%>% round(.,2),
           MAE = mean(abs(p_hat-p_true)) %>% round(.,2),
           Spearman = cor(p_hat,p_true,use = 'complete',method = "spearman")%>% round(.,2),
           Spearman = if_else(is.na(Spearman), 0, Spearman)) %>%
    ungroup() %>%
    group_by(Method) %>% mutate(RMSE_Total = sqrt(mean((p_hat-p_true)^2))%>% round(.,2),
                                MAE_Total = mean(abs(p_hat-p_true)) %>% round(.,2),
                                Spearman_Total = cor(p_hat,p_true,use = 'complete',method = "spearman")%>% round(.,2)))

  df2 <- df1 %>%  group_by(Method) %>%   mutate(RMSE_median = median(unique(RMSE)) %>% round(.,2),
                                               Spearman_median=median(unique(Spearman),na.rm = T) %>% round(.,2),
                                               Spearman_sd=sd(unique(Spearman),na.rm = T) %>% round(.,2),
                                               RMSE_mean = mean(RMSE) %>% round(.,2),
                                               Spearman_mean=mean(Spearman,na.rm = T) %>% round(.,2),
                                               Spearman_max = max(Spearman,na.rm =T)%>% round(.,2),
                                               Spearman_min = min(Spearman,na.rm =T)%>% round(.,2),
                                               Spearman_max_celltype = ifelse(Spearman == Spearman_max,as.character(CellType),NA),
                                               Spearman_min_celltype = ifelse(Spearman == Spearman_min,as.character(CellType),NA))

  df3 <- list()
  for (i in 1:length(res)) {
    met_res <- res[[i]]
    if(any(colSums(met_res) == 0)){
      na_names <- colnames(met_res)[colSums(met_res) == 0]
      df3[[i]] <- na_names
      names(df3)[i] <- names(res)[i]
    }else{
      df3[[i]] <-  NA
      names(df3)[i] <- names(res)[i]
    }
  }

  return(list(com_res = df1, res = df2,mis = df3))
}
