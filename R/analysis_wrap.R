#' This function is used to get ensemble results
#'
#' @param allfold_res all gene deconvolution results
#' @param allgene_res two folds deconvolution results
#' @param df_true Optional true cell type fraction
#' @param trueMet character 
#' 
#' @import ggpubr
#' @import Rsolnp
#'
#' @export

get_res_wrap <- function(allfold_res,allgene_res,df_true,trueMet ,feature_selection = FALSE,lambda = NULL,incl_lasso = FALSE,get_metric = TRUE,r2 = FALSE,r2_ind = NULL,scale = "log",transformation = "none",optimization = TRUE){

    # remove incomplete markers scenarios.
    allfold_res <- remove_na_mrk(allfold_res, cv = T)$keep
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
    
    phat_all <- lapply(phat_all,sum_to_one)

    res_all <- list()

    # extract p_hat and ensemble test_bulk
    for (i in 1:2) {
        res_all[[i]]<-lapply(allfold_res,ex_fold,i)
    }

    res_all_gene <- lapply(allgene_res, ex_fold,1)

    res_all_up <- lapply(res_all, get_all_phatyhat,scale = scale,transformation = transformation,r2 = r2, r2_ind,p_true  = df_true)
    res_all_gene_up <- get_all_phatyhat(res_all_gene,scale = scale,transformation = transformation,r2 = r2, r2_ind,p_true  = df_true,allgene = TRUE)

#     p1 = ggarrange(res_all_up[[1]]$p1,res_all_up[[2]]$p1,common.legend = T)
#     p1 = annotate_figure(p1,
# bottom = text_grob("left small size fold, right large size fold", color = "blue",
# hjust = 1, x = 1, face = "italic", size = 10)
# )
# 
# p2 = ggarrange(res_all_up[[1]]$p2,res_all_up[[2]]$p2,common.legend = T)
# p2 =annotate_figure(p2,
#                 bottom = text_grob("left small size fold, right large size fold, rm", color = "blue",
#                                    hjust = 1, x = 1, face = "italic", size = 10)
# )


    rmse_fold <- lapply(res_all_up, get_rmse)

    rmse_gene <- get_rmse(res_all_gene_up)

    res_all_weight <- lapply(rmse_fold, get_weight_cv,feature_selection = FALSE,lambda = NULL,incl_lasso = FALSE,optimization = TRUE)


    tmp <- nrow(res_all_weight[[1]][[1]])+1
    res_all_weight[[1]][[1]][tmp,] <- 1/ncol(res_all_weight[[1]][[1]])
    rownames(res_all_weight[[1]][[1]])[tmp] <- "mean_sub"
    res_all_weight[[2]][[1]][tmp,] <- 1/ncol(res_all_weight[[2]][[1]])
    rownames(res_all_weight[[2]][[1]])[tmp] <- "mean_sub"
    res_all_weight_allgene <- get_weight_cv(rmse_gene,feature_selection = FALSE,lambda = NULL,incl_lasso = FALSE,optimization = TRUE)
    res_all_weight_allgene[[1]][tmp,] <- 1/ncol(res_all_weight_allgene[[1]])
    rownames(res_all_weight_allgene[[1]])[tmp] <- "mean_sub"

    phat <- get_phat_fin(res_all_weight,phat_all,fold = 2)
    names(phat) <- paste0("2fold","_",names(phat))
    phat_allgene <- get_phat_fin(res_all_weight_allgene,phat_all,fold = 1)
    names(phat_allgene) <- paste0("allgene","_",names(phat_allgene))
    phat <- append(phat,phat_allgene)

    phat_sub2 <- get_phat_fin(res_all_weight[[2]],phat_all,fold = 1)
    names(phat_sub2) <- paste0("1fold_lar","_",names(phat_sub2))

    phat_sub1 <- get_phat_fin(res_all_weight[[1]],phat_all,fold = 1)
    names(phat_sub1) <- paste0("1fold_smal","_",names(phat_sub1))

    phat_sub = append(phat_sub1,phat_sub2)
    phat = append(phat,phat_sub)
    ### rm certain dmeths
    rmse_fold_rm <- lapply(res_all_up, get_rmse,rm = T)

    rmse_gene_rm <- get_rmse(res_all_gene_up,rm = T)

    res_all_weight_rm <- lapply(rmse_fold_rm, get_weight_cv,feature_selection = FALSE,lambda = NULL,incl_lasso = FALSE,optimization = TRUE)


    tmp <- nrow(res_all_weight_rm[[1]][[1]])+1
    res_all_weight_rm[[1]][[1]][tmp,] <- 1/ncol(res_all_weight_rm[[1]][[1]])
    rownames(res_all_weight_rm[[1]][[1]])[tmp] <- "mean_sub"
    res_all_weight_rm[[2]][[1]][tmp,] <- 1/ncol(res_all_weight_rm[[2]][[1]])
    rownames(res_all_weight_rm[[2]][[1]])[tmp] <- "mean_sub"
    res_all_weight_rm_allgene <- get_weight_cv(rmse_gene_rm,feature_selection = FALSE,lambda = NULL,incl_lasso = FALSE,optimization = TRUE)
    res_all_weight_rm_allgene[[1]][tmp,] <- 1/ncol(res_all_weight_rm_allgene[[1]])
    rownames(res_all_weight_rm_allgene[[1]])[tmp] <- "mean_sub"

    phat2 <- get_phat_fin(res_all_weight_rm,phat_all,fold = 2)
    names(phat2) <- paste0("2foldrm","_",names(phat2))
    phat_allgene2 <- get_phat_fin(res_all_weight_rm_allgene,phat_all,fold = 1)
    names(phat_allgene2) <- paste0("allgenerm","_",names(phat_allgene2))

    phat <- append(phat,phat2)
    phat <- append(phat,phat_allgene2)

    phat_sub2 <- get_phat_fin(res_all_weight_rm[[2]],phat_all,fold = 1)
    names(phat_sub2) <- paste0("1foldrm_lar","_",names(phat_sub2))

    phat_sub1 <- get_phat_fin(res_all_weight_rm[[1]],phat_all,fold = 1)
    names(phat_sub1) <- paste0("1foldrm_smal","_",names(phat_sub1))

    phat_sub = append(phat_sub1,phat_sub2)
    phat = append(phat,phat_sub)


    mean_phat <-  Reduce(`+`, phat_all) / length(phat_all)

    # append the ensemble result to original onefold result
    allgene_res[[1]] <- append_ensemblewt(allgene_res[[1]],phat = phat)
    allgene_res[[1]] <- append_ensemblewt(allgene_res[[1]],phat =list(Mean = mean_phat))


    cl <- makeCluster(20)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    re_org_out_all <- foreach(i = 1:length(allgene_res), .export=c("get_metric_wrap", "ExportRes_all_new","getOutput_ensemble_cv","ord_name","sum_to_one"),.packages=c("dplyr","reshape2","tidyverse")) %dopar% {
        ress <- get_metric_wrap(allgene_res[[i]], true = df_true,trueMet =trueMet )
        ress <- ExportRes_all_new(ress)
        return(ress)
    }
    stopCluster(cl)
    re_org_out_all <- bind_rows(re_org_out_all)

    re_org_out_all$NormalizationTC <- paste0(re_org_out_all$TNormalization,"_",re_org_out_all$CNormalization)
    re_org_out_all$ind <- paste0(re_org_out_all$NormalizationTC,"_",re_org_out_all$Scale)

    return(list(re_org_out_all = re_org_out_all, rmse_fold = rmse_fold,rmse_gene=rmse_gene,res_all_gene_up = res_all_gene_up))
}



get_rmse <- function(res,rm = FALSE){
    rmse <- res[[1]][["rmse"]]
    rmse_ind <- names(rmse)


    rmse <- tibble::enframe(rmse)
    rmse$ind <- rmse_ind
    rmse <- rmse %>% mutate(ind = gsub("BayesPrism_GEP","BayesPrism",ind),
                            sce_name = stringr::str_match(ind, '(.*)_(.*)')[,-1][,1],
                            ind = stringr::str_match(ind, '(.*)_(.*)')[,-1][,2]
    )

    rmse_min <-rmse %>%
       # filter(str_detect(name,"none_none|TPM_TPM|CPM_CPM|TPM_CPM|all|TMM_TMM|TMM_CPM")) %>%
        filter(!(str_detect(name,"Mathys|Welch|Vel|Li_a")&(str_detect(ind,"TPM_TPM")))) %>%
        group_by(ind) %>%
        mutate(min_rmse = min(value)) %>%
        filter(value == min_rmse) %>%
        filter(ind != "GEP")

    r2 <- res[[1]][["r2"]]
    r2_ind <- names(r2)

    r2 <- tibble::enframe(r2)
    r2$ind <- r2_ind
    r2 <- r2 %>% mutate(ind = gsub("BayesPrism_GEP","BayesPrism",ind),
                        sce_name = stringr::str_match(ind, '(.*)_(.*)')[,-1][,1],
                        ind = stringr::str_match(ind, '(.*)_(.*)')[,-1][,2]
    )

    r2_max <-r2 %>%
     #   filter(str_detect(name,"none_none|TPM_TPM|CPM_CPM|TPM_CPM|all")) %>%
        filter(!(str_detect(name,"Mathys|Welch|Vel|Li_a")&(str_detect(ind,"TPM_TPM")))) %>%
        group_by(ind) %>%
        mutate(max_r2 = max(value)) %>%
        filter(value == max_r2) %>%
        filter(ind != "GEP")

    if(rm){
        rmse_min <-rmse_min %>%filter(!ind %in% c("TOASTP","LS Fit","DSA","deconf","Q Prog","ssFrobenius","ssKL","linearRegression","logRegression"))}

    rmse_min <- rmse_min[order(rmse_min$ind, rmse_min$value, decreasing=TRUE),]
    rmse_minClean <- rmse_min[!duplicated(rmse_min$ind),]
    r2_max <- rmse_min[order(r2_max$ind, r2_max$value, decreasing=TRUE),]
    r2_maxClean <- r2_max[!duplicated(r2_max$ind),]

    yhat_all <- res[[1]]$yhat_all[rmse_minClean$name]
    # names(yhat_all) <- stringr::str_match(names(yhat_all), '(.*)_(.*)')[,-1][,2]
    test_bulk <- res[[1]]$test_bulk
    return(list(rmse = rmse, rmse_min = rmse_minClean,yhat_all = yhat_all,test_bulk = test_bulk,r2 = r2,r2_max = r2_maxClean))
}
