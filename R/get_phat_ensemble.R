# # Error in nnls(A = as.matrix(Yhat_re), b = Y_re) :
# # long vectors (argument 1) are not supported in .Fortran
# allgenes_ensemble_cv <- function(Yhat_re, Y_re, preprossessed = FALSE, w_lasso = NULL,optimization = FALSE) {
#   
#   ######## Get Yhat ##########
#   if (!preprossessed) {
#     Y_re = c(Y_re)
#     Yhat_re <- do.call(cbind, lapply(Yhat_re, function(x) c(x)))
#   }
#   
#   sse <- function(x, y) {
#     sum((x - y)^2, na.rm = T)
#   }
#   sae <- function(x, y) {
#     sum(abs(x - y), na.rm = T)
#   }
#   rmsd <- function(x, y) {
#     sqrt(mean((x - y)^2, na.rm = T))
#   }
#   
#   
#   ######## SSE #########
#   message("Searching ENSEMBLE weight by Sum of Squared Errors or Sum of Abs Errors or RMSD")
#   sses <- apply(Yhat_re, 2, function(x) {
#     sse(Y_re, x)
#   })
#   sse.wt <- 1/sses/sum(1/sses)
#   # sum of absolute errors
#   saes <- apply(Yhat_re, 2, function(x) {
#     sae(Y_re, x)
#   })
#   sae.wt <- 1/saes/sum(1/saes)
#   
#   # RMSD
#   rmsds <- apply(Yhat_re, 2, function(x) {
#     rmsd(Y_re, x)
#   })
#   rmsd.wt <- 1/rmsds/sum(1/rmsds)
#   
#   if(optimization){
#     ######## NNLS #########
#     w_nnls <- NA
#     message("Searching ENSEMBLE weight by NNLS -- Minimizing MSE of Y measurement")
#     fitnnls <- nnls(A = as.matrix(Yhat_re), b = Y_re)
#     w_nnls <- fitnnls$x
#     w_nnls <- w_nnls/sum(w_nnls)
#     
#     ########
#     message("Searching ENSEMBLE weight by SOLNP")
#     solnp.wt <- mysolnp(Y_re = Y_re,Yhat_re = Yhat_re)
#   }else{
#     w_nnls <- NA
#     solnp.wt <- NA
#   }
#   
#   
#   
#   if (is.null(w_lasso)) {
#     
#     weight.mat <- rbind(sse.wt, sae.wt, rmsd.wt,w_nnls,solnp.wt)
#     rownames(weight.mat) <- c("inverse SSE", "inverse sae", "rmsd","NNLS","solnp")
#   } else {
#     weight.mat <- rbind(sse.wt, sae.wt, rmsd.wt,w_nnls,solnp.wt, w_lasso)
#     rownames(weight.mat) <- c("inverse SSE", "inverse sae", "rmsd","NNLS","solnp", "LASSO")
#   }
#   
#   return(w_table = weight.mat)
# }
# 
# 
# ####################
# get_phat_fin = function(obj, phat_allgene, w_lasso = NULL,fold = 2,ct = FALSE) {
#   if(ct){
#     isse = cal_prop_ct(obj, 1, phat_allgene,fold = fold)
#     nnls = cal_prop_ct(obj, 4, phat_allgene,fold = fold)
#     solnp = cal_prop_ct(obj, 5, phat_allgene,fold = fold)
#     mean = NULL
#   }else{
#     isse = cal_prop(obj, 1, phat_allgene,fold = fold)
#     nnls = cal_prop(obj, 4, phat_allgene,fold = fold)
#     solnp = cal_prop(obj, 5, phat_allgene,fold = fold)
#     mean = cal_prop(obj, 6, phat_allgene,fold = fold)
#   }
#   
#   #
#   if (!is.null(w_lasso)) {
#     res_all_weight_lasso = lapply(obj, function(x) x[[1]][3, ])
#     res_all_weight_lasso = bind_rows(res_all_weight_lasso)
#     res_all_weight_lasso[is.na(res_all_weight_lasso)] = 0
#     k_lasso = nrow(res_all_weight_lasso) + 1
#     res_all_weight_lasso[k_lasso, ] = colMeans(res_all_weight_lasso)
#     res_all_weight_lasso[k_lasso, ] = res_all_weight_lasso[k_lasso,
#     ]/sum(res_all_weight_lasso[k_lasso, ])
#     phat_sub_lasso = phat_allgene[colnames(res_all_weight_lasso)]
#     
#     lasso = 0
#     for (j in 1:length(phat_sub_lasso)) {
#       lasso = phat_sub_lasso[[j]] * as.vector(res_all_weight_lasso[k_lasso, j]) + lasso
#     }
#     res = list(isse = isse, nnls = nnls,
#                lasso = lasso,solnp=solnp )
#   } else {
#     res = list(isse = isse,nnls = nnls,solnp=solnp,mean = mean)
#   }
#   
#   return(res)
# }


cal_prop = function(obj, n, phat_allgene,fold = 2) {
  if(fold<2){
    res_all_weight_obj = obj[[1]][n, ]
  }else{
    res_all_weight_obj_old = lapply(obj, function(x) x[[1]][n, ])
    # res_all_weight_obj = bind_rows(res_all_weight_obj)
    check = lapply(res_all_weight_obj_old, function(x) colnames(x))
    if(length(intersect(check[[1]],check[[2]])) == 0){
      res_all_weight_obj = cbind(res_all_weight_obj_old[[1]],res_all_weight_obj_old[[2]])
      res_all_weight_obj[2,] = NA
    }else{
      df1 = res_all_weight_obj_old[[1]][1,]
      df2 = res_all_weight_obj_old[[2]][1,]
      df1[setdiff(names(df2), names(df1))] <- NA
      df2[setdiff(names(df1), names(df2))] <- NA
      res_all_weight_obj = rbind(df1, df2)
      rm(df1,df2)
      
    }
    
  }
  
  
  
  res_all_weight_obj[is.na(res_all_weight_obj)] = 0
  k_obj = nrow(res_all_weight_obj) + 1
  res_all_weight_obj[k_obj, ] = colMeans(res_all_weight_obj)
  res_all_weight_obj[k_obj, ] = res_all_weight_obj[k_obj, ]/sum(res_all_weight_obj[k_obj,])
  
  # rm not detected in allgene
  sub_ind = intersect(colnames(res_all_weight_obj),names(phat_allgene))
  res_all_weight_obj = res_all_weight_obj[,sub_ind]
  phat_sub_obj = phat_allgene[sub_ind]
  
  ensemble_p = 0
  for (j in 1:length(phat_sub_obj)) {
    ensemble_p = phat_sub_obj[[j]] * as.vector(res_all_weight_obj[k_obj, j]) + ensemble_p
  }
  
  return(ensemble_p)
}


cal_prop_ct = function(obj, n, phat_allgene,fold = 2) {
  if(fold<2){
    res_all_weight_obj =lapply(obj[[1]],function(x) x[n, ])
    res_all_weight_obj = bind_rows(res_all_weight_obj)
  }else{
    res_all_weight_obj_old = lapply(obj, function(x) bind_rows(lapply(x[[1]],function(y) y[n, ])) )
    # res_all_weight_obj = bind_rows(res_all_weight_obj)
    check = lapply(res_all_weight_obj_old, function(x) colnames(x))
    if(length(intersect(check[[1]],check[[2]])) == 0){
      res_all_weight_obj = cbind(res_all_weight_obj_old[[1]],res_all_weight_obj_old[[2]])
      res_all_weight_obj[2,] = NA
    }else{
      df1 = res_all_weight_obj_old[[1]][1,]
      df2 = res_all_weight_obj_old[[2]][1,]
      df1[setdiff(names(df2), names(df1))] <- NA
      df2[setdiff(names(df1), names(df2))] <- NA
      res_all_weight_obj = rbind(df1, df2)
      rm(df1,df2)
      
    }
    
  }
  
  
  
  res_all_weight_obj[is.na(res_all_weight_obj)] = 0
  k_obj = nrow(res_all_weight_obj) + 1
  res_all_weight_obj[k_obj, ] = colMeans(res_all_weight_obj)
  res_all_weight_obj[k_obj, ] = res_all_weight_obj[k_obj, ]/sum(res_all_weight_obj[k_obj,])
  
  # rm not detected in allgene
  sub_ind = intersect(colnames(res_all_weight_obj),names(phat_allgene))
  res_all_weight_obj = res_all_weight_obj[,sub_ind]
  phat_sub_obj = phat_allgene[sub_ind]
  
  ensemble_p = 0
  for (j in 1:length(phat_sub_obj)) {
    ensemble_p = phat_sub_obj[[j]] * as.vector(res_all_weight_obj[k_obj, j]) + ensemble_p
  }
  
  return(ensemble_p)
}

########### get_weight_cv ##################

get_weight_cv = function(x, feature_selection = TRUE, lambda = NULL, incl_lasso = FALSE,optimization = FALSE) {
  set.seed(2020)
  Yhat = x[["yhat_all"]]
  Yhat <- do.call(cbind, lapply(Yhat, function(x) c(x)))
  Y = as.matrix(x[["test_bulk"]])
  Y = c(Y)
  
  if (feature_selection) {
    if (is.null(lambda)) {
      cv <- cv.glmnet(Yhat, Y, alpha = 1, intercept = 0)
      model <- glmnet(Yhat, Y, alpha = 1, lambda = cv$lambda.min,
                      intercept = 0)
    } else {
      model <- glmnet(Yhat, Y, alpha = 1, lambda = lambda, intercept = 0)
    }
    
    # Dsiplay regression coefficients
    w_lasso = data.matrix(coef(model))[-1, ]
    w_lasso = w_lasso[w_lasso > 0]
    ind = as.data.frame(w_lasso)
    Yhat_sub = Yhat[, rownames(ind)]
  } else {
    Yhat_sub = Yhat
  }
  
  if (incl_lasso) {
    allweights = allgenes_ensemble_cv(Yhat_re = Yhat_sub, Y_re = Y,
                                      preprossessed = TRUE, w_lasso = w_lasso,optimization = optimization)
  } else {
    allweights = allgenes_ensemble_cv(Yhat_re = Yhat_sub, Y_re = Y,
                                      preprossessed = TRUE, w_lasso = NULL,optimization = optimization)
  }
  
  return(list(allweights = as.data.frame(allweights), Yhat_re = Yhat_sub, Y_re = Y))
}

get_weight_cv_ct = function(x, feature_selection = TRUE, lambda = NULL, incl_lasso = FALSE,optimization = FALSE) {
  set.seed(2020)
  Yhat = x[["yhat_all"]]
  allweights <- list()
  for (i in 1:length(Yhat)) {
    Yhat_int <- do.call(cbind, lapply(Yhat[[i]], function(x) c(x[[i]])))
    Y = as.matrix(x[["test_bulk"]])
    Y = c(Y)
    
    if (feature_selection) {
      if (is.null(lambda)) {
        cv <- cv.glmnet(Yhat_int, Y, alpha = 1, intercept = 0)
        model <- glmnet(Yhat_int, Y, alpha = 1, lambda = cv$lambda.min,
                        intercept = 0)
      } else {
        model <- glmnet(Yhat_int, Y, alpha = 1, lambda = lambda, intercept = 0)
      }
      
      # Dsiplay regression coefficients
      w_lasso = data.matrix(coef(model))[-1, ]
      w_lasso = w_lasso[w_lasso > 0]
      ind = as.data.frame(w_lasso)
      Yhat_sub = Yhat_int[, rownames(ind)]
    } else {
      Yhat_sub = Yhat_int
    }
    
    if (incl_lasso) {
      allweights[[i]] = as.data.frame(allgenes_ensemble_cv(Yhat_re = Yhat_sub, Y_re = Y,preprossessed = TRUE, w_lasso = w_lasso,optimization = optimization))
    } else {
      allweights[[i]] = as.data.frame(allgenes_ensemble_cv(Yhat_re = Yhat_sub, Y_re = Y,preprossessed = TRUE, w_lasso = NULL,optimization = optimization))
    }
  }
  
  
  return(list(allweights = allweights, Yhat_re = Yhat_sub, Y_re = Y))
}




############ get phat ###############

get_all_phatyhat = function(rds,scale = "log",transformation = "CPM",r2 = r2, r2_ind,allgene = FALSE,p_true,test_bulk) {
  # function to calculate yhat and get phat list
  phat_all = list()
  yhat_all = list()
  ind_rds = names(rds)
  # if(transformation == "none"){
  #   rds_sub = rds[which(str_detect(ind_rds, "none_none"))]
  #   if (length(rds_sub) == 0) {
  #     rds_sub = rds[which(str_detect(ind_rds, "none"))]
  #   }
  # }else if(transformation == "CPM"){
  #   rds_sub = rds[which(str_detect(ind_rds, "CPM_CPM"))]
  #   if (length(rds_sub) == 0) {
  #     rds_sub = rds[which(str_detect(ind_rds, "CPM"))]
  #   }
  # }else{
  #   rds_sub = rds[which(str_detect(ind_rds, "TPM_TPM"))]
  #   if (length(rds_sub) == 0) {
  #     rds_sub = rds[which(str_detect(ind_rds, "TPM"))]
  #   }
  # }
  # 
  # if(scale == "log"){
  #   rds_sub = rds_sub[which(str_detect(names(rds_sub), "log"))]
  # }else{
  #   rds_sub = rds_sub[which(str_detect(names(rds_sub), "linear"))]
  # }
  
  # test_bulk = rds_sub[[1]]$ensemble$test_bulk
  if(allgene){
    if(transformation != "none"){
      test_bulk = Normalization(test_bulk,transformation)
    }
    if(scale == "log"){
      test_bulk = log2(test_bulk+1)
    }
  }
  if(r2){
    rds <- rds[r2_ind]
  }
  # for (i in 1:length(rds)) {
  #   res <- rds[[i]]$p_hat
  #   # res <- lapply(res, ord_name)
  #   phat <- lapply(res, sum_to_one)
  #   p <- rds[[i]][["p"]]
  #   names(rds[[i]][["p_hat"]]) <- paste0(paste0(p, collapse = "_"),
  #                                       "_", names(phat))
  #   cl <- makeCluster(50)
  #   registerDoParallel(cl)
  #   yhat <- foreach(i = 1:length(res)) %dopar% {
  #     ress <-predict(lm(t(data.matrix(test_bulk)) ~ 0 + res[[i]]))
  #     return(ress)
  #   }
  #   stopCluster(cl)
  #
  #   yhat <- lapply(yhat, t)
  #   names(yhat) <- names(rds[[i]][["p_hat"]])
  #   yhat_all <- append(yhat_all, yhat)
  #   phat_all <- append(phat_all, rds[[i]][["p_hat"]])
  # }
  
  rds_phat <- unlist(lapply(rds, function(x) x$p_hat),recursive = F)
  rds_phat <- lapply(rds_phat, sum_to_one)
  #rds_phat <- lapply(rds_phat, ord_name, df_true =p_true)
  names(rds_phat) <- gsub(".","_",names(rds_phat),fixed = T)
  cl <- makeCluster(20)
  registerDoSNOW(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  new_res <- foreach(i = 1:length(rds_phat)) %dopar% {
    res <- lm(t(data.matrix(test_bulk)) ~ rds_phat[[i]])
    res_sum <- summary(res)
    #r2 <- mean(sapply(res_sum, function(x) x$r.squared))
    rmse <- sqrt(mean(res$residuals^2))
    yhat <-predict(res)
    y_pearson_sample = mean(diag(cor(test_bulk,t(yhat),method = "spearman")))
    y_pearson_gene = mean(diag(cor(t(test_bulk),yhat,method = "spearman")))
    y_pearson_total = cor(c(t(test_bulk)),c(yhat),method = "spearman")
    if(!is.null(p_true)){
      p_rmse = sqrt(mean((as.matrix(p_true) - rds_phat[[i]])^2))
      p_cor = cor(c(as.matrix(p_true)),c(rds_phat[[i]]),method = "spearman")
    }else{
      p_rmse = NULL
      p_cor = NULL
    }
    
    return(list(yhat = yhat,r2 = r2,rmse = rmse,p_rmse = p_rmse,p_cor = p_cor,y_pearson_sample = y_pearson_sample,y_pearson_total = y_pearson_total,y_pearson_gene = y_pearson_gene))
    
  }
  
  stopCluster(cl)
  res =list()
  res$phat_all <- rds_phat
  names(new_res)<-names(rds_phat)
  res$yhat_all <- lapply(new_res,function(x) x[[1]])
  res$r2 <- sapply(new_res,function(x) x[[2]])
  res$rmse <- sapply(new_res,function(x) x[[3]])
  res$p_rmse <- sapply(new_res,function(x) x[[4]])
  res$p_cor <- sapply(new_res,function(x) x[[5]])
  res$y_pearson_sample <- sapply(new_res,function(x) x[[6]])
  res$y_pearson_total <- sapply(new_res,function(x) x[[7]])
  res$y_pearson_gene <- sapply(new_res,function(x) x[[8]])
  res$metric = cbind(res$r2,res$rmse,res$p_rmse,res$p_cor,res$y_pearson_sample,res$y_pearson_gene,res$y_pearson_total)
  colnames(res$metric) = c("r2","rmse","p_rmse","p_cor","y_pearson_sample","y_pearson_gene","y_pearson_total")
  res$metric = as.data.frame(res$metric)
  res$metric =res$metric%>% mutate(ind = rownames(res$metric),)
  res$metric =res$metric %>% mutate(ind = gsub("BayesPrism_GEP","BayesPrism",ind),
                                    ,sce_name = stringr::str_match(ind, '(.*)_(.*)')[,-1][,1],
                                    ind = stringr::str_match(ind, '(.*)_(.*)')[,-1][,2]
  )
  p1 = ggscatter(res$metric, x = "p_rmse", y = "rmse",
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "RMSE of cell type proportion", ylab = "RMSE of bulk expresseion",color = "ind",palette = "safe",title = paste(scale,transformation,"allgene:",allgene))
  p2 = res$metric %>%
    filter(!ind %in% c("TOASTP","LS Fit","DSA","deconf","Q Prog","ssFrobenius","ssKL","linearRegression","logRegression")) %>%
    ggscatter( x = "p_rmse", y = "rmse",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "RMSE of cell type proportion", ylab = "RMSE of bulk expresseion",color = "ind",palette = "safe",title = paste(scale,transformation,"allgene:",allgene))
  
  p3 = ggscatter(res$metric, x = "p_cor", y = "y_pearson_sample",
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "Spearman of cell type proportion", ylab = "Spearman of bulk expresseion (sample)",color = "ind",palette = "safe",title = paste(scale,transformation,"allgene:",allgene))
  p4= res$metric %>%
    filter(!ind %in% c("TOASTP","LS Fit","DSA","deconf","Q Prog","ssFrobenius","ssKL","linearRegression","logRegression")) %>%
    ggscatter( x = "p_cor", y = "y_pearson_sample",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Spearman of cell type proportion", ylab = "Spearman of bulk expresseion (sample)",color = "ind",palette = "safe",title = paste(scale,transformation,"allgene:",allgene))
  
  p5 = ggscatter(res$metric, x = "p_cor", y = "y_pearson_gene",
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "Spearman of cell type proportion", ylab = "Spearman of bulk expresseion (gene)",color = "ind",palette = "safe",title = paste(scale,transformation,"allgene:",allgene))
  p6= res$metric %>%
    filter(!ind %in% c("TOASTP","LS Fit","DSA","deconf","Q Prog","ssFrobenius","ssKL","linearRegression","logRegression")) %>%
    ggscatter( x = "p_cor", y = "y_pearson_gene",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Spearman of cell type proportion", ylab = "Spearman of bulk expresseion (gene)",color = "ind",palette = "safe",title = paste(scale,transformation,"allgene:",allgene))
  
  p7 = ggscatter(res$metric, x = "p_cor", y = "y_pearson_total",
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "Spearman of cell type proportion", ylab = "Spearman of bulk expresseion (total)",color = "ind",palette = "safe",title = paste(scale,transformation,"allgene:",allgene))
  p8= res$metric %>%
    filter(!ind %in% c("TOASTP","LS Fit","DSA","deconf","Q Prog","ssFrobenius","ssKL","linearRegression","logRegression")) %>%
    ggscatter( x = "p_cor", y = "y_pearson_total",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Spearman of cell type proportion", ylab = "Spearman of bulk expresseion (total)",color = "ind",palette = "safe",title = paste(scale,transformation,"allgene:",allgene))
  res$test_bulk <- test_bulk
  return(list(res=res,p1 =p1,p2= p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7,p8=p8))
}

get_all_phatyhat_ct = function(rds,scale = "log",transformation = "CPM",r2 = r2, r2_ind,allgene = FALSE,p_true) {
  # function to calculate yhat and get phat list
  phat_all = list()
  yhat_all = list()
  ind_rds = names(rds)
  if(transformation == "none"){
    rds_sub = rds[which(str_detect(ind_rds, "none_none"))]
    if (length(rds_sub) == 0) {
      rds_sub = rds[which(str_detect(ind_rds, "none"))]
    }
  }else if(transformation == "CPM"){
    rds_sub = rds[which(str_detect(ind_rds, "CPM_CPM"))]
    if (length(rds_sub) == 0) {
      rds_sub = rds[which(str_detect(ind_rds, "CPM"))]
    }
  }else{
    rds_sub = rds[which(str_detect(ind_rds, "TPM_TPM"))]
    if (length(rds_sub) == 0) {
      rds_sub = rds[which(str_detect(ind_rds, "TPM"))]
    }
  }
  
  if(scale == "log"){
    rds_sub = rds_sub[which(str_detect(names(rds_sub), "log"))]
  }else{
    rds_sub = rds_sub[which(str_detect(names(rds_sub), "linear"))]
  }
  
  test_bulk = rds_sub[[1]]$ensemble$test_bulk
  if(allgene){
    if(transformation != "none"){
      test_bulk = Normalization(test_bulk,transformation)
    }
    if(scale == "log"){
      test_bulk = log2(test_bulk+1)
    }
  }
  if(r2){
    rds <- rds[r2_ind]
  }
  
  
  rds_phat <- unlist(lapply(rds, function(x) x$p_hat),recursive = F)
  rds_phat <- lapply(rds_phat, sum_to_one)
  rds_phat <- lapply(rds_phat, ord_name, df_true =p_true)
  names(rds_phat) <- gsub(".","_",names(rds_phat),fixed = T)
  cl <- makeCluster(20)
  registerDoSNOW(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  new_res <- foreach(i = 1:length(rds_phat)) %dopar% {
    yhat =  list()
    r2 = rmse = p_rmse=rmse = p_cor =NULL
    for (j in 1:ncol(p_true)) {
      res <- lm(t(data.matrix(test_bulk)) ~ rds_phat[[i]][,j])
      res_sum <- summary(res)
      r2[j] <- mean(sapply(res_sum, function(x) x$r.squared))
      rmse[j] <- sqrt(mean(res$residuals^2))
      yhat[[j]] <-predict(res)
      p_rmse[j] = sqrt(mean((as.matrix(p_true[,j]) - rds_phat[[i]][,j])^2))
      p_cor[j] = cor(c(as.matrix(p_true[,j])),c(rds_phat[[i]][,j]))
      
    }
    return(list(yhat = yhat,r2 = r2,rmse = rmse,p_rmse = p_rmse,p_cor = p_cor))
  }
  
  stopCluster(cl)
  res =list()
  res$phat_all <- rds_phat
  names(new_res)<-names(rds_phat)
  res$yhat_all <- lapply(new_res,function(x) x[[1]])
  res$r2 <- sapply(new_res,function(x) x[[2]])
  res$rmse <- sapply(new_res,function(x) x[[3]])
  res$p_rmse <- sapply(new_res,function(x) x[[4]])
  res$p_cor <- sapply(new_res,function(x) x[[5]])
  res$test_bulk <- test_bulk
  return(res)
}


get_all_phat <- function(rds) {
  phat_all <- list()
  for (i in 1:length(rds)) {
    res <- rds[[i]]$p_hat
    res <- lapply(res, ord_name)
    phat <- lapply(res, sum_to_one)
    p <- rds[[i]][["p"]]
    names(rds[[i]][["p_hat"]]) <- paste0(p$TScaling, "_", p$data_name,
                                         "_", p$Marker.Method, "_", p$Scale, "_", names(phat))
    phat_all <- append(phat_all, rds[[i]][["p_hat"]])
  }
  rds$phat_all <- phat_all
  
  return(rds)
}


fitted_Y <- function(df, test_bulk) {
  res <- predict(lm(t(data.matrix(test_bulk)) ~ 0 + df))
  return(res)
}

# fitted_r2 <- function(df, test_bulk) {
#   res <- summary(lm(t(data.matrix(test_bulk)) ~ 0 + df))$r.squared
#   return(res)
# }

mysolnp <- function(Y_re,Yhat_re){
  fn <- function(x){
    norm((t(Y_re) - t(x) %*% t(Yhat_re)),type = "2")
  }
  #EQUALITY CONSTRAINT
  eqn <- function(x){
    sum(x)
  }
  constraints <- 1
  x0 <- matrix(1/ncol(Yhat_re), ncol(Yhat_re), 1) #STARTING PARAMETER VECTOR (NOT SURE WHAT STARTING VALUES TO PUT HERE)
  res <- solnp(x0, fun=fn, eqfun=eqn, eqB=constraints,LB = rep(0,ncol(Yhat_re)),UB = rep(1,ncol(Yhat_re)))
  return(res$pars)
}

