#This function is intended to prepare data markers for deconvolution
############ analyze_dset
analyze_dset <- function(dset, method, q, n_markers, gamma, dmeths, verb, normalize,scale) {

  p_hat <- list()
  p_hat_bc <- list()
  timing_bc <- list()
  timing <- list()
  mark_list <- list()

  data_type <- dset$annotation$data_type


  if(data_type == "microarray"){
    data_type <- "microarray-gene"
  }else if(data_type == "singlecell-rna"){
    data_type <-"rna-seq"
  }

  # Pure samples & signature matrix
  annotation <- dset$annotation
  pure_samples <- annotation$pure_samples
  pure <- unlist(pure_samples)
  K <- length(pure_samples)

  sig_matrix <- my_GenerateSCReference(count_sc = dset$data$data_c,meta_sc = dset$data$meta_ref,"deconv_clust")

  # Get Marker Genes

  if(method == "limma"){
    design <- 1*dset$annotation$pure
    v <- voom(dset$data$data_c, design=design, plot=FALSE)
    fit <- lmFit(v, design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, trend=TRUE)
    markers <- marker.fc(fit2, log2.threshold = log2(2))
    mrk <- marker_strategies(markers, marker_strategy, data_c)

  }else if(method %in% c("t","binom","wilcox")){

    out <-  findMarkers(dset$data$data_c, groups=dset$data$meta_ref$deconv_clust,direction="up",test=method,pval.type ="all")


    mrks_old <-  lapply(out, function(x) x@rownames[1:n_markers])

    mrks <-  lapply(1:length(mrks_old), function(i){
      x<-  which(rownames(dset$data$data_c) %in% mrks_old[[i]])
      names(x) <- rownames(dset$data$data_c)[which(rownames(dset$data$data_c) %in% mrks_old[[i]])]
      return(x)
    })
    names(mrks) <- names(mrks_old)
    rm(out,mrks_old)

  }else if(method == "combined"){
    out <- multiMarkerStats(
      t=findMarkers(dset$data$data_c, groups=dset$data$meta_ref$deconv_clust,direction="up",test="t",pval.type ="all"),
      wilcox=findMarkers(dset$data$data_c, groups=dset$data$meta_ref$deconv_clust,direction="up",test="wilcox",pval.type ="all"))

    mrks_old <-  lapply(out, function(x) x@rownames[1:n_markers])

    mrks <-  lapply(1:length(mrks_old), function(i){
      x<-  which(rownames(dset$data$data_c) %in% mrks_old[[i]])
      names(x) <- rownames(dset$data$data_c)[which(rownames(dset$data$data_c) %in% mrks_old[[i]])]
      return(x)
    })
    names(mrks) <- names(mrks_old)
    rm(out,mrks_old)

  }else if(method == "TOAST"){
    mrks <- csDeCompress(Y_raw = dset$data$data_t,K = K,nMarker = n_markers*K)$finalsigs
  }else if(method == "none"){
    mrks <- NULL
  }else if(method == "linseed"){
    mrks <- rownames(linCor(yref = dset$data$data_t,n.types = K)[["sig"]])
  }else if(method == "COSG"){
    mrks <- as.list(my_cosg(object =dset$data$data_c,group_info = dset$data$meta_ref$deconv_clust, n_genes_user = n_markers)$names)
  }else if(method %in% c("ratio","regression","diff","p.value") ){
    prc <- process_markers(t(cbind(dset$data$data_c,dset$data$data_t)), pure_samples, n_markers, data_type, gamma, markers = NULL,
                           marker_method= method)
    n_markers <- prc$n_markers
    mrks <- prc$mrkrs
  }else if(method == "roc"){
    if(any(table(dset$data$meta_ref$deconv_clust)== 1)){
      stop("Only 1 cell detected for 1 or more cell types.")
    }else{
      rownames(dset$data$meta_ref) = dset$data$meta_ref$SamplesName
      mydata <- CreateSeuratObject(counts  = dset$data$data_c, min.cells = 1, min.features  = 1, project = "tmp",meta.data = dset$data$meta_re)
      Idents(mydata) <- mydata@meta.data$deconv_clust
      
      all.markers <- FindAllMarkers(object = mydata ,test.use = "roc",min.cells.group = 1)
      if(any(table(all.markers$cluster)<2)){
        stop("Only 1 or 0 marker detected for 1 or more cell types.")
      }else{
        tmp_type = unique(dset$data$meta_ref$deconv_clust)
        
        mrks = lapply(tmp_type, function(x){
          tmp = all.markers %>% dplyr::filter(cluster == x) %>% slice_head(n = n_markers)
          return(tmp$gene)
        })  
        names(mrks) = tmp_type
      }
    }
  }else if(method == "limma"){
    

  }else if(method == "Hedge"){
    mrks = my_hedge(ref_mtx = dset$data$data_c,meta_ref = dset$data$meta_ref,n_mrk = n_markers)
  }
  

    mrks_old = NULL
  

    # Run deconvolution methods
    p_truth <- annotation$mixture[-pure, ]
    n_choose <- lengths(mrks)

    for (mth in dmeths) {
      dcnv_out <-run_deconv_method(method_name = mth,to_deconv= dset$data$data_t,ref_matrix = dset$data$data_c,sig_matrix = sig_matrix,meta_ref = dset$data$meta_ref,pure_samples = pure_samples,markers= mrks,
                                   data_type = data_type , verb = verb, marker_method = method,batchcorrec = F,scale = scale)

      p_hat[[mth]] <- dcnv_out$estimate
      timing[[mth]] <- dcnv_out$time

    }

    return(list(n_choose = n_choose, p_hat = p_hat, p_truth = p_truth, timing = timing, markers=mrks,markers_old = mrks_old))
  


}

#################### analyze #################


analyze <- function(method, q, n_markers,gamma, dmeths = NULL, verb = TRUE, normalize = TRUE,scale = scale,
                    datasets = NULL) {
  sig <- paste(method, q, gamma)
  updt(paste(sig, "Starting."), init = TRUE)

  start_time <- Sys.time()

  output <- lapply(datasets, function(dset) analyze_dset(dset, method =method, q = q,n_markers =n_markers, gamma=gamma,
                                                         dmeths = dmeths, verb = verb, normalize = normalize,scale = scale ))
  end_time <- Sys.time()


    p_hat <- lapply(output, "[[", "p_hat")
    p_truth <- lapply(output, "[[", "p_truth")
    n_choose <- lapply(output, "[[", "n_choose")
    markers <- lapply(output, "[[", "markers")
    markers_old <- lapply(output, "[[", "markers_old")
    timing <- lapply(output, "[[", "timing")

    return(list(p_hat = p_hat, p_truth = p_truth, n = n_choose, time = timing,markers=markers,markers_old=markers_old,all_time =  end_time - start_time))
  


}

