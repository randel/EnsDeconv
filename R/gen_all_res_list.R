#' This function is used when there are multiple references.
#'
#' @param count_bulk Bulk gene expression data.
#'  The count_bulk parameter expects a two-dimensional numeric matrix in a gene-by-sample format. 
#'          It must be convertible using \code{as.matrix}. Optionally, the data can be in its original scale.
#'
#' @param ref_list Reference data list.
#'  The ref_list is a list of lists, where each sublist contains \code{ref_matrix} and \code{meta_ref}. 
#'          The top-level list should be named with a vector of \code{data_name}, indicating the bulk-reference pair.
#'          The sublists should contain:
#'          \itemize{
#'          \item{ref_matrix}{A matrix with rows as genes and columns as samples.}
#'          \item{meta_ref}{Metadata for reference data, including "SamplesName" (column names of ref_matrix) 
#'            and "deconv_clust" (deconvolution clusters, e.g., cell types).}
#'          \item{data_name}{A description of the data, formatted as "bulk data name_reference data name".}
#'          }
#'          
#' @param enableFileSaving Enable Saving of Intermediate Output
#'  (Optional) A boolean flag that controls the saving of intermediate outputs as separate files. 
#'          When set to TRUE, intermediate outputs of the analysis will be saved to files. 
#'          If not explicitly set, this parameter defaults to FALSE, meaning that intermediate 
#'          outputs will not be saved by default.
#'          
#' @param exportRef Enable output of reference generated per scenario 
#'  (Optional) A boolean flag that controls the output of reference generated per scenario . 
#'          When set to TRUE, reference generated per scenario will be output. 
#'          If not explicitly set, this parameter defaults to FALSE, meaning that the results will not contain
#'          references per scenario.
#'         
#'
#' @param parallel_comp Use parallel computing.
#'  (Optional) A logical flag indicating whether to perform computations in parallel. 
#'          Defaults to FALSE.
#'
#' @param ncore Number of cores for parallel execution.
#'  (Optional) Sets the number of cores for parallel processing when \code{parallel_comp} is TRUE. 
#'          Default is 5. Only effective if parallel computing is enabled.
#'
#'
#' @param true_frac True cell type proportions.
#'  (Optional) A two-dimensional numeric matrix indicating the true cell type proportions 
#'          in the samples. The matrix should be formatted with samples as rows and cell types as columns.
#'
#' @param params Ensemble learning parameters.
#'  (Optional) A dataframe specifying parameters for ensemble learning. 
#'          For more details, refer to the \code{get_params} function.
#' @param outpath 
#'
#' @import parallel
#' @importFrom progress progress_bar
#' @importFrom doSNOW registerDoSNOW
#' @importFrom Biobase exprs
#' @importFrom foreach foreach %dopar%
#' @importFrom Matrix t
#' @importFrom xbioc pVar
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom sparseMatrixStats rowVars
#' @importFrom  Seurat FindAllMarkers CreateSeuratObject
#' @importFrom scran findMarkers
#' @import glmnet
#' @export
#'
gen_all_res_list = function(count_bulk,ref_list,enableFileSaving,exportRef = FALSE,
                            outpath = NULL,true_frac = NULL,params = NULL,parallel_comp = FALSE,ncore){
  if(enableFileSaving){
    if(!is.null(outpath)){
      dir.create(outpath,showWarnings = F)
    }
  }
  
  
  # before parallel computing
  if(is.null(params)){
    params <- get_params(data_type = "singlecell-rna", data_name = names(ref_list))
  }
  
  count_bulk = filterzerovar(count_bulk)
  
  ref_list = lapply(ref_list, ref_prep,count_bulk = count_bulk)
  
  res_all = list()
  exclude <- c()
  if(parallel_comp){
    unique_methods <- unique(params$dmeths)
    testparams <- do.call(rbind, lapply(unique_methods, function(method) {
      params[params$dmeths == method, ][1, ]
    }))
    method_times <- c()
    for (i in 1:nrow(testparams)) {
      p <- testparams[i, ]
      Dataset <- get_input_ensemble(count_bulk = count_bulk, ref_matrix = ref_list[[p$data_name]]$ref_matrix, meta_bulk = NULL,
                                    meta_ref = ref_list[[p$data_name]]$meta_ref, true_frac = true_frac, params = p)
      time_taken <- system.time({
        a <- try(analyze(p$Marker.Method, q = p$Quantile, n_markers = p$n_markers, gamma = p$gamma, dmeths = p$dmeths,
                         normalize = p$Normalize, datasets = Dataset, scale = p$Scale, exportRef = exportRef))
      })
      
      if (inherits(a, "try-error") || time_taken[3] > p$time_limit) {
        warning(sprintf("Method %s timed out and is ignored", p$dmeths))
        exclude <- c(exclude, p$dmeths)
      } else {
        method_times <- c(method_times, time_taken[3])
      }
    }
    params <- params[!params$dmeths %in% exclude, ]
    avg_time_per_method <- mean(method_times)
    estimated_total_time <- avg_time_per_method * nrow(params) / ncore
    message(sprintf("Estimated total time for parallel execution: %.2f seconds (%.2f minutes)", 
                    estimated_total_time, estimated_total_time / 60))
    pb <- progress_bar$new(
      format = "Current : :current [:bar] :elapsed | percent: :percent",
      total = nrow(params),
      clear = FALSE,
      force = TRUE,
      width = 60)
    
    progress_letter <- rep(1:10, 10)  # token reported in progress bar
    
    # allowing progress bar to be used in foreach -----------------------------
    progress <- function(n){
      pb$tick(tokens = list(letter = progress_letter[n]))
    }
    
    
    opts <- list(progress = progress)
    
    os = get_os()
    if(os == "windows"){
      cl = makeCluster(ncore, outfile="")
      registerDoSNOW(cl)
      clusterCall(cl, function(x) .libPaths(x), .libPaths())
    }else if(os == "osx"){
      cl = makeCluster(ncore, setup_strategy = "sequential")
      registerDoSNOW(cl)
      clusterCall(cl, function(x) .libPaths(x), .libPaths())
    }else if(os == "Linux"){
      cl <- makeCluster(ncore,type = "SOCK")
      registerDoSNOW(cl)
    }
    
    
    res_all = foreach(i = 1:nrow(params),.options.snow = opts, .errorhandling='pass',
                      .packages = c("nnls","xbioc","Biobase","scran","preprocessCore","glmnet","edgeR","Seurat","dplyr","sparseMatrixStats")) %dopar% {
                        #res_all = foreach(i = 1:nrow(params),.options.snow = opts, .errorhandling='pass') %dopar% {
                        
                        p = params[i,]
                        # logdir <- paste0(outpath,p$data_name,"/Analysis/")
                        # dir.create(logdir, showWarnings = FALSE, recursive = TRUE)
                        # capture.output({
                        p = params[i,]
                        
                        # Prepare data sets
                        Dataset = get_input_ensemble(count_bulk = count_bulk, ref_matrix = ref_list[[p$data_name]]$ref_matrix, meta_bulk = NULL,
                                                     meta_ref = ref_list[[p$data_name]]$meta_ref, true_frac = true_frac,params = p)
                        
                        
                        # analyze data
                        a <- analyze(p$Marker.Method,q =  0,n_markers = p$n_markers, gamma = 1,dmeths = p$dmeths,
                                     normalize = p$Normalize, datasets = Dataset,scale = p$Scale,exportRef = exportRef)
                        if(enableFileSaving){
                          if(!is.null(outpath)){
                            dir.create(paste0(outpath,p$data_name), showWarnings = FALSE, recursive = TRUE)
                            dir.create(paste0(outpath,p$data_name,"/cases/"), showWarnings = FALSE, recursive = TRUE)
                            saveRDS(list(a = a, p = p), file = paste0(outpath,p$data_name,"/cases/", paste0(params[i, ], collapse = "_"),  ".rds"))
                            
                          }
                        }
                        
                        base::message(base::sprintf("Remaining %i ", i),
                                      "scenarios.")
                        # }, file = paste0(logdir, "log", i, ".txt"))
                        gc()
                        res_list= list(a = a, p = p,ensemble = 0)
                        
                        return(res_list)
                      }
    stopCluster(cl)
    for (i in 1:nrow(params)) {
      names(res_all)[i] =  paste0(params[i, ], collapse = "_")
    }
  }
else {
    for (i in 1:nrow(params)) {
      p <- params[i,]
      if (p$dmeths %in% exclude) {
        res_all[[i]] <- NULL  # if demths is in exclude list，skip
        next
      }
      
      Dataset <- get_input_ensemble(count_bulk = count_bulk, ref_matrix = ref_list[[p$data_name]]$ref_matrix, meta_bulk = NULL,
                                    meta_ref = ref_list[[p$data_name]]$meta_ref, true_frac = true_frac, params = p)
      
      # recording time
      time_taken <- system.time({
        a <- try(analyze(p$Marker.Method, q = p$Quantile, n_markers = p$n_markers, gamma = p$gamma, dmeths = p$dmeths,
                         normalize = p$Normalize, datasets = Dataset, scale = p$Scale, exportRef = exportRef))
      })
      
      if (inherits(a, "try-error") || time_taken[3] > p$time_limit) {
        warning(sprintf("Method %s time out and is ignored", p$dmeths))
        exclude <- c(exclude, p$dmeths)  # if this method time out in one scenario，blacklist it
        res_all[[i]] <- NULL  
        next  # skip if time out
      }
      
      gc()
      res_all[[i]] <- list(a = a, p = p)
      
      if (enableFileSaving) {
        if (!is.null(outpath)) {
          saveRDS(list(a = a, p = p), file = paste0(outpath, paste0(params[i, ], collapse = "_"), ".rds"))
        }
      }
    }
  }
  
  # filter NULL
  res_all <- res_all[!sapply(res_all, is.null)]
  
  if (enableFileSaving) {
    if (!is.null(outpath)) {
      saveRDS(res_all, paste0(outpath, "Res_list.rds"))
    }
  }
  
  return(res_all)
}
#####################3 sc_prep###############
#This function is used for prepare reference mtx for gen_all_res_list
ref_prep = function(obj,count_bulk){
  
  # rm_zero(obj$ref_matrix, type = class(obj$ref_matrix)[[1]])
  
  obj$ref_matrix=obj$ref_matrix[,pmatch(obj$meta_ref$SamplesName,colnames(obj$ref_matrix))]
  obj$ref_matrix=filterzerovar(obj$ref_matrix)
  gene = intersect(rownames(count_bulk),rownames(obj$ref_matrix))
  obj$ref_matrix = obj$ref_matrix[pmatch(gene,rownames(obj$ref_matrix)),]
  return(obj)
}