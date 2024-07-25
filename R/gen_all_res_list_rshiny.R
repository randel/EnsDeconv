#' This function is used when parallel computing is not called to generate Rshiny progress bar.
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
#' @param outputPath Destination for Saved Output Files
#'  (Optional) Specifies the ile path where output files should be saved, applicable 
#'          only if \code{enableFileSaving} is set to TRUE. Providing this path directs the function 
#'          to save all intermediate output files to the specified location. 
#'          If \code{enableFileSaving} is FALSE or not set, the value of "outputPath" is ignored.
#'          This parameter should be a valid file system path.
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
#'
#' @import parallel
#' @importFrom progress progress_bar
#' @importFrom doSNOW registerDoSNOW
#' @importFrom Biobase exprs
#' @importFrom e1071 svm
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
gen_all_res_list_rshiny = function(count_bulk,meta_bulk = NULL,ref_list,enableFileSaving = FALSE,exportRef = FALSE,
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
  progress <- shiny::Progress$new()
  on.exit(progress$close())

  progress$set(message = "Running scenarios", value = 0)
  exclude <- c()
    for(i in 1:nrow(params)){
      progress$inc(1/nrow(params), detail = paste("Doing part", i))
      p = params[i,]
      if (p$dmeths %in% exclude) {
        res_all[[i]] <- NULL  # if demths is in exclude list，skip
        warning(sprintf("part %s is ignored due to method time out", i))
        next
      }
      Dataset = get_input_ensemble(count_bulk = count_bulk, ref_matrix = ref_list[[p$data_name]]$ref_matrix, meta_bulk = meta_bulk,
                                   meta_ref = ref_list[[p$data_name]]$meta_ref, true_frac = true_frac,params = p)


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
      res_all[[i]] = list(a = a, p = p)
      names(res_all)[i] =  paste0(params[i, ], collapse = "_")
      if(enableFileSaving){
        if(!is.null(outpath)){
          saveRDS(list(a = a, p = p), file = paste0(outpath, paste0(params[i, ], collapse = "_"),  ".rds"))
        }
      }
      

    }

  # filter NULL
  res_all <- res_all[!sapply(res_all, is.null)]
  
  if(enableFileSaving){
    if(!is.null(outpath)){
      saveRDS(res_all,paste0(outpath,"Res_list.rds"))
    }
  }
  


  return(res_all)

}

