#' This function is used when parallel computing is not called to generate Rshiny progress bar.
#'
#' @param count_bulk Bulk gene expression.
#' (Required)  Two-dimensional numeric. Must be in gene x sample format. Must implemented \code{as.matrix}
#' (Optional) In original scale.
#' @param ref_list List of of list. Must implement as.list.
#'
#'
#' (Required) The i-th element of the top-level list is a list of \code{ref_matrix}, \code{meta_ref}. Names of the top level list should be vector of
#' \code{data_name}, namely : bulk-reference.
#'
#' \itemize{
#' \item{'ref_matrix'}{ Reference matrix.}
#' \item{'meta_ref'}{ Meta data for reference matrix..}
#' \item{'data_name'}{Data description. Character in format "Bulk data name_reference data name"}
#' }
#'
#' @param customed_markers Self-defined markers.
#' (Optional) List of one-dimensional string Names of the list should match the deconv_clust.
#' @param markers_range Specific for markerpen.
#' (Optional)
#' @param true_frac True cell type proportions for bulk gene expresseeion.
#' (Optional) Two-dimensional numeric. Must be in samples by celltype.
#' @param params Parameters dataframe for ensemble learning, more details could refer to \code{get_params}.
#' @param outpath (Optional) Path to save output.
#' @param data_name Data description.
#' (Optional) Only input when you want default params. Character in format "Bulk data name-reference data name"
#' @param parallel_comp Logical.
#' @param  ncore 	The number of cores to use for parallel execution.
#' @param rm.duplicated Logical. Remove duplicated genes after maker gene selection. Default: FALSE.
#' @param mrkpen Logical. Apply markerpen on marker gene list. Default: FALSE.
#'
#' @import parallel
#' @importFrom progress progress_bar
#' @importFrom doSNOW registerDoSNOW
#' @importFrom Biobase exprs
#' @importFrom e1071 svm
#' @importFrom foreach foreach %dopar%
#' @importFrom Matrix t
#' @importFrom dtangle dtangle
#' @importFrom MuSiC music_prop
#' @importFrom EPIC EPIC
#' @importFrom xbioc pVar
#' @importFrom ICeDT ICeDT
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom hspe hspe
#' @importFrom sva ComBat
#' @importFrom RVenn overlap_pairs Venn
#' @importFrom FARDEEP fardeep 
#' @importFrom ComICS dcq
#' @importFrom sparseMatrixStats rowVars
#' @importFrom  Seurat FindAllMarkers CreateSeuratObject
#' @importFrom scran findMarkers
#' @import glmnet
#' @import reticulate
#' @export
#'
gen_all_res_list_rshiny = function(count_bulk,meta_bulk = NULL,ref_list,customed_markers = NULL,markers_range = NULL,true_frac = NULL,params = NULL,
                            outpath = NULL,parallel_comp = FALSE,ncore,rm.duplicated =FALSE,mrkpen = FALSE,dmeths = NULL){
  
  if(!is.null(outpath)){
    dir.create(outpath,showWarnings = F)
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
    for(i in 1:nrow(params)){
      progress$inc(1/nrow(params), detail = paste("Doing part", i))
      p = params[i,]
      Dataset = get_input_ensemble(count_bulk = count_bulk, ref_matrix = ref_list[[p$data_name]]$ref_matrix, meta_bulk = meta_bulk,
                                   meta_ref = ref_list[[p$data_name]]$meta_ref, true_frac = true_frac,params = p)
      
      
      a <- analyze(p$Marker.Method,q =  p$Quantile,n_markers = p$n_markers, gamma = p$gamma,dmeths = p$dmeths,
                   normalize = p$Normalize, datasets = Dataset,scale = p$Scale,
                   customed_markers = customed_markers,batchcorrec = p$batchcorrec,rm.duplicated = rm.duplicated,mrkpen = mrkpen)
      gc()
      res_all[[i]] = list(a = a, p = p)
      names(res_all)[i] =  paste0(params[i, ], collapse = "_")
      if(!is.null(outpath)){
        saveRDS(list(a = a, p = p), file = paste0(outpath, paste0(params[i, ], collapse = "_"),  ".rds"))
      }
      
    }
    
  
  if(!is.null(outpath)){
    saveRDS(res_all,paste0(outpath,"Res_list.rds"))
  }
  
  
  return(res_all)
  
}

