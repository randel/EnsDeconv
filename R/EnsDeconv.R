#' This function is a wrapper for EnsDeconv
#'
#' @param count_bulk Bulk gene expression.
#' (Required)  Two-dimensional numeric. Must be in gene x sample format. Must implemented \code{as.matrix}
#' (Optional) In original scale.
#' @param meta_bulk Meta data for bulk gene expression.
#' (Optional) Dataframe, much include character variable SamplesName.
#' (Required) If reference data is single cell rna seq.
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
#' @param os Operation system. Default is "win". For mac user, please specify as "OS"
#' @param rm.duplicated Logical. Remove duplicated genes after maker gene selection. Default: FALSE.
#' @param mrkpen Logical. Apply markerpen on marker gene list. Default: FALSE.
#' 
#' @importFrom matrixcalc frobenius.norm
#' @return A list containing the output of the EnsDeconv algorithm and output from each scenarios
#' @export

EnsDeconv <- function(count_bulk,meta_bulk = NULL,ref_list,customed_markers = NULL,markers_range = NULL,true_frac = NULL,params = NULL,
                        outpath = NULL,parallel_comp = FALSE,ncore,os = "win",rm.duplicated =FALSE,mrkpen = FALSE,dmeths = NULL,trueMet){


   allgene_res = gen_all_res_list(count_bulk = as.matrix(count_bulk), meta_bulk = meta_bulk, ref_list = ref_list, true_frac =true_frac, outpath =outpath, ncore =ncore, parallel_comp = parallel_comp, params = params,dmeths = dmeths)
  
  # Check available scenarios
  ind = sapply(allgene_res, function(x){
    length(x[["a"]][["p_hat"]][[1]])
  })
  allgene_res = allgene_res[which(ind == 1)]
  
  EnsDeconv_p = adaptive_L1(allgene_res)
  
  return(list(EnsDeconv = EnsDeconv_p,allgene_res = allgene_res))
}
