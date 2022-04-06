#' This function is a wrapper for EnsDeconv
#'
#' @param count_bulk Bulk gene expression.
#' 
#' (Required)  Two-dimensional numeric. Must be in gene x sample format. Must implemented \code{as.matrix}.
#' (Optional) In original scale.
#' @param ref_list List of of list. Must implement as.list.
#'
#' (Required) The i-th element of the top-level list is a list of \code{ref_matrix}, \code{meta_ref}. Names of the top level list should be vector of
#' \code{data_name}, namely : bulk-reference.
#'
#' \itemize{
#' \item{'ref_matrix'}{ Reference matrix, rows are genes, columns are samples.}
#' \item{'meta_ref'}{ Meta data for reference data. Must include two variables: "SamplesName" (column names of ref_matrix), "deconv_clust" (deconvolution cluster, eg. cell types).}
#' \item{'data_name'}{Data description. Character in format "bulk data name_reference data name"}
#' }
#'
#' @param customed_markers (Optional) Self-defined markers.
#' 
#' @param true_frac (Optional) True cell type proportions for bulk gene expression
#' Two-dimensional numeric. Must be in samples by cell type.
#' @param params Parameters Dataframe for ensemble learning, more details could refer to \code{get_params}.
#' @param outpath (Optional) Path to save output.
#' @param parallel_comp Logical. Use parallel computing or not. Default: FALSE.
#' @param  ncore 	The number of cores to use for parallel execution.
#' @param rm.duplicated Logical. Remove duplicated genes after maker gene selection. Default: FALSE.
#' @param mrkpen Logical. Apply markerpen on marker gene list. Default: FALSE.
#' @param markers_range (Optional) Specific for markerpen.
#' @importFrom matrixcalc frobenius.norm
#' @return A list containing the output of the EnsDeconv algorithm (EnsDeconv) and a list output from each scenario (allgene_res). 
#' @export

EnsDeconv <- function(count_bulk,ref_list,customed_markers = NULL,true_frac = NULL,params = NULL,
                        outpath = NULL,parallel_comp = FALSE,ncore,rm.duplicated =FALSE,mrkpen = FALSE,markers_range = NULL,dmeths = NULL,inrshiny = FALSE){

if(inrshiny){
  allgene_res = gen_all_res_list_rshiny(count_bulk = as.matrix(count_bulk), meta_bulk = NULL, ref_list = ref_list, true_frac =true_frac, outpath =outpath, ncore =ncore, parallel_comp = parallel_comp, params = params,dmeths = dmeths)
}else{
   allgene_res = gen_all_res_list(count_bulk = as.matrix(count_bulk), meta_bulk = NULL, ref_list = ref_list, true_frac =true_frac, outpath =outpath, ncore =ncore, parallel_comp = parallel_comp, params = params,dmeths = dmeths)
}
  # Check available scenarios
  ind = sapply(allgene_res, function(x){
    length(x[["a"]][["p_hat"]][[1]])
  })
  allgene_res = allgene_res[which(ind == 1)]
  
  EnsDeconv_p = CTS_EnsDeconv_wrapper(allgene_res)
  
  return(list(EnsDeconv = EnsDeconv_p,allgene_res = allgene_res))
}
