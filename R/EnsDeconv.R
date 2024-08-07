#' EnsDeconv: A Wrapper Function for Ensemble Deconvolution
#'
#' This function serves as a wrapper for the EnsDeconv algorithm, providing an interface for ensemble deconvolution analysis on gene expression data.
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
#' @param outpath 
#' @param inrshiny 
#'
#' @importFrom matrixcalc frobenius.norm
#' @importFrom quadprog solve.QP
#'
#' @return A list containing two elements:
#'         - \code{EnsDeconv}: The output of the EnsDeconv algorithm.
#'         - \code{allgene_res}: A list of results from each scenario analyzed.
#'
#' @export
#'
#' @examples
#' # Example usage
#' # EnsDeconv(count_bulk, ref_list)
#'



EnsDeconv <- function(count_bulk,ref_list,enableFileSaving = FALSE,exportRef = FALSE,outpath = NULL,parallel_comp = FALSE,ncore =5,
                      true_frac = NULL,params = NULL,inrshiny = FALSE){

if(inrshiny){
  allgene_res = gen_all_res_list_rshiny(count_bulk = as.matrix(count_bulk), ref_list = ref_list,enableFileSaving = enableFileSaving,exportRef = exportRef, outpath =outpath, true_frac =true_frac, ncore =ncore, parallel_comp = parallel_comp, params = params)
}else{
   allgene_res = gen_all_res_list(count_bulk = as.matrix(count_bulk), ref_list = ref_list,enableFileSaving = enableFileSaving,exportRef = exportRef,outpath =outpath,  true_frac =true_frac, ncore =ncore, parallel_comp = parallel_comp, params = params)
}
  # Check available scenarios
  ind = sapply(allgene_res, function(x){
    length(x[["a"]][["p_hat"]][[1]])
  })
  allgene_res = allgene_res[which(ind == 1)]
  
  EnsDeconv_p = CTS_EnsDeconv_wrapper(allgene_res)
  
  return(list(EnsDeconv = EnsDeconv_p,allgene_res = allgene_res))
}
