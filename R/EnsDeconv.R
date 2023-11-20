#' EnsDeconv: A Wrapper Function for Ensemble Deconvolution
#'
#' This function serves as a wrapper for the EnsDeconv algorithm, providing an interface for ensemble deconvolution analysis on gene expression data.
#'
#' @param count_bulk Bulk gene expression data.
#' @details The count_bulk parameter expects a two-dimensional numeric matrix in a gene-by-sample format. 
#'          It must be convertible using `as.matrix`. Optionally, the data can be in its original scale.
#'
#' @param ref_list Reference data list.
#' @details The ref_list is a list of lists, where each sublist contains `ref_matrix` and `meta_ref`. 
#'          The top-level list should be named with a vector of `data_name`, indicating the bulk-reference pair.
#'          The sublists should contain:
#'          - `ref_matrix`: A matrix with rows as genes and columns as samples.
#'          - `meta_ref`: Metadata for reference data, including "SamplesName" (column names of ref_matrix) 
#'            and "deconv_clust" (deconvolution clusters, e.g., cell types).
#'          - `data_name`: A description of the data, formatted as "bulk data name_reference data name".
#'          
#' @param enableFileSaving Enable Saving of Intermediate Output
#' @details (Optional) A boolean flag that controls the saving of intermediate outputs as separate files. 
#'          When set to TRUE, intermediate outputs of the analysis will be saved to files. 
#'          If not explicitly set, this parameter defaults to FALSE, meaning that intermediate 
#'          outputs will not be saved by default.
#'
#' @param outputPath Destination for Saved Output Files
#' @details (Optional) Specifies the ile path where output files should be saved, applicable 
#'          only if "enableFileSaving" is set to TRUE. Providing this path directs the function 
#'          to save all intermediate output files to the specified location. 
#'          If "enableFileSaving" is FALSE or not set, the value of "outputPath" is ignored.
#'          This parameter should be a valid file system path.
#'          
#' @param parallel_comp Use parallel computing.
#' @details (Optional) A logical flag indicating whether to perform computations in parallel. 
#'          Defaults to FALSE.
#'
#' @param ncore Number of cores for parallel execution.
#' @details (Optional) Sets the number of cores for parallel processing when "parallel_comp" is TRUE. 
#'          Default is 5. Only effective if parallel computing is enabled.
#'
#'
#' @param true_frac True cell type proportions.
#' @details (Optional) A two-dimensional numeric matrix indicating the true cell type proportions 
#'          in the samples. The matrix should be formatted with samples as rows and cell types as columns.
#'
#' @param params Ensemble learning parameters.
#' @details (Optional) A dataframe specifying parameters for ensemble learning. 
#'          For more details, refer to the `get_params` function.
#'
#'
#' @importFrom matrixcalc frobenius.norm
#' @importFrom quadprog solve.QP
#'
#' @return A list containing two elements:
#'         - `EnsDeconv`: The output of the EnsDeconv algorithm.
#'         - `allgene_res`: A list of results from each scenario analyzed.
#'
#' @export
#'
#' @examples
#' # Example usage
#' # EnsDeconv(count_bulk, ref_list)
#'



EnsDeconv <- function(count_bulk,ref_list,enableFileSaving = FALSE,outpath = NULL,parallel_comp = FALSE,ncore =5,
                      true_frac = NULL,params = NULL,inrshiny = FALSE){

if(inrshiny){
  allgene_res = gen_all_res_list_rshiny(count_bulk = as.matrix(count_bulk), ref_list = ref_list,enableFileSaving = enableFileSaving, outpath =outpath, true_frac =true_frac, ncore =ncore, parallel_comp = parallel_comp, params = params)
}else{
   allgene_res = gen_all_res_list(count_bulk = as.matrix(count_bulk), ref_list = ref_list,enableFileSaving = enableFileSaving,outpath =outpath,  true_frac =true_frac, ncore =ncore, parallel_comp = parallel_comp, params = params)
}
  # Check available scenarios
  ind = sapply(allgene_res, function(x){
    length(x[["a"]][["p_hat"]][[1]])
  })
  allgene_res = allgene_res[which(ind == 1)]
  
  EnsDeconv_p = CTS_EnsDeconv_wrapper(allgene_res)
  
  return(list(EnsDeconv = EnsDeconv_p,allgene_res = allgene_res))
}
