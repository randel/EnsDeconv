#' This function allows you generate parameters for ensemble deconvolution.
#'
#' @param Normalization normalization methods for bulk gene expression and reference gene expression
#'
#' (Option) Vector of string. Defaults are "CPM", "TPM" and "none"
#'
#' @param Scale Scaling methods for gene expression data.
#'
#' (Option) Vector of string. Defaults are "log" and "linear"
#'
#' @param data_type Reference data data type
#'
#' (Required) One-dimensional string. Options are
#'
#' \itemize{
#' \item{'singlecell-rna'}{ Single cell RNA seq reference.}
#' \item{'microarray'}{ Microarray reference.} }
#'
#' @param data_name Name for bulk_reference
#'
#' (Options) One-dimensional string. Can be in format : bulk_reference
#'
#' @param n_markers Number of markers
#'
#' (Option) One-dimensional numeric, default is 150
#'
#' How many markers genes to use for deconvolution. Can either be a single integer, vector of integers. All cell types use the same number of markers.
#'
#' @param dmeths Deconvolution methods to used. Default is NULL, EnsDeconv will run all possible deconvolution methods for specific reference data sets.
#' For RNA seq deconvolution, EnsDeconv will run "dtangle", "hspe","CIBERSORT","EPIC","MuSiC","BisqueRNA","GEDIT", "ICeDT","DeconRNASeq","FARDEEP","DCQ".
#' For microarray deconvolution, EnsDeconv will not run "MuSiC" and "BisqueRNA".
#'
#'
#' @param Marker.Method Method used to choose marker genes
#' \itemize{
#' \item{'ratio'}{ selects and ranks markers by the ratio of the mean expression of each gene in each cell type to the mean of that gene in all other cell types.}
#' \item{'regression '}{ selects and ranks markers by estimated regression coefficients in a series of regressions with single covariate that is indicator of
#' each type.}
#' \item{'diff'}{ selects and ranks markers based upon the difference, for each cell type, between the median expression of a gene by each cell type and the
#' median expression of that gene by the second most highly expressed cell type.}
#' \item{'p.value'}{ selects and ranks markers based upon the p-value of a t-test between the median expression of a gene by each cell type and the median
#' expression of that gene by the second most highly expressed cell type.}
#' \item{'t'}{ Perform pairwise Welch t-tests between groups of cells, possibly after blocking on uninteresting factors of variation.}
#' \item{'wilcox'}{ Perform pairwise Wilcoxon rank sum tests between groups of cells, possibly after blocking on uninteresting factors of variation.}
#' \item{'binom'}{ Perform pairwise binomial tests between groups of cells, possibly after blocking on uninteresting factors of variation.}
#' \item{'TOAST'}{ Iteratively searches for cell type-specific features }
#' \item{'Hedge'}{ Iteratively searches for cell type-specific features }
#' \item{'roc'}{ Iteratively searches for cell type-specific features }
#' }
#'
#' @param teqc Logical. Use same normalization between bulk data and reference data. Default: TRUE.
#' @param batchcorrec Logical. Apply batch correction or not. Default: FALSE
#' @return  a data.frame that each row corresponding to specific scenario
#' @export

get_params = function(TNormalization = c("CPM","none","TPM","TMM","QN"),CNormalization = c("CPM","none","TPM","TMM","QN"), Scale = c("log","linear"), data_type, data_name, n_markers = 50, Marker.Method = c("t","wilcox","combined","none","p.value","regression"),dmeths = NULL,teqc = TRUE,batchcorrec =FALSE){
 
   if(is.null(dmeths)){
    if(data_type == "singlecell-rna"){
      dmeths <- c("dtangle", "hspe","CIBERSORT","EPIC","MuSiC","BisqueRNA","GEDIT", "ICeDT","DeconRNASeq","FARDEEP","DCQ")
    }else{
      dmeths <- c("dtangle", "hspe","CIBERSORT","EPIC","GEDIT", "ICeDT","DeconRNASeq","FARDEEP","DCQ")
    }
  }
  norm_params <-  list(TNormalization =TNormalization, CNormalization =CNormalization , data_type = data_type , data_name=data_name, 
                     n_markers = n_markers, Marker.Method=Marker.Method ,Scale=Scale, all_markers = TRUE,dmeths = dmeths)
  params <- expand.grid(norm_params, stringsAsFactors = FALSE)
  if(data_type %in% c("singlecell-rna","rna-seq")){
    # params<-params[!(params$Scale=="linear" & params$Marker.Method=="t"),]
    # params<-params[!(params$Scale=="log" & params$Marker.Method=="none"),]
    params<-params[!(params$CNormalization=="CPM" & params$Marker.Method=="none"),]
    params<-params[!(params$CNormalization=="TPM" & params$Marker.Method=="none"),]
    params<-params[!(params$CNormalization=="QN"),]
  }else if(data_type %in% c("microarray","microarray-gene")){
    params<-params[!(params$CNormalization=="TPM"),]
    params<-params[!(params$CNormalization=="CPM"),]
  }
  if(teqc){
    params = params %>% dplyr::filter(TNormalization == CNormalization) %>%unique()
  }else{
    params = params   %>%unique() %>%
      dplyr::filter(TNormalization != CNormalization & TNormalization !="TMM" &CNormalization !="TMM")
  }

  #params$batchcorrec = batchcorrec

  
  params =  params%>% 
    filter(!(Scale == "linear" & dmeths %in% c("dtangle","hspe"))) %>%
    filter(!(Scale == "log" & dmeths %in% c("ICeDT","BisqueRNA","MuSiC"))) %>% 
    filter(!(TNormalization == "QN" & dmeths == "GEDIT"))
  
  
  return(params)
}
