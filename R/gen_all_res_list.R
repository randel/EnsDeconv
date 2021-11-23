#' This function is used when there are multiple references.
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
#'
#' @export
#' @import parallel
#' @importFrom progress progress_bar
#' @importFrom doSNOW registerDoSNOW
#' @importFrom Biobase exprs
#' @importFrom limma voom lmFit contrasts.fit eBayes
#' @importFrom CellMix MarkerList
#' @importFrom deconvSeq getdge getb0.rnaseq getx1.rnaseq
#' @importFrom e1071 svm
#' @importFrom foreach foreach %dopar%
#' @importFrom Matrix t
#' @importFrom dtangle dtangle
#' @importFrom MuSiC music_prop
#' @importFrom EPIC EPIC
#' @importFrom CellMix ged
#' @importFrom xbioc pVar
#' @importFrom TED run.Ted
#' @importFrom biomaRt useMart getBM
#' @importFrom ICeDT ICeDT
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom hspe hspe
#' @importFrom sva ComBat
#' @importFrom RVenn overlap_pairs Venn
#' @import FARDEEP
#' @import ComICS
#' @import Seurat
#' @import effsize
#' @examples
#'
#' data(testdata)
#'
#' params  = get_params("none","none","log","singlecell-rna","test_test",20,"p.value")
#'
#' gen_all_res_list(count_bulk = testdata$count_bulk,meta_bulk = testdata$meta_bulk,ref_list = testdata$ref_list, true_frac = testdata$true_frac,ncv_input =2,
#' outpath ="D:/ensemble deconvolution/test/" ,ncore = 12,parallel_comp = TRUE, params = params)
#'
#'
gen_all_res_list = function(count_bulk,meta_bulk = NULL,ref_list,customed_markers = NULL,markers_range = NULL,true_frac = NULL,params = NULL,
                            outpath,parallel_comp = FALSE,ncore,os = "win",rm.duplicated =FALSE,mrkpen = FALSE,dmeths = NULL){

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

  if(parallel_comp){

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

    if(os == "win"){
      cl = makeCluster(ncore, outfile="")
    }else{
      cl = makeCluster(ncore, setup_strategy = "sequential")
    }

  registerDoSNOW(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
    res_all = foreach(i = 1:nrow(params),.options.snow = opts, .errorhandling='pass',
                      .packages = c("SCDC","MuSiC","Bisque","TED","TOAST","EPIC", "dtangle","CellMix","nnls","xbioc","Biobase","scran","DeCompress","ICeDT","DeconRNASeq","preprocessCore","glmnet","sva","hspe","RVenn","markerpen","FARDEEP","ComICS","edgeR","BayICE","reticulate","effsize","Seurat","dplyr")) %dopar% {
      p = params[i,]
      logdir <- paste0(outpath,p$data_name,"/Analysis/")
      dir.create(logdir, showWarnings = FALSE, recursive = TRUE)
      capture.output({
        p = params[i,]

        # Prepare data sets
        Dataset = get_input_ensemble(count_bulk = count_bulk, ref_matrix = ref_list[[p$data_name]]$ref_matrix, meta_bulk = meta_bulk,
                                     meta_ref = ref_list[[p$data_name]]$meta_ref, true_frac = true_frac,params = p)
        # ensemble_dat = lapply(Dataset,function(x) x[["ensemble"]])
# 
#         if(is.null(dmeths)){# Get deconvolution methods
#         dmeths_ori <- c("dtangle", "hspe","deconf","ssFrobenius","ssKL","DSA","Q Prog","LS Fit","CIBERSORT","logRegression","linearRegression","EPIC","TOASTP","MuSiC","Bisque","GEDIT", "ICeDT","DeconRNASeq", "BayesPrism")
# 
#         if (p$Scale == "log"){
#           # dmeths = c("hspe","dtangle", "LS Fit","CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog","logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
#           rm_dmeths <- c("Bisque","MuSiC","logRegression")
#           dmeths <- setdiff(dmeths_ori,rm_dmeths)
#         }else if(p$TNormalization == "none"& p$CNormalization == "none" & p$data_type == "singlecell-rna"){
#           dmeths <- dmeths_ori
#         }else if(p$TNormalization == "none"& p$CNormalization == "none" & p$data_type != "singlecell-rna"){
#           # dmeths = c("LS Fit", "CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq","GEDIT")
#           rm_dmeths <- c("Bisque","MuSiC")
#           dmeths <- setdiff(dmeths_ori,rm_dmeths)
#         }else if(p$data_type == "singlecell-rna" & p$Scale == "linear"){
#           # dmeths = c("MuSiC","LS Fit","CIBERSORT", "EPIC","BayesPrism_GEP","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
#           rm_dmeths <- c("Bisque","dtangle","hspe","GEDIT")
#           dmeths <- setdiff(dmeths_ori,rm_dmeths)
#         }else{
#           # dmeths = c("LS Fit","CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
#           rm_dmeths <- c("Bisque","MuSiC","dtangle","hspe","GEDIT")
#           dmeths <- setdiff(dmeths_ori,rm_dmeths)
#         }
# 
#         if(p$Marker.Method %in% c("linseed","TOAST")){
#           rm_dmeths = c("hspe","dtangle", "deconf", "ssFrobenius", "ssKL", "DSA","TOASTP","MuSiC","Bisque")
#           dmeths = setdiff(dmeths,rm_dmeths)
#         }
# 
#         if(p$Marker.Method == "none" & p$data_type == "singlecell-rna"){
#           dmeths = c("Bisque","MuSiC")
#         }}

        # analyze data
        a <- analyze(p$Marker.Method,q =  p$Quantile,n_markers = p$n_markers, gamma = p$gamma,dmeths = p$dmeths,
                         normalize = p$Normalize, datasets = Dataset,scale = p$Scale,
                         customed_markers = customed_markers,markers_range = markers_range,pval.type = p$pval.type,feature.select = feature.select,batchcorrec = p$batchcorrec,rm.duplicated = rm.duplicated,mrkpen = mrkpen)
        if(!is.null(outpath)){
          dir.create(paste0(outpath,p$data_name), showWarnings = FALSE, recursive = TRUE)
          dir.create(paste0(outpath,p$data_name,"/cases/"), showWarnings = FALSE, recursive = TRUE)
          saveRDS(list(a = a, p = p,ensemble = ensemble_dat), file = paste0(outpath,p$data_name,"/cases/", paste0(params[i, ], collapse = "_"),  ".rds"))
          
        }
          base::message(base::sprintf("Remaining %i ", i),
                      "scenarios.")
      }, file = paste0(logdir, "log", i, ".txt"))
      gc()
        res_list= list(a = a, p = p,ensemble = 0)
      
      return(res_list)
    }
    stopCluster(cl)
    for (i in 1:nrow(params)) {
      names(res_all)[i] =  paste0(params[i, ], collapse = "_")
    }

  }else{

    for(i in 1:nrow(params)){

      p = params[i,]
      require(SCDC)
      Dataset = get_input_ensemble(count_bulk = count_bulk, ref_matrix = ref_list[[p$data_name]]$ref_matrix, meta_bulk = meta_bulk,
                                   meta_ref = ref_list[[p$data_name]]$meta_ref, true_frac = true_frac,params = p)
      ensemble_dat = NULL
      if(ncv_input >1){
        ensemble_dat = lapply(Dataset,function(x) x[["ensemble"]])
      }

      # Get deconvolution methods
      # dmeths_ori <- c("dtangle", "hspe","deconf","ssFrobenius","ssKL","DSA","Q Prog","LS Fit","CIBERSORT","logRegression","linearRegression","EPIC","TOASTP","MuSiC","Bisque","GEDIT", "ICeDT","DeconRNASeq", "BayesPrism")
      # 
      # if (p$Scale == "log"){
      #   # dmeths = c("hspe","dtangle", "LS Fit","CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog","logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
      #   rm_dmeths <- c("Bisque","MuSiC","BayesPrism","logRegression")
      #   dmeths <- setdiff(dmeths_ori,rm_dmeths)
      # }else if(p$TNormalization == "none"& p$CNormalization == "none" & p$data_type == "singlecell-rna"){
      #   dmeths <- dmeths_ori
      # }else if(p$TNormalization == "none"& p$CNormalization == "none" & p$data_type != "singlecell-rna"){
      #   # dmeths = c("LS Fit", "CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq","GEDIT")
      #   rm_dmeths <- c("Bisque","MuSiC","BayesPrism")
      #   dmeths <- setdiff(dmeths_ori,rm_dmeths)
      # }else if(p$data_type == "singlecell-rna" & p$Scale == "linear"){
      #   # dmeths = c("MuSiC","LS Fit","CIBERSORT", "EPIC","BayesPrism_GEP","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
      #   rm_dmeths <- c("Bisque","BayesPrism","dtangle","hspe","GEDIT")
      #   dmeths <- setdiff(dmeths_ori,rm_dmeths)
      # }else{
      #   # dmeths = c("LS Fit","CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
      #   rm_dmeths <- c("Bisque","MuSiC","BayesPrism","dtangle","hspe","GEDIT")
      #   dmeths <- setdiff(dmeths_ori,rm_dmeths)
      # }
      # 
      # if(p$Marker.Method %in% c("linseed","TOAST")){
      #   rm_dmeths = c("hspe","dtangle", "deconf", "ssFrobenius", "ssKL", "DSA","TOASTP","MuSiC","Bisque")
      #   dmeths = setdiff(dmeths,rm_dmeths)
      # }
      # 
      # if(p$Marker.Method == "none" & p$data_type == "singlecell-rna"){
      #     dmeths = c("Bisque","MuSiC")
      # }

        a <- analyze(p$Marker.Method,q =  p$Quantile,n_markers = p$n_markers, gamma = p$gamma,dmeths = p$dmeths,
                         normalize = p$Normalize, datasets = Dataset,scale = p$Scale,
                         customed_markers = customed_markers,markers_range = markers_range,pval.type = p$pval.type,feature.select = feature.select,batchcorrec = p$batchcorrec,rm.duplicated = rm.duplicated,mrkpen = mrkpen)
      gc()
      res_all[[i]] = list(a = a, p = p,ensemble = 0)
      names(res_all)[i] =  paste0(params[i, ], collapse = "_")
      if(!is.null(outpath)){
        saveRDS(list(a = a, p = p,ensemble = 0), file = paste0(outpath, paste0(params[i, ], collapse = "_"),  ".rds"))
      }
      
    }

  }
  if(!is.null(outpath)){
    saveRDS(res_all,paste0(outpath,"Res",ncv_input,"_list.rds"))
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

