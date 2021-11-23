#' A parent Function
#'
#' This function allows you generate results for cross validation using different deconvolution methods.
#' @param count_bulk Bulk gene expression.
#'
#' (Required)  Two-dimensional numeric. Must be in gene x sample format. Must implemented \code{as.matrix}
#'
#' (Optional) In original scale.
#'
#' @param meta_bulk Meta data for bulk gene expression.
#'
#' (Optional) Dataframe, much include character variable SamplesName.
#'
#' (Required) If reference data is single cell rna seq.
#'
#' @param ref_matrix Reference matrix.
#'
#' (Required) Two-dimensional numeric.  Could be signature matrix (gene x celltype), non-single cell RNA seq reference matrix (gene x sample),
#' single cell RNA seq reference matrix (gene x cell). Could be  \code{as(ref_matrix, "dgCMatrix")} or \code{as.matrix}
#'
#' @param meta_ref Meta data for reference matrix.
#'
#' (Optional) Dataframe, much include character variable SamplesName, and deconv_clust.
#'
#' (Required)  If reference data is single cell rna seq.
#'
#' deconv_clust : character strings for "cell types".
#'
#' @param customed_markers Self-defined markers.
#'
#' (Optional) List of one-dimensional string Names of the list should match the deconv_clust.
#'
#' @param markers_range Specific for markerpen.
#'
#' (Optional)
#'
#' @param true_frac True cell type proportions for bulk gene expresseeion.
#'
#' (Optional) Two-dimensional numeric. Must be in samples by celltype.
#'
#' @param params Parameters dataframe for ensemble learning, more details could refer to \code{get_params}.
#'
#' @param ncv_input Number of folds.
#'
#' (Optional) Default is 2
#'
#' @param outpath (Required) Path to save output.
#'
#' @param data_name Data description.
#'
#' (Optional) Only input when you want default params. Character in format "Bulk data name_reference data name"
#'
#' @param parallel_comp Logical.
#'
#' @param  ncore 	The number of cores to use for parallel execution.
#'
#' @param os Operation system. Default is "win". For mac user, please specify as "OS"
#'
#'
#'
#' @examples
#'
#' data(testdata)
#'
#' params  = get_params("none","none","log","singlecell-rna","test_test",20,"p.value")
#'
#' gen_all_res(count_bulk =testdata$count_bulk ,meta_bulk = testdata$meta_bulk,ref_matrix = testdata$ref_list$Nowakowski$ref_matrix,
#' meta_ref =  testdata$ref_list$Nowakowski$meta_ref,true_frac = testdata$true_frac, params = params, ncv_input =1,
#' outpath = "D:/ensemble deconvolution/test/",data_name = "test-test",parallel_comp = TRUE, ncore = 2)


gen_all_res <- function(count_bulk,meta_bulk = NULL,ref_matrix,meta_ref = NULL,customed_markers = NULL,markers_range = NULL,true_frac = NULL,params = NULL,
                       ncv_input =2,outpath,data_name = NULL,parallel_comp = FALSE,ncore = NULL,os = "win",myseed =2020,feature.select = FALSE,clustering = TRUE,
                       rm.duplicated = TRUE,mrkpen = FALSE,dmeths = NULL){

  dir.create(outpath,showWarnings = FALSE, )
  # before parallel computing
  if(is.null(params)){
    norm_params <- list(TNormalization = c("none", "CPM"), CNormalization = c("none", "CPM"), data_type = "singlecell-rna", data_name = data_name,
                       Quantile = 0.9,  Marker.Method =c("p.value","regression","ration","diff"),gamma = 1, Scale = c("linear","log"), all_markers = TRUE)
    params <- expand.grid(norm_params, stringsAsFactors = FALSE)
    params <- params %>% filter(TNormalization == CNormalization)
  }

    # ref_matrix = rm_zero(ref_matrix, type = class(ref_matrix)[[1]])

    # count_bulk = rm_zero(count_bulk)
    # filter zero variance gene
    ref_matrix <- filterzerovar(ref_matrix)
    count_bulk <- filterzerovar(count_bulk)

    gene <- intersect(rownames(ref_matrix), rownames(count_bulk))

    count_bulk <- count_bulk[gene,]
    ref_matrix <- ref_matrix[pmatch(gene,rownames(ref_matrix)),]

    if(!is.null(meta_ref)){
      ref_matrix<-ref_matrix[,pmatch(meta_ref$SamplesName,colnames(ref_matrix))]
    }

    res_all <- list()

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
        cl <- parallel::makeCluster(ncore)
      }else{
        cl <- makeCluster(ncore, setup_strategy = "sequential")
      }

      # clusterExport(cl, c("get_input_ensemble", "Normalization","cv_set","analyze","analyze_dset","run_deconv_method","deconv_method_switch",
      # "process_markers","updt","find_markers","combine_Y_refs","get_marker_list","CIBERSORT","CoreAlg","doPerm","get_set_basis","sum_to_one"))
     # registerDoSNOW(cl)
      registerDoParallel(cl)
     clusterCall(cl, function(x) .libPaths(x), .libPaths()) # .options.snow = opts,
      res_all <- foreach(i = 1:nrow(params), .errorhandling='pass', .packages = c("SCDC","MuSiC","BisqueRNA","TED","TOAST","EPIC",
                                                                                                     "dtangle","CellMix","nnls","xbioc","Biobase","scran",
                                                                                                     "BiocParallel","stats","DeCompress","ICeDT","DeconRNASeq",
                                                                                                     "preprocessCore","hspe","glmnet","sva","markerpen","RVenn","FARDEEP","ComICS","edgeR","BayICE")) %dopar% {
        p <- params[i,]
        logdir <- paste0(outpath,p$data_name,"/Analysis/")
        dir.create(logdir, showWarnings = FALSE, recursive = TRUE)
        capture.output({
          p <- params[i,]
          Dataset <- get_input_ensemble(count_bulk = count_bulk, meta_bulk = meta_bulk, ref_matrix = ref_matrix,  meta_ref = meta_ref, true_frac = true_frac,
                                       params = p,ncv_input = ncv_input,myseed = myseed,clustering = clustering)
          ensemble_dat <- lapply(Dataset,function(x) x[["ensemble"]])

          if(is.null(dmeths)){
            # Get deconvolution methods
            dmeths_ori <- c("dtangle", "hspe","deconf","ssFrobenius","ssKL","DSA","Q Prog","LS Fit","CIBERSORT","logRegression","linearRegression","EPIC","TOASTP",
                            "MuSiC","BisqueRNA","GEDIT", "ICeDT","DeconRNASeq", "BayesPrism","DCQ","FARDEEP")

            if (p$Scale == "log"){
              # dmeths = c("hspe","dtangle", "LS Fit","CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog","logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
              rm_dmeths <- c("BisqueRNA","MuSiC","logRegression")
              dmeths <- setdiff(dmeths_ori,rm_dmeths)
            }else if(p$TNormalization == "none"& p$CNormalization == "none" & p$data_type == "singlecell-rna"){
              dmeths <- dmeths_ori
            }else if(p$TNormalization == "none"& p$CNormalization == "none" & p$data_type != "singlecell-rna"){
              # dmeths = c("LS Fit", "CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq","GEDIT")
              rm_dmeths <- c("BisqueRNA","MuSiC")
              dmeths <- setdiff(dmeths_ori,rm_dmeths)
            }else if(p$data_type == "singlecell-rna" & p$Scale == "linear"){
              # dmeths = c("MuSiC","LS Fit","CIBERSORT", "EPIC","BayesPrism_GEP","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
              rm_dmeths <- c("BisqueRNA","dtangle","hspe","GEDIT")
              dmeths <- setdiff(dmeths_ori,rm_dmeths)
            }else{
              # dmeths = c("LS Fit","CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
              rm_dmeths <- c("BisqueRNA","MuSiC","dtangle","hspe","GEDIT")
              dmeths <- setdiff(dmeths_ori,rm_dmeths)
            }

            if(p$Marker.Method %in% c("linseed","TOAST")){
              rm_dmeths = c("hspe","dtangle", "deconf", "ssFrobenius", "ssKL", "DSA","TOASTP","MuSiC","BisqueRNA")
              dmeths = setdiff(dmeths,rm_dmeths)
            }

            if(p$Marker.Method == "none" & p$data_type == "singlecell-rna"){
              dmeths = c("BisqueRNA","MuSiC")
            }
          }


          dir.create(paste0(outpath,p$data_name), showWarnings = FALSE, recursive = TRUE)
          # dir.create(paste0(outpath,p$data_name,"/",ncv_input,"/"), showWarnings = FALSE, recursive = TRUE)


            a <- analyze(p$Marker.Method,q =  p$Quantile,n_markers = p$n_markers, gamma = p$gamma,dmeths = dmeths,
                         normalize = p$Normalize,scale = p$Scale, datasets = Dataset,
                         customed_markers = customed_markers,markers_range = markers_range,pval.type = p$pval.type,feature.select = feature.select,
                         batchcorrec = p$batchcorrec,rm.duplicated = rm.duplicated,mrkpen= mrkpen)
          # saveRDS(list(a = a, p = p,ensemble = ensemble_dat), file = paste0(outpath,p$data_name,"/",ncv_input,"/", paste0(params[i, ], collapse = "_"), ".rds"))
        }, file = paste0(logdir, "log", i, ".txt"))

        gc()
        res_list<- list(a = a, p = p,ensemble = ensemble_dat)
        return(res_list)
      }
      parallel::stopCluster(cl)
      for (i in 1:nrow(params)) {
        names(res_all)[i] <-  paste0(params[i, ], collapse = "_")
      }

    }else{

      for(i in 1:nrow(params)){

        p <- params[i,]
        require(SCDC)
        Dataset <- get_input_ensemble(count_bulk = count_bulk, ref_matrix = ref_matrix, meta_bulk = meta_bulk, meta_ref = meta_ref, true_frac = true_frac,
                                     params = p,ncv_input = ncv_input,clustering = clustering)

        ensemble_dat <- NULL
        if(ncv_input >1){
          ensemble_dat <- lapply(Dataset,function(x) x[["ensemble"]])
        }

        # Get deconvolution methods
        dmeths_ori <- c("dtangle", "hspe","deconf","ssFrobenius","ssKL","DSA","Q Prog","LS Fit","CIBERSORT","logRegression","linearRegression","EPIC","TOASTP",
                        "MuSiC","BisqueRNA","GEDIT", "ICeDT","DeconRNASeq", "BayesPrism")

        if (p$Scale == "log"){
          # dmeths = c("hspe","dtangle", "LS Fit","CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog","logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
          rm_dmeths <- c("BisqueRNA","MuSiC","BayesPrism","logRegression")
          dmeths <- setdiff(dmeths_ori,rm_dmeths)
        }else if(p$TNormalization == "none"& p$CNormalization == "none" & p$data_type == "singlecell-rna"){
          dmeths <- dmeths_ori
        }else if(p$TNormalization == "none"& p$CNormalization == "none" & p$data_type != "singlecell-rna"){
          # dmeths = c("LS Fit", "CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq","GEDIT")
          rm_dmeths <- c("BisqueRNA","MuSiC","BayesPrism")
          dmeths <- setdiff(dmeths_ori,rm_dmeths)
        }else if(p$data_type == "singlecell-rna" & p$Scale == "linear"){
          # dmeths = c("MuSiC","LS Fit","CIBERSORT", "EPIC","BayesPrism_GEP","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
          rm_dmeths <- c("BisqueRNA","BayesPrism","dtangle","hspe","GEDIT")
          dmeths <- setdiff(dmeths_ori,rm_dmeths)
        }else{
          # dmeths = c("LS Fit","CIBERSORT", "EPIC","deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", "logRegression","linearRegression","TOASTP","ICeDT","DeconRNASeq")
          rm_dmeths <- c("BisqueRNA","MuSiC","BayesPrism","dtangle","hspe","GEDIT")
          dmeths <- setdiff(dmeths_ori,rm_dmeths)
        }

        if(p$Marker.Method %in% c("linseed","TOAST")){
          rm_dmeths = c("hspe","dtangle", "deconf", "ssFrobenius", "ssKL", "DSA","TOASTP","MuSiC","BisqueRNA")
          dmeths = setdiff(dmeths,rm_dmeths)
        }

        if(p$Marker.Method == "none" & p$data_type == "singlecell-rna"){
            dmeths = c("BisqueRNA","MuSiC")
        }
        dir.create(paste0(outpath,p$data_name), showWarnings = FALSE, recursive = TRUE)
        dir.create(paste0(outpath,p$data_name,"/",ncv_input,"/"), showWarnings = FALSE, recursive = TRUE)


        a <- analyze(p$Marker.Method,q =  p$Quantile,n_markers = p$n_markers, gamma = p$gamma,dmeths = dmeths,
                     normalize = p$Normalize, datasets = Dataset,scale = p$Scale,
                     customed_markers = customed_markers,markers_range = markers_range,pval.type = p$pval.type,feature.select = feature.select,
                     batchcorrec = p$batchcorrec,rm.duplicated = rm.duplicated,mrkpen = mrkpen)

        gc()
        res_all[[i]] <- list(a = a, p = p,ensemble = ensemble_dat)
        names(res_all)[i] <-  paste0(params[i, ], collapse = "_")
      }

    }
########################################################################

  saveRDS(res_all,paste0(outpath,data_name,"_",ncv_input,"_res_list.rds"))

  return(res_all)
}

############################################################################


