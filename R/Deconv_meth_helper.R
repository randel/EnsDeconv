############## run_deconv_method ############
run_deconv_method <- function(method_name,sig_matrix, to_deconv,ref_matrix= NULL,meta_ref = NULL, pure_samples,pure, markers,data_type, verb,data_prepossessed =FALSE,marker_method,batchcorrec = FALSE,scale) {
  NA_array <- array(NA, c(ncol(to_deconv), length(markers)))
  out <- tryCatch({
    deconv_method_switch(method_name = method_name, to_deconv = to_deconv,ref_matrix = ref_matrix,meta_ref = meta_ref,  pure_samples = pure_samples, markers = markers,verb = verb,data_prepossessed = data_prepossessed,marker_method = marker_method,batchcorrec = batchcorrec,scale = scale,sig_matrix = sig_matrix)
  }, error = function(cond) {
    message(paste0(">>>Failed ", method_name, "."))
    cat(paste("Caught", cond))
    return(list(out = NA_array, time = NA))
  })

  return(list(estimate = out$phat, time = out$tme))

}


############# deconv_method_switch ###########################

deconv_method_switch <- function(method_name, to_deconv, ref_matrix, meta_ref,pure_samples, markers, verb, sc.eset.list, data_prepossessed, marker_method,batchcorrec,scale,sig_matrix) {


  pure <- unlist(pure_samples)
  K <- length(pure_samples)

  if(marker_method == "customed"){
    markers <- lapply(markers,intersect,rownames(to_deconv))
  }

  if(method_name %in% c("deconf","ssFrobenius","ssKL","DSA")){
    mrk_venn <- Venn(markers)
    dup_mrk <- unlist(overlap_pairs(mrk_venn))
    markers = lapply(markers, function(x) setdiff(x,dup_mrk))
  }

    #if (!(method_name %in% c("dtangle_scRNA", "deconf", "ssFrobenius", "ssKL", "DSA","dtangle_GEP","hspe","DSA_Patrick"))) {
      if(marker_method %in% c("linseed","TOAST")){
        to_deconv <- to_deconv[markers, ]
        sig_matrix <- sig_matrix[markers, ]
        ref_matrix <- ref_matrix[markers, ]
      }else if(marker_method == "none"){
        markers <- rownames(to_deconv)
      }else if(!method_name %in% c("deconf","ssFrobenius","ssKL","DSA")){
        method_name_new <-ifelse(method_name %in% c("dtangle_scRNA","dtangle_GEP","hspe"),"dtangle",method_name)
        markers <- switch(method_name_new, dtangle = lapply(markers, function(x){
          if(any(markers[[1]]%in% rownames(to_deconv))){
            x
          }else{
           return(rownames(to_deconv)[x])
          }
        } ), markers)
        to_deconv <- to_deconv[unlist(markers), ]
        sig_matrix <- sig_matrix[unlist(markers), ]
        ref_matrix <- ref_matrix[unlist(markers), ]
        markers_toast <- markers
        

        if(class(markers[[1]][1]) == "character"){
          markers_toast <- markers
        }else{
          test_mq <- function(t) rownames(to_deconv)[t]
          markers_toast <- lapply(markers_toast,test_mq)
        }
      }

    #}




  # if(method_name == "dtangle_GEP"){
  #   # pure sample indexes
  #   pure_samples_GEP <- as.list(1:ncol(sig_matrix))
  #   names(pure_samples_GEP) <- colnames(sig_matrix) # cell type names
  #   pure_GEP <- unlist(pure_samples_GEP)
  # }

  if(batchcorrec){
    to_deconv <- B_mode(mix = to_deconv,ref = sig_matrix,scale = scale,phat = output)}

  tme <- system.time(output <- switch(method_name,
                                      dtangle_GEP = dtangle::dtangle(Y = t(cbind(sig_matrix, to_deconv)), pure_samples = pure_samples_GEP,
                                                            markers = markers[names(pure_samples_GEP)])$estimates[-pure_GEP, ],
                                      dtangle = dtangle::dtangle(Y = t(cbind(ref_matrix, to_deconv)), pure_samples = pure_samples,
                                                           markers = markers[names(pure_samples)])$estimates[-unlist(pure_samples), ],
                                      hspe = hspe::hspe(Y = t(cbind(ref_matrix, to_deconv)), pure_samples = pure_samples, markers = markers)$estimates,
                                      deconf = t(ged(object = to_deconv, x = markers, method = "deconf", verbose = verb)@fit@H),
                                      ssFrobenius = t(ged(object = to_deconv, x = markers, method = "ssFrobenius", verbose = verb)@fit@H),
                                      ssKL = t(ged(object = to_deconv, x = markers, method = "ssKL", verbose = verb)@fit@H),
                                      DSA = t(ged(object = to_deconv, x = markers, method = "DSA", verbose = verb)@fit@H),
                                      `Q Prog` = t(ged(object = to_deconv, x = sig_matrix, method = "qprog", verbose = verb)@fit@H),
                                      `LS Fit` = t(ged(object = to_deconv, x = sig_matrix, method = "lsfit", verbose = verb)@fit@H),
                                      CIBERSORT = CIBERSORT(sig_matrix, to_deconv)[, 1:K, drop = FALSE],
                                      logRegression =  mylog(to_deconv = to_deconv,sig_matrix = sig_matrix),
                                      linearRegression = mylinear(to_deconv = to_deconv,sig_matrix = sig_matrix),
                                      EPIC = EPIC::EPIC(bulk = to_deconv, reference = list(refProfiles = sig_matrix,
                                                                                     sigGenes = rownames(sig_matrix)))$cellFractions[, 1:K],
                                      TOASTP = t(MDeconv(to_deconv, markers_toast,epsilon = 1e-3, verbose = FALSE)$H),
                                      MuSiC = mymusic_prop(count_bulk = to_deconv, count_sc =  ref_matrix, meta_sc = meta_ref,clusters = 'deconv_clust', markers = NULL, normalize = FALSE,
                                                           samples = 'SubjectName', verbose = F)$Est.prop.weighted,
                                      Bisque =t(my_bisque(to_deconv, count_sc = ref_matrix,meta_sc = meta_ref,cell.types = "deconv_clust",
                                                             subject.names = "SubjectName")$bulk.props),
                                      GEDIT = GEDIT_wrap(ref = sig_matrix,mix = to_deconv),
                                      SCDC = SCDC_prop(bulk.eset = T.eset, sc.eset = C.eset, ct.varname = "deconv_clust", sample = "SubjectName",
                                                       ct.sub = unique(as.character(meta_ref$deconv_clust)), iter.max = 200)$prop.est.mvw,
                                      SCDC_EN = SCDC_ENSEMBLE(bulk.eset = T.eset, sc.eset.list =  sc.eset.list, ct.varname = "deconv_clust",
                                                              sample = "SubjectName", ct.sub =  ct.sub, search.length = 0.01),
                                      deconvSeq = deconvSeq_wrap(bulkdata = to_deconv, singlecelldata = ref_matrix,
                                                                 celltypes.sc = as.character(meta_ref$deconv_clust)),
                                      ICeDT = t(ICeDT::ICeDT(Y = to_deconv, Z = sig_matrix, tumorPurity = NULL, refVar = NULL)$rho[-1,]),
                                      DeconRNASeq = DeconRNASeq_wrap(to_deconv = as.data.frame(to_deconv),sig_matrix = as.data.frame(sig_matrix)),
                                      BayesPrism = run.Ted(ref.dat =t(sig_matrix), X = t(to_deconv), input.type="GEP")$res$final.gibbs.theta,
                                      FARDEEP = FARDEEP::fardeep(sig_matrix, to_deconv, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta,
                                      DCQ = ComICS::dcq(reference_data = sig_matrix, mix_data = to_deconv, marker_set = as.data.frame(rownames(sig_matrix)), alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 10)$average,
                                      BayICE = t(BayICE(ref.set=ref_matrix, mix.set=to_deconv,ref.id=meta_ref$deconv_clust,iter=2000,burnin = 0.6,
                                             thinning = 3)[["w"]]),
                                      DWLS = solveDampenedWLS(sig_matrix,to_deconv),
                                      stop("Method not found.")))

  output <- output/rowSums(output) #
  output[output  < 0] <- 0

  return(list(phat = output, tme = tme))
}


############# util ##############
sze <- function(unit = "Mb", envir = globalenv()) {
  os <- function(x) {
    object.size(get(x))
  }
  os_pretty <- function(x) {
    format(object.size(get(x)), units = unit)
  }
  ord <- order(sapply(ls(envir), os), decreasing = TRUE)
}


updt <- function(msg, init = FALSE) {
  VERBOSE <- TRUE
  TIME <- Sys.time()
  if (VERBOSE) {
    t2 <- Sys.time()
    if (init) {
      TIME <<- t2
      FIRST <<- t2
    }
    MSG <- paste(msg, " (", format(t2 - FIRST, digits = 1), " ellapsed)", sep = "")
    TIME <<- t2
    # cat(blue(paste0("[", t2, "]")), MSG, "\n")
    return(MSG)
  }
  return("")
}

########## CIBERSORT ############

# CIBERSORT R script v1.03 (last updated 07-10-2015) Note: Signature matrix
# construction is not currently available; use java version for full
# functionality.  Author: Aaron M. Newman, Stanford University
# (amnewman@stanford.edu) Requirements: R v3.0 or later. (dependencies below
# might not work properly with earlier versions) install.packages('e1071')
# install.pacakges('parallel') install.packages('preprocessCore') if
# preprocessCore is not available in the repositories you have selected, run the
# following: source('http://bioconductor.org/biocLite.R')
# biocLite('preprocessCore') Windows users using the R GUI may need to Run as
# Administrator to install or update packages.  This script uses 3 parallel
# processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.  Usage: Navigate to directory containing R script
# In R: source('CIBERSORT.R') results <-
# CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN) Options: i) perm
# = No. permutations; set to >=100 to calculate p-values (default = 0) ii) QN =
# Quantile normalization of input mixture (default = TRUE) Input: signature
# matrix and mixture file, formatted as specified at
# http://cibersort.stanford.edu/tutorial.php Output: matrix object containing all
# results and tabular data written to disk 'CIBERSORT-Results.txt' License:
# http://cibersort.stanford.edu/CIBERSORT_License.txt


# dependencies


# Core algorithm
CoreAlg <- function(X, y) {

  # try different values of nu
  svn_itor <- 3

  res <- function(i) {
    if (i == 1) {
      nus <- 0.25
    }
    if (i == 2) {
      nus <- 0.5
    }
    if (i == 3) {
      nus <- 0.75
    }
    model <- svm(X, y, type = "nu-regression", kernel = "linear", nu = nus, scale = F)
    model
  }

  if (Sys.info()["sysname"] == "Windows")
    out <- mclapply(1:svn_itor, res, mc.cores = 1) else out <- mclapply(1:svn_itor, res, mc.cores = 1)

  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)

  # do cibersort
  t <- 1
  while (t <= svn_itor) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    if (any(weights > 0)) {
      w <- weights/sum(weights)
    } else {
      w <- rep(0, length(weights))
    }
    u <- sweep(X, MARGIN = 2, w, "*")
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  # pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  # get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q < 0)] <- 0
  if (any(q > 0)) {
    w <- (q/sum(q))
  } else {
    w <- rep(0, length(q))
  }

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list(w = w, mix_rmse = mix_rmse, mix_r = mix_r)

}

# do permutations
doPerm <- function(perm, X, Y) {
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while (itor <= perm) {
    # print(itor)

    # random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])

    # standardize mixture
    yr <- (yr - mean(yr))/sd(yr)

    # run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    # store correlation
    if (itor == 1) {
      dist <- mix_r
    } else {
      dist <- rbind(dist, mix_r)
    }

    itor <- itor + 1
  }
  newList <- list(dist = dist)
}

# main function: Disable quantile normalization (recommended for RNA-Seq data) https://cibersort.stanford.edu/runcibersort.php
CIBERSORT <- function(sig_matrix, mixture_file, perm = 0, QN = F) {

  # read in data
  X <- sig_matrix
  Y <- mixture_file

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  # order
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]

  P <- perm  #number of permutations

  # anti-log if max < 50 in mixture file if(max(Y) < 50) {Y <- 2^Y}

  # quantile normalization of mixture file
  if (QN == TRUE) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  # intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, ]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY, ]

  # standardize sig matrix
  X <- (X - mean(X))/sd(as.vector(X))

  # empirical null distribution of correlation coefficients
  if (P > 0) {
    nulldist <- sort(doPerm(P, X, Y)$dist)
  }

  # print(nulldist)

  header <- c("Mixture", colnames(X), "P-value", "Correlation", "RMSE")
  # print(header)
  # print("Finished initialisation.")

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  # iterate through mixtures
  while (itor <= mixtures) {

    # print(paste0("Iteration: ", itor))

    y <- Y[, itor]

    # standardize mixture
    y <- (y - mean(y))/sd(y)

    # run SVR core algorithm
    result <- CoreAlg(X, y)

    # get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    # calculate p-value
    if (P > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r))/length(nulldist))
    }

    # print output
    out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
    if (itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }

    itor <- itor + 1

  }

  # save results write.table(rbind(header,output), file='CIBERSORT-Results.txt',
  # sep='\t', row.names=F, col.names=F, quote=F)

  # return matrix object containing all results
  obj <- rbind(header, output)
  obj <- obj[, -1]
  obj <- obj[-1, ]
  obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
  obj
}


################## my linear / my log #########
mylinear <- function(to_deconv, sig_matrix){
  res <- t(lm(to_deconv ~ 0 + sig_matrix)$coefficients)
  colnames(res) <- gsub("sig_matrix","",colnames(res))
  return(res)
}

mylog <- function(to_deconv, sig_matrix){
  res <-t(lm(log(to_deconv) ~ 0 + log(sig_matrix))$coefficients)
  colnames(res) <- gsub("log(sig_matrix)","",colnames(res),fixed = T)
  return(res)
}

############# deconvSeq ################
# Wrapper function for deconvSeq
#
#  bulkdata Bulk gene expression.
# (Required)  Two-dimensional numeric. Must be in gene x sample format. Must implemented \code{as.matrix}
#  singlecelldata  Single cell reference matrix.
#  celltypes.sc character strings for cell type.
#
#
deconvSeq_wrap <- function(bulkdata,singlecelldata,celltypes.sc){
  #To avoid "Design matrix not of full rank" when removing 1 CT

  design.singlecell <- model.matrix(~ -1 + as.factor(celltypes.sc))
  colnames(design.singlecell) <- levels(as.factor(celltypes.sc))
  rownames(design.singlecell) <- colnames(singlecelldata)

  dge.singlecell <- getdge(singlecelldata,design.singlecell, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")
  b0.singlecell <- getb0.rnaseq(dge.singlecell, design.singlecell, ncpm.min =1, nsamp.min = 4)
  dge.tissue <- getdge(bulkdata, NULL, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")

  RESULTS <- getx1.rnaseq(NB0 = "top_fdr",b0.singlecell, dge.tissue)$x1
  return(RESULTS)
}

############# DeconRNASeq ################

# Wrapper function for DeconRNASeq
DeconRNASeq_wrap <- function(to_deconv,sig_matrix){
  x <- DeconRNASeq(datasets = to_deconv, signatures = sig_matrix, use.scale = TRUE, fig = FALSE)$out.all
  rownames(x) <- colnames(to_deconv)
  return(x)
}

################ GEDIT ##################
GEDIT_wrap = function(ref,mix){
  mix_names = colnames(mix)
  write.table(ref, file="./tmp/ref.tsv", quote=FALSE, sep='\t', col.names = NA)

  write.table(mix, file="./tmp/mix.tsv", quote=FALSE, sep='\t', col.names = NA)
  x ='python ./scripts/GEDIT.py -mix ./tmp/mix.tsv -ref ./tmp/ref.tsv'
  res1 = system(x, intern = TRUE)
  if(length(res1) == 3){
    mix = gsub("\\\\", "/",res1[[2]])
    ref = gsub("\\\\", "/",res1[[3]])
  }else{
    mix = gsub("\\\\", "/",res1[[1]])
    ref = gsub("\\\\", "/",res1[[2]])
  }

  res = GEDIT_R(mix = mix, ref = ref)
  rownames(res) =mix_names

  return(res)
}

GEDIT_R = function(mix,ref){
  A_GEDITDecon = function(MixMat, RefMat){
    Predictions = C_GLMRegressionLoop(MixMat, RefMat, alpha, lambda, intercept)
    if (dim(Predictions)[2] == 1){
      Predictions[,2] = Predictions[,1]
    }
    Predictions = Predictions[2:dim(Predictions)[1],]
    Adjusted = DD_Adjust(Predictions)
    Adjusted = round(Adjusted[,],4)
    Adjusted = t(Adjusted)
    return(Adjusted)}

  C_GLMRegressionLoop = function(Samples, Refmat, Alph, Lamb, Inter)
  {
    Predictions = as.data.frame(matrix(0, ncol = dim(Samples)[2], nrow = 1+dim(Refmat)[2]))
    for (i in 1:dim(Samples)[2]){
      OutMat = as.matrix(coef(glmnet(as.matrix(Refmat), as.matrix(Samples[,i]),
                                     lower.limits = 0.0, alpha = Alph, lambda = Lamb, intercept = Inter)))
      LastColumn = dim(OutMat)[2]
      Predictions[,i] = OutMat[,LastColumn]
    }
    names(Predictions) = names(Samples)
    row.names(Predictions) = c("intercept",names(Refmat))
    return(Predictions)
  }

  D_AdjustLoop = function(Predictions){
    Adjusted = Predictions
    for (i in 1:length(Predictions)){
      Adjusted[[i]] = DD_Adjust(Predictions[[i]])}
    return(Adjusted)}

  DD_Adjust = function(Predictions){
    Adjusted = Predictions
    numPreds = dim(Predictions)[2]
    numCTs   = dim(Predictions)[1]
    for (i in 1:numPreds){
      Adjusted[i] = round(Predictions[i]/sum(Predictions[1:numCTs,i]),5)
      if (sum(Predictions[1:numCTs,i]) == 0.0){
        Adjusted[i] = round(1.0/numCTs,5)
      }
    }
    return(Adjusted)
  }

  alpha = 0 #as.numeric(args[5])
  lambda = 0 #as.numeric(args[6])
  intercept = FALSE #as.logical(args[7])
  mixture = read.table(mix, header = TRUE, row.names = 1, sep = "\t")
  ref = read.table(ref, header = TRUE, row.names = 1, sep = "\t")

  predictions = as.matrix(A_GEDITDecon(mixture, ref))

  return(predictions)
}

