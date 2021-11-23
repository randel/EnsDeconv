########### process_markers ################
# Determines number of markers n_markers, marker list mrkrs, and gamma}.
# Y Expression matrix.
#
# (Required) Two-dimensional numeric. Must implement as.matrix.
#
# Each row contains expression measurements for a particular sample. Each columm contains the measurements of the same gene over all individuals. Can either contain just the mixture samples to be deconvolved or both the mixture samples and the reference samples. See pure_samples} and references} for more details.
# pure_samples The pure sample indicies.
#
# (Optional) List of one-dimensional integer. Must implement as.list}.
#
# The i-th element of the top-level list is a vector of indicies (rows of Y} or references}) that are pure samples of type i. If references} is not specified then this argument identifies which rows of Y} correspond to pure reference samples of which cell-types. If references} is specified then this makes same idenficiation but for the references} matrix instead.
# data_type Type of expression measurements.
#
# (Optional) One-dimensional string.
#
# An optional string indicating the type of the expression measurements. This is used to set gamma to a pre-determined value based upon the data type. Valid values are for probe-level microarray as ``microarray-probe'', gene-level microarray as ``microarray-gene'' or rna-seq as ``rna-seq''. Alternatively can set gamma} directly.
# n_markers Number of marker genes.
#
#(Optional) One-dimensional numeric.
#
#How many markers genes to use for deconvolution. Can either be a single integer, vector of integers (one for each cell type), or single or vector of percentages (numeric in 0 to 1). If a single integer then all cell types use that number of markers. If a vector then the i-th element determines how many marker genes are used for the i-th cell type. If single percentage (in 0 to 1) then that percentage of markers are used for all types. If vector of percentages then that percentage used for each type, respectively. If not specified then top 10\% of genes are used.
#gamma Expression adjustment term.
#
#(Optional) One-dimensional positive numeric.
#
#If provided as a single positive number then that value will be used for gamma} and over-ride the value of gamma chosen by the data_type} argument. If neither gamma} nor data_type} are specified then gamma} will be set to one.
#markers Marker gene indices.
#
#(Optional) List of one-dimensional integer.
#
#Top-level list should be same length as pure_samples}, i.e. one element for each cell type. Each element of the top-level list is a vector of indicies (columns of Y}) that will be considered markers of that particular type. If not supplied then dtangle} finds markers internally using find_markers}. Alternatively, one can supply the output of find_markers} to the markers argument.
# marker_method Method used to rank marker genes.
#
process_markers <- function(Y, pure_samples, n_markers, data_type, gamma, markers, marker_method) {

    K <- length(pure_samples)

    if (is.null(gamma))
        gamma <- get_gamma(data_type)

    if (is.null(markers)) {
        markers <- find_markers(Y = Y, pure_samples = pure_samples, data_type = data_type,
            gamma = gamma, marker_method = marker_method)
        if (is.null(n_markers)) {
            n_markers <- sapply(floor(0.1 * lengths(markers$L)), min, ncol(Y)/K)
        }
    }
    markers <- get_marker_list(markers)

    if (is.null(n_markers)) {
        n_markers <- lengths(markers)
    } else {
        if (length(n_markers) == 1)
            n_markers <- rep(n_markers, K)

        wq_markers <- which(n_markers < 1)

        n_markers[wq_markers] <- floor(n_markers[wq_markers] * lengths(markers)[wq_markers])

    }

    n_markers <- sapply(n_markers, max, 1)

    mrkrs <- lapply(1:K, function(i) {
        markers[[i]][1:n_markers[i]]
    })
    names(mrkrs) <- names(pure_samples)

    return(list(n_markers = n_markers, mrkrs = mrkrs, gamma = gamma))
}

################# find markers ######################
#Find marker genes for each cell type.
#return List with four elements. ``L'' is respective ranked markers for each cell type and ``V'' is the corresponding values of the ranking method (higher are better) used to determine markers and sort them, ``M'' is the matrix used to create the other two arguments after sorting and subsetting, and ``sM'' is a sorted version of M.
#Deconvolve cell type mixing proportions from gene expression data.

find_markers <- function(Y, references = NULL, pure_samples = NULL, data_type = NULL, gamma = NULL, marker_method = marker_method) {

    cmbd <- combine_Y_refs(Y, references, pure_samples)
    Y <- cmbd$Y
    pure_samples <- cmbd$pure_samples

    if (any(lengths(pure_samples) == 1) & marker_method == "p.value") {
        message("Can't use p.value method.")
        marker_method <- "diff"
    }
    if (is.null(gamma))
        gamma <- get_gamma(data_type)
    K <- length(pure_samples)
    N <- dim(Y)[2]
    pure <- unlist(pure_samples)
    C <- array(0, c(K, N))
    colnames(C) <- colnames(Y)
    if (marker_method == "ratio") {
        avg_exp_fn <- function(x) {
            colMeans(2^(Y[x, , drop = FALSE]))/gamma
        }
        eta_hats <- t(sapply(pure_samples, avg_exp_fn))
        C <- t(sapply(1:K, function(i) {
            eta_hats[i, ]/apply(eta_hats[-i, , drop = FALSE], 2, sum)
        }))
    } else if (marker_method == "regression") {
        for (i in 1:K) {
            X <- as.numeric(pure %in% pure_samples[[i]])
            Yp <- as.matrix(Y[pure, ])
            m <- lm(Yp ~ 1 + X)
            cfdf <- data.frame(t(coef(m)))
            C[i, ] <- cfdf$X
        }
    } else if (marker_method == "diff") {
        for (i in 1:K) {
            C[i, ] <- apply(Y[pure_samples[[i]], , drop = FALSE], 2, median)
        }
        less_second <- function(x) {
            x - sort(x, decreasing = TRUE)[2]
        }
        C <- apply(C, 2, less_second)
    } else if (marker_method == "p.value") {
        for (i in 1:K) {
            C[i, ] <- apply(Y[pure_samples[[i]], , drop = FALSE], 2, mean)
        }
        calc_pvals <- function(i) {
            x <- C[, i]
            top <- which(x == max(x))[1]
            second <- order(x, decreasing = TRUE)[2]
            pvs <- rep(NA, length(x))
            for (j in 1:length(x)) {
                pvs[j] <- tryCatch({
                    x1 <- Y[pure_samples[[j]], i]
                    x2 <- Y[pure_samples[[second]], i]
                    n1 <- length(x1)
                    n2 <- length(x2)
                    sd1 <- sd(x1)
                    sd2 <- sd(x2)
                    sp <- sqrt(((n1 - 1) * sd1 + (n2 - 1) * sd2)/(n1 + n2 - 2))
                    t.value <- (mean(x1) - mean(x2))/(sp * sqrt((1/n1) + (1/n2)))
                    tmp <- pt(abs(t.value), df = n1 + n2 - 2)
                    tmp
                })
            }
            pvs[-top] <- 0
            return(pvs)
        }
        C <- sapply(1:ncol(C), calc_pvals)
    } else {
        stop("Marker method not found.")
    }
    pick_top <- function(x) {
        m <- which(x == max(x, na.rm = TRUE))
        if (length(m) != 1)
            return(c(NA, NaN))
        return(c(m, x[m]))
    }
    M <- apply(C, 2, pick_top)
    M <- data.frame(t(M))
    colnames(M) <- c("top", "value")
    M$rn <- 1:N
    rownames(M) <- colnames(Y)
    M$Cell.Type <- names(pure_samples)[M$top]
    if (marker_method == "p.value") {
        diffmm <- find_markers(Y = Y, pure_samples = pure_samples, data_type = data_type,
                               gamma = gamma, marker_method = "diff")$M
        M$diff <- diffmm$value
        iM <- M[complete.cases(M), ]
        sM <- iM[order(iM$top, -iM$value, -iM$diff), ]
    } else {
        iM <- M[complete.cases(M), ]
        sM <- iM[order(iM$top, -iM$value), ]
    }
    L <- lapply(1:K, function(i) {
        vals <- sM[sM$top == i, "rn"]
        names(vals) <- rownames(sM[sM$top == i, ])
        return(vals)
    })
    V <- lapply(1:K, function(i) {
        vals <- sM[sM$top == i, "value"]
        names(vals) <- rownames(sM[sM$top == i, ])
        return(vals)
    })
    names(L) <- names(pure_samples)
    names(V) <- names(pure_samples)
    return(list(L = L, V = V, M = M, sM = sM))
}

################# get marker list ###############
get_marker_list <- function(value) {
    typ <- typeof(value[[1]])
    if (typ == "list")
        return(value$L)
    return(value)
}

################ get gamma ##############
get_gamma <- function(data_type) {

    if (is.null(data_type))
        return(1)

    if (data_type == "microarray-probe") {
        gamma <- dtangle:::gma$ma_probe
    } else if (data_type == "microarray-gene") {
        gamma <- dtangle:::gma$ma_gene
    } else if (data_type == "rna-seq") {
        gamma <- dtangle:::gma$rna_seq
    }
    return(gamma)
}

################ combine_Y_refs############
combine_Y_refs <- function(Y, references, pure_samples) {

    if (is.null(pure_samples)) {
        pure_samples <- lapply(1:nrow(references), identity)
        names(pure_samples) <- rownames(references)
    }

    if (is.null(colnames(Y)) & !is.null(colnames(references))) {
        colnames(Y) <- colnames(references)
    }

    if (!is.null(references) & is.null(colnames(references)))
        colnames(references) <- colnames(Y)

    if (!is.null(references))
        Y <- as.matrix(rbind(as.matrix(references), as.matrix(Y)))

    if (is.null(colnames(Y)))
        colnames(Y) <- 1:ncol(Y)

    return(list(Y = Y, pure_samples = pure_samples))
}


################ Get markers through limma ########################
marker.fc <- function(fit2, log2.threshold = 1, output_name = "markers"){

    topTable_RESULTS = topTable(fit2, coef = 1:ncol(cont.matrix), number = Inf, adjust.method = "BH", p.value = 0.05, lfc = log2.threshold)
    AveExpr_pval <- topTable_RESULTS[,(ncol(topTable_RESULTS)-3):ncol(topTable_RESULTS)]
    topTable_RESULTS <- topTable_RESULTS[,1:(ncol(topTable_RESULTS)-4)]

    if(length(grep("ERCC-",topTable_RESULTS$gene)) > 0){ topTable_RESULTS <- topTable_RESULTS[-grep("ERCC-",topTable_RESULTS$gene),] }

    markers <- apply(topTable_RESULTS,1,function(x){
        temp = sort(x)
        ((temp[ncol(topTable_RESULTS)] - temp[ncol(topTable_RESULTS)-1]) >= log2.threshold) | (abs(temp[1] - temp[2]) >= log2.threshold)

    })

    topTable_RESULTS = topTable_RESULTS[markers,]

    markers <- cbind.data.frame(rownames(topTable_RESULTS),
                                t(apply(topTable_RESULTS, 1, function(x){
                                    temp = max(x)
                                    if(temp < log2.threshold){
                                        temp = c(min(x),colnames(topTable_RESULTS)[which.min(x)])
                                    } else {
                                        temp = c(max(x),colnames(topTable_RESULTS)[which.max(x)])
                                    }
                                    temp
                                })))

    colnames(markers) <- c("gene","log2FC","CT")
    markers$log2FC = as.numeric(as.character(markers$log2FC))
    markers <- markers %>% arrange(CT,desc(log2FC))

    markers$AveExpr <- AveExpr_pval$AveExpr[match(markers$gene,rownames(AveExpr_pval))]
    markers$gene <- as.character(markers$gene)
    markers$CT <- as.character(markers$CT)

    #write.table(markers, file = output_name, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

    return(markers)

}



############# get_mrkpen_bulk ##############
# Only provide version for Astrocytes, oligodendrocytes, microglia, endothelial, neurons
get_mrkpen_bulk = function(bulk,mrk_old,markers_range,scale){
    #bulk row: obs, col: gene
    bulk = t(bulk)
    if(scale == "log"){
        bulk = 2^(bulk)-1
    }
    ast_re = refine_markers(bulk, markers_range$astrocytes, names(mrk_old$Astro),
                            lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
    # Remove selected markers from the expression matrix
    mat_rest = bulk[, setdiff(colnames(bulk), ast_re$markers)]

    # Markers for oligodendrocytes
    oli_re = refine_markers(mat_rest, markers_range$oligodendrocytes, names(mrk_old$Oligo),
                            lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
    mat_rest = mat_rest[, setdiff(colnames(mat_rest), oli_re$markers)]

    # Markers for microglia
    mic_re = refine_markers(mat_rest, markers_range$microglia, names(mrk_old$Micro),
                            lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
    mat_rest = mat_rest[, setdiff(colnames(mat_rest), mic_re$markers)]

    # Markers for endothelial
    end_re = refine_markers(mat_rest, markers_range$endothelial, names(mrk_old$Endo),
                            lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
    mat_rest = mat_rest[, setdiff(colnames(mat_rest), end_re$markers)]

    # Markers for neurons
    neu_re = refine_markers(mat_rest, markers_range$neurons, names(mrk_old$Neuro),
                            lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)

    # Refined markers
    markers_re = list(Astro       = ast_re$markers,
                      Oligo = oli_re$markers,
                      Micro        = mic_re$markers,
                      Endo      = end_re$markers,
                      Neuro          = neu_re$markers)

    return(markers_re)
}


my_hedge = function(ref_mtx,meta_ref,n_mrk){
    meta_ref$ind = 1:nrow(meta_ref)
    tmp_type = unique(meta_ref$deconv_clust)
    mrks = lapply(tmp_type, function(c_type){
        tmp <- apply(ref_mtx,1, function(i){
            tmp_hed = sapply(1:(length(tmp_type)-1), function(x){
                tmp_type_sub = tmp_type[-which(tmp_type ==c_type )]
                tmp_meta = meta_ref %>% dplyr::filter(deconv_clust %in% c(c_type,tmp_type_sub[x]))
                cohen.d(d = i[tmp_meta$ind],f = tmp_meta$deconv_clust,hedges.correction=TRUE)$estimate
            })
            return(min(tmp_hed))
        },simplify = TRUE)
        tmp = sort(tmp,decreasing = T)
        return(names(tmp)[1:n_mrk])
    })
    names(mrks) = tmp_type
    return(mrks)
}

