############ Construct Design Matrix, Library Size, and Subject-level Variation of Relative abundance for MuSiC ############
## These functions are for cell type specific mean expression, cross-subject variance and mean library size for MuSiC deconvolution

# Cross-subject Mean of Relative Abudance
#
# This function is for calculating the cross-subject mean of relative abundance for selected cell types.
#
#  x ExpressionSet, single cell dataset
#  non.zero logical, if true, remove all gene with zero expression
#  markers vector or list of gene names
#  clusters character, the phenoData used as clusters
#  samples character,the phenoData used as samples
#  select.ct vector of cell types included, default as \code{NULL}. If \code{NULL}, include all cell types in \code{x}
# @return gene by cell type matrix of average relative abundance
#

music_M.theta = function(x, non.zero, markers, clusters, samples, select.ct){
  if(!is.null(select.ct)){
    s.ct = sampleNames(x)[as.character(pVar(x, clusters)) %in% select.ct]
    x <- x[, s.ct, drop = FALSE]
  }
  if(non.zero){  ## eliminate non expressed genes
    nz.gene = rownames(x)[( rowSums(exprs(x)) != 0 )]
    x <- x[nz.gene, , drop = FALSE]
  }

  clusters <- as.character(pVar(x, clusters))
  samples <- as.character(pVar(x, samples))
  M.theta <- sapply(unique(clusters), function(ct){
    my.rowMeans(sapply(unique(samples), function(sid){
      y = exprs(x)[,clusters %in% ct & samples %in% sid, drop = FALSE]
      rowSums(y)/sum(y)
    }), na.rm = TRUE)
  })

  if(!is.null(select.ct)){
    m.ct = match(select.ct, colnames(M.theta))
    M.theta = M.theta[, m.ct]
  }

  if (!is.null(markers)){
    ids <- intersect(unlist(markers), rownames(x))
    m.ids = match(ids, rownames(x))
    M.theta <- M.theta[m.ids, ]
  }
  return(M.theta)
}

# Subject and cell type specific relative abudance
#
# This function is for calculating the subject and cell type specific relative abundance for selected cell types.
#
#  x ExpressionSet, single cell dataset
#  non.zero logical, defualt as F. If true, remove all gene with zero expression
#  markers vector or list of gene names
#  clusters character, the phenoData used as clusters
#  samples character,the phenoData used as samples
#  select.ct vector of cell types included, default as \code{NULL}. If \code{NULL}, include all cell types in \code{x}
# @return gene*subject by cell type matrix of relative abundance
#
# @export
music_Theta <- function(x, non.zero = FALSE, clusters, samples, select.ct = NULL){
  if(!is.null(select.ct)){
    s.ct = sampleNames(x)[as.character(pVar(x, clusters)) %in% select.ct]
    x <- x[, s.ct, drop = FALSE]
  }
  if(non.zero){
    nz.gene = rownames(x)[(rowSums(exprs(x)) != 0)]
    x <- x[nz.gene, , drop = FALSE]
  }
  nGenes = nrow(x);

  clusters <- as.character(pVar(x, clusters))
  samples <- as.character(pVar(x, samples))
  Theta <- sapply(unique(clusters), function(ct){
    sapply(unique(samples), function(sid){
      y = exprs(x)[,clusters %in% ct & samples %in% sid, drop = FALSE]
      rowSums(y)/sum(y)
    })
  })
  n.ct = length(unique(clusters));
  if(!is.null(select.ct)){
    m.ct = match(select.ct, colnames(Theta))
    Theta = Theta[, m.ct]
    n.ct = length(select.ct)
  }

  return(Theta = Theta)
}

# Cross-subject Corvriance of Relative Abudance
#
# This function is for calculating the cross-subject covariance of relative abundance for selected cell types.
#
#  x ExpressionSet, single cell dataset
#  non.zero logical, if true, remove all gene with zero expression
#  markers vector or list of gene names
#  clusters character, the phenoData used as clusters
#  samples character,the phenoData used as samples
#  select.ct vector of cell types included, default as \code{NULL}. If \code{NULL}, include all cell types in \code{x}
# @return celltype^2 by gene matrix of covariance
#
# @export
music_Sigma.ct = function(x, non.zero, markers, clusters, samples, select.ct){
  if(!is.null(select.ct)){
    s.ct = sampleNames(x)[as.character(pVar(x, clusters)) %in% select.ct]
    x <- x[, s.ct, drop = FALSE]
  }
  if(non.zero){  ## eliminate non expressed genes
    nz.gene = rownames(x)[( rowSums(exprs(x)) != 0 )]
    x <- x[nz.gene, , drop = FALSE]
  }
  nGenes = nrow(x);

  clusters <- as.character(pVar(x, clusters))
  samples <- as.character(pVar(x, samples))
  Sigma <- sapply(unique(clusters), function(ct){
    sapply(unique(samples), function(sid){
      y = exprs(x)[,clusters %in% ct & samples %in% sid, drop = FALSE]
      rowSums(y)/sum(y)
    })
  })
  n.sub = length(unique(samples));
  if(!is.null(select.ct)){
    m.ct = match(select.ct, colnames(Sigma))
    Sigma = Sigma[, m.ct]
    n.ct = length(select.ct)
  }
  Sigma.ct = sapply(1:nGenes, function(g){cov(Sigma[nGenes*(0:(n.sub-1)) + g, ])})
  if (!is.null(markers)){
    ids <- intersect(unlist(markers), rownames(x))
    m.ids = match(ids, rownames(x))
    Sigma.ct <- Sigma.ct[ , m.ids]
  }
  return(Sigma.ct = Sigma.ct)
}

# Cross-subject Varirance of Relative Abudance
#
# This function is for calculating the cross-subject variance of relative abundance for selected cell types.
#
#  x ExpressionSet, single cell dataset
#  non.zero logical, if true, remove all gene with zero expression
#  markers vector or list of gene names
#  clusters character, the phenoData used as clusters
#  samples character,the phenoData used as samples
#  select.ct vector of cell types included, default as \code{NULL}. If \code{NULL}, include all cell types in \code{x}
# @return gene by cell type matrix of variance
#
# @export
music_Sigma = function(x, non.zero, markers, clusters, samples, select.ct){
  if(!is.null(select.ct)){
    s.ct = sampleNames(x)[as.character(pVar(x, clusters)) %in% select.ct]
    x <- x[, s.ct, drop = FALSE]
  }
  if(non.zero){  ## eliminate non expressed genes
    nz.gene = rownames(x)[( rowSums(exprs(x)) != 0 )]
    x <- x[nz.gene, , drop = FALSE]
  }

  clusters <- as.character(pVar(x, clusters))
  samples <- as.character(pVar(x, samples))
  Sigma <- sapply(unique(clusters), function(ct){
    apply(sapply(unique(samples), function(sid){
      y = exprs(x)[,clusters %in% ct & samples %in% sid, drop = FALSE]
      rowSums(y)/sum(y)
    }), 1, var, na.rm = TRUE)
  })

  if(!is.null(select.ct)){
    m.ct = match(select.ct, colnames(Sigma))
    Sigma = Sigma[, m.ct]
  }

  if (!is.null(markers)){
    ids <- intersect(unlist(markers), rownames(x))
    m.ids = match(ids, rownames(x))
    Sigma <- Sigma[m.ids, ]
  }
  return(Sigma = Sigma)
}

# Cell type specific library size
#
# This function is for calculating the cell type specific library size for selected cell types.
#
#  x ExpressionSet, single cell dataset
#  non.zero logical, if true, remove all gene with zero expression
#  clusters character, the phenoData used as clusters
#  samples character,the phenoData used as samples
#  select.ct vector of cell types included, default as \code{NULL}. If \code{NULL}, include all cell types in \code{x}
# @return subject by cell type matrix of library
#
# @export
music_S = function(x, non.zero, clusters, samples, select.ct){
  if(!is.null(select.ct)){
    s.ct = sampleNames(x)[as.character(pVar(x, clusters)) %in% select.ct]
    x <- x[, s.ct, drop = FALSE]
  }
  if(non.zero){  ## eliminate non expressed genes
    nz.gene = rownames(x)[( rowSums(exprs(x)) != 0 )]
    x <- x[nz.gene, , drop = FALSE]
  }

  clusters <- as.character(pVar(x, clusters))
  samples <- as.character(pVar(x, samples))

  S <- sapply(unique(clusters), function(ct){
    my.rowMeans(sapply(unique(samples), function(sid){
      y = exprs(x)[, clusters %in% ct & samples %in% sid, drop = FALSE]
      sum(y)/ncol(y)
    }), na.rm = TRUE)
  })
  S[S == 0] = NA
  M.S = colMeans(S, na.rm = TRUE)

  if(!is.null(select.ct)){
    m.ct = match(select.ct, colnames(S))
    S = S[, m.ct]
  }
  return(S = S)
}

# Cell type specific library size
#
# This function is for calculating the cell type specific library size for selected cell types.
#
# @inheritParams music_S
# @inheritParams music_M.theta
#  x ExpressionSet, single cell dataset
#  non.zero logical, if true, remove all gene with zero expression
#  clusters character, the phenoData used as clusters
#  samples character,the phenoData used as samples
#  select.ct vector of cell types included, default as \code{NULL}. If \code{NULL}, include all cell types in \code{x}
# @return subject by cell type matrix of library
#
# @export
# @seealso
# \code{\link{music_S}}, \code{\link{music_M.theta}},
music_Design.matrix = function(x, non.zero, markers, clusters, samples, select.ct){
  S = music_S(x = x, non.zero = non.zero, clusters = clusters, samples = samples, select.ct = select.ct)
  M.theta = music_M.theta(x = x, non.zero = non.zero, markers = markers, clusters = clusters, samples = samples,
                           select.ct = select.ct)
  S[S == 0] = NA
  M.S = colMeans(S, na.rm = TRUE)
  D <- t(t(M.theta)*M.S)
  return(D)
}

# Prepare Design matrix and Cross-subject Variance for MuSiC Deconvolution
#
# This function is used for generating cell type specific cross-subject mean and variance for each gene. Cell type specific library size is also calcualted.
#
#  x ExpressionSet, single cell dataset
#  non.zero logical, default as TRUE. If true, remove all gene with zero expression.
#  markers vector or list of gene names. Default as NULL. If NULL, then use all genes provided.
#  clusters character, the phenoData used as clusters;
#  samples character,the phenoData used as samples;
#  select.ct vector of cell types. Default as NULL. If NULL, then use all cell types provided.
#  cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data.
#  ct.cov logical. If TRUE, use the covariance across cell types.
#  verbose logical, default as TRUE.
# @return a list of
#     * gene by cell type matrix of Design matrix
#     * subject by celltype matrix of Library size
#     * vector of average library size for each cell type
#     * gene by celltype matrix of average relative abudance
#     * gene by celltype matrix of cross-subject variation
#
# @export
mymusic_basis = function(count_sc,meta_sc, non.zero = TRUE, markers = NULL, clusters, samples, select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE){
  if(!is.null(select.ct)){
    s.ct = colnames(count_sc)[as.character(meta_sc[[clusters]]) %in% select.ct]
    count_sc <- count_sc[, s.ct, drop = FALSE]
    meta_sc = meta_sc[as.character(meta_sc[[clusters]]) %in% select.ct,]
  }
  if(non.zero){  ## eliminate non expressed genes
    nz.gene = rownames(count_sc)[( rowSums2(count_sc) != 0 )]
    count_sc <- count_sc[nz.gene, , drop = FALSE]
  }

  clusters <- as.character(meta_sc[[clusters]])
  samples <- as.character(meta_sc[[samples]])

  M.theta <- sapply(unique(clusters), function(ct){
    my.rowMeans(sapply(unique(samples), function(sid){
      y = count_sc[,clusters %in% ct & samples %in% sid, drop = FALSE]
      rowSums2(y)/sum(y)
    }), na.rm = TRUE)
  })
  rownames(M.theta) = rownames(count_sc)
  if(verbose){message("Creating Relative Abudance Matrix...")}
  if(ct.cov){
    nGenes = nrow(count_sc);
    n.ct = length(unique(clusters));
    nSubs = length(unique(samples))

    Theta <- sapply(unique(clusters), function(ct){
      sapply(unique(samples), function(sid){
        y = count_sc[,clusters %in% ct & samples %in% sid, drop = FALSE]
        return( rowSums2(y)/sum(y) )
      })
    })
    if(!is.null(select.ct)){
      m.ct = match(select.ct, colnames(Theta))
      Theta = Theta[, m.ct]
    }

    Sigma.ct = sapply(1:nGenes, function(g){
      sigma.temp = Theta[nGenes*(0:(nSubs - 1)) + g, ];
      Cov.temp = cov(sigma.temp)
      Cov.temp1 = cov(sigma.temp[rowSums(is.na(Theta[nGenes*(0:(nSubs - 1)) + 1, ])) == 0, ])
      Cov.temp[which(colSums(is.na(sigma.temp))>0), ] = Cov.temp1[which(colSums(is.na(sigma.temp))>0), ]
      Cov.temp[, which(colSums(is.na(sigma.temp))>0)] = Cov.temp1[, which(colSums(is.na(sigma.temp))>0)]
      return(Cov.temp)
    })
    colnames(Sigma.ct) = rownames(count_sc);

    if (!is.null(markers)){
      ids <- intersect(unlist(markers), rownames(count_sc))
      m.ids = match(ids, rownames(count_sc))
      Sigma.ct <- Sigma.ct[ , m.ids]
    }
    if(verbose){message("Creating Covariance Matrix...")}
  }else{
    Sigma <- sapply(unique(clusters), function(ct){
      apply(sapply(unique(samples), function(sid){
        y = count_sc[,clusters %in% ct & samples %in% sid, drop = FALSE]
        rowSums2(y)/sum(y)
      }), 1, var, na.rm = TRUE)
    })
    rownames(Sigma) = rownames(count_sc)
    if(!is.null(select.ct)){
      m.ct = match(select.ct, colnames(Sigma))
      Sigma = Sigma[, m.ct]
    }

    if (!is.null(markers)){
      ids <- intersect(unlist(markers), rownames(count_sc))
      m.ids = match(ids, rownames(count_sc))
      Sigma <- Sigma[m.ids, ]
    }
    if(verbose){message("Creating Variance Matrix...")}
  }

  S <- sapply(unique(clusters), function(ct){
    my.rowMeans(sapply(unique(samples), function(sid){
      y = count_sc[, clusters %in% ct & samples %in% sid, drop = FALSE]
      sum(y)/ncol(y)
    }), na.rm = TRUE)
  })
  if(verbose){message("Creating Library Size Matrix...")}

  S[S == 0] = NA
  M.S = colMeans(S, na.rm = TRUE)
  #S.ra = relative.ab(S, by.col = FALSE)
  #S.ra[S.ra == 0] = NA
  #S[S == 0] = NA
  #M.S = mean(S, na.rm = TRUE)*ncol(S)*colMeans(S.ra, na.rm = TRUE)

  if(!is.null(cell_size)){
    if(!is.data.frame(cell_size)){
      stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
    }else if(sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))){
      stop("Cell type names in cell_size must match clusters")
    }else if (any(is.na(as.numeric(cell_size[, 2])))){
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]),]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }

  D <- t(t(M.theta)*M.S)

  if(!is.null(select.ct)){
    m.ct = match(select.ct, colnames(D))
    D = D[, m.ct]
    S = S[, m.ct]
    M.S = M.S[m.ct]
    M.theta = M.theta[, m.ct]
  }

  if (!is.null(markers)){
    ids <- intersect(unlist(markers), rownames(count_sc))
    m.ids = match(ids, rownames(count_sc))
    D <- D[m.ids, ]
    M.theta <- M.theta[m.ids, ]
  }

  if(ct.cov){
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, Sigma.ct = Sigma.ct))
  }else{
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, Sigma = Sigma))
  }
}

