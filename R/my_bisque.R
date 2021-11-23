#Generate reference profile for cell types identified in single-cell data
#returnsc.ref Matrix. Reference profile with number of gene rows by number
#  of cell types columns.
#
my_GenerateSCReference <- function(count_sc,meta_sc, cell.types) {
  cell.labels <- factor(meta_sc[[cell.types]])
  all.cell.types <- levels(cell.labels)
  aggr.fn <- function(cell.type) {
    rowMeans2(count_sc[,cell.labels == cell.type, drop=F])
  }
  template <- numeric(nrow(count_sc))
  sc.ref <- vapply(all.cell.types, aggr.fn, template)
  rownames(sc.ref) = rownames(count_sc)
  return(sc.ref)
}

my_GenerateSCReference_MED <- function(count_sc,meta_sc, cell.types) {
  cell.labels <- factor(meta_sc[[cell.types]])
  all.cell.types <- levels(cell.labels)
  aggr.fn <- function(cell.type) {
    rowMedians(count_sc[,cell.labels == cell.type, drop=F])
  }
  template <- numeric(nrow(count_sc))
  sc.ref <- vapply(all.cell.types, aggr.fn, template)
  rownames(sc.ref) = rownames(count_sc)
  return(sc.ref)
}

#Calculate cell proportions based on single-cell data
#'
#Returns proportion of each cell type out of total cells for each individual
#in the single-cell Expression
#param subject.names A character string. Name of phenoData attribute in
#  sc.eset that indicates individual ID.
#param cell.types A character string. Name of phenoData attribute in sc.eset
#  that indicates cell type
#returnsc.props Matrix. Cell proportions with number of cell types rows
#  by number of individuals columns
#
my_CalculateSCCellProportions <- function(count_sc,meta_sc, subject.names, cell.types) {
  individual.labels <- factor(meta_sc[[subject.names]])
  individuals <- levels(individual.labels)
  cell.labels <- as.factor(meta_sc[[cell.types]])
  aggr.fn <- function(individual) {
    table(cell.labels[individual.labels == individual]) /
      length(cell.labels[individual.labels == individual])
  }
  sc.props <- sapply(individuals, aggr.fn)
  return(sc.props)
}

SemisupervisedTransformBulk <- function(gene, Y.train, X.pred) {
  # Learns linear transformation of observed bulk to match distribution of
  # weighted sum of reference
  #
  # Used with vapply, processes one gene
  Y.train.scaled <- scale(Y.train[gene,,drop=T])
  Y.center <- attr(Y.train.scaled, "scaled:center")
  Y.scale <- attr(Y.train.scaled, "scaled:scale")
  n <- length(Y.train.scaled)
  # Shrinkage estimator that minimizes MSE for scaling factor
  shrink.scale <- sqrt(sum((Y.train[gene,,drop=T]-Y.center)^2)/n+1)
  X.pred.scaled <- scale(X.pred[gene,,drop=T])
  Y.pred <- matrix((X.pred.scaled * shrink.scale) + Y.center,
                         dimnames=list(colnames(X.pred), gene))
  return(Y.pred)
}

my_bisque <- function(                  count_bulk,
                                        count_sc,
                                        meta_sc,
                                        markers=NULL,
                                        cell.types="cellType",
                                        subject.names="SubjectName",
                                        verbose=TRUE) {
  if (! cell.types %in% colnames(meta_sc)) {
    stop(sprintf("Cell type label \"%s\" ", cell.types),
         "not found in single-cell meta data.")
  }
  else if (! subject.names %in% colnames(meta_sc)) {
    stop(sprintf("Individual label \"%s\"", subject.names),
         " not found in single-cell meta data.")
  }
  n.sc.individuals <- length(unique(meta_sc[[subject.names]]))

  if (n.sc.individuals == 1) {
    stop("Only one individual detected in single-cell data. At least ",
         "two subjects are needed (three or more recommended).")
  }
  else if (n.sc.individuals == 2) {
    warning("Only two individuals detected in single-cell data. While ",
            "Bisque will run, we recommend at least three subjects for",
            " reliable performance.")
  }
  n.cell.types <- length(unique(meta_sc[[cell.types]]))

  if (n.cell.types == 1) {
    stop("Single-cell pheno data indicates only one cell type",
         " present. No need for decomposition.")
  }
  if (verbose) {
    message(sprintf("Decomposing into %i cell types.",
                    n.cell.types))
  }
  markers <- rownames(count_sc)

  genes <- intersect(rownames(count_sc),rownames(count_bulk))

  sc.ref <- my_GenerateSCReference(count_sc,meta_sc,cell.types)[genes,,drop=F]
  sc.props <- my_CalculateSCCellProportions(count_sc,meta_sc, subject.names, cell.types)
  sc.props <- sc.props[colnames(sc.ref),,drop=F]

  if (verbose) {
    message("Inferring bulk transformation from single-cell alone.")
  }
  Y.train <- sc.ref %*% sc.props
  # X.pred is the bulk for the remaining samples to be decomposed.
  X.pred <- count_bulk[genes,,drop=F]
  sample.names <- colnames(count_bulk)
  template <- numeric(length(sample.names))
  names(template) <- sample.names
  if (verbose) {
    message("Applying transformation to bulk samples and decomposing.")
  }
  # Y.pred is the transformed bulk for samples to be decomposed.
  Y.pred <- matrix(vapply(X=genes,
                          FUN=SemisupervisedTransformBulk,
                          FUN.VALUE=template,
                          Y.train, X.pred,
                          USE.NAMES=TRUE),
                   nrow=length(sample.names))

  # Columns in Y.pred with NaN indicate transformation could not be learned
  #   for that gene.
  indices <- apply(Y.pred, MARGIN=2,
                   FUN=function(column) {anyNA(column)})
  if (any(indices)) {
    if (verbose) {
      n.dropped <- sum(indices)
      message(sprintf("Dropped an additional %i genes", n.dropped),
              " for which a transformation could not be learned.")
    }
    if (sum(!indices) == 0) {
      stop("Zero genes left for decomposition.")
    }
    Y.pred <- Y.pred[,!indices,drop=F]
    sc.ref <- sc.ref[!indices,,drop=F]
  }
  # limsolve nnls matrices and vectors
  E <- matrix(1,nrow=n.cell.types, ncol=n.cell.types)
  f <- rep(1, n.cell.types)
  G <- diag(n.cell.types)
  h <- rep(0, n.cell.types)
  results <- as.matrix(apply(Y.pred, 1,function(b) {
    sol <- limSolve::lsei(sc.ref, b,
                          E, f, G, h)
    sol.p <- sol$X
    sol.r <- sqrt(sol$solutionNorm)
    return(append(sol.p, sol.r))
  }))
  rownames(results) <- append(colnames(sc.ref), "rnorm")
  colnames(results) <- sample.names
  rnorm <- results["rnorm",,drop=T]
  names(rnorm) <- sample.names
  Y.pred <- t(Y.pred)
  rownames(Y.pred) <- rownames(sc.ref)
  colnames(Y.pred) <- sample.names
  results <- list(bulk.props=results[colnames(sc.ref),,drop=F],
                  sc.props=sc.props,
                  rnorm=rnorm,
                  genes.used=rownames(sc.ref),
                  transformed.bulk=Y.pred)
  return(results)
}
