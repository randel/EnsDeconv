# EnsDencov/R/prepare.R

#' Prepare data for analysis
#'
#' This function prepares single-cell and bulk RNA sequencing data for further analysis.
#'
#' @param dat A Seurat object containing single-cell RNA data.
#' @param bulk_tissue A matrix containing bulk RNA sequencing data.
#' @return A list containing processed bulk and single-cell data.
#' @export
prepare <- function(dat, bulk_tissue) {
  # Get sample names; if counts data is missing, use data instead
  if (length(dat@assays[["RNA"]]@counts@Dimnames[[2]]) > 0) {
    sample_names <- dat@assays[["RNA"]]@counts@Dimnames[[2]]
  } else {
    sample_names <- dat@assays[["RNA"]]@data@Dimnames[[2]]
  }
  
  # Create meta_list containing sample names and cell types
  meta_list <- list(sample_names, as.character(dat@meta.data[["cell_type"]]))
  names(meta_list) <- c("SamplesName", "deconv_clust")
  
  # Get the gene names from single-cell and bulk data
  genes_sc <- rownames(dat)
  genes_bulk<-rownames(bulk_tissue)
  both_genes <- intersect(genes_sc, genes_bulk)
  
  # Subset bulk tissue data, keeping only shared genes
  bulk_tissue_new <- bulk_tissue[both_genes, , drop = FALSE]
  
  # Subset single-cell data, keeping only shared genes
  dat <- dat[both_genes, ]
  
  # Extract the expression matrix from single-cell data
  if (length(dat@assays[["RNA"]]@counts@Dimnames[[2]]) > 0) {
    ref_matrix <- dat@assays[["RNA"]]@counts
  } else {
    ref_matrix <- dat@assays[["RNA"]]@data
  }
  
  # Create reference data list
  ref_list <- list(ref_matrix, meta_list)
  names(ref_list) <- c("ref_matrix", "meta_ref")
  
  # Put the reference data list into a named list
  ref_list <- list(ref_list)
  names(ref_list) <- "xx"
  
  # Create the final data list, including processed bulk data and the reference data list
  data <- list(bulk_tissue_new, ref_list)
  names(data) <- c("bulk", "ref_list")
  
  return(data)
  
}

