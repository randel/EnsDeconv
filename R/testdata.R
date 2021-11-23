#' A data example
#'
#' A data list for demonstration.
#'
#'
#' @name testdata
#' @docType data
#' @return A list containing
#' \item{count_bulk}{}
#' \item{meta_bulk}{}
#' \item{ref_list}{}
#' \item{true_frac}{}
#' @importFrom parallel stopCluster
#' @importFrom sparseMatrixStats colSums2 rowSums2 rowVars
#' @importFrom edgeR cpm
#' @importFrom caret createFolds
#' @importFrom dplyr select arrange filter arrange top_n
#' @importFrom reshape2 melt
#' @importFrom nnls nnls
#' @importFrom stats lm coef median sd pt complete.cases
#' @importFrom limma topTable
#' @importFrom ggridges stat_density_ridges theme_ridges
#' @importFrom ggpubr ggarrange annotate_figure text_grob get_palette
#' @importFrom scran findMarkers
#'
#' @examples
#'
#' data(testdata)
#'
NULL
