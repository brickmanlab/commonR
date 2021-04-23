#' @title scRNAseq marker database
#' @description This is a marker database of all external scRNAseq experiments in the scViz app. LFC > 1 and padj < 0.05
#' @format A data frame with 62 rows and 3 variables:
#' \describe{
#'   \item{\code{publication}}{factor Study from where the scRNAseq data was taken from}
#'   \item{\code{cluster}}{factor Cell type or cluster}
#'   \item{\code{genes}}{factor Comma separated vector containing gene markers} 
#'}
#' @details The markers in this database are calculated per study, that is, the data was not integrated.
"marker_database"