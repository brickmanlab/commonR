#' Manipulation with markers
#'
#' @importFrom dplyr "%>%"
#' @param seu Seurat object
#' @param filename When specified, save into a file
#' @param positive Only return positive markers (Default: True)
#' @param min_n Number of markers per identity (Default: 10)
#' @param alpha Threshhold for adjusted p-value (Default: 0.05)
#' @param logfold Minimum required logfold change (Default: 1)
#' @export
get_all_markers <- function(seu, filename=NULL, positive=T, min_n=10, alpha=0.05, logfold=1) {
  markers <- Seurat::FindAllMarkers(seu, only.pos = positive, min.pct = 0.25)
  markers_filt <- markers[markers$p_val_adj < alpha, ] %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = min_n, wt = avg_logFC)

  if (!is.null(file)) {
    print(paste0("Saving markers into ", filename))
    utils::write.csv(markers_filt, file = filename)
  }

  return(markers_filt)
}
