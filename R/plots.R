#' Various plots
#'
#' @importFrom dplyr "%>%"
#' @param seu Seurat object
#' @param gene Gene of interest
#' @param conditions Specify which conditions to compare, otherwise use 
#' ones specified with `Idents()` (Default: NULL)
#' @param asterisk Replace p-values with asterisks (Default: False)
#' 
#' @import ggplot2 ggsignif
#' 
#' @export
plot_violin_sig <- function(seu, gene, conditions=NULL, asterisk=F) {
  
  df <- NULL
  if (is.null(conditions)) {
    conditions <- as.vector(unique(Seurat::Idents(seu)))
  }
  for (con in conditions) {
    cond <- seu[gene, Seurat::Idents(seu) == con]
    cond_df <- data.frame(
      gene = as.vector(Seurat::GetAssayData(cond, slot = 'data')), 
      group = con, 
      split = con)
    rownames(cond_df) <- colnames(cond)
    df <- rbind(df, cond_df)
  }
  
  p <- ggplot(data=df, aes(x=factor(group, levels = conditions), y=gene)) + 
    geom_violin(aes(fill=group)) +
    geom_jitter(height = 0, width = 0.1) +
    geom_signif(comparisons = utils::combn(conditions, 2, simplify = F),
                map_signif_level = asterisk, textsize=5, test = 'wilcox.test',
                step_increase = .1) +
    labs(x = gene, y = '')
  
  return (p)
}