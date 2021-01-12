#' Visualize Violin plot with significance
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

#' Visualize gene expression using boxplot
#' 
#' @importFrom dplyr "%>%" filter
#' @param seu Seurat object
#' @param features Genes of interest
#' @param slot Which matrix to use (Default: 'data')
#' @import ggplot2 reshape2
#' 
#' @export
plot_boxplot <- function(seu, features = NULL, slot = 'data') {
  mat <- Seurat::GetAssayData(seu, slot = slot) %>% as.data.frame
  clusters <- as.vector(unique(Seurat::Idents(seu)))
  genes <- intersect(rownames(seu), features)
  df <- NULL
  for (cluster in clusters) {
    cells <- rownames(data.frame(by = Seurat::Idents(seu)) %>% filter(by == cluster))
    for (gene in genes) {
      gene_exp <- reshape2::melt(mat[gene, cells])
      gene_exp$variable <- gene
      gene_exp$cluster <- cluster
      df <- rbind(df, gene_exp)
    }
  }
  ggplot(df, aes(x=variable, y=value, fill=cluster)) +
    geom_boxplot() +
    xlab("Genes") + ylab("Gene Expression")
}

#' Visualize cellular proportion 
#' 
#' @importFrom dplyr "%>%" count group_by
#' @param seu Seurat object
#' @param group.by Column name from metadata
#' @import ggplot2
#' 
#' @export
plot_proportion_bar <- function(seu, group.by = NULL) {
  df <- data.frame(
    cluster = as.vector(Seurat::Idents(seu)), 
    group = as.factor(seu@meta.data[[group.by]])
  )
  df <- df %>% group_by(cluster, group) %>% count(name = 'value')
  
  plt <- ggplot(data=df, aes(x=cluster, y=value, fill=group)) +
    geom_bar(stat='identity', position="fill") +
    scale_y_continuous(name = 'Cellular proportion (%)', labels = seq(0, 100, 25)) +
    xlab('') +
    guides(fill=guide_legend(group.by))
  
  return(plt)
}
