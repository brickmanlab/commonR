% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_significant}
\alias{plot_significant}
\title{Visualize Violin plot with significance}
\usage{
plot_significant(
  seu,
  gene,
  conditions = NULL,
  test = "wilcox.test",
  asterisk = F,
  pos_only = F,
  plot_type = "violin",
  add_jitter = F
)
}
\arguments{
\item{seu}{Seurat object}

\item{gene}{Gene of interest}

\item{conditions}{Specify which conditions to compare, otherwise use 
ones specified with `Idents()` (Default: NULL)}

\item{test}{Test (wilcox.test, t.test or custom func) (Default: wilcox.test)}

\item{asterisk}{Replace p-values with asterisks (Default: False)}

\item{pos_only}{Display only significant results (Default: False)}

\item{plot_type}{Type of plot: violin, boxplot or jitter (Default: violin)}

\item{add_jitter}{Display points (Default: False)}
}
\description{
Visualize Violin plot with significance
}
