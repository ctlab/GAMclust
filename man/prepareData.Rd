% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{prepareData}
\alias{prepareData}
\title{Prepare gene expression matrix}
\usage{
prepareData(
  E,
  gene.id.type = NULL,
  keep.top.genes = 12000,
  use.PCA = TRUE,
  use.PCA.n = 50,
  repeats = seq_len(ncol(E)),
  network.annotation
)
}
\arguments{
\item{E}{Expression matrix with rownames as gene symbols.}

\item{gene.id.type}{Gene ID type.}

\item{keep.top.genes}{Which top of the most expressed genes to keep for the further analysis.}

\item{use.PCA}{Whether to reduce matrix dimentionality by PCA or not.}

\item{repeats}{Here you may collapse biological replicas by providing vector with repeated sample names}

\item{network.annotation}{Metabolic network annotation.}
}
\value{
Expression matrix prepared for the analysis.
}
\description{
Prepare gene expression matrix
}
