\name{plot.cartography}
\alias{plot.cartography}
\title{
The function to plot gene expression cartography over the SpaceFold axis.
}
\description{
The function to plot gene expression cartography over the SpaceFold axis. If pdfjam is installed, it will output one pdf showing an array  of panels of the dimension #of cell types x #number of genes. Otherwise, the output will be one pdfs with each page containing an individual panel.
}
\usage{
plot.cartography (sf.obj,
						raw.or.norm=c("raw","norm"),
						selected.genes,
						selected.cell.types,
						denoise=TRUE,
						bin.by=c("equal.size","equal.step"),
						n.bins=20,
						span=0.75,
						show.raw =FALSE,
						overlay.hist=TRUE,
						col.hist=adjustcolor("grey",0.5),
						pdf.prefix,
						...)
}
\arguments{
  \item{sf.obj}{A SpaceFold object}
  \item{raw.or.norm}{A character variable to denote whether plot raw expression (deconvovled reads) or   expression normalized by total expression of each cell type.}
  \item{selected.genes}{A character vector of gene symbols or ENSEMBL IDs. (automatically determines based on sf.obj@feature) }
  \item{selected.cell.types}{A character vector of selected cell types or names of grouped cell types.}
  \item{denoise}{A logical variable to denote whether to plot cartography denoised by markov diffusion (SpaceFold v2.0) or the sliding window average (SpaceFold v1.0, used by the paper). Default=TRUE.}
  \item{bin.by}{How to bin spots. Used when denoise=FALSE. "equal.size": each bin has the same number of spots. "equal.step": each bin has the same length on the SpaceFold axis. Default="equal.size".}
  \item{span}{The parameter which controls the degree of smoothing by loess. Default=0.75. If denoise=TRUE, we recommend to use 0.25.}
  \item{show.raw}{A logical variable to denote whether to overlay the raw deconvovled expression in each spot. Default=FALSE.}
  \item{overlay.hist}{A logical variable to denote whether to overlay histogram of cell type distribution. Default=FALSE.}
  \item{col.hist}{Color of histogram.}
  \item{pdf.prefix}{prefix of pdf output}
  \item{return.raw}{A logical variable to denote whether to return the data used to generate the plot. Default=FALSE.}
}
\details{


}
\value{
}

