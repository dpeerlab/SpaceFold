\name{plot.beeswarm}
\alias{plot.beeswarm}
\title{
Make beeswarm plot
}
\description{
Plot the distribution of spots from the top quantile of each cell type
}
\usage{
plot.beeswarm (sf.obj,
						 palette=colorRampPalette(brewer.pal(12, "Paired")),
						 q.cut=0.95,
						 use.background=FALSE,
						 pdf.prefix="output",
						 height=5,
						 width=12,
						 mar=c(8,3,1,1),
						 cell.type.order=NULL)
}
\arguments{
  \item{sf.obj}{A S4 object of the class "SpaceFold". The @SpaceFold.axis slot should be computed.}
  \item{palette}{A function to generate the color palette of each cell type. Default uses colorRampPalette(brewer.pal(12, "Paired")).}
  \item{q.cut} {A numeric number between 0 and 1 to denote the upper quantile of spots to plot. Default=0.95, i.e., spots containing the top 5 percent of each cell type.}
  \item{use.background} {A logical vector to denote if subsetting spots above the defined background. Default=FALSE.}
  \item{pdf.prefix}{The prefix of pdf filename.}
  \item{height, width, mar}{Height width and margin of pdf.}
  \item{cell.type.order}{A character vector speficying the order of cell types for plotting. Default=NULL, which uses the colnames of theta matrix.}  
}


