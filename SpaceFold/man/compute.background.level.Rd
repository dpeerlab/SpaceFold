\name{compute.background.level}
\alias{compute.background.level}
\title{
Compute background level for each cell type
}
\description{
Compute background level for each cell type by fitting a mixture model on cell type fraction theta and amount of reads of each cell type Znk.
}
\usage{
compute.background.level (sf.obj,
									        theta.cutoffs.user=NULL,
									        Znk.cutoffs.user=NULL)
}
\arguments{
  \item{sf.obj}{A S4 object of the class "SpaceFold"}
  \item{theta.cutoffs.user}{A user defined cutoff for cell type fraction. If is not null, background level will be set by taking the maximum of theta.cutoffs.user and the one determined by the mixture model. Default=NULL.}
  \item{Znk.cutoffs.user}{A user defined cutoff for amount of reads of each cell type. If is not null, background level will be set by taking the maximum of Znk.cutoffs.user and the one determined by the mixture model. Default=NULL.}
}
\details{


}
\value{
returns an S4 object of the class "SpaceFold", with @selected.spot.matrix modified in the @Data entry. 
}

