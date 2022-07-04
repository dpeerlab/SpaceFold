\name{new.sf}
\alias{new.sf}
\title{
Constructor of SpaceFold S4 object
}
\description{
Create a new SpaceFold object from a BayesPrism object.
}
\usage{
new.sf (bp.obj,
				feature=NULL,
				meta=NULL)
}
\arguments{
  \item{bp.obj}{A S4 object of the class "BayesPrism"}
  \item{feature}{A dataframe to annotate genes. rownames to be matched the gene names used by BayesPrism, i.e., colnames(bp.obj$prism$mixture). It usually contains the following columns: ENSEMBL ID, gene symbol, etc... depending on the type of input genes used for plotting later. Recommend to have ENSEMBL IDs and gene symbols as columns.}
  \item{meta}{A dataframe to annotate spatial spot (optional). Required if wants to subset the SpaceFold object.}
}
\details{


}
\value{
returns an S4 object of the class "SpaceFold"
}

