\name{merge.cell.type}
\alias{merge.cell.type}
\title{
Merge cell types
}
\description{
Sum over multiple cell types for theta, Z and Znk, perform normalization to get Z.normed, and take the union of @data@selected.spot.matrix.
}
\usage{
merge.cell.type(sf.obj,
					      grouping.list)
}
\arguments{
  \item{sf.obj}{A S4 object of the class "SpaceFold".}
  \item{grouping.list}{A named list containing the (character vector of) cell types to be summed. The names(grouping.list) corresponds to the names of merged cell type. When define names for merged cell types, avoid the use of names already appeared in cell type names.}
}
\details{


}
\value{
returns an S4 object of the class "SpaceFold". The cell type dimension of theta, Z, Z.normed, Znk and selected.spot.matrix will expand corrspondingly. 
}

