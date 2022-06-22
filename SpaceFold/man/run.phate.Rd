\name{run.phate}
\alias{run.phate}
\title{
Generate SpaceFold axis using PHATE
}
\description{
Generate SpaceFold axis using PHATE
}
\usage{
run.phate (sf.obj,
					 anchorCellTypes ="all",
					 renorm.to.one = TRUE,
					 center=TRUE,
					 scale=TRUE,
					 if.pseudo.axis=TRUE,
					 if.invert=FALSE,
					 mds.solver="smacof",
					 n.jobs=20,
					 knn=10,
					 ...)
}
\arguments{
  \item{sf.obj}{A S4 object of the class "SpaceFold"}
  \item{anchorCellTypes}{A character vector to denote the names of cell types that are informative of the physical location of the underlying axis. Dafault="all". }
  \item{renorm.to.one}{A logical variable to denote if to nomalize the seleceted cell type fractions to sum to one. Default TRUE. rescale the cordinate between 0 and 1. Default=TRUE.}
  \item{center}{A logical variable to denote if to center the mean of cell type fraction to 0. Default=TRUE.}
  \item{scale}{A logical variable to denote if to scale the seleceted cell type fractions across spatial spots to sd=1. Default TRUE.}
  \item{if.pseudo.axis}{A logical variable to denote if to return a scaled cordinate with min=zero and max=one.}
  \item{if.invert}{ A logical variable to denote whehter to flip the spacefold cordinate. Default=FALSE}   \item{mds.solver}{MDS solver used by PHATE. Default="smacof". See ?phate for more details.} 
  \item{n.jobs}{ Number of threads used by PHATE. Default=20. } 
  \item{knn}{ The knn paramter used by PHATE. Default=10. See ?phate for more details. } 
  \item{...}{ additional paramters passed to phate.} 
}
\details{


}
\value{
returns an S4 object of the class "SpaceFold"
}

