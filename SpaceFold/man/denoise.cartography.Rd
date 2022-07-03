\name{denoise.cartography}
\alias{denoise.cartography}
\title{
Denoise the cartography
}
\description{
Apply data diffusion with a adaptive Gaussion kernel to denoise gene expression cartography. This is a feature of SpaceFold2.0.
}
\usage{
denoise.cartography (sf.obj, 
                     ka=5, 
                     tansition=1000，
                     force.recompute=FALSE)
}
\arguments{
  \item{sf.obj}{A S4 object of the class "SpaceFold".}
  \item{ka}{The kth nearest neighbor spot to define the sigma in the adaptive Gaussion kernel.}
  \item{tansition}{The diffusion time. Larger value yields higher smoothness.}
  \item{force.recompute}{A logical vector to control if recompute denoised data for all cell types. Default=FALSE (only recompute the newly added cell type (groups)). }
}
\details{


}
\value{
returns an S4 object of the class "SpaceFold" with the slot denoisedCartography
}

