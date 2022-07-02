\name{select.spot}
\alias{select.spot}
\title{
Select spot manually
}
\usage{
select.spot (sf.obj,
						cell.type,
						selected.idx,
						op=c("and","or","new"))
}
\arguments{
  \item{sf.obj}{A S4 object of the class "SpaceFold".}
  \item{cell.type}{A character variable to denote the cell type from which to select slots}
  \item{selected.idx}{A logical vector of the same length as the number of spots to denote the spots to be selected.}
  \item{cell.type}{A character variable to denote the operation with the existing selection."and"=select from existing subset. "or"=take the union of the exsting and new selectiong. "new"=only use the new selection.}
  
}
\details{


}
\value{
returns an S4 object of the class "SpaceFold", with selected.spot.matrix modified in place.
}

