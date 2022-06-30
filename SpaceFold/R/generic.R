#define generic functions for S4 classes

#' generic show function for an S4 object SpaceFoldData
#' @exportMethod
setMethod("show", "SpaceFoldData",
	function(object) {
		
		cat("Number of spatial spots: ",
			nrow(object@theta), "\n")
		
		cat("Number of genes: ",
			dim(object@Z)[[2]], "\n")

		cat("Fraction of cell types: \n")
		theta.summary <- apply(object@theta, 2, summary)	
		print(round(theta.summary,3))
	
		cat("Total reads of each cell type : \n")
		Znk.summary <- apply(object@Znk, 2, summary)	
		print(round(Znk.summary,3))
		
		if(!is.na(object@selected.spot.matrix)){
			cat("Number of spots above background for each cell type: \n")
			print(apply(object@selected.spot.matrix,2,sum))
		}
				 
   }
)





#' generic show function for an S4 object SpaceFold
#' @exportMethod
setMethod("show", "SpaceFold",
	function(object) {
		cat("Data info: \n")
		show(object@data)
		cat("\n")
		
		if(!is.na(object@SpaceFold.axis))
			cat("SpaceFold.axis computed.")
		 
		if(length(object@denoised.cartography@Z.denoised)){
			cat("Cartography denoised. \n")
			cat("Denoising parameter. \n")
			print(object@denoised.cartography@control_param)
		}
		
		cat("SpaceFold parameter. \n")
		print(object@control_param)	
   }
)



#generic $ function for S4 object SpaceFold
#' @exportMethod
setMethod("$", "SpaceFold",
    function(x, name)
{
    x@meta[,name]
})



#' subset SpaceFold object by spot. Slot denoised.cartography will not be subsutted. As denoising depends on the spot.
#' @param sf.obj a SpaceFold output object
#' @param sub.idx a logical/numeric vector 
subset.data <- function(sf.obj, 
					   sub.idx,
					   drop){
					   		
	#subset data slot
	data <- new("SpaceFoldData")
	
	if(length(dim(sf.obj@data@Z))>1) data@Z <- sf.obj@data@Z[sub.idx,,,drop= drop]	
	if(sum(dim(sf.obj@data@theta))>2) data@theta <- sf.obj@data@theta[sub.idx,,drop= drop]
	if(sum(dim(sf.obj@data@Znk))>2) data@Znk <- sf.obj@data@Znk[sub.idx,,drop= drop]
	if(length(dim(sf.obj@data@Z.normed))>1) data@Z.normed <- sf.obj@data@Z.normed[sub.idx,,,drop= drop]
	if(sum(dim(sf.obj@data@selected.spot.matrix))>2) data@selected.spot.matrix <- sf.obj@data@selected.spot.matrix[sub.idx,,drop= drop]
	sf.obj@data <- data
	
	#subset SpaceFold.axis slot
	if(!is.na(sf.obj@SpaceFold.axis)) 
		sf.obj@SpaceFold.axis<- sf.obj@SpaceFold.axis[sub.idx,,drop= drop]
	
	if(nrow(sf.obj@meta)>0) 
		sf.obj@meta<- sf.obj@meta[sub.idx,,drop= drop]
		
	sf.obj
}


#generic [ function for S4 object SpaceFold
#' @exportMethod
setMethod("[", c("SpaceFold", "integer", "missing", "ANY"),
    ## we won't support subsetting on j; dispatching on 'drop' doesn't
    ## make sense (to me), so in rebellion we'll quietly ignore it.
    function(x, i, j, ..., drop)
{
    subset.data(sf.obj=x, sub.idx=i, drop=drop)
})

#generic [ function for S4 object SpaceFold
#' @exportMethod
setMethod("[", c("SpaceFold", "logical", "missing", "ANY"),
    ## we won't support subsetting on j; dispatching on 'drop' doesn't
    ## make sense (to me), so in rebellion we'll quietly ignore it.
    function(x, i, j, ..., drop)
{
    subset.data(sf.obj=x, sub.idx=i, drop=drop)
})











