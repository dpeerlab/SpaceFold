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
		
		if(!is.na(object@selected.spot.matrix[1,1])){
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
		
		if(!is.na(object@SpaceFold.axis[1,1]))
			cat("SpaceFold.axis computed. \n")
		 
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



#generic [ function for S4 object SpaceFold
#' @exportMethod
setMethod("[", c("SpaceFold", "integer", "missing", "ANY"),
    ## we won't support subsetting on j; dispatching on 'drop' doesn't
    ## make sense (to me), so in rebellion we'll quietly ignore it.
    function(x, i, j, ..., drop)
{
    subset.spot(sf.obj=x, sub.idx=i, drop=drop)
})

#generic [ function for S4 object SpaceFold
#' @exportMethod
setMethod("[", c("SpaceFold", "logical", "missing", "ANY"),
    ## we won't support subsetting on j; dispatching on 'drop' doesn't
    ## make sense (to me), so in rebellion we'll quietly ignore it.
    function(x, i, j, ..., drop)
{
    subset.spot(sf.obj=x, sub.idx=i, drop=drop)
})






