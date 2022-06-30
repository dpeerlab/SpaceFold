setClass("SpaceFoldData",
         slots = c(
           Z = "array",
           theta = "matrix",
           Znk = "matrix",
           Z.normed = "array",
           selected.spot.matrix = "matrix"
         ),
         prototype = list(
           Z = array(),
           theta = matrix(),
           Znk = matrix(),
           Z.normed = array(),
           selected.spot.matrix= matrix()
         )
)

#validator to be implemented


setClass("denoisedCartography",
         slots = c(
           Z.denoised = "list",
           Z.normed.denoised = "list",
           control_param = "list"
         ),
         prototype = list(
           Z.denoised = list(),
           Z.normed.denoised = list(),
           control_param = list(ka = NA_real_,
           						tansition = NA_real_)
         )
)


setClass("SpaceFold",
         slots = c(
           data = "SpaceFoldData",
           SpaceFold.axis = "matrix",
           denoised.cartography = "denoisedCartography",
           feature = "list",
           meta = "list",
           control_param = "list"
         ),
         prototype = list(
           data = new("SpaceFoldData"),
           SpaceFold.axis = matrix(),
           denoised.cartography = new("denoisedCartography"),
           feature = data.frame(),
           meta = data.frame(),
           control_param = list(theta.cutoffs = NA_real_,
           					 	Znk.cutoffs = NA_real_)
         )
)





new.sf <- function(bp.obj,
				   feature=NULL,
				   meta=NULL){
	
	dat.raw <- new("SpaceFoldData", 
					theta = get.theta(bp.obj),
				    Z = get.Z(bp.obj))
	
	dat.raw <- compute.Znk(dat.raw)
	dat.raw <- norm.by.sf(dat.raw)
	
	sf.obj <- new("SpaceFold", 
				   data = dat.raw)
	
	sf.obj <- add.feature (sf.obj, feature)
	sf.obj <- add.meta (sf.obj, meta)

	sf.obj
}






