# denoise cartography

get.dist.to.knn <- function(dist.mat, ka){
	apply(dist.mat, 1, function(dist.i){
		dist.i[rank(dist.i,ties.method="first")==ka]
	})
}


denoise.cartography <- function(sf.obj, 
								ka=5, 
								tansition=1000, 
								force.recompute=FALSE){
	
	cord <- sf.obj@SpaceFold.axis[,1]
	stopifnot(!is.null(cord))
	
	
	#extract exsiting cell types in denoised cartography
	
	cell.type.old <- names(sf.obj@denoised.cartography@Z.denoised)
	
	distance.matrix <- as.matrix(dist(cord))
	
	all.cell.types <- colnames(sf.obj@data@theta)
	
	if(!force.recompute) all.cell.types <- all.cell.types[! all.cell.types %in% cell.type.old]
	
	denoise.dat.list <- lapply(all.cell.types, function(ct.i){
		
		selected.spot <- sf.obj@data@selected.spot.matrix[, ct.i]
		
		if(sum(selected.spot) < ka*2){
			#do not perform denoise if too few spots are above background
			cat("cell type ", ct.i, " has less than ", ka*2, " spots above background. It will not be denoised.\n")
			Z.denoise.i <- sf.obj@data@Z[selected.spot,,ct.i]
			Z.normed.denoise.i <- sf.obj@data@Z.normed[selected.spot,,ct.i]
		}
		else{
			distance.matrix.i <- distance.matrix[selected.spot, selected.spot]
		
			Z.i <- sf.obj@data@Z[selected.spot,,ct.i]
			Z.normed.i <- sf.obj@data@Z.normed[selected.spot,,ct.i]
		
			sigma <- get.dist.to.knn(dist.mat= distance.matrix.i, ka= ka)
		
			A.mat <- (distance.matrix.i/sigma)^2
			A.mat <- (exp(-(A.mat)) + exp(-t(A.mat)))/2
		
			A.mat <- A.mat/rowSums(A.mat)
			A.mat <- A.mat %^% tansition
		
			Z.denoise.i <- A.mat %*% Z.i
			Z.normed.denoise.i <- A.mat %*% Z.normed.i
		}
		
		list(Z.denoise = Z.denoise.i, Z.normed.denoise = Z.normed.denoise.i)
		
	})
	
	names(denoise.dat.list) <- all.cell.types
	
	
	if(!force.recompute){
		sf.obj@denoised.cartography <- new("denoisedCartography", 
										Z.denoised=c(sf.obj@denoised.cartography@Z.denoised,
													 lapply(denoise.dat.list, "[[", "Z.denoise")),
										Z.normed.denoised=c(sf.obj@denoised.cartography@Z.normed.denoised,
															lapply(denoise.dat.list, "[[", "Z.normed.denoise")),
										control_param = list(ka=c(sf.obj@denoised.cartography@control_param$ka,ka), 
															 tansition=c(sf.obj@denoised.cartography@control_param$tansition, tansition)))
	}
	else{
		sf.obj@denoised.cartography <- new("denoisedCartography", 
										Z.denoised=lapply(denoise.dat.list, "[[", "Z.denoise"),
										Z.normed.denoised =lapply(denoise.dat.list, "[[", "Z.normed.denoise"),
										control_param = list(ka=ka, tansition=tansition))
	}
	
	sf.obj	
}

