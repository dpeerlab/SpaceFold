
get.theta <- function(bp.obj){
	
	if(is(bp.obj,"BayesPrism")) {
		theta <- bp.obj@posterior.theta_f@theta
		#theta <- bp.obj@posterior.initial.cellType@theta
	}
	if(is(bp.obj,"BayesPrismST")) {
		theta <- bp.obj@posterior.cellType@theta	
	}
	
	theta

}


get.Z <- function(bp.obj){
	
	if(is(bp.obj,"BayesPrism")) {
		Z <- bp.obj@posterior.initial.cellType@Z
	}
	if(is(bp.obj,"BayesPrismST")) {
		Z <- bp.obj@posterior.cellType@Z	
	}
	
	Z

}

compute.Znk <- function(sfd.obj){
	
	sfd.obj@Znk <- apply(sfd.obj@Z,c(1,3),sum)
	
	sfd.obj 
}



#utility functions to add annotation / subset BayesPrism output


#' add annotation dataframe (meta) to a BayesPrism output
#' @param bp.obj a BayesPrism output object
#' @param meta a dataframe that annotates the visium spots, usually with the following columns:
	#barcode sample tissue row col imagerow  imagecol Cluster height width sum_umi sum_gene
	#add.meta will create an entry meta in bp.obj$para by matching barcode to rownames(bp.obj$para$X)
add.meta <- function(sf.obj, meta=NULL){
	
	bp.barcode <- dimnames(sf.obj@data@Z)[[1]]
	
	
	if(!is.null(meta)) {
		stopifnot(sum(bp.barcode %in% meta$barcode)==length(bp.barcode))
		meta.matched <- meta[match(bp.barcode, meta $barcode),]
	}
	else meta.matched <- data.frame(barcode= bp.barcode)
	
	sf.obj@meta <- meta.matched
	sf.obj
}


#' add gene annotation dataframe (feature) to a BayesPrism output
#' @param bp.obj a BayesPrism output object
#' @param feature a dataframe that annotates the genes (colnames of X), with rownames to be matched the colnames(bp.obj$para$X),
	#usually with the following columns: ESEMBLE ID, gene symbol, etc...
	#add.feature will create an entry feature in bp.obj$para by matching rownames of feature to colnames(bp.obj$para$X)
add.feature <- function(sf.obj, feature=NULL){
	
	bp.genes <- dimnames(sf.obj@data@Z)[[2]]
		
	if(!is.null(feature)) {
		stopifnot(sum(bp.genes %in% rownames(feature))==length(bp.genes))
		feature.matched <- feature[match(bp.genes, rownames(feature)),]
	}
	else feature.matched <- data.frame(feature = bp.genes)
	
	sf.obj@feature <- feature.matched
	sf.obj
}


#' sum over multiple cell states / cell types for theta, Znkg and Znk(if exist)
#' @param bp.obj a BayesPrism output object
#' @param grouping.list a named list containing the cell states / cell types to be summed, with name denoting the corresponding cell type 
#' @param return.merge.only a logical variable denoting if only returns the regrouped results or the full bp.obj with a new entry "XXX.regrouped" added to the $res entry
merge.cell.type <- function(sf.obj,
					   		grouping.list,
					   		get.total=TRUE){
	
	#summing over theta (or Znk)
	sum.theta <- function(theta, grouping.list){
		
		theta.sum <- matrix(NA,
							nrow=nrow(theta),
							ncol=length(grouping.list),
							dimnames=list(rownames(theta),names(grouping.list)))
		
		for(i in 1:length(grouping.list)){
			cell.types.i <- grouping.list[[i]]
			if(length(cell.types.i)==1) theta.sum[,i] <- theta[, cell.types.i]
			else theta.sum[,i] <- rowSums(theta[, cell.types.i])
		}
		
		#replace the previous group definition
		theta <- theta[,!colnames(theta) %in% colnames(theta.sum)]
				
		cbind(theta, theta.sum)
	}
	
	#summing over Znkg
	sum.Zngk <- function(Zngk, grouping.list){
				
		uniq.ct <- dimnames(Zngk)[[3]] [!dimnames(Zngk)[[3]] %in% names(grouping.list)]
		
		Zngk.all <- array(NA,
						  dim=c(dim(Zngk)[1], 
						  		dim(Zngk)[2], 
						  		length(uniq.ct)+length(grouping.list)),
						  dimnames = list(dimnames(Zngk)[[1]],
						   				  dimnames(Zngk)[[2]], 
	 			                          c(uniq.ct,names(grouping.list))))
	 	
	 	Zngk.all[,,1:length(uniq.ct)] <- Zngk[,, uniq.ct]
	 			                   
		for(i in 1:length(grouping.list)){
			cell.types.i <- grouping.list[[i]]
			Zngk.i <- Zngk[,,cell.types.i]
			if(length(dim(Zngk.i))==2) Zngk.all[,,length(uniq.ct)+i] <- Zngk.i
			else Zngk.all[,,length(uniq.ct)+i] <- rowSums(Zngk.i[,,cell.types.i, drop=F], dims=2)			
		}
		Zngk.all
	}
	
	merge.seletion.matrix <- function(selected.spot.matrix, grouping.list){
		
		selected.spot.matrix.merged <- do.call(cbind, 
												  lapply(grouping.list, 
												  		function(i) apply(selected.spot.matrix[,i,drop=F],1,sum)>0))
		colnames(selected.spot.matrix.merged) <- names(grouping.list)
		
		#replace the previous group definition
		selected.spot.matrix <- selected.spot.matrix[,!colnames(selected.spot.matrix) %in% colnames(selected.spot.matrix.merged)]
		
		selected.spot.matrix <- cbind(selected.spot.matrix , selected.spot.matrix.merged)
		selected.spot.matrix[,"total"] <- TRUE
		selected.spot.matrix
	}
	
	theta <- sf.obj@data@theta
	Z <- sf.obj@data@Z
	Znk <- sf.obj@data@Znk
	selected.spot.matrix <- sf.obj@data@selected.spot.matrix	
	
	all.cell.types <- colnames(theta)
	
	if(get.total) grouping.list <- c(grouping.list, list(total = all.cell.types))
	
	stopifnot(!is.null(names(grouping.list)))
	stopifnot(all(unlist(grouping.list) %in% all.cell.types))
	
	dat.merged <- new("SpaceFoldData",
					   Z = sum.Zngk(Z, grouping.list),
					   theta = sum.theta(theta, grouping.list),
					   Znk = sum.theta(Znk, grouping.list))
	
	dat.merged <- norm.by.sf(dat.merged)
	
	if(!is.null(selected.spot.matrix))
		dat.merged@selected.spot.matrix <- merge.seletion.matrix (selected.spot.matrix, grouping.list)
	
	
	sf.obj@data <- dat.merged
		
	sf.obj
}


#' function to compute normalized expression per cell type per spot, Znkg, by the size factor of each cell type in each spot, Znk
#' @param bp.obj a BayesPrism output object
norm.by.sf <- function(sfd.obj){
	
	Z <- sfd.obj@Z
	Znk <- sfd.obj@Znk
	
	if(is.null(Znk)) {
		sf.obj <- compute.Znk(sf.obj)
		Znk <- sfd.obj@Znk
	}
	
	Z.normed <- Z
	for(k in 1:dim(Z.normed)[3]) Z.normed[,,k] <- Z[,,k] / Znk[,k]
	
	sfd.obj@Z.normed <- Z.normed
	
	sfd.obj
}



#' subset SpaceFold object by spot. Slot denoised.cartography will not be subsutted. As denoising depends on the spot.
#' @param sf.obj a SpaceFold output object
#' @param sub.idx a logical/numeric vector 
subset.spot<- function(sf.obj, 
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


#' subset SpaceFold object by gene.
#' @param sf.obj a SpaceFold output object
#' @param sub.idx a logical/numeric vector
subset.gene <- function(sf.obj,
							gene.name,
							return.all=FALSE,
							drop=FALSE){
								
	gene.idx <- get.gene.idx (sf.obj, gene.name)
	
	sf.obj.sub <- new("SpaceFold")
	
	#subset data slot
	data <- new("SpaceFoldData")
	if(length(dim(sf.obj@data@Z))>1) data@Z <- sf.obj@data@Z[,gene.idx,,drop= drop]	
	if(length(dim(sf.obj@data@Z.normed))>1) data@Z.normed <- sf.obj@data@Z.normed[,gene.idx,,drop= drop]	
	
	if(return.all){
		if(sum(dim(sf.obj@data@theta))>2) data@theta <- sf.obj@data@theta
		if(sum(dim(sf.obj@data@Znk))>2) data@Znk <- sf.obj@data@Znk
		if(sum(dim(sf.obj@data@selected.spot.matrix))>2) data@selected.spot.matrix <- sf.obj@data@selected.spot.matrix
	}
	
	sf.obj.sub@data <- data
	
	#subset denoised.cartography slot
	denoised.data <- new("denoisedCartography")
	if(length(sf.obj@denoised.cartography@Z.denoised)>0) {
		denoised.data@Z.denoised <- lapply(sf.obj@denoised.cartography@Z.denoised, function(i)i[,gene.idx,drop=drop])
		names(denoised.data@Z.denoised) <- names(sf.obj@denoised.cartography@Z.denoised)
	}
	if(length(sf.obj@denoised.cartography@Z.normed.denoised)>0){
		denoised.data@ Z.normed.denoised <- lapply(sf.obj@denoised.cartography@ Z.normed.denoised, function(i)i[,gene.idx,drop=drop])
		names(denoised.data@ Z.normed.denoised) <- names(sf.obj@denoised.cartography@ Z.normed.denoised)
	}
	
	if(return.all)
		denoised.data@control_param <- sf.obj@denoised.cartography@control_param
	
	sf.obj.sub@denoised.cartography <- denoised.data
	
	denoised.data@control_param <- sf.obj@denoised.cartography@control_param
	
	if(nrow(sf.obj@feature)>0) 
		sf.obj.sub@feature <- sf.obj@feature[gene.idx,,drop= drop]
	
	if(return.all){
		sf.obj.sub@meta <- sf.obj@ meta
		sf.obj.sub@ SpaceFold.axis <- sf.obj@ SpaceFold.axis
		sf.obj.sub@ control_param <- sf.obj@ control_param
	}
		
	sf.obj.sub

}




combine.genes <- function(sf.obj.list, 
						  theta, 
						  Znk, 
						  selected.spot.matrix, 
						  meta, 
						  SpaceFold.axis){
	
	all.identical <- function(l) all(mapply(identical, head(l, 1), tail(l, -1)))
	
	stopifnot(all.identical(lapply(sf.obj.list, function(i)dimnames(i@data@Z)[[1]])))
	stopifnot(all.identical(lapply(sf.obj.list, function(i)dimnames(i@data@Z)[[3]])))
	
	all.genes <- unlist(lapply(sf.obj.list, function(i)dimnames(i@data@Z)[[2]]))
	all.spots <- dimnames(sf.obj.list[[1]]@data@Z)[[1]]
	all.cell.types <- dimnames(sf.obj.list[[1]]@data@Z)[[3]]
	
		
	sf.merged <- new("SpaceFold")
	sf.merged@meta <- meta
	sf.merged@ SpaceFold.axis <- SpaceFold.axis

	sf.merged@data@theta <- theta
	sf.merged@data@Znk <- Znk
	sf.merged@data@selected.spot.matrix <- selected.spot.matrix

	sf.merged@data@Z <- array(dim=c(length(all.spots), length(all.genes),length(all.cell.types)),
								dimnames=list(all.spots, all.genes, all.cell.types))
	sf.merged@data@Z.normed <- sf.merged@data@Z
	
	for(i in 1:length(sf.obj.list)){
		sf.merged@data@Z[,i,] <- sf.obj.list[[i]]@data@Z
		sf.merged@data@Z.normed[,i,] <- sf.obj.list[[i]]@data@Z.normed
	}
	
	
	sf.merged@denoised.cartography@Z.denoised <- lapply(all.cell.types, function(ct.i){
											 do.call(cbind, lapply(sf.obj.list, function(sf.i) sf.i@denoised.cartography@Z.denoised[[ct.i]]))
											 })
	names(sf.merged@denoised.cartography@Z.denoised) <- all.cell.types
	
	sf.merged@denoised.cartography@Z.normed.denoised <- lapply(all.cell.types, function(ct.i){
											 do.call(cbind, lapply(sf.obj.list, function(sf.i) sf.i@denoised.cartography@Z.normed.denoised[[ct.i]]))
											 })
	names(sf.merged@denoised.cartography@Z.normed.denoised) <- all.cell.types
	
	sf.merged@feature <- do.call(rbind, lapply(sf.obj.list, function(sf.i) sf.i@feature))	
	
	sf.merged
	
}









