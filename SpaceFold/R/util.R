
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




# # 
# #' add labels, e.g. created by selecting spots in the Loupe browser, to the meta entry of the BayesPrism output
# #' @param bp.obj a BayesPrism output object
# #' @label.df a dataframe (generated by selecting spots from loupe browser). First column is barcode, second column is the label. 
# add.region.labels <- function(bp.obj, label.df){
	# stopifnot(!is.null(bp.obj$para$meta))
	
	# colnames(label.df)[colnames(label.df)=="Barcode"] <- "barcode"
	
	# meta.new <- merge(bp.obj$para$meta, 
									  # label.df, 
									  # by = "barcode",
									  # all.x=TRUE)
	
	# meta.new[,ncol(meta.new)][is.na(meta.new[,ncol(meta.new)])] <- "unassigned"
	
	# bp.obj$para$meta <- meta.new
	
	# bp.obj
# }










# #' subset BayesPrism object based on region
#' @param sf.obj a SpaceFold output object
#' @param col.name the column name of the meta dataframe
#' @param val the value on the right hand side comparison operator
#' @param operator the character denoting the comaprison, e.g. "==", ">", "<". default is "==".
subset.sf <- function(sf.obj, 
					  col.name, 
					  val,
					  operator="=="){
	
	if(is.character(val)) val <- paste("\'", val, "\'",sep="")
	
	exp.text <- paste( "sf.obj@meta$", col.name, " ", operator, " ",  val,sep="")
	sub.idx <- eval ( parse(text = exp.text)  )
	
	subset.data(sf.obj, sub.idx)
	
}


# # 
# #' update cell type names in the BayesPrism object
# #' @param bp.obj a BayesPrism output object
# #' @param old.name.vec a character vector containing the old names of cell types to be substituted
# #' @param new.name.vec a character vector containing the new names of cell types used to replace the old names
# update.name <- function(bp.obj, old.name.vec, new.name.vec){
	# stopifnot(length(old.name.vec)==length(new.name.vec))
	
	# current.cell.types.all  <- colnames(bp.obj$res$first.gibbs.res$theta.merged)
	# new.name.vec <- new.name.vec[old.name.vec %in% current.cell.types.all]	
	# old.name.vec <- old.name.vec[old.name.vec %in% current.cell.types.all]
	
	# if(length(old.name.vec)==0) return(bp.obj) # nothing to replace
	# update_name_each <- function(input){
		# # cell types are at the second dimension for both theta and Znkg 
		# dimnames(input)[[2]][match(old.name.vec , dimnames(input)[[2]])] <- new.name.vec
		# return(input)
	# }
	
	# bp.obj$res$first.gibbs.res$theta.merged <- update_name_each(bp.obj$res$first.gibbs.res$theta.merged)
	# bp.obj$res$first.gibbs.res$Znkg.merged <- update_name_each(bp.obj$res$first.gibbs.res$Znkg.merged)
	# bp.obj$res$final.gibbs.theta <- update_name_each(bp.obj$res$final.gibbs.theta)
	
	# bp.obj
# }













