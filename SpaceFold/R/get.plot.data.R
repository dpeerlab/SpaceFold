#functions to project gene expression onto the spacefold axis

#' function that obtains denoised expression data, cordinates, binned data, and fitted curves to plot (after applying the filter)
#' returns a list (over genes) of list (over cell types) of exp.curv object (outputed by fit.exp.curv)
get.plot.dat.denoised <- function(sf.obj,
								  raw.or.norm=c("raw","norm"),
								  selected.genes,
								  selected.cell.types,
								  span){
	
	gene.idx <- get.gene.idx(sf.obj = sf.obj, selected.genes= selected.genes)
	
	if(raw.or.norm=="raw") {
		exp.dat.raw <- sf.obj@data@Z
		exp.dat.denoised <- sf.obj@denoised.cartography@Z.denoised
	}
	if(raw.or.norm=="norm") {
		exp.dat.raw <- sf.obj@data@Z.normed
		exp.dat.denoised <- sf.obj@denoised.cartography@Z.normed.denoised
	}
	
	dat.list <- list()
	for(g in 1:length(selected.genes)){
		ct.list<- list() 
		for(k in 1:length(selected.cell.types)){
			spot.idx <- sf.obj@data@selected.spot.matrix[, selected.cell.types[k]]
			cord.vec <- sf.obj@SpaceFold.axis[spot.idx,1]			
			
			exp.vec.raw <- exp.dat.raw[spot.idx,  gene.idx[g], selected.cell.types[k]]	
			exp.vec.denoised <- exp.dat.denoised[[selected.cell.types[k]]][,  gene.idx[g]]
			
			
			cord.ordered <- cord.vec[order(cord.vec)]
			exp.ordered <- exp.vec.denoised[order(cord.vec)]
			
			if(length(cord.ordered)==0) {
				ct.list[[k]]<-NULL
				next
			}
			
			loess.fit <- loess.sd(y= exp.ordered, x= cord.ordered, nsigma = 2, span=span)
			
			ct.list[[k]] <- list(exp.vec = exp.vec.raw, 
				 				 cord.vec = cord.vec,
				 				 grid.point = cord.ordered,
				 				 exp.bin.mean = exp.vec.denoised,
				 				 cord.bin.mid = cord.vec,
				 				 loess_fit.mean =loess.fit$y,
				 				 loess_fit.up = loess.fit$upper, 
				 				 loess_fit.low = loess.fit$lower)
			
		}
		names(ct.list) <- selected.cell.types
		ct.list[sapply(ct.list, is.null)] <- NULL
		dat.list[[g]] <- ct.list
	} 
	names(dat.list) <- selected.genes
	
	dat.list
}



#' function that obtains raw expression data, cordinates, binned data, and fitted curves to plot (after applying the filter)
#' returns a list (over genes) of list (over cell types) of exp.curv object (outputed by fit.exp.curv)
get.plot.dat <- function(sf.obj,
						raw.or.norm=c("raw","norm"),
						selected.genes,
						selected.cell.types,
						bin.by,
						n.bins,
						span){
	
	gene.idx <- get.gene.idx(sf.obj = sf.obj, selected.genes= selected.genes)
	
	if(raw.or.norm=="raw") exp.dat <- sf.obj@data@Z
	if(raw.or.norm=="norm") exp.dat <- sf.obj@data@Z.normed
	
	dat.list <- list()
	for(g in 1:length(selected.genes)){
		ct.list<- list() 
		for(k in 1:length(selected.cell.types)){
			spot.idx <- sf.obj@data@selected.spot.matrix[, selected.cell.types[k]]
			exp.vec <- exp.dat[spot.idx,  gene.idx[g], selected.cell.types[k]]
			cord.vec <- sf.obj@SpaceFold.axis[spot.idx,1]
			
			if(length(cord.vec)==0) {
				ct.list[[k]]<-NULL
				next
			}
			
			ct.list[[k]] <- fit.exp.curv (exp.vec= exp.vec, 
						 								 cord.vec= cord.vec, 
						 								 n.bins= n.bins,
						 								 bin.by= bin.by,
						 								 span = span)
			
		}
		
		names(ct.list) <- selected.cell.types
		ct.list[sapply(ct.list, is.null)] <- NULL
		dat.list[[g]] <- ct.list
	} 
	names(dat.list) <- selected.genes
	
	dat.list
}






