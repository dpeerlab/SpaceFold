#functions to project gene expression onto the spacefold axis


#set background level using user provided values
#' @param theta.cut cutoff for theta of background spots, a numerical number or a named numerical vector of the same length as ncol(theta).
#' @param Znk.cut cutoff for Znk of background spots, a numerical number or a named numerical vector of the same length as ncol(theta).
set.background.level <-  function(sf.obj, 
								  theta.cutoffs, 
								  Znk.cutoffs){
		
	Znk <- sf.obj@data@Znk
	theta <- sf.obj@data@theta
	
	if(is.null(Znk)) 
		sf.obj <- compute.Znk(sf.obj)
	
	if(length(theta.cutoffs)==1) {
		theta.cutoffs <- rep(theta.cutoffs, ncol(theta))
		names(theta.cutoffs) <- colnames(theta)
	}
	if(length(Znk.cutoffs)==1) {
		Znk.cutoffs <- rep(Znk.cutoffs, ncol(theta))
		names(Znk.cutoffs) <- colnames(theta)
	}

	
	stopifnot(identical(names(theta.cutoffs),colnames(theta)))
	stopifnot(identical(names(Znk.cutoffs),colnames(Znk)))
	
	
	selected.spot.matrix <- do.call(cbind, 
									lapply(1:ncol(theta), 
											function(col.idx) {
														theta[, col.idx] > theta.cutoffs[col.idx] & 
														Znk[, col.idx] > Znk.cutoffs[col.idx] }))
	colnames(selected.spot.matrix) <- colnames(theta)
	
	sf.obj@data@selected.spot.matrix <- selected.spot.matrix
	
	sf.obj@control_param$theta.cutoffs <- theta.cutoffs
	sf.obj@control_param$Znk.cutoffs <- Znk.cutoffs
	
	sf.obj
}




#set background level by fitting a mixture model on theta and Znk (for each cell type)
#' @param sf.obj a BayesPrism output object
#' @param which.theta use first(initial) or final(updated) theta 
compute.background.level <-  function(sf.obj,
									  theta.cutoffs.user=NULL,
									  Znk.cutoffs.user=NULL){
		
	Znk <- sf.obj@data@Znk
	theta <- sf.obj@data@theta

	if(is.null(Znk)) 
		sf.obj <- compute.Znk(sf.obj)
	Znk <- sf.obj@data@Znk
	
	theta.cutoffs <-c()
	print("fitting mixture models on theta...")
	for(i in 1:ncol(theta)){
		print(colnames(theta)[i])

		fit=gammamixEM(theta[,i], k=2,maxit=10000,maxrestarts=100)
		cls <- apply(fit$posterior,1,which.max)
		if(length(unique(cls))==1) {
			print("refit")
			fit=gammamixEM(theta[,i], k=2,maxit=10000,maxrestarts=100, mom.start=F)
			cls <- apply(fit$posterior,1,which.max)
		}
		print(table(cls))
		theta.cutoffs <- c(theta.cutoffs, median(c(range(theta[cls==1,i]),range(theta[cls==2,i]))))
	}
	names(theta.cutoffs) <- colnames(theta)

	if(!is.null(theta.cutoffs.user)) theta.cutoffs[theta.cutoffs<theta.cutoffs.user] <- theta.cutoffs.user

	Znk.cutoffs <-c()
	print("fitting mixture models on Znk...")
	for(i in 1:ncol(Znk)){
		print(colnames(Znk)[i])
	
		fit= Mclust(Znk[,i], G=2)
		cls <- fit$classification
		print(table(cls))
		Znk.cutoffs <- c(Znk.cutoffs, median(c(range(Znk[cls==1,i]),range(Znk[cls==2,i]))))
	}
	names(Znk.cutoffs) <- colnames(Znk)

	if(!is.null(Znk.cutoffs)) Znk.cutoffs[Znk.cutoffs<Znk.cutoffs.user] <- Znk.cutoffs.user

	set.background.level (sf.obj = sf.obj, 
						  theta.cutoffs = theta.cutoffs,
						  Znk.cutoffs = Znk.cutoffs)	

}




