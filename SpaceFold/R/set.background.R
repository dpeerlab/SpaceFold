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
									  posterior.cutoff=0.7,
									  theta.cutoffs.user=NULL,
									  Znk.cutoffs.user=NULL,
									  seed=NULL){
	
	if(is.null(seed)) set.seed(seed)
		
	Znk <- sf.obj@data@Znk
	theta <- sf.obj@data@theta

	if(is.null(Znk)) 
		sf.obj <- compute.Znk(sf.obj)
	Znk <- sf.obj@data@Znk
	
	theta.cutoffs <-c()
	cat("fitting mixture models on theta... \n")
	cat("Current cell type: ")
	for(i in 1:ncol(theta)){
		cat(colnames(theta)[i], " ")

		capture.output({fit=gammamixEM(theta[,i], k=2,maxit=10000,maxrestarts=100, mom.start=F)})		
		cls <- apply(fit$posterior,1,which.max)
		posterior <- fit$posterior
		
		#detemin which cluster has higher mean
		cls.mean <- fit$gamma.pars["alpha",] * fit$gamma.pars["beta",]
		
		if(cls.mean[1] > cls.mean[2]) max.cls <- 1
		else max.cls <- 2
		
		high.cls.idx <- posterior[,max.cls] > posterior.cutoff
		
		#loop until there are non-zero spots with high posterior in the cluster with higher mean
		while(length(high.cls.idx)==0){
			posterior.cutoff <- posterior.cutoff - 0.05
			high.cls.idx <- posterior[,max.cls] > posterior.cutoff
		}
		
		theta.cutoff.i <- min(theta[high.cls.idx,i])
		
		theta.cutoffs <- c(theta.cutoffs, theta.cutoff.i)
	}
	names(theta.cutoffs) <- colnames(theta)

	cat("\n")

	if(!is.null(theta.cutoffs.user)){
		if(length(theta.cutoffs.user)==1) {
			theta.cutoffs.user <- rep(theta.cutoffs.user, length(theta.cutoffs))
			names(theta.cutoffs.user) <- names(theta.cutoffs)
		}
		
		stopifnot(!is.null(names(theta.cutoffs.user)))
		
		auto.cutoff.selected <- theta.cutoffs[names(theta.cutoffs.user)]
		
		user.higher.idx <- names(auto.cutoff.selected)[auto.cutoff.selected < theta.cutoffs.user]
		
		if(length(user.higher.idx)>0)
			theta.cutoffs[user.higher.idx] <- theta.cutoffs.user[user.higher.idx]
	} 


	Znk.cutoffs <-c()
	cat("fitting mixture models on Znk... \n")
	cat("Current cell type: ")
	for(i in 1:ncol(Znk)){
		cat(colnames(Znk)[i], " ")
		
		capture.output({fit <- Mclust(Znk[,i], G=2)})
		posterior <- fit$z
		
		#detemin which cluster has higher mean
		cls.mean <- fit$parameters$mean
		
		if(cls.mean[1] > cls.mean[2]) max.cls <- 1
		else max.cls <- 2
		
		high.cls.idx <- posterior[,max.cls] > posterior.cutoff
		
		#loop until there are non-zero spots with high posterior in the cluster with higher mean
		while(length(high.cls.idx)==0){
			posterior.cutoff <- posterior.cutoff - 0.05
			high.cls.idx <- posterior[,max.cls] > posterior.cutoff
		}
		
		Znk.cutoff.i <- min(Znk[high.cls.idx,i])
		
		Znk.cutoffs <- c(Znk.cutoffs, Znk.cutoff.i)
	}
	names(Znk.cutoffs) <- colnames(Znk)
	cat("\n")	
	
	if(!is.null(Znk.cutoffs.user)){
		if(length(Znk.cutoffs.user)==1) {
			Znk.cutoffs.user <- rep(Znk.cutoffs.user, length(Znk.cutoffs))
			names(Znk.cutoffs.user) <- names(Znk.cutoffs)
		}
		
		stopifnot(!is.null(names(Znk.cutoffs.user)))
		
		auto.cutoff.selected <- Znk.cutoffs[names(Znk.cutoffs.user)]
		
		user.higher.idx <- names(auto.cutoff.selected)[auto.cutoff.selected < Znk.cutoffs.user]
		
		if(length(user.higher.idx)>0)
			Znk.cutoffs[user.higher.idx] <- Znk.cutoffs.user[user.higher.idx]
	} 

	sf.obj@control_param$posterior.cutoff <- posterior.cutoff

	set.background.level (sf.obj = sf.obj, 
						  theta.cutoffs = theta.cutoffs,
						  Znk.cutoffs = Znk.cutoffs)	

}


#selected.idx a logical vector

select.spot <- function(sf.obj,
						cell.type,
						selected.idx,
						op){

	old.selection <- sf.obj@data@selected.spot.matrix[, cell.type] 
	
	stopifnot(length(selected.idx)==length(old.selection))
	
	if(op=="and") 
		new.selection <- old.selection & selected.idx
	
	if(op=="or") 
		new.selection <- old.selection | selected.idx
	
	if(op=="new") 
		new.selection <- selected.idx
	
	sf.obj@data@selected.spot.matrix[, cell.type]  <- new.selection
	
	return(sf.obj)						
}







