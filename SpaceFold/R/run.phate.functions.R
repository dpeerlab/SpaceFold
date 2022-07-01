

#' functions to run PHATE on BayesPrism output
#' default parameters were those used by the spacefold paper
#' @param bp.obj: a (merged) BayesPrism output  
#' @param anchorCellTypes: a character vector denoting the cell types selected for computing embeddings. 
#'			only cell types show common the spatial pattern across all bp.obj in the bp.list should be selected.
#'			anchorCellTypes should be informative about the physical cordinate of the spot, 
#'			i.e., cell types that are known or suspected to vary with the cordinate.
#'			Spacefold is usually robust to selection of cell types.
#'			#default="all', i.e. uses all cell types
#' @param renorm.to.one, nomalize the seleceted cell type fractions to sum to one. Default TRUE.
#' @param center, scale the seleceted cell type fractions across spatial spots to mean=0. Default TRUE.
#' @param scale, scale the seleceted cell type fractions across spatial spots to sd=1. Default TRUE.
#' @param if.pseudo.axis, return a scaled cordinate with min=zero and max=one.
#' @param if.invert whehter flip the spacefold cordinate. Default is FALSE 
#' @param mds.solver solver used by PHATE
#' @param n.jobs number of threads used by PHATE
#' @param knn the knn paramter used by PHATE
#' @param ... additional paramters passed to PHATE
run.phate <- function(sf.obj,
					 anchorCellTypes ="all",
					 renorm.to.one = TRUE,
					 center=TRUE,
					 scale=TRUE,
					 if.pseudo.axis=TRUE,
					 if.uniform=FALSE,
					 if.invert=FALSE,
					 mds.solver="smacof",
					 n.jobs=1,
					 knn=10,
					 ...){
	
	#other assertion arugments to be written...
	
	theta <- sf.obj@data@theta
	
	if(length(anchorCellTypes)==1 & anchorCellTypes=="all") anchorCellTypes <- colnames(theta) 
	
	#subset using relevant cell types
	theta <- theta[, anchorCellTypes]
	
	if(renorm.to.one) theta <- theta / rowSums(theta)
	
	theta <- scale(theta, center = center, scale = scale)
	
	#can add addtional PCA steps...
	
	phate.res <- phate(data= theta, 
					   ndim= 1, 
					   mds.solver= mds.solver,
					   n.jobs= n.jobs,
					   knn= knn,
					   ...)
	
	cord <- phate.res$embedding
	
	if(if.invert) cord <- -cord
	
	if(if.pseudo.axis) {
		cord <- apply(cord,2,function(cord.i){
								cord.i <- cord.i - min(cord.i)
								cord.i <- cord.i / max(cord.i)
								cord.i
					})
	}
	
	if(if.uniform)
		cord <- rank(cord) / length(cord)
	
	#store spacefold parameters
	sf.obj@control_param <- list(anchorCellTypes= anchorCellTypes,
									renorm.to.one= renorm.to.one,
									center= center,
									scale= scale,
									if.pseudo.axis= if.pseudo.axis,
									if.uniform= if.uniform,
									ndim=1,
					 				mds.solver= mds.solver,
					 				n.jobs= n.jobs,
					 				knn= knn,
					 				list(...))
	
	sf.obj@SpaceFold.axis <- as.matrix(cord)
		
	sf.obj
}





#' function to plot beeswarm
#' @param palette a function that generates the palette 
#' @param q.cut upper quantile of spots to plot, default is 0.95 
#'			(=spots containing the top 5% of each cell type) 
#' @param pdf.prefix the prefix of pdf filename.
#' @param cell.type.order a character vector speficying the order of plotting. 
#'			Default is NULL, which uses the colnames of theta matrix
plot.beeswarm <- function(sf.obj,
						 palette=colorRampPalette(brewer.pal(12, "Paired")),
						 q.cut=0.95,
						 use.background=FALSE,
						 pdf.prefix="output",
						 height=5,
						 width=12,
						 mar=c(8,3,1,1),
						 cell.type.order=NULL){
	
	theta <- sf.obj@data@theta
	Znk <- sf.obj@data@Znk
	
	cord <- as.vector(sf.obj@SpaceFold.axis[,1])
	stopifnot(!is.null(cord))
		
	if(use.background) {
		background <- sf.obj@background
		if(is.null(background$theta.cutoffs) | is.null(background$Znk.cutoff))
			stop("please set background!")
	}
									
	
	my.palette <-  palette(ncol(theta))

	cell.type.uniq <-colnames(theta)

	plot.df <- do.call(rbind.data.frame, 
							lapply(1:length(cell.type.uniq), function(idx){
	
							ct.idx <- cell.type.uniq[idx]
							theta.vec <- theta[, ct.idx]
							Znk.vec <- Znk[, ct.idx]
							
							theta.vec.high <- theta.vec
							Znk.vec.high <- Znk.vec
							if(use.background) {
								theta.vec.high <- theta.vec[background$selected.spot.matrix[, ct.idx]]
								Znk.vec.high <- Znk.vec[background$selected.spot.matrix[, ct.idx]]
							}
							
							theta.cutoff <- quantile(theta.vec.high, prob= q.cut)
							Znk.cutoff <- quantile(Znk.vec.high,prob= q.cut)
							
							data.frame(layout="PHATE",
									   cord= cord, 
									   cell.type = ct.idx,
									   theta = theta.vec,           
									   if.above.background =  theta.vec > theta.cutoff & Znk.vec > Znk.cutoff,
									   col.solid = my.palette[idx],
									   col.theta = sapply(  theta.vec, adjustcolor, col=my.palette[idx]),
									   stringsAsFactors=FALSE)
	}))

	if(!is.null(cell.type.order)) plot.df $cell.type <- factor(plot.df $cell.type,levels= cell.type.order)

	#dev.new(height= height, width= width, file=paste(pdf.prefix,"_spacefold_beeswarm.pdf",sep=""))
	pdf(height= height, width= width, file=paste(pdf.prefix,"_spacefold_beeswarm.pdf",sep=""))
	par(mar= mar, xpd=TRUE)

	beeswarm(cord ~ cell.type, data = plot.df[plot.df $if.above.background,],
         pch = 16, pwcol = as.character(col.solid),
         xlab = "",
         ylab = "projected coordinate,colored by selected theta",
         cex=0.5,las=2,cex.axis=0.5,
         corral="wrap",
         bty="n",
         main=paste("theta above", q.cut, "quantile"))

	beeswarm(cord ~ cell.type, data = plot.df,
         pch = 16, pwcol = as.character(col.theta),
         xlab = "",
         ylab = "projected coordinate",
         cex=0.5,las=2,cex.axis=0.5,
         corral="wrap",
         bty="n",
         main=paste("color by theta"))
         
	dev.off()

}




