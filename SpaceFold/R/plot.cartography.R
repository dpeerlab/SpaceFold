#functions to project gene expression onto the spacefold axis



get.gene.idx <- function(sf.obj, selected.genes){
	#determine if selected.genes is symbol or ID, 
	#and return a vector gene.idx denoting the index of selected genes.
	#search by maximum overlapping column
	which.col <- which.max(apply(sf.obj@feature,2,
								 function(col.i) length(intersect(as.character(col.i), selected.genes))))
	gene.idx <- match(selected.genes, sf.obj@feature[, which.col]) 	
	if(sum(is.na(gene.idx))>0) stop(paste(c("selected.genes", selected.genes[is.na(gene.idx)], "were not found"), collapse=" "))
	gene.idx	
}



#' function that fits the (binned) expression value against the (binned) cordinate
#' @param exp.vec a numeric vector of expression value in selected spots of one gene in one cell type
#' @param cord.vec a numeric vector of spacefold cordinate in selected spots. same length as exp.vec
#' @param n.bins number of bins, default=20
#' @param bin.by method to bin the data, deviding them into each bin of each sample size (equal.size), 
		# or by the same distance over spacefold axis (equal.step)
#' @param span parameter used by loess curve fitting, default=0.75
fit.exp.curv <- function(exp.vec, 
						 cord.vec, 
						 n.bins=20,
						 bin.by="equal.size",
						 span = 0.75){
	
	quantile.cut.points <- seq(0,1,length.out= n.bins+1 )

	if (bin.by=="equal.size"){
		breaks= quantile(cord.vec,prob= quantile.cut.points)
	}
	if (bin.by=="equal.step"){
		breaks= quantile(cord.vec,prob= quantile.cut.points) 
	}
	
	cord.cuts <- cut(cord.vec,
							breaks= breaks,
							labels=paste("cord",quantile.cut.points[-1],sep="-"),include.lowest=TRUE)
	exp.bin.mean <- by(exp.vec, INDICES= cord.cuts, mean)
	exp.bin.se <- by(exp.vec, INDICES= cord.cuts, function(i) sd(i)/sqrt(length(i)) ) 
	exp.bin.mean.up <- exp.bin.mean +2* exp.bin.se
	exp.bin.mean.low <- exp.bin.mean -2* exp.bin.se

	cord.bin.mid <- by(cord.vec, INDICES= cord.cuts, mean)	
	
	#spline fit
	grid.point <- seq(min(cord.bin.mid,na.rm=T),max(cord.bin.mid,na.rm=T),length.out=1000)

	loess_fit.mean <- predict(loess(exp.bin.mean ~ cord.bin.mid, span= span,na.action="na.omit"), grid.point)
	loess_fit.up <- predict(loess(exp.bin.mean.up ~ cord.bin.mid, span= span,na.action="na.omit"), grid.point)
	loess_fit.low <- predict(loess(exp.bin.mean.low ~ cord.bin.mid, span= span,na.action="na.omit"), grid.point)
	
	return(list(exp.bin.mean = exp.bin.mean,
			    cord.bin.mid = cord.bin.mid, 
			    loess_fit.mean= loess_fit.mean, 
			    loess_fit.up= loess_fit.up, 
			    loess_fit.low= loess_fit.low,
			    grid.point= grid.point,
			    exp.vec= exp.vec,
			    cord.vec= cord.vec))
}





get.ylim <- function(dat.list,
					 gene.idx,
					 show.raw){
	
	dat.list.gene <- dat.list[[gene.idx]]
	
	fit.dat.vec <- c( unlist(lapply(dat.list.gene, '[[', 'loess_fit.up')),
				  unlist(lapply(dat.list.gene, '[[', 'loess_fit.low'))	 
				 )
					 	
	if(show.raw) range(c(fit.dat.vec, unlist(lapply(dat.list.gene, '[[', 'exp.vec'))))
	else range(fit.dat.vec)				 	
	
}



#' function that plots expression curve
#' @param ct.name character variable of cell type name
#' @param gene.name character variable of gene name
#' @param show.raw logical variable of whether overlay the unbinned Znkg value
#' @param color the color of the fitted curves
#' @param color.raw the color of the raw points
#' @param pch.raw pch of raw points
#' @param ylim limit of y axis, default=NULL (automatically determined)
#' @param xlim limit of x axis, dafault=c(0,1),
#' @param ... other arguments passed to plot 
plot.exp.curv <- function(exp.curv,
						  ct.name=NULL, 
						  gene.name=NULL,
						  show.raw=FALSE,
						  color="black",
						  color.raw= adjustcolor("black",0.75),
						  cex.raw=0.7,
						  pch.raw=16,
						  ylim=NULL,
						  xlim=c(0,1),
						  xlab=NULL,
						  ylab=NULL,
						  axis.pos="left", #c("left","right")
						  ...){
	
	if(is.null(ylim)) ylim <- range(exp.curv$exp.bin.mean)

	if(axis.pos=="left"){
		plot(exp.curv$exp.bin.mean ~ exp.curv$cord.bin.mid,
		 	type="n",
		 	xlim= xlim, 
		 	ylim= ylim,
		 	bty="n",
		 	main=paste(gene.name, "in", ct.name, sep=" "),
		 	xlab=xlab,
		 	ylab= ylab,
		 	...)
	}
	else{
		plot(exp.curv$exp.bin.mean ~ exp.curv$cord.bin.mid,
		 type="n",
		 xlim= xlim, 
		 ylim= ylim,
		 bty="n",
		 main=paste(gene.name, "in", ct.name, sep=" "),
		 xlab=xlab,
		 ylab="",
		 axes=FALSE,
		 ...)
		 mtext(ylab,side=4,line=4,cex=2) 
		axis(4, ylim= ylim)
	}

	
	polygon(c(exp.curv$grid.point,rev(exp.curv$grid.point)),
		c(exp.curv$loess_fit.low,rev(exp.curv$loess_fit.up)),
		border=NA,col=adjustcolor(color,0.25))
	lines(exp.curv$grid.point, exp.curv$loess_fit.mean,lwd=1,col= color)
	points(x= exp.curv$cord.bin.mid, y=exp.curv$exp.bin.mean,pch=16,col= color)
	if(show.raw) points(x= exp.curv$cord.vec, y=exp.curv$exp.vec,pch= pch.raw,col= color.raw,cex= cex.raw)
	NULL
}



get.pdf.name <- function(pdf.idx){
	if(pdf.idx>= 10 & pdf.idx < 100) pdf.name <- paste("tmp.0", pdf.idx,".pdf",sep="")
	if(pdf.idx < 10) pdf.name <- paste("tmp.00", pdf.idx,".pdf",sep="")
	if(pdf.idx >= 100) pdf.name <- paste("tmp.", pdf.idx,".pdf",sep="")
	pdf.name
}




#visualize gene expression level along the projected axis
# use this function only in pdfjam is installed. Otherwise, use plot.cartography.nopdfjam.
#' @param sf.obj a BayesPrism output object
#' @param raw.or.norm plotting raw value or value normalized by size factor
#' @param selected.genes a character vector specifiying the gene name / symbols to be plotted, which is automatically determined by matching columns of  sf.obj$para$feature 
#' @param selected.cell.types a character vector specifiying the cell types / cell type groups to be plotted, which is automatically determined by matching colnames of sf.obj$res$first.gibbs.res$theta.merged and sf.obj$res$res.regrouped$theta0
#' @param bin.by method to bin the data, deviding them into each bin of each sample size (equal.size), 
		# or by the same distance over spacefold axis (equal.step)
#' @param n.bins number of bins, default=20
#' @param span parameter used by loess curve fitting, default=0.75
#' @param show.raw whether plot the unbinned Znkg/Znkg.norm
#' @param overlay.hist whether plot the histogram of the cordiantes containing each cell type separately or on top of the expression curve.
#' @param col.hist color of the histogram
#' @param pdf.prefix the prefix of the pdf output
plot.cartography.pdfjam <- function(sf.obj,
						raw.or.norm=c("raw","norm"),
						denoise =TRUE,
						selected.genes,
						selected.cell.types,
						bin.by=c("equal.size","equal.step"),
						n.bins=20,
						span=0.75,
						show.raw =FALSE,
						overlay.hist=TRUE,
						col.hist=adjustcolor("grey",0.5),
						pdf.prefix,
						return.raw,
						...){
	#assertions...
	if(denoise) plot.dat <- get.plot.dat.denoised  (sf.obj, raw.or.norm, selected.genes, selected.cell.types, span)
	else plot.dat <- get.plot.dat  (sf.obj, raw.or.norm, selected.genes, selected.cell.types, bin.by, n.bins, span)
	
	color.num <- length(selected.cell.types)
	if(color.num<=8) my.palette <- brewer.pal(color.num,"Dark2")
	if(color.num>8 && color.num<=12) my.palette <- brewer.pal(color.num,"Paired")
	if(color.num > 12)  my.palette <-  colorRampPalette(brewer.pal(12, "Paired"))(color.num)
	

	ylim.all <- lapply(selected.genes, FUN= get.ylim, dat.list= plot.dat, show.raw= show.raw)
	
	
	pdf.idx <- 1
	
	if(!overlay.hist) {
		for (ct.idx in 1:length(selected.cell.types)){
			pdf(get.pdf.name(pdf.idx),useDingbats=F,pointsize=8)
			par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
			par(mar=c(5.1, 5.1, 4.1, 4.1))
			hist(plot.dat[[1]][[ct.idx]]$cord.vec,
				main= selected.cell.types[ct.idx], 
				xlim= c(0,1), 
				breaks=50, 
				col= adjustcolor(my.palette[ct.idx],0.5),
				xlab="SpaceFold cordinate",
				ylab="Cell Type Frequency")
			dev.off()
			pdf.idx <- pdf.idx + 1
			
			for (gene.idx in 1:length(selected.genes)){
				pdf(get.pdf.name(pdf.idx),useDingbats=F,pointsize=8)
				par(mar=c(5.1, 5.1, 4.1, 4.1))
				par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
				plot.exp.curv (plot.dat[[gene.idx]][[ct.idx]],
						  ct.name= selected.cell.types[ct.idx], 
						  gene.name= selected.genes[gene.idx],
						  show.raw= show.raw,
						  ylim= ylim.all[[gene.idx]],
						  color= my.palette[ct.idx],
						  xlab="SpaceFold cordinate",
						  ylab=paste(raw.or.norm, "expression value"),
						  axis.pos= "left")
				dev.off()
				pdf.idx <- pdf.idx + 1
			}
		}
	}
	else{
		for (ct.idx in 1:length(selected.cell.types)){
			for (gene.idx in 1:length(selected.genes)){
				pdf(get.pdf.name(pdf.idx),useDingbats=F,pointsize=8)
				par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
				par(mar=c(6.1, 5.1, 4.1, 6.1))
				hist(plot.dat[[1]][[ct.idx]]$cord.vec,
				 	main= NULL, 
				 	xlim= c(0,1), 
					breaks=50, 
					col= col.hist,
					xlab="SpaceFold cordinate",
					ylab="Cell Type Frequency")
				 par(new=TRUE) 
				 
				 #par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
				 #par(mar=c(5.1, 5.1, 4.1, 2.1))
				 plot.exp.curv (plot.dat[[gene.idx]][[ct.idx]],
						  	ct.name= selected.cell.types[ct.idx], 
						  	gene.name= selected.genes[gene.idx],
						  	show.raw= show.raw,
						  	ylim= ylim.all[[gene.idx]],
						  	color= my.palette[ct.idx],
						  	xlab="SpaceFold cordinate",
						  	ylab=paste(raw.or.norm, "expression value"),
						  	axis.pos= "right")
				 dev.off()
				 pdf.idx <- pdf.idx + 1
			}
		}
	}	

	#merge pdfs
	ncol=length(selected.genes)
	nrow = length(selected.cell.types)
	if(!overlay.hist) ncol <- ncol+1
	system(paste("pdfjam", " ./tmp*.pdf ",  "--nup ", ncol, "x",  nrow , " --landscape  --outfile ", pdf.prefix,".pdf",sep="")) # needs to install pdfjam
	system("rm ./tmp*.pdf")	

	if(return.raw) return(plot.dat)
}


#visualize gene expression level along the projected axis
# make individual plot and output in a single pdf
#' @param sf.obj a BayesPrism output object
#' @param raw.or.norm plotting raw value or value normalized by size factor
#' @param selected.genes a character vector specifiying the gene name / symbols to be plotted, which is automatically determined by matching columns of  sf.obj$para$feature 
#' @param selected.cell.types a character vector specifiying the cell types / cell type groups to be plotted, which is automatically determined by matching colnames of sf.obj$res$first.gibbs.res$theta.merged and sf.obj$res$res.regrouped$theta0
#' @param bin.by method to bin the data, deviding them into each bin of each sample size (equal.size), 
		# or by the same distance over spacefold axis (equal.step)
#' @param n.bins number of bins, default=20
#' @param span parameter used by loess curve fitting, default=0.75
#' @param show.raw whether plot the unbinned Znkg/Znkg.norm
#' @param overlay.hist whether plot the histogram of the cordiantes containing each cell type separately or on top of the expression curve.
#' @param col.hist color of the histogram
#' @param pdf.prefix the prefix of the pdf output
plot.cartography.nopdfjam <- function(sf.obj,
							 raw.or.norm=c("raw","norm"),
							 denoise=TRUE,
							 selected.genes,
							 selected.cell.types,
							 bin.by=c("equal.size","equal.step"),
							 n.bins=20,
							 span=0.75,
							 show.raw =FALSE,
							 overlay.hist=TRUE,
							 col.hist=adjustcolor("grey",0.5),
							 pdf.prefix,
							 return.raw,
							 ...){
	#assertions...
	if(denoise) plot.dat <- get.plot.dat.denoised  (sf.obj, raw.or.norm, selected.genes, selected.cell.types, span)
	else plot.dat <- get.plot.dat  (sf.obj, raw.or.norm, selected.genes, selected.cell.types, bin.by, n.bins, span)
	
	color.num <- length(selected.cell.types)
	if(color.num<=8) my.palette <- brewer.pal(color.num,"Dark2")
	if(color.num>8 && color.num<=12) my.palette <- brewer.pal(color.num,"Paired")
	if(color.num > 12)  my.palette <-  colorRampPalette(brewer.pal(12, "Paired"))(color.num)
	
	pdf(paste(pdf.prefix,".pdf",sep=""),useDingbats=F,pointsize=8)

	if(!overlay.hist){
		for (ct.idx in 1:length(selected.cell.types)){	
			par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
			par(mar=c(5.1, 5.1, 4.1, 2.1))
			hist(plot.dat[[1]][[ct.idx]]$cord.vec,
				 main= selected.cell.types[ct.idx], 
				 xlim= c(0,1), 
				 breaks=50, 
				 col= adjustcolor(my.palette[ct.idx],0.5),
				 xlab="SpaceFold cordinate",
				 ylab="Cell Type Frequency")
		}
		axis.pos <- "left"
	}
	
	#making plots
	for (gene.idx in 1:length(selected.genes)){
		#get the range
		ylim.gene <- get.ylim (plot.dat, gene.idx, show.raw)
		for (ct.idx in 1:length(selected.cell.types)){
			par(mar=c(5.1, 5.1, 4.1, 5.1))
			par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
				
			if(overlay.hist){
				hist(plot.dat[[1]][[ct.idx]]$cord.vec,
				 main= "", 
				 xlim= c(0,1), 
				 breaks=50, 
				 col= col.hist,
				 xlab="SpaceFold cordinate",
				 ylab="Cell Type Frequency")
				 par(new=TRUE) 
			}
			
			par(mar=c(5.1, 5.1, 4.1, 5.1))
			par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
			plot.exp.curv (plot.dat[[gene.idx]][[ct.idx]],
						  ct.name= selected.cell.types[ct.idx], 
						  gene.name= selected.genes[gene.idx],
						  show.raw= show.raw,
						  ylim= ylim.gene,
						  color= my.palette[ct.idx],
						  xlab="SpaceFold cordinate",
						  ylab=paste(raw.or.norm, "expression value"),
						  axis.pos= "right")
		}
	}
	
	dev.off()
	
	if(return.raw) return(plot.dat)
}



#visualize gene expression level along the projected axis
# If pdfjam is installed use it. Otherwise,  make individual plot and output in a single pdf.
#' @param sf.obj a BayesPrism output object
#' @param raw.or.norm plotting raw value or value normalized by size factor
#' @param selected.genes a character vector specifiying the gene name / symbols to be plotted, which is automatically determined by matching columns of  sf.obj$para$feature 
#' @param selected.cell.types a character vector specifiying the cell types / cell type groups to be plotted, which is automatically determined by matching colnames of sf.obj$res$first.gibbs.res$theta.merged and sf.obj$res$res.regrouped$theta0
#' @param bin.by method to bin the data, deviding them into each bin of each sample size (equal.size), 
		# or by the same distance over spacefold axis (equal.step)
#' @param n.bins number of bins, default=20
#' @param span parameter used by loess curve fitting, default=0.75
#' @param show.raw whether plot the unbinned Znkg/Znkg.norm
#' @param overlay.hist whether plot the histogram of the cordiantes containing each cell type separately or on top of the expression curve.
#' @param col.hist color of the histogram
#' @param pdf.prefix the prefix of the pdf output
plot.cartography <- function(sf.obj,
						raw.or.norm=c("raw","norm"),
						denoise =TRUE,
						selected.genes,
						selected.cell.types,
						bin.by=c("equal.size","equal.step"),
						n.bins=20,
						span=0.75,
						show.raw =FALSE,
						overlay.hist=TRUE,
						col.hist=adjustcolor("grey",0.5),
						pdf.prefix,
						return.raw=FALSE,
						...){
	#check if pdfjam is installed
	check.pdfjam <- tryCatch( system("pdfjam -V"), error=function(e) {return (-1)},
											warning=function(w)return(-1))
		
	if(check.pdfjam != -1){
		plot.cartography.pdfjam (sf.obj,
						raw.or.norm,
						denoise,
						selected.genes,
						selected.cell.types,
						bin.by,
						n.bins,
						span,
						show.raw,
						overlay.hist,
						col.hist,
						pdf.prefix,
						return.raw,
						...)
	}
	else{
		plot.cartography.nopdfjam (sf.obj,
						raw.or.norm,
						denoise,
						selected.genes,
						selected.cell.types,
						bin.by,
						n.bins,
						span,
						show.raw,
						overlay.hist,
						col.hist,
						pdf.prefix, 
						return.raw,
						...)
	}
}



