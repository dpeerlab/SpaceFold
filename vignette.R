library(SpaceFold)


load("/data/peer/tinyi/RU_data/visium/dat/bp_dat/from_reprocessed_anndata/SI/merge_finer_states_061122_fullgene_package/SI.bp.sig.rdata")


gene.tab <- read.table("/lila/data/peer/tinyi/RU_data/visium/dat/visium_dat/SI/WTSI2_all_reads_D1_outs/raw_feature_bc_matrix/features.tsv.gz",sep="\t",check.names=F,header=F)[,1:2]
colnames(gene.tab) <- c("id",'symbol')
rownames(gene.tab) <- gene.tab[,"id"]


sf.obj <- new.sf(bp.obj = SI.bp, feature= gene.tab)

sf.obj <- compute.background.level(sf.obj, theta.cutoffs.user=0.001, Znk.cutoffs.user=20)

sf.obj <- run.phate(sf.obj, if.invert=TRUE)


myorder <- c("myofibroblast","glial","stromal 3","stromal 1/2","blood endothelium","lymphatic",
			 "Paneth","Lgr5+ stem","Lgr5+ progenitor","TA","bottom zone enterocyte","mid zone enterocyte","top zone enterocyte",
			 "secretory progenitor","goblet cycling","goblet 1","goblet 2","tuft","neurons/enteroendocrine",
			 "resting B cell","cycling/GC B cell","plasma","T cell","cDC/monocyte","macrophage","pDC")

myorder <- c("myofibroblast","glial","stromal 3","stromal 1/2","blood endothelium","lymphatic",
			 "paneth","Lgr5+ stem","Lgr5+ progenitor","TA","bottom zone enterocyte","mid zone enterocyte","top zone enterocyte",
			 "secretory progenitor","goblet cycling","goblet 1","goblet 2","tuft","neurons/enteroendocrine",
			 "resting B cell","cycling/GC B cell","plasma","T cell","cDC/monocyte","macrophage","pDC")

plot.beeswarm(sf.obj, cell.type.order= myorder)




group.list <- list(enterocytes=c("bottom zone enterocyte","mid zone enterocyte","top zone enterocyte"),
					goblet=c("goblet 1","goblet 2","secretory progenitor", "goblet cycling"))


sf.obj <- merge.cell.type (sf.obj= sf.obj,
					   grouping.list = group.list)


sf.obj <- denoise.cartography (sf.obj= sf.obj, ka=5, tansition=1000)



secreted_factors=c('Ntn1','Il33','Wnt2','Ccl21a','Rspo3','Reln')
plot.cartography (sf.obj= sf.obj,
							 raw.or.norm="norm",
							 denoise=T,
							 selected.genes= secreted_factors,
							 selected.cell.types=c('lymphatic','blood endothelium'),
							 show.raw =FALSE,
							 overlay.hist=TRUE,
							 span=0.25,
							 pdf.prefix="SI.secreted_factors.v2")

plot.cartography (sf.obj= sf.obj,
							 raw.or.norm="norm",
							 denoise=F,
							 bin.by="equal.size",
							 selected.genes= secreted_factors,
							 selected.cell.types=c('lymphatic','blood endothelium'),
							 show.raw =FALSE,
							 overlay.hist=TRUE,
							 span=0.75,
							 pdf.prefix="SI.secreted_factors.v1")


