createModuleSummary <- function(ObjectList, modmeth="VBSR", numclus=1, supertype="refined") {
	alllabels = ObjectList$datasets[[1]]$responder[ObjectList$datasets[[1]]$keepsamps]
	samps2pheno = alllabels
	samps2pheno[which(alllabels == ObjectList$phenotype_class_vals_label[2])] = ObjectList$phenotype_class_vals[2]
	samps2pheno[which(alllabels == ObjectList$phenotype_class_vals_label[1])] = ObjectList$phenotype_class_vals[1]

	nonrespond_idxs = names(samps2pheno)[which(samps2pheno == ObjectList$phenotype_class_vals[1])]
	responder_idxs = names(samps2pheno)[which(samps2pheno == ObjectList$phenotype_class_vals[2])]

	# include modmeth
	modsumm_name = paste0(ObjectList$outdir,"/supermod_rewiring/supermod.",modmeth,".",numclus,"/",supertype,"summ.rds")
	if (file.exists(modsumm_name)){
		modsumm <- readRDS(modsumm_name)
	} else {
        stop("Rewiring file does not exist, please run runrewiring() before creating module summary.")
	}

	# set up output html page
	codedir <- paste0(system.file("extdata",package="TraRe"),"/RewiringReport/")
	dir.create(file.path(paste0(ObjectList$outdir,"/module_summ.",modmeth,".",numclus,".",supertype,"/")), showWarnings = FALSE)
	htmlinfo <- create_index_page(outdir = paste0(ObjectList$outdir,"/module_summ.",modmeth,".",numclus,".",supertype), runtag = "", codedir = codedir)

	orderobj <- geneOrder(modsumm, ObjectList$datasets[[1]]$keepsamps, ObjectList$datasets[[1]]$keeplabels, ObjectList$datasets[[1]]$norm_expr_mat_keep)
	createLegendPlot(htmlinfo)

	# Different Sections in the html summary
	superModuleStatistics(orderobj$modregs, orderobj$modtargs, orderobj$mat, ObjectList$datasets[[1]]$keeplabels, htmlinfo)
	correlationOfModuleGene(supertype, orderobj$regorder, orderobj$targetorder, orderobj$mat, orderobj$cormats, ObjectList$keeplabels, htmlinfo, ObjectList$phenotype_class_vals)
	expressionTableOfModuleGenes(supertype, modsumm$nodesumm, htmlinfo)
	expressionPlotsOfModuleGenes(supertype, orderobj$regorder, orderobj$targetorder, orderobj$mat, samps2pheno, ObjectList$phenotype_class_vals, htmlinfo)
	bipartiteGraphsSumm(numclus, modsumm, modmeth, htmlinfo)
	nullDistributionOfRewiringStatistic(orderobj$mat, ObjectList$datasets[[1]]$keeplabels, modmeth, supertype, htmlinfo)
	rankdf <- violinPlots(ObjectList$datasets[[1]]$norm_expr_mat_keep, ObjectList$datasets[[1]]$keepsamps, ObjectList$datasets[[1]]$keeplabels, modsumm$nodesumm, modsumm$fulledgesumm, modsumm$appendmat, htmlinfo)
	regulatorSummaryAndRank(rankdf, htmlinfo)
}