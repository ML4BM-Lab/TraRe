createModuleSummary <- function(ObjectList,run_mode=2,cluster_num=1) {
	modmeths <- names(ObjectList$datasets[[cluster_num]]$rundata$modules)

	alllabels = ObjectList$datasets[[cluster_num]]$responder[ObjectList$datasets[[cluster_num]]$keepsamps]
	samps2pheno = alllabels
	samps2pheno[which(alllabels == ObjectList$phenotype_class_vals_label[2])] = ObjectList$phenotype_class_vals[2]
	samps2pheno[which(alllabels == ObjectList$phenotype_class_vals_label[1])] = ObjectList$phenotype_class_vals[1]

	nonrespond_idxs = names(samps2pheno)[which(samps2pheno == ObjectList$phenotype_class_vals[1])]
	responder_idxs = names(samps2pheno)[which(samps2pheno == ObjectList$phenotype_class_vals[2])]

	mymod <- ''
	modsumm <- NULL
	runmoddata <- NULL
	if (run_mode == 1){
	  mymod <- 'supermod'
	  modsumm <- readRDS(paste0(ObjectList$outdir,"/supermodule_",cluster_num,"/rawsumm.rds"))
	  runmoddata <- readRDS(paste0(ObjectList$outdir,"/supermodule_",cluster_num,"/rawrunmoddata.rds"))
	}else if (run_mode == 2) {
	  mymod <- 'refinedsupermod'
	  modsumm <- readRDS(paste0(ObjectList$outdir,"/supermodule_",cluster_num,"/refinedsumm.rds"))
	  runmoddata <- readRDS(paste0(ObjectList$outdir,"/supermodule_",cluster_num,"/refinedrunmoddata.rds"))
	} 

	# set up output html page
	codedir <- paste0(system.file("extdata",package="TraRe"),"/RewiringReport/")
	dir.create(file.path(paste0(ObjectList$outdir,"/mod_summ_html/")), showWarnings = FALSE)
	htmlinfo <- create_index_page(outdir = paste0(ObjectList$outdir,"/mod_summ_html"), runtag = "", codedir = codedir)
	# imgdir <- paste0(htmlinfo$htmldir, htmlinfo$imgstr)

	orderobj <- geneOrder(modsumm, ObjectList$datasets[[cluster_num]]$keepsamps, ObjectList$datasets[[cluster_num]]$keeplabels, ObjectList$datasets[[cluster_num]]$norm_expr_mat_keep)
	createLegendPlot(htmlinfo)

	# Different Sections in the html summary
	superModuleStatistics(orderobj$modregs, orderobj$modtargs, orderobj$mat, ObjectList$datasets[[cluster_num]]$keeplabels, htmlinfo)
	correlationOfModuleGene(mymod, orderobj$regorder, orderobj$targetorder, orderobj$mat, orderobj$cormats, ObjectList$keeplabels, htmlinfo, ObjectList$phenotype_class_vals)
	expressionTableOfModuleGenes(mymod, modsumm$nodesumm, htmlinfo)
	expressionPlotsOfModuleGenes(mymod, orderobj$regorder, orderobj$targetorder, orderobj$mat, samps2pheno, ObjectList$phenotype_class_vals, htmlinfo)
	bipartiteGraphsSumm(cluster_num, modsumm, runmoddata, ObjectList$datasets[[cluster_num]]$norm_expr_mat_keep, nonrespond_idxs, responder_idxs, modmeths[cluster_num], htmlinfo)
	nullDistributionOfRewiringStatistic(orderobj$mat, ObjectList$datasets[[cluster_num]]$keeplabels, modmeths[cluster_num], mymod, htmlinfo)
	rankdf <- violinPlots(ObjectList$datasets[[cluster_num]]$norm_expr_mat_keep, ObjectList$datasets[[cluster_num]]$keepsamps, ObjectList$datasets[[cluster_num]]$keeplabels, modsumm$nodesumm, modsumm$fulledgesumm, modsumm$appendmat, htmlinfo)
	regulatorSummaryAndRank(rankdf, htmlinfo)
}

