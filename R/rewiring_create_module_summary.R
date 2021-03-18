#' Create module summary
#'
#' Create module summary for supermodule after the rewiring step.
#'
#' @param ObjectList Output from `preparerewiring()`containing some required parameters.
#' @param modmeth Method that uses to generate the supermodule (default:"VBSR").
#' @param numclus Integer specifying the cluster detected from 'runrewiring()' (default:1).
#' @param supertype Secification of the type of supermodule, "refined" or "raw", (default:"refined")
#' @param numdataset Integer specifying the dataset used, 1 when there is only one dataset (default:1).
#'
#' @return htmlinfo the index page to the module summary.
#'
#' @examples
#'
#' ## Lets assume that we have already generated the ObjectList, we will load it from
#' ## the folder containing the examples files. After running the 'runrewiring()', running this 
#' ## function will generate a module summary based on the ouput from 'runrewiring()'.
#'
#' objectlist <- readRDS(file=paste0(system.file('extdata',package='TraRe'),
#'                       '/prepared_rewiring_example.rds'))
#'
#'
#' ## We are going to create the folder containing
#' ## the graphs, reports, etc, and then we are deleting it.
#' ## If you want to keep it, do not run the last line.
#'
#' ## We are modifying output directory for this example.
#' objectlist$outdir <- paste(getwd(),'examplefolder',sep='/')
#'
#' runrewiring(ObjectList = objectlist)
#' createModuleSummary(ObjectList = objectlist, modmeth = "VBSR", numclus = 1, supertype = "refined", numdataset = 1)
#' unlink(objectlist$outdir,recursive = TRUE)
#'
#'
#' @export
createModuleSummary <- function(ObjectList, modmeth = "VBSR", numclus = 1, supertype = "refined", numdataset = 1) {
    
    alllabels <- ObjectList$datasets[[numdataset]]$responder[ObjectList$datasets[[numdataset]]$keepsamps]
    samps2pheno <- alllabels
    samps2pheno[which(alllabels == ObjectList$phenotype_class_vals_label[2])] <- ObjectList$phenotype_class_vals[2]
    samps2pheno[which(alllabels == ObjectList$phenotype_class_vals_label[1])] <- ObjectList$phenotype_class_vals[1]
    
    nonrespond_idxs <- names(samps2pheno)[which(samps2pheno == ObjectList$phenotype_class_vals[1])]
    responder_idxs <- names(samps2pheno)[which(samps2pheno == ObjectList$phenotype_class_vals[2])]
    
    # include modmeth
    modsumm_name <- paste0(ObjectList$outdir, "/supermod_rewiring/supermodule", numdataset,".", modmeth, ".", numclus, "/", supertype, "summ.rds")
    if (file.exists(modsumm_name)) {
        modsumm <- readRDS(modsumm_name)
    } else {
        stop(paste0("Rewiring file ", modsumm_name," does not exist, please run runrewiring() before creating module summary."))
    }
    
    # set up output html page
    dir.create(paste0(ObjectList$outdir, "/rewiring_module_summary"))
    codedir <- paste0(system.file("extdata", package = "TraRe"), "/RewiringReport/")
    dir.create(file.path(paste0(ObjectList$outdir, "/rewiring_module_summary/dataset", numdataset, ".", modmeth, ".cluster", numclus, ".", supertype, 
        "/")), showWarnings = FALSE)
    htmlinfo <- create_index_page(outdir = paste0(ObjectList$outdir, "/rewiring_module_summary/dataset", numdataset, ".", modmeth, ".cluster", numclus, 
        ".", supertype), runtag = "", codedir = codedir)
    
    orderobj <- geneOrder(modsumm, ObjectList$datasets[[numdataset]]$keepsamps, ObjectList$datasets[[numdataset]]$keeplabels, ObjectList$datasets[[numdataset]]$norm_expr_mat_keep)
    createLegendPlot(htmlinfo)
    
    # Different Sections in the html summary
    superModuleStatistics(orderobj$modregs, orderobj$modtargs, orderobj$mat, ObjectList$datasets[[numdataset]]$keeplabels, htmlinfo)
    correlationOfModuleGene(supertype, orderobj$regorder, orderobj$targetorder, orderobj$mat, orderobj$cormats, ObjectList$keeplabels, 
        htmlinfo, ObjectList$phenotype_class_vals)
    expressionTableOfModuleGenes(supertype, modsumm$nodesumm, htmlinfo)
    expressionPlotsOfModuleGenes(supertype, orderobj$regorder, orderobj$targetorder, orderobj$mat, samps2pheno, ObjectList$phenotype_class_vals, 
        htmlinfo)
    bipartiteGraphsSumm(numclus, modsumm, numdataset, modmeth, htmlinfo)
    nullDistributionOfRewiringStatistic(orderobj$mat, ObjectList$datasets[[numdataset]]$keeplabels, modmeth, supertype, htmlinfo)
    rankdf <- violinPlots(ObjectList$datasets[[numdataset]]$norm_expr_mat_keep, ObjectList$datasets[[numdataset]]$keepsamps, ObjectList$datasets[[numdataset]]$keeplabels, 
        modsumm$nodesumm, modsumm$fulledgesumm, modsumm$appendmat, htmlinfo)
    regulatorSummaryAndRank(rankdf, htmlinfo)

    return(htmlinfo)
}
