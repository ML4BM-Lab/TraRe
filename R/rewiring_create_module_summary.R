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
    dir_prefix <- paste0("/supermod_rewiring/supermodule", numdataset,".", modmeth, ".", numclus)
    modsumm_name <- paste0(ObjectList$outdir, dir_prefix, "/", supertype, "summ.rds")
    if (file.exists(modsumm_name)) {
        modsumm <- readRDS(modsumm_name)
    } else {
        stop(paste0("Rewiring file ", modsumm_name," does not exist, please run runrewiring() before creating module summary."))
    }
    
    # set up output html page
    dir.create(paste0(ObjectList$outdir, "/supermod_rewiring/rewiring_module_summary"))
    codedir <- paste0(system.file("extdata", package = "TraRe"), "/RewiringReport/")
    dir.create(file.path(paste0(ObjectList$outdir, "/supermod_rewiring/rewiring_module_summary/dataset", numdataset, ".", modmeth, ".cluster", numclus, ".", supertype, 
        "/")), showWarnings = FALSE)
    htmlinfo <- create_index_page(outdir = paste0(ObjectList$outdir, "/supermod_rewiring/rewiring_module_summary/dataset", numdataset, ".", modmeth, ".cluster", numclus, 
        ".", supertype), runtag = "", codedir = codedir)
    
    orderobj <- geneOrder(modsumm, ObjectList, numdataset)
    createLegendPlot(htmlinfo)
    ref_cluster_index <- paste0("<a href = '../../supermodule", numdataset,".", modmeth, ".", numclus, "/index.html'>Return to Cluster Summary</a><br>")
    write(ref_cluster_index, paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)

    ref_curr_index <- paste0("<a href = '../rewiring_module_summary/dataset", numdataset, ".", modmeth, ".cluster", numclus, ".", supertype, "/index.html'>Complete Rewiring Module Summary ",supertype,"</a><br>")
    clustersumm_dir <- paste0(ObjectList$outdir, dir_prefix, "/index.html")
    write(ref_curr_index, clustersumm_dir, append = TRUE)

    alllabels <- ObjectList$datasets[[numdataset]]$responder[ObjectList$datasets[[numdataset]]$keepsamps]
    samps2pheno <- alllabels
    samps2pheno[which(alllabels == ObjectList$phenotype_class_vals_label[2])] <- ObjectList$phenotype_class_vals[2]
    samps2pheno[which(alllabels == ObjectList$phenotype_class_vals_label[1])] <- ObjectList$phenotype_class_vals[1]
    
    nonrespond_idxs <- names(samps2pheno)[which(samps2pheno == ObjectList$phenotype_class_vals[1])]
    responder_idxs <- names(samps2pheno)[which(samps2pheno == ObjectList$phenotype_class_vals[2])]

    # obj_regs = intersect(modsumm$runmoddata$regulators,ObjectList$datasets[[numdataset]]$allregs)
    # obj_targs = intersect(modsumm$runmoddata$target,ObjectList$datasets[[numdataset]]$alltargs)
    obj_runmoddata <- list(regulators = orderobj$modregs, target_genes = orderobj$modtargs)
    norm_expr_mat_keep <- ObjectList$datasets[[numdataset]]$norm_expr_mat_keep
    name2idx <- ObjectList$datasets[[numdataset]]$name2idx
    # obj_modsumm <- summarize_module(norm_expr_mat_keep, obj_runmoddata, ObjectList$datasets[[numdataset]]$name2idx, nonrespond_idxs, responder_idxs) 
    
    obj_nodesumm <- module_node_summary(norm_expr_mat_keep, obj_runmoddata, name2idx, nonrespond_idxs, responder_idxs)
    obj_edgesumm <- module_edge_summary(norm_expr_mat_keep, obj_runmoddata, name2idx, nonrespond_idxs, responder_idxs)

    # Different Sections in the html summary
    superModuleStatistics(orderobj$modregs, orderobj$modtargs, orderobj$mat, ObjectList$datasets[[numdataset]]$keeplabels, htmlinfo)
    correlationOfModuleGene(supertype, orderobj$regorder, orderobj$targetorder, orderobj$mat, orderobj$cormats, ObjectList$keeplabels, 
        htmlinfo, ObjectList$phenotype_class_vals)
    expressionTableOfModuleGenes(supertype, obj_nodesumm, htmlinfo)
    expressionPlotsOfModuleGenes(supertype, orderobj$regorder, orderobj$targetorder, orderobj$mat, samps2pheno, ObjectList$phenotype_class_vals, 
        htmlinfo)
    # bipartiteGraphsSumm(ObjectList, numclus, modsumm, numdataset, modmeth, supertype, htmlinfo)
    bipartiteGraphsSumm(numclus, obj_nodesumm, obj_edgesumm, numdataset, modmeth, htmlinfo)
    nullDistributionOfRewiringStatistic(orderobj$mat, ObjectList$datasets[[numdataset]]$keeplabels, modmeth, supertype, htmlinfo)
    rankdf <- violinPlots(ObjectList$datasets[[numdataset]]$norm_expr_mat_keep, ObjectList$datasets[[numdataset]]$keepsamps, ObjectList$datasets[[numdataset]]$keeplabels, 
        obj_nodesumm, modsumm$fulledgesumm, orderobj$modtargs, htmlinfo)
    regulatorSummaryAndRank(rankdf, htmlinfo)

    return(htmlinfo)
}
