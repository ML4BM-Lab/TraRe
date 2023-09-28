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
#' ## To create html summary on previously run supermodel in new data, we will have
#' ## generated a new ObjectList using the new data by running ‘preparerewiring()’.
#' ##  Running the #' ## ‘createModuleSummary()’ with the new ObjectList will create
#' ## a module summary based on the new dataset but with the gene module from the
#' ## 'runrewiring()' results.
#'
#' objectlist <- readRDS(file=paste0(system.file('extdata',package='TraRe'),
#'                       '/prepared_rewiring_example.rds'))
#'
#' ## We are going to create the folder containing
#' ## the graphs, reports, etc, and then we are deleting it.
#' ## If you want to keep it, do not run the last line.
#'
#' ## We are modifying output directory for this example.
#' objectlist$outdir <- paste(getwd(),'examplefolder',sep='/')
#'
#' ## runrewiring(ObjectList = objectlist)
#' ##createModuleSummary(ObjectList = objectlist, modmeth = "VBSR", numclus = 1,
#' ##                    supertype = "refined", numdataset = 1)
#' ##unlink(objectlist$outdir,recursive = TRUE)
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
  createLegendPlot(htmlinfo,1)
  
  write(paste0("<table style='width:100%' bgcolor='gray'><tr><td><h1>", "Complete Rewiring Module ",numclus, " Summary", "</h1></td></tr></table><br>\n"),
        paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
  # ref_cluster_index <- paste0("<a href = '../../supermodule", numdataset,".", modmeth, ".", numclus, "/index.html'>Return to Cluster Summary</a><br>")
  ref_cluster_index <- paste0("<button onclick=\"window.location.href='../../supermodule",numdataset,".", modmeth, ".", numclus, "/index.html'\";>\n\tReturn to Cluster Summary\n</button><br><br>")

  write(ref_cluster_index, paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
  # ref_curr_index <- paste0("<a href = '../rewiring_module_summary/dataset", numdataset, ".", modmeth, ".cluster", numclus, ".", supertype, "/index.html'> Complete Rewiring Module",numclus," Summary ",supertype,"</a><br>")
  ref_curr_index <- paste0("<button onclick=\"window.location.href='../rewiring_module_summary/dataset",numdataset,".", modmeth, ".cluster", numclus,".",supertype,"/index.html'\";>\n\tComplete Rewiring Module",numclus," Summary ",supertype,"\n</button><br><br>")
  clustersumm_dir <- paste0(ObjectList$outdir, dir_prefix, "/index.html")
  write(ref_curr_index, clustersumm_dir, append = TRUE)
  
  pheno <- ObjectList$datasets[[numdataset]]$pheno
  phenosamples <- ObjectList$datasets[[numdataset]]$phenosamples
  
  samps2pheno <- ObjectList$datasets[[numdataset]]$phenosamples
  # Decide if add extra parameter
  phenotype_class_vals <- c("pheno0","pheno1")
  samps2pheno[pheno] <- phenotype_class_vals[2]
  samps2pheno[!pheno] <- phenotype_class_vals[1]
  
  nonrespond_idxs <- phenosamples[!pheno]
  responder_idxs <- phenosamples[pheno]
  
  obj_runmoddata <- list(regulators = orderobj$modregs, target_genes = orderobj$modtargs)
  norm_expr_mat_keep <- ObjectList$datasets[[numdataset]]$lognorm_est_counts
  name2idx <- ObjectList$datasets[[numdataset]]$name2idx
  
  obj_nodesumm <- module_node_summary(norm_expr_mat_keep, obj_runmoddata, name2idx, nonrespond_idxs, responder_idxs)
  obj_edgesumm <- module_edge_summary(norm_expr_mat_keep, obj_runmoddata, name2idx, nonrespond_idxs, responder_idxs)
  
  # Different Sections in the html summary
  superModuleStatistics(orderobj$modregs, orderobj$modtargs, orderobj$mat, as.numeric(pheno), htmlinfo)
  correlationOfModuleGene(supertype, orderobj$regorder, orderobj$targetorder, orderobj$mat, orderobj$cormats,
                          htmlinfo,phenotype_class_vals )
  expressionTableOfModuleGenes(supertype, obj_nodesumm, htmlinfo)
  plotzlim<- createLegendPlot(htmlinfo,2,orderobj$mat)
  expressionPlotsOfModuleGenes(supertype, orderobj$regorder, orderobj$targetorder, orderobj$mat, samps2pheno, phenotype_class_vals,
                               htmlinfo,plotzlim)
  bipartiteGraphsSumm(numclus, obj_nodesumm, obj_edgesumm, numdataset, modmeth, htmlinfo)
  
  nullDistributionOfRewiringStatistic(orderobj$mat, pheno, modmeth, supertype, htmlinfo)
  rankdf <- violinPlots(norm_expr_mat_keep, phenosamples, pheno,
                        obj_nodesumm, modsumm$fulledgesumm, orderobj$modtargs, htmlinfo)
  regulatorSummaryAndRank(rankdf, htmlinfo)
  write(ref_cluster_index, paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
  return(htmlinfo)
}
