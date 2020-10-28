#' GRN inference via selected model.
#'
#' Gene Regulatory Network inference via model selection.
#'
#' @param lognorm_est_counts Matrix of log-normalized estimated counts of the gene expression data (Nr Genes x Nr samples)
#' @param target_filtered_idx Index of the target genes on the lognorm_est_counts matrix.
#' @param regulator_filtered_idx Index of the regulatory genes on the lognorm_est_counts matrix.
#' @param link_mode Chosen method(s) to link module eigengenes to regulators. The available options are "VBSR", "LASSOmin", "LASSO1se" and "LM". By default, all methods are chosen.
#' @param graph_mode Chosen method(s) to generate the edges in the bipartite graph. The available options are "VBSR", "LASSOmin", "LASSO1se" and "LM". By default, all methods are chosen.
#' @param module_rep Method selected for use. Default set to LINKER.
#' @param NrModules Number of modules that are a priori to be found (note that the final number of modules discovered may differ from this value). By default, 100 modules.
#' @param corrClustNrIter output from preparedata(). By default, 100.
#' @param Nr_bootstraps Number of bootstrap of Phase I. By default, 10.
#' @param FDR The False Discovery Rate correction used for the enrichment analysis. By default, 0.05.
#' @param NrCores Nr of computer cores for the parallel parts of the method. Note that the parallelization is NOT initialized in any of the functions. By default, 30.
#'
#'
#' @return List containing the GRN raw results, GRN modules, GRN graphs and GEA results.
#'
#' @examples
#'    \dontrun{
#'    #(we generally will load it. We will assume it contains 1000 genes and 43 samples, i.e 1000x43)
#'    lognorm_est_counts <- readRDS(file="../geneexpressiondata.rds")
#'
#'    #suppose lognorm_est_counts first 10 genes are drivers
#'    regulator_filtered_idx <- c(1,2,3,4,5,6,7,8,9,10)
#'
#'    #and the rest are targets
#'    target_filtered_idx <- c(11:1000)
#'    Gene_set_collections <- readRDS(file=".../gsc.rds")
#'
#'    #The rest have default values.
#'    foo <- LINKER_run(lognorm_est_counts,regulator_filtered_idx,
#'                       target_filtered_idx,Gene_set_collections)
#'    }
#' @export
LINKER_run<-function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx,
                     link_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                     graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                     module_rep="MEAN",
                     NrModules=100,
                     corrClustNrIter=100,
                     Nr_bootstraps=10,
                     FDR=0.05,
                     NrCores=30)

{
  res<-list()
  modules<-list()
  graphs<-list()
  for(i in 1:length(link_mode)){
    res[[ link_mode[i] ]]<-LINKER_runPhase1(lognorm_est_counts, target_filtered_idx,  regulator_filtered_idx,
                                            NrModules,NrCores=NrCores,
                                            mode=link_mode[i], used_method=module_rep,
                                            corrClustNrIter=corrClustNrIter,Nr_bootstraps=Nr_bootstraps)

    modules[[ link_mode[i] ]]<-LINKER_extract_modules(res[[ link_mode[i] ]])
    print(paste0("Link mode ",link_mode[i]," completed!"))

    graphs[[ link_mode[i] ]]<-list()
    for(j in 1:length(graph_mode)){
      graphs[[ link_mode[i] ]][[ graph_mode[j] ]] <- LINKER_compute_modules_graph(modules[[ link_mode[i] ]], lognorm_est_counts, mode=graph_mode[j])
      print(paste0("Graphs for (",link_mode[i],",",graph_mode[j], ") computed!"))
    }

  }

  return(list(raw_results=res,modules=modules,graphs=graphs))


}
