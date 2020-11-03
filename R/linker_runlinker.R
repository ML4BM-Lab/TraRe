#' GRN inference via selected model.
#'
#' Gene Regulatory Network inference via model selection. Consists of two phases,
#' `LINKER_runPhase1()` and `LINKER_runPhase2()`. Help them for more information.
#'
#' @param lognorm_est_counts Matrix of log-normalized estimated counts of the gene expression
#' data (Nr Genes x Nr samples)
#' @param target_filtered_idx Index of the target genes on the lognorm_est_counts matrix.
#' @param regulator_filtered_idx Index of the regulatory genes on the lognorm_est_counts matrix.
#' @param link_mode Chosen method(s) to link module eigengenes to regulators. The available options are
#' "VBSR", "LASSOmin", "LASSO1se" and "LM". By default, all methods are chosen.
#' @param graph_mode Chosen method(s) to generate the edges in the bipartite graph. The available options
#' are "VBSR", "LASSOmin", "LASSO1se" and "LM". By default, all methods are chosen.
#' @param module_rep Method selected for use. Default set to MEAN.
#' @param NrModules Number of modules that are a priori to be found (note that the final number of modules
#' discovered may differ from this value). By default, 100 modules.
#' @param corrClustNrIter output from preparedata(). By default, 100.
#' @param Nr_bootstraps Number of bootstrap of Phase I. By default, 10.
#' @param FDR The False Discovery Rate correction used for the enrichment analysis. By default, 0.05.
#' @param NrCores Nr of computer cores for the parallel parts of the method. Note that the parallelization
#' is NOT initialized in any of the functions. By default, 3.
#'
#'
#' @return List containing the GRN raw results, GRN modules and GRN graphs.
#'
#' @examples
#'    ## For this example, we are going to join 15 drivers and 100 targets from the example folder.
#'
#'    drivers <- readRDS(paste0(system.file("extdata",package="TraRe"),'/tfs_cliques_example.rds'))
#'    targets <- readRDS(paste0(system.file("extdata",package="TraRe"),'/targets_linker_example.rds'))
#'
#'    lognorm_est_counts <- rbind(drivers[seq_len(15),],targets[seq_len(100),])
#'
#'    ## We create the index for drivers and targets in the log-normalized gene expression matrix.
#'    L <- 15
#'    regulator_filtered_idx <- seq_len(L)
#'    target_filtered_idx <- L+c(seq_len(100))
#'
#'
#'    ## We recommend VBSR.
#'    \dontrun{
#'    linkeroutput <- LINKER_run(lognorm_est_counts,target_filtered_idx,regulator_filtered_idx,
#'                               link_mode="VBSR",graph_mode="VBSR",NrModules=3,Nr_bootstraps=1)
#'    }
#'
#'
#' @export
LINKER_run<-function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx,
                     link_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                     graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                     module_rep="MEAN",
                     NrModules=100,
                     corrClustNrIter=100,
                     Nr_bootstraps=10,
                     FDR=0.05,
                     NrCores=3)

{
  res<-list()
  modules<-list()
  graphs<-list()
  for(i in seq_along(link_mode)){
    res[[ link_mode[i] ]]<-LINKER_runPhase1(lognorm_est_counts, target_filtered_idx,  regulator_filtered_idx,
                                            NrModules,NrCores=NrCores,
                                            mode=link_mode[i], used_method=module_rep,
                                            corrClustNrIter=corrClustNrIter,Nr_bootstraps=Nr_bootstraps)

    modules[[ link_mode[i] ]]<-LINKER_extract_modules(res[[ link_mode[i] ]])
    print(paste0("Link mode ",link_mode[i]," completed!"))

    graphs[[ link_mode[i] ]]<-list()
    for(j in seq_along(graph_mode)){
      graphs[[ link_mode[i] ]][[ graph_mode[j] ]] <- LINKER_runPhase2(modules[[ link_mode[i] ]],
                                                                                  lognorm_est_counts, mode=graph_mode[j])
      print(paste0("Graphs for (",link_mode[i],",",graph_mode[j], ") computed!"))
    }

  }

  return(list(raw_results=res,modules=modules,graphs=graphs))


}
