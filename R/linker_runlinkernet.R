#' Single gene network approaches.
#'
#' `NET_run()` generate a single GRN. `NET_compute_graph_all_LASSO1se()` defines
#' the statistics of drivers and targets to be the Lasso method, choosing 1 standard error from the minimum RSS.
#' `NET_compute_graph_all_LASSOmin()` uses Lasso method, choosing the minimum RSS point. `NET_compute_graph_all_LM()`
#' uses a linear model and `NET_compute_graph_all_VBSR()` uses a Variational Bayes Spike Regression.
#'
#' @param lognorm_est_counts Matrix of log-normalized estimated counts of the gene expression
#' data (Nr Genes x Nr samples)
#' @param target_filtered_idx Index of the target genes on the lognorm_est_counts matrix.
#' @param regulator_filtered_idx Index of the regulatory genes on the lognorm_est_counts matrix.
#' @param graph_mode Chosen method(s) to generate the edges in the bipartite graph. The available options are "VBSR",
#' "LASSOmin", "LASSO1se" and "LM". By default, all methods are chosen.
#' @param FDR The False Discovery Rate correction used for the enrichment analysis. By default, 0.05.
#' @param NrCores Nr of computer cores for the parallel parts of the method. Note that
#' the parallelization is NOT initialized in any of the functions. By default, 3.
#'
#'
#' @return List containing the GRN graphs.
#'
#' @examples
#'
#'    ## Assume we have run the rewiring method and we have discovered a rewired module.
#'    ## After we have selected the drivers and targets from that modules, we can build
#'    ## a single GRN to study it separately.
#'
#'
#'    ## Imagine our rewired module consists of 5 driver genes and 40 target genes.
#'
#'    drivers <- readRDS(paste0(system.file("extdata",package="TraRe"),'/tfs_cliques_example.rds'))
#'    targets <- readRDS(paste0(system.file("extdata",package="TraRe"),'/targets_linker_example.rds'))
#'
#'    lognorm_est_counts <- rbind(drivers[1:5,],targets[1:30,])
#'    regulator_filtered_idx <- 1:5
#'    target_filtered_idx <- 5+c(1:30)
#'
#'    ## We recommend VBSR (rest of parameters are set by default)
#'    graph <- NET_run(lognorm_est_counts,target_filtered_idx,
#'                      regulator_filtered_idx,graph_mode="VBSR")
#' @export NET_run
NET_run<-function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx,
                  graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                  FDR=0.05,
                  NrCores=3)
{
  graphs<-list()
  for(j in seq_along(graph_mode)){
    graphs[[ graph_mode[j] ]] <- switch( graph_mode[j],
                                         "VBSR" = NET_compute_graph_all_VBSR(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx),
                                         "LASSOmin" = NET_compute_graph_all_LASSOmin(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx),
                                         "LASSO1se" = NET_compute_graph_all_LASSO1se(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx),
                                         "LM" = NET_compute_graph_all_LM(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx)
    )

    print(paste0("Graphs for (" ,graph_mode[j], ") computed!"))
  }

  return(list(graphs=graphs))
}
#' @export
#' @rdname NET_run
#' @param alpha feature selection parameter in case of a LASSO model to be chosen.
NET_compute_graph_all_LASSO1se<-function(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx, alpha=1-1e-06)
{

  X<-lognorm_est_counts[regulator_filtered_idx,]

  driverMat<-matrix(data = NA, nrow = length(target_filtered_idx), ncol = length(regulator_filtered_idx))
  `%dopar%` <- foreach::`%dopar%`
  #compute the LASSO1se
  idx_gene<-NULL
  driverMat<-foreach::foreach(idx_gene=seq_along(target_filtered_idx), .combine = rbind, .packages="glmnet") %dopar%
    {
      y<-lognorm_est_counts[target_filtered_idx[idx_gene],]
      fit = glmnet::cv.glmnet(t(X), y, alpha = alpha)

      b_o = stats::coef(fit,s = fit$lambda.1se)
      b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
      b_opt
    }

  rownames(driverMat)<-rownames(lognorm_est_counts)[target_filtered_idx]
  colnames(driverMat)<-rownames(lognorm_est_counts)[regulator_filtered_idx]

  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]

  g<-igraph::graph_from_incidence_matrix(driverMat)

  return(g)

}
#' @export
#' @rdname NET_run
NET_compute_graph_all_LASSOmin<-function(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx, alpha=1-1e-06)
{

  X<-lognorm_est_counts[regulator_filtered_idx,]

  driverMat<-matrix(data = NA, nrow = length(target_filtered_idx), ncol = length(regulator_filtered_idx))

  `%dopar%` <- foreach::`%dopar%`
  #compute the LASSOmin
  idx_gene<-NULL
  driverMat<-foreach::foreach(idx_gene=seq_along(target_filtered_idx), .combine = rbind, .packages="glmnet") %dopar%
    {
      y<-lognorm_est_counts[target_filtered_idx[idx_gene],]
      fit = glmnet::cv.glmnet(t(X), y, alpha = alpha)

      b_o = stats::coef(fit,s = fit$lambda.min)
      b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
      b_opt
    }

  rownames(driverMat)<-rownames(lognorm_est_counts)[target_filtered_idx]
  colnames(driverMat)<-rownames(lognorm_est_counts)[regulator_filtered_idx]

  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]

  g<-igraph::graph_from_incidence_matrix(driverMat)

  return(g)

}
#' @export
#' @rdname NET_run
NET_compute_graph_all_LM<-function(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx)
{
  GEA_per_regulator<-list()
  i<-1

  X<-t(lognorm_est_counts[regulator_filtered_idx,])

  #compute the LM
  NrTotalEdges<-length(target_filtered_idx)*length(regulator_filtered_idx)
  Pthre<-0.05/(NrTotalEdges)

  `%dopar%` <- foreach::`%dopar%`
  idx_gene<-NULL
  driverMat<-foreach::foreach(idx_gene=seq_len(target_filtered_idx), .combine=rbind) %dopar%

    {
      y<-lognorm_est_counts[target_filtered_idx[idx_gene],]
      driverVec<-numeric(length=ncol(X))
      for(i in seq_len(ncol(X)))
      {
        x<-X[,i]
        fit = stats::lm(y~x)
        s<-summary(fit)
        driverVec[i]<-(s$coefficients[2,"Pr(>|t|)"] < Pthre)
      }
      driverVec
    }

  target_genes<-rownames(lognorm_est_counts)[target_filtered_idx]
  regulators<-rownames(lognorm_est_counts)[regulator_filtered_idx]

  rownames(driverMat)<-target_genes
  colnames(driverMat)<-regulatory_genes

  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]

  g<-igraph::graph_from_incidence_matrix(driverMat)

  return(g)

}
#' @export
#' @rdname NET_run
NET_compute_graph_all_VBSR<-function(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx)
{

  X<-lognorm_est_counts[regulator_filtered_idx,]

  driverMat<-matrix(data = NA, nrow = length(target_filtered_idx), ncol = length(regulator_filtered_idx))

  `%dopar%` <- foreach::`%dopar%`
  #compute the VBSR
  idx_gene<-NULL
  driverMat<-foreach::foreach(idx_gene=seq_along(target_filtered_idx), .combine = rbind, .packages = "vbsr") %dopar%
    {
      y<-lognorm_est_counts[target_filtered_idx[idx_gene],]
      res<-vbsr::vbsr(y,t(X),n_orderings = 15,family='normal')
      betas<-res$beta
      betas[res$pval > 0.05/(length(target_filtered_idx)*length(regulator_filtered_idx))]<-0
      betas
    }
  rownames(driverMat)<-rownames(lognorm_est_counts)[target_filtered_idx]
  colnames(driverMat)<-rownames(lognorm_est_counts)[regulator_filtered_idx]

  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]

  g<-igraph::graph_from_incidence_matrix(driverMat)

  return(g)

}
