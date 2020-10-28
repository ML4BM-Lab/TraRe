#' GRN inference via selected model.
#'
#' Gene Regulatory Network inference via model selection. `NET_compute_graph_all_LASSO1se()` defines
#' the statistics of drivers and targets to be the Lasso method, choosing 1 standard error from the minimum RSS.
#' `NET_compute_graph_all_LASSOmin()` uses Lasso method, choosing the minimum RSS point. `NET_compute_graph_all_LM()`
#' uses a linear model and `NET_compute_graph_all_VBSR()` uses a Variational Bayes Spike Regression.
#'
#'
#' @param lognorm_est_counts Matrix of log-normalized estimated counts of the gene expression data (Nr Genes x Nr samples)
#' @param target_filtered_idx Index of the target genes on the lognorm_est_counts matrix.
#' @param regulator_filtered_idx Index of the regulatory genes on the lognorm_est_counts matrix.
#' @param graph_mode Chosen method(s) to generate the edges in the bipartite graph. The available options are "VBSR", "LASSOmin", "LASSO1se" and "LM". By default, all methods are chosen.
#' @param FDR The False Discovery Rate correction used for the enrichment analysis. By default, 0.05.
#' @param NrCores Nr of computer cores for the parallel parts of the method. Note that the parallelization is NOT initialized in any of the functions. By default, 30.
#'
#'
#' @return List containing the GRN graphs and GEA results.
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
#'                      target_filtered_idx,Gene_set_collections)
#'    }
#' @export NET_run
NET_run<-function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx,
                  graph_mode=c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                  FDR=0.05,
                  NrCores=30)
{
  graphs<-list()
  for(j in 1:length(graph_mode)){
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
#' @param Data Matrix of log-normalized estimated counts of the gene expression data (Nr Genes x Nr samples).
#' @param regulators_idx Index of the regulatory genes on the lognorm_est_counts matrix.
#' @param target_idx Index of the target genes on the lognorm_est_counts matrix.
#' @param alpha feature selection parameter.
NET_compute_graph_all_LASSO1se<-function(Data, regulators_idx, target_idx, alpha=1-1e-06)
{

  X<-Data[regulators_idx,]

  driverMat<-matrix(data = NA, nrow = length(target_idx), ncol = length(regulators_idx))
  `%dopar%` <- foreach::`%dopar%`
  #compute the LASSO1se
  driverMat<-foreach::foreach(idx_gene=1:length(target_idx), .combine = rbind, .packages="glmnet") %dopar%
    {
      y<-Data[target_idx[idx_gene],]
      fit = glmnet::cv.glmnet(t(X), y, alpha = alpha)

      b_o = stats::coef(fit,s = fit$lambda.1se)
      b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
      b_opt
    }

  rownames(driverMat)<-rownames(Data)[target_idx]
  colnames(driverMat)<-rownames(Data)[regulators_idx]

  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]

  g<-igraph::graph_from_incidence_matrix(driverMat)

  return(g)

}
#' @export
#' @rdname NET_run
NET_compute_graph_all_LASSOmin<-function(Data, regulators_idx, target_idx, alpha=1-1e-06)
{

  X<-Data[regulators_idx,]

  driverMat<-matrix(data = NA, nrow = length(target_idx), ncol = length(regulators_idx))

  `%dopar%` <- foreach::`%dopar%`
  #compute the LASSOmin
  driverMat<-foreach::foreach(idx_gene=1:length(target_idx), .combine = rbind, .packages="glmnet") %dopar%
    {
      y<-Data[target_idx[idx_gene],]
      fit = glmnet::cv.glmnet(t(X), y, alpha = alpha)

      b_o = stats::coef(fit,s = fit$lambda.min)
      b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
      b_opt
    }

  rownames(driverMat)<-rownames(Data)[target_idx]
  colnames(driverMat)<-rownames(Data)[regulators_idx]

  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]

  g<-igraph::graph_from_incidence_matrix(driverMat)

  return(g)

}
#' @export
#' @rdname NET_run
NET_compute_graph_all_LM<-function(Data, regulators_idx, target_idx)
{
  GEA_per_regulator<-list()
  i<-1

  X<-t(Data[regulators_idx,])

  #compute the LM
  NrTotalEdges<-length(target_idx)*length(regulators_idx)
  Pthre<-0.05/(NrTotalEdges)

  `%dopar%` <- foreach::`%dopar%`
  driverMat<-foreach::foreach(idx_gene=1:length(target_idx), .combine=rbind) %dopar%

    {
      y<-Data[target_idx[idx_gene],]
      driverVec<-numeric(length=ncol(X))
      for(i in 1:ncol(X))
      {
        x<-X[,i]
        fit = stats::lm(y~x)
        s<-summary(fit)
        driverVec[i]<-(s$coefficients[2,"Pr(>|t|)"] < Pthre)
      }
      driverVec
    }

  target_genes<-rownames(Data)[target_idx]
  regulators<-rownames(Data)[regulators_idx]

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
NET_compute_graph_all_VBSR<-function(Data, regulators_idx, target_idx)
{

  X<-Data[regulators_idx,]

  driverMat<-matrix(data = NA, nrow = length(target_idx), ncol = length(regulators_idx))

  `%dopar%` <- foreach::`%dopar%`
  #compute the VBSR
  driverMat<-foreach::foreach(idx_gene=1:length(target_idx), .combine = rbind, .packages = "vbsr") %dopar%
    {
      y<-Data[target_idx[idx_gene],]
      res<-vbsr::vbsr(y,t(X),n_orderings = 15,family='normal')
      betas<-res$beta
      betas[res$pval > 0.05/(length(target_idx)*length(regulators_idx))]<-0
      betas
    }
  rownames(driverMat)<-rownames(Data)[target_idx]
  colnames(driverMat)<-rownames(Data)[regulators_idx]

  regulated_genes<-which(rowSums(abs(driverMat))!=0)
  regulatory_genes<-which(colSums(abs(driverMat))!=0)
  driverMat<-driverMat[regulated_genes,]
  driverMat<-driverMat[,regulatory_genes]

  g<-igraph::graph_from_incidence_matrix(driverMat)

  return(g)

}
