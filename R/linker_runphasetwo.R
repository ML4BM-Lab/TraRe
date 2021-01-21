#' Phase II : bipartitive graphs generation
#'
#' Run second phase of the linker method where a bipartitive graph is generated from the phase I output.
#' A bipartite graph is a set of graph nodes decomposed into two disjoint sets such that no
#' two graph nodes within the same set are adjacent.
#'
#' @param Data Matrix of log-normalized estimated counts of the gene expression
#' data (Nr Genes x Nr samples)
#' @param mode Chosen method(s) to link module eigengenes to regulators. The available options are
#' "VBSR", "LASSOmin", "LASSO1se" and "LM". By default, all methods are chosen.
#' @param NrCores Nr of computer cores for the parallel parts of the method. Note that the parallelization
#' is NOT initialized in any of the functions. By default, 2.
#' @param modules Modules obtained from the phase I linker output.
#' @param alpha alpha parameter if a LASSO model is chosen.
#'
#'
#' @examples
#'
#'    ## We are going to proceed in the same manner as in the `liner_runphaseone()` example
#'    ## but we start at the end of it, loading the output from the example folder.
#'
#'    phaseone <- readRDS(paste0(system.file("extdata",package="TraRe"),
#'                               '/linker_phaseone_example.rds'))
#'
#'    ## We are loading drivers and targets to generate lognorm_est_counts, as we need it
#'    ## for the phase 2.
#'
#'    drivers <- readRDS(paste0(system.file("extdata",package="TraRe"),
#'                                          '/tfs_linker_example.rds'))
#'    targets <- readRDS(paste0(system.file("extdata",package="TraRe"),
#'                                          '/targets_linker_example.rds'))
#'
#'    lognorm_est_counts <- as.matrix(rbind(drivers,targets))
#'
#'    ## Now we proceed to call `LINKER_runPhase2()`.
#'    ## We first, we need to extract modules from the `LINKER_runPhase1()` output.
#'
#'    modules_phaseone<-LINKER_extract_modules(phaseone)
#'
#'    ## Now we generate the bipartitive graph from the extracted modules
#'
#'
#'    graph <- LINKER_runPhase2(modules=modules_phaseone,Data=lognorm_est_counts,
#'                              NrCores=1,mode="LM")
#'
#'
#' @return igraph object containing the related drivers and targets in the form of a bipartitive graph.
#' @export

LINKER_runPhase2<-function(modules,Data,NrCores, mode="VBSR",alpha=1-1e-06)
{

  bp_g<-list()
  i<-1

  # this will register nr of cores/threads, keep this here so
  # the user can decide how many cores based on their hardware.
  cl <- parallel::makeCluster(NrCores)
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`

  mod_idx<-NULL
  bp_g<-foreach::foreach(mod_idx=seq_along(modules), .packages = c("vbsr","glmnet","igraph"))%dopar%
  {
    targetgenes<-unlist(modules[[mod_idx]]$target_genes)
    regulators<-unlist(modules[[mod_idx]]$regulators)
    X<-Data[regulators,]

    #We need to handle the special case where only one regulator regulates a module/community
    if(length(regulators)<2){

      non_zero_beta<-modules[[mod_idx]]$regulatory_program[which(modules[[mod_idx]]$regulatory_program != 0)]
      if(length(non_zero_beta) != 1){
        warning("NON_ZERO_BETA != 1")
      }

      driverMat<-matrix(data = non_zero_beta, nrow = length(targetgenes), ncol = length(regulators))

      rownames(driverMat)<-targetgenes
      colnames(driverMat)<-regulators
    }
    else{

      driverMat<-matrix(data = NA, nrow = length(targetgenes), ncol = length(regulators))

      for(idx_gene in seq_along(targetgenes))
      {
        y<-Data[targetgenes[idx_gene],]

        if(mode=="VBSR")
        {
          res<-vbsr::vbsr(y,t(X),n_orderings = 15,family='normal')
          betas<-res$beta
          betas[res$pval > 0.05/(length(regulators)*length(targetgenes))]<-0
          driverMat[idx_gene,]<-betas
        }
        else if(mode=="LASSOmin")
        {
          fit <- glmnet::cv.glmnet(t(X), y, alpha = alpha)

          b_o <- stats::coef(fit,s = fit$lambda.min)
          b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
          driverMat[idx_gene,]<-b_opt
        }
        else if(mode=="LASSO1se")
        {
          fit <- glmnet::cv.glmnet(t(X), y, alpha = alpha)

          b_o <- stats::coef(fit,s = fit$lambda.1se)
          b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
          driverMat[idx_gene,]<-b_opt
        }
        else if(mode=="LM")
        {
          for(idx_regs in seq_along(regulators))
          {
            x<-t(X)[,idx_regs]
            fit <- stats::lm(y~x)
            s<-summary(fit)
            driverMat[idx_gene,idx_regs]<-s$coefficients[2,"Pr(>|t|)"]<0.05/(length(targetgenes)*length(regulators))
          }
        }
        else
        {
          warning("MODE NOT RECOGNIZED")
        }

      }

      rownames(driverMat)<-targetgenes
      colnames(driverMat)<-regulators

      regulated_genes<-which(rowSums(abs(driverMat))!=0)
      regulatory_genes<-which(colSums(abs(driverMat))!=0)

      # We need to treat the special cases independently
      if(length(regulated_genes)<2){
        driverMat<-driverMat[regulated_genes,]
        driverMat<-driverMat[regulatory_genes]
      }
      else if(length(regulatory_genes)<2){
        driverMat<-driverMat[,regulatory_genes]
        driverMat<-driverMat[regulated_genes]
      }
      else{
        driverMat<-driverMat[,regulatory_genes]
        driverMat<-driverMat[regulated_genes,]
      }

    }

    igraph::graph_from_incidence_matrix(driverMat)
  }
  parallel::stopCluster(cl)

  return(bp_g)
}
