#' Preprocess prior to LINKER_run
#'
#' Preprocess data prior to GRN inference process.  Help them for more information.
#'
#' @param data_matrix Matrix of log-normalized estimated counts of the gene expression
#' data (Nr Genes x Nr samples) or SummarizedExperiment object.
#' @param geneinfo array of the regulator gene names that are in the specified data matrix.
#' @param nassay If SummarizedExperiment object is passed as input to lognorm_est_counts, name of the
#' assay containing the desired matrix. Default: 1 (first item in assay's list).
#' @param low_var_genes_th Perform a filtering to drop out low variance (<th) genes across samples. Default: 0.25
#' @param low_var_samples_th Perform a filtering to drop out low variance (<th) samples across genes. Default: 1
#' @param verbose whether to show or not messages during function run.

#' @return TraReObj containing preprocessed input matrix
#'
#' @examples
#'    ##
#'
#' @export
linker_preprocessing <- function(data_matrix, geneinfo, nassay = 1, low_var_genes_th = 0.25, low_var_samples_th = 1, verbose = TRUE){

  if (inherits(data_matrix, "SummarizedExperiment")) {
    # Get lognorm_est_counts expression matrix.
    lognorm_est_counts <- SummarizedExperiment::assays(data_matrix)[[nassay]]
  }else{
    # Get lognorm_est_counts expression matrix.
    lognorm_est_counts <- data_matrix
  }

  #get shape of our data matrix prior to dropping out low variance genes/samples
  data_shape <- dim(lognorm_est_counts)

  if (low_var_genes_th){
    genes_var_cond <- which(apply(lognorm_est_counts,1,stats::var) >= low_var_genes_th)
    lognorm_est_counts <- lognorm_est_counts[genes_var_cond,]
    if (verbose){
    methods::show(paste0(data_shape[1] - nrow(lognorm_est_counts) ,
                         ' genes have been dropped out according to variance (across samples) threshold: ',
                         low_var_genes_th, ', from: ', data_shape[1], ' to: ', nrow(lognorm_est_counts)))
      }
  }

  if (low_var_samples_th){
    samples_var_cond <- which(apply(lognorm_est_counts,2,stats::var) >= low_var_samples_th)
    lognorm_est_counts <- lognorm_est_counts[,samples_var_cond]
    if (verbose){
    methods::show(paste0(data_shape[2] - ncol(lognorm_est_counts),
                         ' samples have been dropped out according to variance (across genes) threshold: ',
                         low_var_samples_th, ', from: ', data_shape[2], ' to: ', ncol(lognorm_est_counts)))
      }
  }

  #get genenames
  genenames <- rownames(lognorm_est_counts)
  #get regulator and target idx
  regulator_filtered_idx <- which(genenames %in% geneinfo)
  target_filtered_idx <- which(!genenames %in% geneinfo)

  ##Create LinkerObj and save lognorm, target and regulator
  methods::setClass("LinkerClass", slots=list(lognorm_counts="matrix",
                                              target_idx="numeric", regulator_idx="numeric"))
  LinkerObj <- methods::new("LinkerClass", lognorm_counts = lognorm_est_counts,
                                           target_idx = target_filtered_idx,
                                           regulator_idx = regulator_filtered_idx)

  return(LinkerObj)
}
