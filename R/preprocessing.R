#' Preprocess prior to LINKER_run and run_rewiring method
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
#' @export trare_preprocessing
trare_preprocessing <- function(data_matrix, geneinfo = NULL, nassay = 1, low_var_genes_th = 0.25, low_var_samples_th = 1, verbose = TRUE){

  if (inherits(data_matrix, "SummarizedExperiment")) {

    # Get lognorm_est_counts expression matrix.
    lognorm_est_counts <- SummarizedExperiment::assays(data_matrix)[[nassay]]
    if (is.null(geneinfo)){
      geneinfo <- rownames(lognorm_est_counts)[which(SummarizedExperiment::rowData(data_matrix)[['gene_info']] == 1)]
    }

  }else{

    # Get lognorm_est_counts expression matrix.
    lognorm_est_counts <- data_matrix
    if (is.null(geneinfo)){
      stop('geneinfo matrix has to be specified!')
    }
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

  ## Include phenotype in case is detected
  phenotype <- 0
  if (inherits(data_matrix, "SummarizedExperiment")) {
    if (dim(SummarizedExperiment::colData(data_matrix))[2]){
        phenotype_df <- SummarizedExperiment::colData(data_matrix)
        if ('phenotype'%in%colnames(phenotype_df)){
            message('Also including phenotype for rewiring future uses')
          if (inherits(phenotype_df[,'phenotype'], "factor")){
            phenotype <- as.numeric(phenotype_df[colnames(lognorm_est_counts), 'phenotype'],labels = c(0,1))) - 1
          }else{
            phenotype <- as.numeric(factor(phenotype_df[colnames(lognorm_est_counts), 'phenotype'], labels = c(0,1))) - 1
        }
    }
  }

  ##Create TraReObj and save lognorm, target and regulator
  methods::setClass("TraReClass", slots=list(lognorm_counts="matrix", target_idx="numeric", 
                                              regulator_idx="numeric", pheno = "numeric"))

  TraReObj <- methods::new("TraReClass", lognorm_counts = lognorm_est_counts,
                                           target_idx = target_filtered_idx,
                                           regulator_idx = regulator_filtered_idx,
                                           pheno = phenotype)

  return(TraReObj)
}

#' @export
#' @rdname trare_preprocessing
#' @param TraReObj TraReObj from trare_preprocesing
#' @param phenotype_f file containing samples as rownames and a column call 'phenotype' containing the labels
rewiring_add_phenotype <- function(TraReObj, phenotype_f){

  TraReObj@lognorm_counts <- TraReObj@lognorm_counts[,rownames(phenotype_f)]
  TraReObj@pheno <- as.numeric(factor(phenotype_f[,'phenotype'],labels = c(0,1))) - 1
  
  return(TraReObj)
}