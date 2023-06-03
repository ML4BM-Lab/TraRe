#' Preprocess prior to LINKER_run and run_rewiring method
#' 
#' Preprocess data prior to GRN inference process.  Help them for more information.
#'
#' @param data_matrix Matrix of log-normalized estimated counts of the gene expression
#' data (Nr Genes x Nr samples) or SummarizedExperiment object (with or without the sample phenotype information).
#' @param geneinfo array of the regulator gene names that are in the specified data matrix.
#' @param nassay If SummarizedExperiment object is passed as input to lognorm_est_counts, name of the
#' assay containing the desired matrix. Default: 1 (first item in assay's list).
#' @param low_var_genes_th Perform a filtering to drop out low variance (<th) genes across samples. Default: 0.25
#' @param low_var_samples_th Perform a filtering to drop out low variance (<th) samples across genes. Default: 1
#' @param verbose whether to show or not messages during function run.
#' @return TraReObj containing preprocessed input matrix
#'
#' @examples
#' #For this example, we are going to load a example matrix
#' lognorm_est_counts_p <- paste0(system.file('extdata',package='TraRe'),
#'                                  '/expression_rewiring_example.txt')
#' lognorm_est_counts <- as.matrix(read.delim(lognorm_est_counts_p, header=TRUE,row.names=1))
#' 
#' # Load gene info, its an array of regulators' names.
#' gene_info_p <- paste0(system.file('extdata',package='TraRe'),
#'                          '/geneinfo_rewiring_example.txt')
#' gene_info <- read.delim(gene_info_p,header=TRUE)
#' geneinfo <- gene_info[gene_info[,'regulator'] == 1,'uniq_isos']
#'
#' #TraReObj <- trare_preprocessing(data_matrix = lognorm_est_counts,
#'   ##                                  geneinfo = geneinfo, verbose = FALSE)
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
            message('Also including phenotype for rewiring future use')
            phenotype <- as.numeric(factor(phenotype_df[colnames(lognorm_est_counts), 'phenotype'], labels = c(0,1))) - 1
        }
    }
  }

  ##Create TraReObj and save lognorm, target and regulator
  # TraReClass <- methods::setClass("TraReClass", slots=list(lognorm_counts="matrix", target_idx="numeric", 
  #                                             regulator_idx="numeric", pheno = "numeric"))

  # TraReObj <- methods::new("TraReClass", lognorm_counts = lognorm_est_counts,
  #                                          target_idx = target_filtered_idx,
  #                                          regulator_idx = regulator_filtered_idx,
  #                                          pheno = phenotype)
  
  TraReObj <- TraReClass(lognorm_counts = lognorm_est_counts,
                           target_idx = target_filtered_idx,
                           regulator_idx = regulator_filtered_idx,
                           pheno = phenotype)

  return(TraReObj)
}
#' @export TraReClass
#' @importFrom methods new
#' @rdname trare_preprocessing
#' @slot lognorm_counts Matrix of log-normalized estimated counts of the gene expression data (Nr Genes x Nr samples),
#' @slot target_idx Index array of the target genes on the lognorm_counts matrix.
#' @slot regulator_idx Index array of the regulatory genes on the lognorm_counts matrix.
#' @slot pheno When provided dataframe containing samples as rownames and a column named 'phenotype' containing the binary phenotype labels (character of numeric). By default 0.
#' 
#' @include preprocessing.R
#' @return class generator function for class "TraReClass"
 
##Create TraReObj and save lognorm, target and regulator
TraReClass <- methods::setClass("TraReClass", 
                                slots=list(lognorm_counts="matrix", 
                                           target_idx="numeric", 
                                           regulator_idx="numeric", 
                                           pheno = "numeric"))


#' @export
#' @rdname trare_preprocessing
#' @param TraReObj TraReObj from trare_preprocesing
#' @param phenotype_f dataframe containing samples as rownames and a column named 'phenotype' containing the binary phenotype labels (character of numeric)
rewiring_add_phenotype <- function(TraReObj, phenotype_f){

  TraReObj@lognorm_counts <- TraReObj@lognorm_counts[,rownames(phenotype_f)]
  TraReObj@pheno <- as.numeric(factor(phenotype_f[,'phenotype'],labels = c(0,1))) - 1
  
  return(TraReObj)
}
