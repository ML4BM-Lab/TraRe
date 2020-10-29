#' Prepare rewiring data for running the method.
#'
#' Prepare all the files neccessary for running the gene regulatory network rewiring method.
#'
#' @param name Desired name of the folder which is generated.
#' @param linker_saved_file Output file from linker function path.
#' @param expr_matrix_file Lognorm counts of the gene expression matrix path.
#' @param gene_info_file path of a two-column file containing genes and 'regulator' boolean variable.
#' @param phenotype_file path of a two-column file containing used samples and Responder or No Responder 'Class' (NR,R).
#' @param final_signif_thresh Significance threshold for the rewiring method. The lower the threshold, the restrictive the method.
#' @param regulator_info_col_name Column name of the gene_info_file. By default, 'regulator'.
#' @param phenotype_col_name Column name of the phenotype_file. By default, 'Class'.
#' @param phenotype_class_vals_string Boolean terms of the phenotype_file at the 'Class' column. By default, (NR,R).
#'
#'
#' @return Return a list containing linker output file path, expression matrix path, geneinfo file path, phenotype file path,
#' regulator info column name, the phenotype column name, the phenotype class strings (NR,R by default), initial permutations
#' parameters, the retest threshold, the retest permutations and the significance threshold.
#'
#' @examples
#' \dontrun{
#' wd <- "mypath"
#' linker_saved_file<-c(paste(wd,'rds','su2c.rds',sep="/"))
#' expr_matrix_file<-c(paste(wd,'expression','su2c_expression_ARSI.txt',sep="/"))
#' gene_info_file<-c(paste(wd,'geneinfo','su2c_geneinfo.txt',sep="/"))
#' phenotype_file<-c(paste(wd,'clinical','su2c_clinical_ARSI.txt',sep="/"))
#' final_signif_thresh<-c(0.051)
#'
#' rewiringinput<-preparerewiring(name="example",linker_rds_v,expr_mat_v,gene_info_v,
#'                                phenotype_file_v,final_signif_thresh)
#' }
#' @export
preparerewiring<- function(name,linker_saved_file,
                           expr_matrix_file,
                           gene_info_file,
                           phenotype_file,
                           final_signif_thresh,
                           regulator_info_col_name='regulator',
                           phenotype_col_name='Class',
                           phenotype_class_vals_string='NR,R'){

  newdir<-paste(name,paste(final_signif_thresh,collapse="_"),sep="_")
  outdir<-paste(getwd(),newdir,sep="/")

  dir.create(file.path(outdir), showWarnings = FALSE)

  rewobjects<-list()

  for (i in 1:length(linker_saved_file)){

    rewobject<-list()

    phenotype_class_vals <- unlist(strsplit(phenotype_class_vals_string, ","))

    # read in linker output
    rundata <- readRDS(linker_saved_file[i]); #used outside
    methods::show('RDS done')

    # read in expression matrix
    input_expr_mat <- as.matrix(utils::read.table(expr_matrix_file[i], header = T,
                                           row.names = 1, sep = "\t", quote = ""))
    methods::show(paste(c("Expression Matrix Size", dim(input_expr_mat))))

    # read in gene info
    gene_info_df <- utils::read.table(gene_info_file[i], header = TRUE, sep = "\t", quote = "")
    rownames(gene_info_df) <- gene_info_df[, 1]
    methods::show(paste(c("Gene Info Table Size", dim(gene_info_df))))

    # find intersection of genes
    keepgenes <- intersect(rownames(input_expr_mat), rownames(gene_info_df))
    methods::show(paste(c("NumGenes Kept", length(keepgenes))))

    norm_expr_mat_keep <- input_expr_mat[keepgenes, ] #used outside
    #show(norm_expr_mat_keep)
    gene_info_df_keep <- gene_info_df[keepgenes, ]
    name2idx <- 1:dim(norm_expr_mat_keep)[1]
    names(name2idx) <- rownames(norm_expr_mat_keep)

    # divide genes into regulators and targets
    allregs <- keepgenes[which(gene_info_df_keep[, regulator_info_col_name] == 1)]
    alltargs <- keepgenes[which(gene_info_df_keep[, regulator_info_col_name] == 0)]
    methods::show(paste(c("NumRegs and NumTargs", length(allregs), length(alltargs))))

    # read in phenotype file
    pheno_df <- utils::read.table(phenotype_file[i], header = T, row.names = 1,
                           sep = "\t", quote = "", stringsAsFactors = F);
    methods::show(paste(c("Phenotype Table Size", dim(pheno_df))))

    # clean up phenotype column
    responder <- pheno_df[, phenotype_col_name]
    responder[which(pheno_df[, phenotype_col_name] == phenotype_class_vals[1])] <- 0
    responder[which(pheno_df[, phenotype_col_name] == phenotype_class_vals[2])] <- 1
    names(responder) <- make.names(rownames(pheno_df))
    pheno_df <- cbind(pheno_df, responder)
    #show(responder)
    methods::show(paste('responder: ',responder,sep=''))

    # find intesection of sample ids, keepsamps
    keepsamps <- intersect(colnames(norm_expr_mat_keep),
                           names(responder)[which(responder == 0 | responder == 1 )]) #usedoutside
    methods::show(keepsamps)
    # keeplabels is numeric class id, 0 or 1, in keep samps order
    keeplabels <- as.numeric(responder[keepsamps]) #used outside
    methods::show(paste('keeplabels: ',keeplabels,sep=""))

    class_counts <- as.numeric(table(keeplabels)) #used outside
    methods::show(paste(c("Class Per Counts", class_counts)))

    rewobject$'rundata'<-rundata
    rewobject$'norm_expr_mat_keep'<-norm_expr_mat_keep
    rewobject$'keepsamps'<-keepsamps
    rewobject$'keeplabels'<-keeplabels
    rewobject$'class_counts'<-class_counts
    rewobject$'final_signif_thresh'<-final_signif_thresh
    rewobject$'outdir'<-outdir

    rewobjects[[i]]<-rewobject
  }

  return (rewobjects)

}


