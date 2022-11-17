#' Prepare rewiring data for running the method.
#'
#' Prepare neccessary files for running `runrewiring()`
#'
#' @param name Desired name of the folder which is generated. The chosen
#' threshold will be `paste()` to the folder's name.
#' @param linker_output_p Output file from linker function path. RDS format is required.
#' @param lognorm_est_counts_p Lognorm counts of the gene expression matrix path.
#' @param SEObject_p SummarizedExperiment objects path.
#' @param gene_info_p Path of a two-column file containing genes and 'regulator' boolean variable.
#' @param phenotype_p Path of a two-column file containing used samples and Responder or No Responder 'Class' (NR,R).
#' @param nassays Name of assays in case SummarizedObject is provided.
#' @param final_signif_thresh Significance threshold for the rewiring method. The lower the threshold, the restrictive the method.
#' @param regulator_info_col_name Column name of the gene_info_p. By default, 'regulator'.
#' @param phenotype_col_name Column name of the phenotype_p. By default, 'Class'.
#' @param phenotype_class_vals_string Boolean terms of the phenotype_p at the 'Class' column. By default, (NR,R).
#' @param phenotype_class_vals_string_label Boolean terms of the phenotype_p values at the 'Class' column. By default (0,1).
#' @param orig_test_perms Initial permutations for first test (default: 100) .
#' @param retest_thresh Threshold if a second test is performed (default: 0.08) .
#' @param retest_perms Permutations if a second test is performed (default: 1000) .
#' @param use_graphs Boolean specifying the use of graphs as linkeroutput modules. (TRUE by default).
#' @param outdir Directory for the output folder to be located (default: tempdir())
#' @param nrcores Number of cores to run the parallelization within the rewiring test (default: 3).
#' @param last_cluster Boolean specifying whether to include the last_cluster in the rewiring or not. (default: FALSE)
#'
#'
#' @return Return a list containing: LINKER's output, expression matrix, boolean array from phenotype file,
#' array containing number of c(R,NR) samples, significance threshold and output directory.
#'
#' @examples
#'
#' ## We are going to prepare 4 files that we have in the example folder: the output from LINKER, the
#' ## gene expression matrix, the phenotype file and the gene info file. Note that the LINKER
#' ## output is generated from the gene expression matrix. Note: if rewiring across more than 1 dataset
#' ## is desired, paths will be given as arrays. (i.e. linker_output <- c(path1,path2))
#'
#'
#' linker_output_p <- paste0(system.file('extdata',package='TraRe'),'/linker_rewiring_example.rds')
#'
#' lognorm_est_counts_p <- paste0(system.file('extdata',package='TraRe'),
#'                         '/expression_rewiring_example.txt')
#'
#' gene_info_p <- paste0(system.file('extdata',package='TraRe'),'/geneinfo_rewiring_example.txt')
#'
#' phenotype_p <- paste0(system.file('extdata',package='TraRe'),'/phenotype_rewiring_example.txt')
#'
#' outdir <- system.file('extdata',package='TraRe')
#'
#' prepared <- preparerewiring(name='example',linker_output_p=linker_output_p,
#'                             lognorm_est_counts_p=lognorm_est_counts_p,gene_info_p=gene_info_p,
#'                             phenotype_p=phenotype_p,final_signif_thresh=0.05,
#'                             nrcores=1,outdir=outdir)
#'
#'
#' @export
preparerewiring <- function(name = "defaultname", linker_output_p, lognorm_est_counts_p = NULL, SEObject_p = NULL, gene_info_p = NULL,
    phenotype_p = NULL, nassays = 1, final_signif_thresh = 0.001, regulator_info_col_name = "regulator", phenotype_col_name = "Class",
    phenotype_class_vals_string = "NR,R", phenotype_class_vals_string_label = "0,1", orig_test_perms = 100, retest_thresh = 0.08, retest_perms = 1000,
    use_graphs = TRUE, outdir = tempdir(), nrcores = 3, last_cluster = FALSE) {

    # checks

    if (is.null(linker_output_p)) {
        stop("linker_output_p field required")
    }

    # check if SummarizedObject has been provided.
    se <- TRUE
    if (is.null(SEObject_p)) {

        se <- FALSE

        if (is.null(lognorm_est_counts_p)) {
            stop("lognorm_est_counts_p field required")
        }

        if (is.null(gene_info_p)) {
            stop("gene_info_p field required")
        }
        if (is.null(phenotype_p)) {
            stop("phenotype_p field required")
        }

    }

    # Check for comparison mode

    if (length(linker_output_p) > 1) {
        warning("Data comparison mode selected, only heatmap will be generated.")
    }

    # Create folder name
    foldername <- paste(name, paste(final_signif_thresh, collapse = "_"), sep = "_")

    # Concatenate with outdir path
    outdir <- paste(outdir, foldername, sep = "/")

    rewobjects <- list()
    rewobjects$datasets <- list()

    for (i in seq_along(linker_output_p)) {

        rewobject <- list()

        phenotype_class_vals <- unlist(strsplit(phenotype_class_vals_string, ","))
        phenotype_class_vals_label <- unlist(strsplit(phenotype_class_vals_string_label, ","))


        if (se) {

            # check for SummarizedExperiment
            seobject <- readRDS(SEObject_p[i])
            input_expr_mat <- SummarizedExperiment::assays(seobject)[[nassays[i]]]

        } else {

            input_expr_mat <- as.matrix(utils::read.delim(lognorm_est_counts_p[i],header=TRUE,row.names=1))
        }

        ExpMatSize <- c("Expression Matrix Size: ", "(", paste(dim(input_expr_mat), collapse = ","), ")")
        message(ExpMatSize)

        # read in gene info
        if (se) {

            gene_info_df <- SummarizedExperiment::rowData(seobject)

        } else {

            gene_info_df <- utils::read.delim(gene_info_p[i])

        }

        rownames(gene_info_df) <- gene_info_df[, 1]

        GeneInfoTableSize <- c("Gene Info Table Size: ", "(", paste(dim(gene_info_df), collapse = ","), ")")
        message(GeneInfoTableSize)

        # read in linker output
        rundata <- readRDS(linker_output_p[i])

        if (use_graphs){
            #Format linkeroutput modules to graphs
            message("\nFrom here on, graphs will be taken into account for Rewiring and some of them may be dropped out.\nHence, it is recommended to take preparerewiring object's run data to proceed with further analysis\n")
            rundata <- graph_to_modules(rundata,gene_info_df,regulator_info_col_name)
        }

        # find intersection of genes
        keepgenes <- intersect(rownames(input_expr_mat), rownames(gene_info_df))

        NumGenesKept <- c("NumGenes Kept: ", length(keepgenes))
        message(NumGenesKept)


        norm_expr_mat_keep <- input_expr_mat[keepgenes, ]
        gene_info_df_keep <- gene_info_df[keepgenes, ]
        name2idx <- seq_len(nrow(norm_expr_mat_keep))
        names(name2idx) <- rownames(norm_expr_mat_keep)

        # divide genes into regulators and targets
        allregs <- keepgenes[which(gene_info_df_keep[, regulator_info_col_name] == 1)]
        alltargs <- keepgenes[which(gene_info_df_keep[, regulator_info_col_name] == 0)]

        NumRegsTargs <- c("NumRegs and NumTargs: [", length(allregs), ",", length(alltargs), "]")
        message(NumRegsTargs)

        # read in phenotype file
        if (se) {

            pheno_df <- SummarizedExperiment::colData(seobject)

        } else {

            pheno_df <- utils::read.delim(phenotype_p[i], row.names = 1)

        }

        PhenTableSize <- c("Phenotype Table Size: ", "(", paste(dim(pheno_df), collapse = ","), ")")
        message(PhenTableSize)

        # clean up phenotype column
        responder <- pheno_df[, phenotype_col_name]
        responder[which(pheno_df[, phenotype_col_name] == phenotype_class_vals[1])] <- 0
        responder[which(pheno_df[, phenotype_col_name] == phenotype_class_vals[2])] <- 1
        names(responder) <- make.names(rownames(pheno_df))
        pheno_df <- cbind(pheno_df, responder)

        # find intersection of sample ids, keepsamps
        keepsamps <- intersect(colnames(norm_expr_mat_keep), names(responder)[which(responder == 0 | responder == 1)])

        # format for the log file.
        keepsamps_s <- split(keepsamps, cut(seq_along(keepsamps), max(length(keepsamps)/3, 2), labels = FALSE))

        keepsamps_s <- paste0(vapply(keepsamps_s, FUN = paste0, collapse = "|", FUN.VALUE = ""), collapse = "\n")

        SampleNames <- c("\nSample Names: [", keepsamps_s, "]\n\nNumber of samples: ", length(keepsamps))
        message(SampleNames)

        # find intersection of expression matrix and phenotype
        norm_expr_mat_keep <- norm_expr_mat_keep[, keepsamps]
        FiltExpMat <- c("Filtered exp matrix: ", "(", paste(dim(norm_expr_mat_keep), collapse = ","), ")")
        message(FiltExpMat)

        # keeplabels is numeric class id, 0 or 1, in keep samps order
        keeplabels <- as.numeric(responder[keepsamps])  #used outside

        # check if NR/R proportions are similar to ensure property functioning of the method.
        klzero <- sum(keeplabels == 0)
        klone <- sum(keeplabels == 1)

        if (min(klone, klzero)/max(klone, klzero) < 0.8) {
            warning(paste0("phenotype samples proportions imbalance ", toString(c(klzero, klone)), " (<80%)."))
        }

        class_counts <- as.numeric(table(keeplabels))

        ClassPerCounts <- c("Class Per Counts: ", "(", paste(class_counts, collapse = ","), ")")
        message(ClassPerCounts)


        rewobject$rundata <- rundata
        rewobject$norm_expr_mat_keep <- norm_expr_mat_keep
        rewobject$keepsamps <- keepsamps
        rewobject$keeplabels <- keeplabels
        rewobject$class_counts <- class_counts
        rewobject$final_signif_thresh <- final_signif_thresh
        rewobject$responder <- responder
        rewobject$gene_info_df_keep <- gene_info_df_keep
        rewobject$name2idx <- name2idx
        rewobject$allregs <- allregs
        rewobject$alltargs <- alltargs

        rewobjects$datasets[[i]] <- rewobject
    }

    rewobjects$regulator_info_col_name <- regulator_info_col_name
    rewobjects$phenotype_class_vals <- phenotype_class_vals
    rewobjects$phenotype_class_vals_label <- phenotype_class_vals_label
    rewobjects$outdir <- outdir
    rewobjects$orig_test_perms <- orig_test_perms
    rewobjects$retest_thresh <- retest_thresh
    rewobjects$retest_perms <- retest_perms
    rewobjects$NrCores <- nrcores
    rewobjects$last_cluster <- last_cluster

    # Create logfile

    logfile <- list(ExpMatSize, GeneInfoTableSize, NumGenesKept, NumRegsTargs, PhenTableSize, SampleNames, FiltExpMat, ClassPerCounts)
    logfile <- vapply(logfile, FUN = paste0, collapse = "", FUN.VALUE = "")

    rewobjects$logfile <- logfile

    return(rewobjects)

}

# Helper function
graph_to_modules <- function(linkeroutput, geneinfo, regulator_info_col_name){

    ## The structure we want to get is
    ## linkeroutput$modules[[link_mode]][[graph_mode]][[num_module]]$(target or regs)

    #load geneinfo_TFs
    gene_info_TFs <- rownames(geneinfo)[geneinfo[,regulator_info_col_name]==1]

    #names of linkeroutput
    linkeroutput_names <- names(linkeroutput$modules)

    linkeroutput <- lapply(linkeroutput_names,function(x){

        graph_modes <- names(linkeroutput$graphs[[x]])

        selected <- graph_modes[1]

        ## Select VBSR if more than 1 is available

        if (length(graph_modes)>1 & 'VBSR'%in%graph_modes){
            selected <- 'VBSR'
        }

        module_list <- lapply(seq_along(linkeroutput$graphs[[x]][[selected]]),function(y){

            graph <- linkeroutput$graphs[[x]][[selected]][[y]]

            totgenes <- unique(names(igraph::V(graph)))

            regulators <- intersect(totgenes,gene_info_TFs)

            if (identical(regulators,character(0))){

                message('Module number ',y,' has been deleted')
                return(NULL)

            }

            list(regulators = regulators,
                 target_genes = intersect(totgenes, rownames(geneinfo)[geneinfo[,regulator_info_col_name]==0]),
                 bootstrap_idx = linkeroutput$modules[[x]][[y]]$bootstrap_idx)

        })

        #generate old index
        orig_index <- unlist(sapply(seq_along(module_list),function(x){

            if (!is.null(module_list[[x]])){return(x)}

        }))


        module_list <- Filter(Negate(function(X) {

            length(X) == 0

        }),module_list)

        return(list(modules=module_list,graphs=linkeroutput$graphs[[x]][[selected]][orig_index]))

    })

    names(linkeroutput) <- linkeroutput_names

    #initialize the final linkeroutut
    new_linkeroutput <- list(modules=list(),graphs=list())

    for (method in names(linkeroutput)){

        #assign the modules
        new_linkeroutput[['modules']][[method]] <- linkeroutput[[method]][['modules']]
        #assign the graphs
        new_linkeroutput[['graphs']][[method]] <- linkeroutput[[method]][['graphs']]
    }

    return(new_linkeroutput)

}
