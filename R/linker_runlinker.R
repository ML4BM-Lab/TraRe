#' GRN inference via selected model.
#'
#' Gene Regulatory Network inference via model selection. Consists of two phases,
#' `LINKER_runPhase1()` and `LINKER_runPhase2()`. Help them for more information.
#'
#' @param lognorm_est_counts Matrix of log-normalized estimated counts of the gene expression
#' data (Nr Genes x Nr samples) or SummarizedExperiment object.
#' @param target_filtered_idx Index array of the target genes on the lognorm_est_counts matrix if
#' SummarizedExperiment object is not provided.
#' @param regulator_filtered_idx Index array of the regulatory genes on the lognorm_est_counts matrix if
#' SummarizedExperiment object is not provided.
#' @param nassay if SummarizedExperiment object is passed as input to lognorm_est_counts, name of the
#' assay containing the desired matrix. Default: 1 (first item in assay's list).
#' @param regulator if SummarizedExperiment object is passed as input to lognorm_est_counts, name of the
#' rowData() variable to build target_filtered_idx and regulator_filtered_idx. This variable must be one
#' for driver genes and zero for target genes. Default: 'regulator'
#' @param link_mode Chosen method(s) to link module eigengenes to regulators. The available options are
#' 'VBSR', 'LASSOmin', 'LASSO1se' and 'LM'. By default, all methods are chosen.
#' @param graph_mode Chosen method(s) to generate the edges in the bipartite graph. The available options
#' are 'VBSR', 'LASSOmin', 'LASSO1se' and 'LM'. By default, all methods are chosen.
#' @param module_rep Method selected for use. Default set to MEAN.
#' @param NrModules Number of modules that are a priori to be found (note that the final number of modules
#' discovered may differ from this value). By default, 100 modules.
#' @param corrClustNrIter output from preparedata(). By default, 100.
#' @param Nr_bootstraps Number of bootstrap of Phase I. By default, 10.
#' @param FDR The False Discovery Rate correction used for the modules and graphs GRN uncovering. By default, 0.05.
#' @param Lambda Lambda variable for Lasso models.
#' @param NrCores Nr of computer cores for the parallel parts of the method. Note that the parallelization
#' is NOT initialized in any of the functions. By default, 2.
#' @param onlymods Whether to infer only modules or modules and graphs. Default: FALSE
#' @param only_train whether to use only training samples within LINKER run. Default: FALSE
#'
#'
#' @return List containing the GRN raw results, GRN modules and GRN graphs.
#'
#' @examples
#'    ## For this example, we are going to join 60 drivers and
#'    ## 200 targets genes from the example folder.
#'
#'    drivers <- readRDS(paste0(system.file('extdata',package='TraRe'),'/tfs_linker_example.rds'))
#'    targets <- readRDS(paste0(system.file('extdata',package='TraRe'),'/targets_linker_example.rds'))
#'
#'    lognorm_est_counts <- as.matrix(rbind(drivers,targets))
#'
#'    ## We create the index for drivers and targets in the log-normalized gene expression matrix.
#'
#'    R<-60
#'    T<-200
#'
#'    regulator_filtered_idx <- seq_len(R)
#'    target_filtered_idx <- R+c(seq_len(T))
#'
#'
#'    ## We recommend to use the default values of the function.
#'    ## For the sake of time, we will select faster (and worse) ones.
#'
#'    linkeroutput <- LINKER_run(lognorm_est_counts,target_filtered_idx=target_filtered_idx,
#'                               regulator_filtered_idx=regulator_filtered_idx,
#'                               link_mode='LASSOmin',graph_mode='LM',NrModules=5,Nr_bootstraps=1,
#'                                NrCores=2,corrClustNrIter=10)
#'
#'
#'
#' @export
LINKER_run <- function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, nassay = 1, regulator = "regulator", link_mode = c("VBSR",
    "LASSOmin", "LASSO1se", "LM"), graph_mode = c("VBSR", "LASSOmin", "LASSO1se", "LM"), module_rep = "MEAN", NrModules = 100, corrClustNrIter = 100,
    Nr_bootstraps = 10, FDR = 0.05, Lambda = 5, NrCores = 1, onlymods = FALSE, only_train=FALSE) {

    # Check for SummarizedExperiment Object

    if (inherits(lognorm_est_counts, "SummarizedExperiment")) {

        # Generate target and regulator indexes
        genenames <- rownames(lognorm_est_counts)
        geneinfo <- SummarizedExperiment::rowData(lognorm_est_counts)

        target_filtered_idx <- which(genenames %in% rownames(geneinfo)[geneinfo$regulator == 0])
        regulator_filtered_idx <- which(genenames %in% rownames(geneinfo)[geneinfo$regulator == 1])

        # Get lognorm_est_counts expression matrix.
        lognorm_est_counts <- SummarizedExperiment::assays(lognorm_est_counts)[[nassay]]

    }

    # checks for lognorm_est_counts

    if (is.null(lognorm_est_counts)) {
        stop("lognorm_est_counts field empty")
    }

    if (!(is.matrix(lognorm_est_counts))) {
        stop("matrix class is required for input dataset")
    }

    if (!is.numeric(lognorm_est_counts[1, 1]) & !is.integer(lognorm_est_counts[1, 1])) {
        stop("non-numeric values inside lognorm_est_counts variable")
    }

    if (is.null(rownames(lognorm_est_counts)) | is.null(colnames(lognorm_est_counts))) {
        stop("null field detected in row names or column names, check lognorm_est_counts matrix")
    }

    # checks for target and regulator filtered index

    if (is.null(target_filtered_idx)) {
        stop("target_filtered_idx field empty")
    }

    if (is.null(regulator_filtered_idx)) {
        stop("regulator_filtered_idx field empty")
    }

    if (length(target_filtered_idx) + length(regulator_filtered_idx) != nrow(lognorm_est_counts)) {
        stop("the total number of genes is not equal to the sum of target_filtered_idx and regulatory_filtered_idx lengths")
    }

    if (!(is.numeric(target_filtered_idx) & is.numeric(regulator_filtered_idx))) {
        stop("targets and regulators index arrays must be numeric")
    }

    # Link modes

    res <- lapply(seq_along(link_mode), function(x) {
        LINKER_runPhase1(lognorm_est_counts = lognorm_est_counts, target_filtered_idx = target_filtered_idx, regulator_filtered_idx = regulator_filtered_idx,
            NrModules = NrModules, NrCores = NrCores, mode = link_mode[x], used_method = module_rep, corrClustNrIter = corrClustNrIter,
            Nr_bootstraps = Nr_bootstraps, FDR = FDR, Lambda = Lambda, only_train=only_train)
    })

    names(res) <- link_mode

    modules <- lapply(seq_along(link_mode), function(x) {
        LINKER_extract_modules(res[[x]])
    })

    names(modules) <- link_mode

    message("Link modes completed!")

    # Graphs
    if (!onlymods){
        graphs <- lapply(link_mode, function(x) {

            lapply(graph_mode, function(y) {

                LINKER_runPhase2(modules[[x]], lognorm_est_counts, mode = y, NrCores = NrCores, FDR = FDR)

            })
        })

        names(graphs) <- link_mode  #Set names for link_mode
        graphs <- lapply(graphs, function(x) stats::setNames(x, graph_mode))  #Set names for graph_mode

        message("Graphs computed!")

        return(list(raw_results = res, modules = modules, graphs = graphs))
    }

    return(list(raw_results = res, modules = modules))
}
