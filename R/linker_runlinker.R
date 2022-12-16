#' GRN inference via selected model.
#'
#' Gene Regulatory Network inference via model selection. Consists of two phases,
#' `LINKER_runPhase1()` and `LINKER_runPhase2()`. Help them for more information.
#'
#' @param TraReObj TraReObj containing preprocessed input matrix, linker_preprocessing output.
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
#' @param train_size Fraction of samples selected for the train samples. Default: 0.8.
#' @param onlymods Whether to infer only modules or modules and graphs. Default: FALSE
#' @param only_train whether to use only training samples within LINKER run. Default: FALSE
#'
#'
#' @return List containing the GRN raw results, GRN modules and GRN graphs.
#'
#' @examples
#' 
#'   ## For this example, we are going to load a example matrix
#'   lognorm_est_counts_p <- paste0(system.file('extdata',package='TraRe'),
#'                                  '/expression_rewiring_example.txt')
#'   lognorm_est_counts <- as.matrix(read.delim(lognorm_est_counts_p, header=TRUE,row.names=1))
#' 
#'   ## Load gene info, its an array of regulators' names.
#'   gene_info_p <- paste0(system.file('extdata',package='TraRe'),
#'                          '/geneinfo_rewiring_example.txt')
#'   gene_info <- read.delim(gene_info_p,header=TRUE)
#'   geneinfo <- gene_info[gene_info[,'regulator'] == 1,'uniq_isos']
#'
#'    ##TraReObj <- trare_preprocessing(data_matrix = lognorm_est_counts,
#'    ##                                  geneinfo = geneinfo, verbose = FALSE)
#'
#'    ## linker_output <- LINKER_run(TraReObj = TraReObj, link_mode='VBSR',
#'    ##                        graph_mode='VBSR',NrModules=100,Nr_bootstraps=10,
#'    ##                        corrClustNrIter=100)
#'
#' @export
LINKER_run <- function(TraReObj, link_mode = c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                       graph_mode = c("VBSR", "LASSOmin", "LASSO1se", "LM"),
                       module_rep = "MEAN", NrModules = 100, corrClustNrIter = 100,
                       Nr_bootstraps = 10, FDR = 0.05, Lambda = 5, train_size = 0.8,
                       onlymods = FALSE, only_train=FALSE) {

    # get lognorm_est_counts, target_filtered_idx, regulator_filtered_idx
    lognorm_est_counts <- TraReObj@lognorm_counts
    target_filtered_idx <- TraReObj@target_idx
    regulator_filtered_idx <- TraReObj@regulator_idx

    # checks

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
            NrModules = NrModules, NrCores = NrCores, train_size = train_size, mode = link_mode[x], used_method = module_rep, corrClustNrIter = corrClustNrIter,
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
