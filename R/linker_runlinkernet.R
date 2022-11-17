#' Single gene network approaches.
#'
#' `NET_run()` generate a single GRN. `NET_compute_graph_all_LASSO1se()` defines
#' the statistics of drivers and targets to be the Lasso method, choosing 1 standard error from the minimum RSS.
#' `NET_compute_graph_all_LASSOmin()` uses Lasso method, choosing the minimum RSS point. `NET_compute_graph_all_LM()`
#' uses a linear model and `NET_compute_graph_all_VBSR()` uses a Variational Bayes Spike Regression.
#'
#' @param lognorm_est_counts Matrix of log-normalized estimated counts of the gene expression
#' data (Nr Genes x Nr samples)
#' @param target_filtered_idx Index array of the target genes on the lognorm_est_counts matrix if
#' SummarizedExperiment object is not provided.
#' @param regulator_filtered_idx Index array of the regulatory genes on the lognorm_est_counts matrix if
#' SummarizedExperiment object is not provided.
#' @param nassay if SummarizedExperiment object is passed as input to lognorm_est_counts, name of the
#' assay containing the desired matrix. Default: 1 (first item in assay's list).
#' @param regulator if SummarizedExperiment object is passed as input to lognorm_est_counts, name of the
#' rowData() variable to build target_filtered_idx and regulator_filtered_idx. This variable must be one
#' for driver genes and zero for target genes. Default: 'regulator'
#' @param graph_mode Chosen method(s) to generate the edges in the bipartite graph. The available options are 'VBSR',
#' 'LASSOmin', 'LASSO1se' and 'LM'. By default, all methods are chosen.
#' @param FDR The False Discovery Rate correction used for the enrichment analysis. By default, 0.05.
#' @param NrCores Nr of computer cores for the parallel parts of the method. Note that
#' the parallelization is NOT initialized in any of the functions. By default, 3.
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
#'    ## We recommend VBSR (rest of parameters are set by default)
#'    graph <- NET_run(lognorm_est_counts,target_filtered_idx=target_filtered_idx,
#'                      regulator_filtered_idx=regulator_filtered_idx,graph_mode='VBSR')
#'
#' @export NET_run
NET_run <- function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, nassay = 1, regulator = "regulator", graph_mode = c("VBSR",
    "LASSOmin", "LASSO1se", "LM"), FDR = 0.05, NrCores = 1) {
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

    if (!is.numeric(lognorm_est_counts[1, 1]) & !is.numeric(lognorm_est_counts[1, 1])) {
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

    graphs <- list()
    for (j in seq_along(graph_mode)) {
        graphs[[graph_mode[j]]] <- switch(graph_mode[j], VBSR = NET_compute_graph_all_VBSR(lognorm_est_counts, regulator_filtered_idx,
            target_filtered_idx, NrCores), LASSOmin = NET_compute_graph_all_LASSOmin(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx,
            NrCores), LASSO1se = NET_compute_graph_all_LASSO1se(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx, NrCores),
            LM = NET_compute_graph_all_LM(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx, NrCores))

        message("Graphs for (", graph_mode[j], ") computed!")
    }

    return(list(graphs = graphs))
}
#' @export
#' @rdname NET_run
#' @param alpha feature selection parameter in case of a LASSO model to be chosen.
NET_compute_graph_all_LASSO1se <- function(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx, alpha = 1 - 1e-06, NrCores = 1) {

    X <- lognorm_est_counts[regulator_filtered_idx, ]

    driverMat <- matrix(data = NA, nrow = length(target_filtered_idx), ncol = length(regulator_filtered_idx))

    # This will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.

    parallClass <- BiocParallel::bpparam()
    parallClass$workers <- NrCores

    # compute the LASSO1se
    Lasso1se <- function(idx, lognorm_est_counts, alpha) {

        y <- lognorm_est_counts[idx, ]
        fit <- glmnet::cv.glmnet(t(X), y, alpha = alpha)

        b_o <- stats::coef(fit, s = fit$lambda.1se)
        b_opt <- c(b_o[2:length(b_o)])  # removing the intercept.

        b_opt

    }
    driverMat <- BiocParallel::bplapply(target_filtered_idx, Lasso1se, lognorm_est_counts, alpha, BPPARAM = parallClass)

    # Transform list of bettas into matrix, as .combine=rbind in foreach
    driverMat <- do.call(rbind, driverMat)

    rownames(driverMat) <- rownames(lognorm_est_counts)[target_filtered_idx]
    colnames(driverMat) <- rownames(lognorm_est_counts)[regulator_filtered_idx]

    regulated_genes <- which(rowSums(abs(driverMat)) != 0)
    regulatory_genes <- which(colSums(abs(driverMat)) != 0)
    driverMat <- driverMat[regulated_genes, ]
    driverMat <- driverMat[, regulatory_genes]

    g <- igraph::graph_from_incidence_matrix(driverMat)

    return(g)

}
#' @export
#' @rdname NET_run
NET_compute_graph_all_LASSOmin <- function(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx, alpha = 1 - 1e-06, NrCores = 1) {

    X <- lognorm_est_counts[regulator_filtered_idx, ]

    driverMat <- matrix(data = NA, nrow = length(target_filtered_idx), ncol = length(regulator_filtered_idx))

    # This will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.
    parallClass <- BiocParallel::bpparam()
    parallClass$workers <- NrCores

    # compute the LASSOmin
    Lassomin <- function(idx, lognorm_est_counts, alpha) {

        y <- lognorm_est_counts[idx, ]
        fit <- glmnet::cv.glmnet(t(X), y, alpha = alpha)

        b_o <- stats::coef(fit, s = fit$lambda.min)
        b_opt <- c(b_o[2:length(b_o)])  # removing the intercept.

        b_opt

    }
    driverMat <- BiocParallel::bplapply(target_filtered_idx, Lassomin, lognorm_est_counts, alpha, BPPARAM = parallClass)

    # Transform list of bettas into matrix, as .combine=rbind in foreach
    driverMat <- do.call(rbind, driverMat)

    rownames(driverMat) <- rownames(lognorm_est_counts)[target_filtered_idx]
    colnames(driverMat) <- rownames(lognorm_est_counts)[regulator_filtered_idx]

    regulated_genes <- which(rowSums(abs(driverMat)) != 0)
    regulatory_genes <- which(colSums(abs(driverMat)) != 0)
    driverMat <- driverMat[regulated_genes, ]
    driverMat <- driverMat[, regulatory_genes]

    g <- igraph::graph_from_incidence_matrix(driverMat)

    return(g)

}
#' @export
#' @rdname NET_run
NET_compute_graph_all_LM <- function(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx, NrCores = 1) {
    GEA_per_regulator <- list()
    i <- 1

    X <- t(lognorm_est_counts[regulator_filtered_idx, ])

    # compute the LM
    NrTotalEdges <- length(target_filtered_idx) * length(regulator_filtered_idx)
    Pthre <- 0.05/(NrTotalEdges)

    # This will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.
    parallClass <- BiocParallel::bpparam()
    parallClass$workers <- NrCores

    LM <- function(idx, lognorm_est_counts, Pthre) {

        y <- lognorm_est_counts[idx, ]
        driverVec <- numeric(length = ncol(X))

        for (i in seq_len(ncol(X))) {
            x <- X[, i]
            fit <- stats::lm(y ~ x)
            s <- summary(fit)
            driverVec[i] <- (s$coefficients[2, "Pr(>|t|)"] < Pthre)
        }
        driverVec
    }
    driverMat <- BiocParallel::bplapply(target_filtered_idx, LM, lognorm_est_counts, Pthre, BPPARAM = parallClass)

    # Transform list of bettas into matrix, as .combine=rbind in foreach
    driverMat <- do.call(rbind, driverMat)

    target_genes <- rownames(lognorm_est_counts)[target_filtered_idx]
    regulators <- rownames(lognorm_est_counts)[regulator_filtered_idx]

    rownames(driverMat) <- target_genes
    colnames(driverMat) <- regulators

    regulated_genes <- which(rowSums(abs(driverMat)) != 0)
    regulatory_genes <- which(colSums(abs(driverMat)) != 0)
    driverMat <- driverMat[regulated_genes, ]
    driverMat <- driverMat[, regulatory_genes]

    g <- igraph::graph_from_incidence_matrix(driverMat)

    return(g)

}
#' @export
#' @rdname NET_run
NET_compute_graph_all_VBSR <- function(lognorm_est_counts, regulator_filtered_idx, target_filtered_idx, NrCores = 1) {

    X <- lognorm_est_counts[regulator_filtered_idx, , drop = FALSE]

    driverMat <- matrix(data = NA, nrow = length(target_filtered_idx), ncol = length(regulator_filtered_idx))

    # This will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.

    parallClass <- BiocParallel::bpparam()
    parallClass$workers <- NrCores

    # compute the VBSR
    VBSR <- function(idx, lognorm_est_counts, X, regulator_filtered_idx, target_filtered_idx) {

        y <- lognorm_est_counts[idx, , drop = FALSE]
        if (nrow(y)) {
            y <- t(y)
        }
        if (length(unique(y)) == 1) {
            rep(0, length(regulator_filtered_idx))

        } else {
            res <- try(vbsr(y, t(X), n_orderings = 15, family = "normal"), silent = TRUE)

            if (inherits(res, "try-error")) {
                rep(0, length(regulator_filtered_idx))

            } else {
                betas <- res$beta
                betas[res$pval > 0.05/(length(target_filtered_idx) * length(regulator_filtered_idx))] <- 0
                betas
            }
        }
    }

    driverMat <- BiocParallel::bplapply(target_filtered_idx, VBSR, lognorm_est_counts, X, regulator_filtered_idx, target_filtered_idx,
        BPPARAM = parallClass)

    # Transform list of bettas into matrix, as .combine=rbind in foreach
    driverMat <- do.call(rbind, driverMat)

    rownames(driverMat) <- rownames(lognorm_est_counts)[target_filtered_idx]
    colnames(driverMat) <- rownames(lognorm_est_counts)[regulator_filtered_idx]

    regulated_genes <- which(rowSums(abs(driverMat)) != 0)
    regulatory_genes <- which(colSums(abs(driverMat)) != 0)

    driverMat <- driverMat[regulated_genes, ]
    driverMat <- driverMat[, regulatory_genes]

    g <- igraph::graph_from_incidence_matrix(driverMat, weighted = TRUE)

    return(g)

}
