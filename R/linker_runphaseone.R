#' Phase I : module generation
#'
#' Run first phase of the linker method where K modules of similarly expressed target genes and
#' relate them to a linear combination of very few regulators, according to the selected model. `LINKER_init()`
#' evaluate kmeans on a train set to generate a initial set of clusters containing drivers and target genes.
#' `LINKER_ReassignGenesToClusters()` reassigning genes based on closed match to new regulatory programs.
#' This functions takes place inside the linkerrun function, so it is not recommended to run it on its own.
#' `LINKER_corrClust()` go through two steps within a loop, learning regulatory program of modules and reassigning
#' genes. `LINKER_extract_modules()` extract all the modules, genes and relevant information. `LINKER_EvaluateTestSet()`
#' fits the selected model with the test data. `LINKER_LearnRegulatoryPrograms()` learns the regulatory program of the modules.
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
#' @param mode Chosen method(s) to link module eigengenes to regulators. The available options are
#' 'VBSR', 'LASSOmin', 'LASSO1se', 'LASSOparam' and 'LM'. Default set to 'VBSR'
#' @param used_method Method selected for use. Default set to MEAN.
#' @param Nr_bootstraps Number of bootstrap of Phase I. By default, 1.
#' @param NrCores Nr of computer cores for the parallel parts of the method. Note that the parallelization
#' is NOT initialized in any of the functions. By default, 2.
#' @param Lambda Lambda variable for Lasso models.
#' @param alpha Alpha variable for Lasso models.
#' @param pmax Maximum numbers of regulators that we want.
#' @param corrClustNrIter Number of iteration for the phase I part of the method.
#' @param FDR The False Discovery Rate correction used for the modules and graphs GRN uncovering. By default, 0.05.
#'
#' @return igraph object containing the modules containing the related drivers and targets within bootstraps.
#'
#' @examples
#'
#'    ## This example is very similar to the `LINKER_run()` function.
#'    ## Again, we are going to join drivers and targets genes to create the working dataset.
#'
#'    drivers <- readRDS(paste0(system.file('extdata',package='TraRe'),'/tfs_linker_example.rds'))
#'    targets <- readRDS(paste0(system.file('extdata',package='TraRe'),'/targets_linker_example.rds'))
#'
#'    lognorm_est_counts <- rbind(drivers,targets)
#'    ## We create the index for drivers and targets in the log-normalized gene expression matrix.
#'
#'    R<-60
#'    T<-200
#'
#'    regulator_filtered_idx <- seq_len(R)
#'    target_filtered_idx <- R+c(seq_len(T))
#'
#'    ## We recommend to use the default values of the function.
#'    ## For the sake of time, we will select faster (and worse) ones.
#'
#'    linkeroutput <- LINKER_runPhase1(lognorm_est_counts,target_filtered_idx=target_filtered_idx,
#'                                     regulator_filtered_idx=regulator_filtered_idx, NrModules=2,
#'                                     mode='LASSOmin',NrCores=2, corrClustNrIter=10,Nr_bootstraps=1)
#'
#' @export LINKER_runPhase1

LINKER_runPhase1 <- function(lognorm_est_counts, target_filtered_idx, regulator_filtered_idx, nassay = 1, regulator = "regulator",
    NrModules, Lambda = 5, alpha = 1 - 1e-06, pmax = 10, mode = "VBSR", used_method = "MEAN", NrCores = 1, corrClustNrIter = 100,
    Nr_bootstraps = 1, FDR = 0.05) {

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

    # Creating the parameters structure
    Parameters <- list(Lambda = Lambda, pmax = pmax, alpha = alpha, mode = mode, used_method = used_method)
    sample_size <- dim(lognorm_est_counts)[2]
    train_size <- round(0.8 * sample_size)
    EvaluateTestSet <- list()
    bootstrap_modules <- list()
    bootstrap_results <- list()


    for (boost_idx in seq_len(Nr_bootstraps)) {

        train_samples <- sample(seq_len(sample_size), train_size, replace = FALSE)
        validation_samples <- setdiff(seq_len(sample_size), train_samples)

        # Scale by driver genes
        Regulator_data_train <- t(scale(t(lognorm_est_counts[regulator_filtered_idx, train_samples])))
        Regulator_data_validation <- t(scale(t(lognorm_est_counts[regulator_filtered_idx, validation_samples])))

        # Scale by samples within targets
        MA_matrix_Var_train <- t(scale(t(lognorm_est_counts[target_filtered_idx, train_samples])))
        MA_matrix_Var_validation <- t(scale(t(lognorm_est_counts[target_filtered_idx, validation_samples])))

        LINKERinit <- LINKER_init(MA_matrix_Var = MA_matrix_Var_train, RegulatorData = Regulator_data_train, NrModules = NrModules,
            NrCores = NrCores, corrClustNrIter = corrClustNrIter, Parameters = Parameters, FDR = FDR)

        tmp <- LINKER_corrClust(LINKERinit)
        bootstrap_results[[boost_idx]] <- tmp

        EvaluateTestSet[[boost_idx]] <- LINKER_EvaluateTestSet(bootstrap_results[[boost_idx]], MA_matrix_Var_validation, Regulator_data_validation,
            used_method = Parameters$used_method)

        R.utils::printf("Bootstrap %d, NrModules %d:\n", boost_idx, bootstrap_results[[boost_idx]]$NrModules)

        print(matrixStats::colMeans2(EvaluateTestSet[[boost_idx]]))
    }
    return(list(bootstrapResults = bootstrap_results, bootstrapTestStats = EvaluateTestSet))

}
#' @export
#' @rdname LINKER_runPhase1
#' @param MA_matrix_Var Matrix of log-normalized estimated counts of the gene expression data, centered and scaled, containing
#' only the train samples.
#' @param RegulatorData Expression matrix containing only the regulators of the train samples.
#' @param Parameters List of parameters containig lambda, pmax, alpha, mode and used method.
#' @param NrModules Number of modules that are a priori to be found (note that the final number of modules
#' discovered may differ from this value). By default, 100 modules.
LINKER_init <- function(MA_matrix_Var, RegulatorData, NrModules, NrCores = 3, corrClustNrIter = 21, Parameters, FDR) {

    if (nrow(MA_matrix_Var) > NrModules) {

        # K-means++-style initialization
        rnd_cent <- matrix()

        # First one at random
        rnd_cent <- (MA_matrix_Var[sample(seq_len(nrow(MA_matrix_Var)), 1), ])
        if (NrModules > 0) {

            if (NrModules > 1) {

                # the second value is the least correlated from the first one
                corr_dists <- apply(MA_matrix_Var, 1, function(x) ((x %*% rnd_cent)/(length(x) - 1))^2)
                rnd_cent <- cbind(as.matrix(rnd_cent), as.matrix(MA_matrix_Var[which(corr_dists == min(corr_dists))[1], ]))
                rnd_cent <- t(rnd_cent)
                if (NrModules > 2) {

                  # compute the rest of the centers
                  Px <- numeric(length = nrow(MA_matrix_Var))
                  for (center_idx in 3:NrModules) {

                    for (i in seq_len(nrow(MA_matrix_Var))) {
                      gene <- MA_matrix_Var[i, ]
                      corr_dists <- apply(rnd_cent, 1, function(x) ((x %*% gene)/(length(x) - 1))^2)
                      Px[i] <- 2 - max(corr_dists)
                    }
                    if (sum(is.finite(Px)) != length(Px)) {
                      message("asdf")
                    }
                    rnd_cent <- rbind(rnd_cent, MA_matrix_Var[sample(seq_len(nrow(MA_matrix_Var)), 1, prob = Px), ])
                  }
                }
            }

        }


        ModuleVectors <- rnd_cent

        Data <- MA_matrix_Var
        Clusters <- numeric()
        for (jj in seq_len(5)) {
            for (i in seq_len(nrow(MA_matrix_Var))) {
                CurrentGeneVector <- Data[i, , drop = FALSE]
                Correlations <- (stats::cor(t(CurrentGeneVector), t(ModuleVectors)))
                corr <- data.matrix(Correlations, rownames.force = NA)
                MaxCorrelation <- max(corr, na.rm = TRUE)
                MaxPosition <- which(signif(corr, digits = 7) == signif(MaxCorrelation, digits = 7))
                MaxPosition <- MaxPosition[1]  # this is new, to avoid two different reassignements

                Clusters[i] <- MaxPosition
            }
            ClusterIDs <- unique(Clusters)
            for (idx in seq_along(ClusterIDs)) {
                genesInModule <- which(Clusters == ClusterIDs[idx])
                cx <- MA_matrix_Var[genesInModule, ]
                if (length(genesInModule) > 1) {
                  if (Parameters$used_method == "LINKER") {
                    clusterSVD <- svd(cx)
                    y <- clusterSVD$v[, 1]
                    if (stats::cor(matrixStats::colMeans2(cx), y) < 0) {
                      y <- -y
                    }

                  } else {
                    ModuleVectors[idx, ] <- matrixStats::colMeans2(cx)
                  }
                } else {
                  ModuleVectors[idx, ] <- cx
                }

            }
        }

    } else {
        stop("The number of modules is too large compared to the total number of genes.")
    }
    ModuleMembership <- as.numeric(Clusters)
    names(ModuleMembership) <- rownames(MA_matrix_Var)

    return(list(MA_matrix_Var = MA_matrix_Var, RegulatorData = RegulatorData, ModuleMembership = ModuleMembership, Parameters = Parameters,
        NrCores = NrCores, corrClustNrIter = corrClustNrIter, FDR = FDR))

}
#' @export
#' @rdname LINKER_runPhase1
#' @param Data Data Matrix of log-normalized estimated counts of the gene expression data, centered and scaled, containing
#' only the train samples.
#' @param Beta Coefficient on which the decision of reassigning genes is based.
#' @param Clusters Number of modules that are a priori to be found (note that the final number of modules discovered may differ from this value).
LINKER_ReassignGenesToClusters <- function(Data, RegulatorData, Beta, Clusters, NrCores = 1) {

    MIN_NUM_GENES_PER_MODULE <- 2

    RegulatorData_rownames <- rownames(RegulatorData)
    Data_rownames <- rownames(Data)

    NrGenes <- nrow(Data)
    NrSamples <- ncol(Data)
    NrReassignGenes <- 0

    ## reassigning genes based on the Beta getting the predictor data
    X <- RegulatorData
    # creating the cluster 'centroids'
    X1 <- data.matrix(X)
    ModuleVectors <- Beta %*% X1
    GeneNames <- rownames(Data)


    # This will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.
    parallClass <- BiocParallel::bpparam()
    parallClass$workers <- NrCores

    # reassigning genes and generate new clusters:
    ReassignGenes <- function(i, Data) {

        OldModule <- Clusters[i]
        CurrentGeneVector <- Data[i, , drop = FALSE]
        Correlations <- (stats::cor(t(CurrentGeneVector), t(ModuleVectors)))

        corr <- data.matrix(Correlations, rownames.force = NA)
        MaxCorrelation <- max(corr, na.rm = TRUE)
        MaxPosition <- which(signif(corr, digits = 7) == signif(MaxCorrelation, digits = 7))
        MaxPosition <- MaxPosition[1]  # this is new, to avoid two different reassignements

        if (MaxPosition != OldModule) {
            NrReassignGenes <- NrReassignGenes + 1
        }
        NewClusters <- MaxPosition

        return(NewClusters)

    }
    nc <- BiocParallel::bplapply(seq_len(NrGenes), ReassignGenes, Data, BPPARAM = parallClass)

    # Transform list of bettas into matrix, as .combine=c in foreach
    nc <- do.call(c, nc)


    # Remove cluster with too few genes. Avoids singularities. Could be solved imposing priors. future work
    for (i in unique(nc)) {
        NrGenesInCluster <- sum(nc == i)

        if (NrGenesInCluster < MIN_NUM_GENES_PER_MODULE) {
            # I need to reassign these genes
            genesInModule <- which(nc == i)
            # remove the cluster
            ModuleVectors[i, ] <- 0

            for (j in genesInModule) {
                CurrentGeneVector <- Data[j, , drop = FALSE]
                # Correlations = abs(cor(t(CurrentGeneVector),t(ModuleVectors)))
                Correlations <- (stats::cor(t(CurrentGeneVector), t(ModuleVectors)))

                corr <- data.matrix(Correlations, rownames.force = NA)
                MaxCorrelation <- max(corr, na.rm = TRUE)
                MaxPosition <- which(signif(corr, digits = 7) == signif(MaxCorrelation, digits = 7))
                MaxPosition <- MaxPosition[1]  # this is new, to avoid two different reassignements

                nc[j] <- MaxPosition

            }
        }
    }

    NrReassignGenes <- length(which(nc != Clusters))
    result <- list(NrReassignGenes = NrReassignGenes, Clusters = nc)
    return(result)
}
#' @export
#' @rdname LINKER_runPhase1
#' @param LINKERinit Initialization object obtained from `LINKER_init()`.
LINKER_corrClust <- function(LINKERinit) {


    NrIterations <- LINKERinit$corrClustNrIter


    if (nrow(LINKERinit$RegulatorData) == 1) {
        stop("Only one driver is detected. More than one driver is needed.\n")
    }

    Data <- LINKERinit$MA_matrix_Var
    Clusters <- LINKERinit$ModuleMembership
    RegulatorData <- LINKERinit$RegulatorData
    Parameters <- LINKERinit$Parameters
    NrCores <- LINKERinit$NrCores

    RegulatorData_rownames <- rownames(RegulatorData)
    Data_rownames <- rownames(Data)

    # main loop We want to end with the regulatory programs, hence we loop for Reassign, regulatory

    # STEP 1: learning the regulatory program for each cluster
    regulatoryPrograms <- LINKER_LearnRegulatoryPrograms(Data, Clusters, RegulatorData, Lambda = Parameters$Lambda, alpha = Parameters$alpha,
        pmax = Parameters$pmax, mode = Parameters$mode, used_method = Parameters$used_method, NrCores = NrCores, FDR = LINKERinit$FDR)

    jj <- 1
    while (jj < NrIterations) {

        # STEP 2: reassigning genes based on closed match to new regulatory programs
        ReassignGenesToClusters <- LINKER_ReassignGenesToClusters(Data, RegulatorData, regulatoryPrograms$Beta, Clusters, NrCores = NrCores)
        jj <- jj + 1

        NrReassignGenes <- ReassignGenesToClusters$NrReassignGenes
        Clusters <- ReassignGenesToClusters$Clusters

        # STEP 1: learning the regulatory program for each cluster
        regulatoryPrograms <- LINKER_LearnRegulatoryPrograms(Data, Clusters, RegulatorData, Lambda = Parameters$Lambda, alpha = Parameters$alpha,
            pmax = Parameters$pmax, mode = Parameters$mode, used_method = Parameters$used_method, NrCores = NrCores, FDR = LINKERinit$FDR)

    }

    # update results structure
    ModuleMembership <- as.matrix(Clusters)
    rownames(ModuleMembership) <- rownames(Data)
    colnames(ModuleMembership) <- c("ModuleNr")

    result <- list(NrModules = length(unique(Clusters)), RegulatoryPrograms = regulatoryPrograms$Beta, AllRegulators = rownames(RegulatorData),
        AllGenes = rownames(Data), ModuleMembership = ModuleMembership)

    training_stats <- LINKER_EvaluateTestSet(result, Data, RegulatorData, used_method = LINKERinit$Parameters$used_method)

    result <- list(NrModules = length(unique(Clusters)), RegulatoryPrograms = regulatoryPrograms$Beta, AllRegulators = rownames(RegulatorData),
        AllGenes = rownames(Data), ModuleMembership = ModuleMembership, trainingStats = training_stats)

    return(result)
}
#' @export
#' @rdname LINKER_runPhase1
#' @param results Matrix of log-normalized estimated counts of the gene expression data (Nr Genes x Nr samples).
LINKER_extract_modules <- function(results) {

    modules <- list()

    enriched_idx <- 1
    NrBootstraps <- length(results$bootstrapResult)
    for (idx_bootstrap in seq_len(NrBootstraps)) {

        NrModules <- results$bootstrapResult[[idx_bootstrap]]$NrModules
        boot_results <- results$bootstrapResult[[idx_bootstrap]]
        boot_idx <- sort(unique(boot_results$ModuleMembership[, ]))

        for (Module_number in seq_len(NrModules)) {


            Module_target_genes_full_name <- boot_results$AllGenes[which(boot_results$ModuleMembership[, ] == boot_idx[Module_number])]
            Module_target_gene_list <- vapply(Module_target_genes_full_name, function(x) strsplit(x, "\\|"), FUN.VALUE = list(""))
            if (length(Module_target_gene_list[[1]]) == 1) {
                # NOT FULL GENECODE NAME, ONLY ONE NAME PER GENE!
                Module_target_genes <- Module_target_genes_full_name
            } else {
                # FULL GENECODE ANNOTATION!
                Module_target_genes <- vapply(Module_target_gene_list, function(x) x[[6]], FUN.VALUE = "")
                Module_target_genes <- unname(Module_target_genes)
            }


            Modules_regulators_full_name <- names(which(boot_results$RegulatoryPrograms[Module_number, ] != 0))
            if (length(Modules_regulators_full_name) == 0) {
                next
            }
            Modules_regulators_list <- vapply(Modules_regulators_full_name, function(x) strsplit(x, "\\|"), FUN.VALUE = list(""))
            if (length(Modules_regulators_list[[1]]) == 1) {
                # NOT FULL GENECODE NAME, ONLY ONE NAME PER GENE!
                Modules_regulators <- Modules_regulators_full_name
            } else {
                # FULL GENECODE ANNOTATION!
                Modules_regulators <- vapply(Modules_regulators_list, function(x) x[[6]], FUN.VALUE = "")
                Modules_regulators <- unname(Modules_regulators)
            }



            modules[[enriched_idx]] <- list(target_genes = Module_target_genes_full_name, regulators = Modules_regulators_full_name,
                regulatory_program = boot_results$RegulatoryPrograms[Module_number, ], training_stats = boot_results$trainingStats[Module_number,
                  ], test_stats = results$bootstrapTestStats[[idx_bootstrap]][Module_number], assigned_genes = which(boot_results$ModuleMembership[,
                  ] == Module_number), bootstrap_idx = idx_bootstrap)
            enriched_idx <- enriched_idx + 1
        }
    }
    return(modules)

}
#' @export
#' @rdname LINKER_runPhase1
#' @param LINKERresults List containing the number of clusters, regulatoryprogram, name of regulators and all genes and module membership.
#' @param MA_Data_TestSet Matrix of log-normalized estimated counts of the gene expression data, centered and scaled, containing
#' only the test samples.
#' @param RegulatorData_TestSet Expression matrix containing only the regulators of the test samples.
LINKER_EvaluateTestSet <- function(LINKERresults, MA_Data_TestSet, RegulatorData_TestSet, used_method = "MEAN") {
    nrSamples <- ncol(MA_Data_TestSet)
    RegulatorNames <- rownames(RegulatorData_TestSet)

    # Iterating over the Modules
    stats <- mat.or.vec(LINKERresults$NrModules, 7)
    Rsquare <- mat.or.vec(LINKERresults$NrModules, 1)
    RsquareAjusted <- mat.or.vec(LINKERresults$NrModules, 1)
    modules <- list()

    for (i in seq_len(LINKERresults$NrModules)) {
        # check regulator presence
        currentRegulators <- RegulatorNames[which(LINKERresults$RegulatoryPrograms[i, ] != 0)]
        stats[i, 1] <- length(currentRegulators)

        # checking the presence of the clusters
        currentClusterGenes <- LINKERresults$AllGenes[which(LINKERresults$ModuleMembership[, 1] == i)]
        stats[i, 2] <- length(currentClusterGenes)

        # predict cluster expression in test set, always calculate but report the totel percentage weight that is represented
        currentWeights <- LINKERresults$RegulatoryPrograms[i, which(LINKERresults$RegulatoryPrograms[i, ] != 0)]

        modules[[i]] <- currentClusterGenes[currentClusterGenes %in% rownames(MA_Data_TestSet)]

        # drop=FALSE, this solves the problem when you have only one regulator, so the previous version is not needed.
        predictions <- (t(RegulatorData_TestSet[currentRegulators, , drop = FALSE])) %*% (currentWeights)  # need to make sure that the first argument remains a matrix.
        predictions <- data.matrix(predictions)
        if (length(modules[[i]]) != 0) {
            if (length(currentClusterGenes) > 1) {
                cx <- MA_Data_TestSet[currentClusterGenes, ]
                module_SVD <- svd(cx)

                if (used_method == "LINKER") {
                  outcome <- module_SVD$v[, 1]
                  if (stats::cor(predictions, outcome) < 0) {
                    outcome <- -outcome
                  }
                } else {
                  outcome <- matrixStats::colMeans2(cx)
                }

                varEx <- module_SVD$d[1]^2/sum(module_SVD$d^2)
            } else {
                outcome <- MA_Data_TestSet[currentClusterGenes, ]
                varEx <- 0
            }

            module_data <- MA_Data_TestSet[currentClusterGenes, ]
            if (nrow(t(module_data)) == 1) {
                module_data <- t(module_data)
            }
            inmodule_corr <- abs(stats::cor(t(module_data), outcome))
            meanInModuleCor <- mean(inmodule_corr)

            homogeneity <- abs(stats::cor(t(module_data), t(module_data)))
            homogeneity <- (sum(homogeneity) - dim(homogeneity)[1])/((dim(homogeneity)[1] - 1) * (dim(homogeneity)[1] - 1))

            # using explained variance as metric, since mean square error is not enough, no baseline interpretation possible

            SStot <- sum((outcome - mean(outcome))^2)
            SSres <- sum((predictions - outcome)^2)
            Rsquare <- 1 - (SSres/SStot)
            RsquareAjusted <- 1 - (1 - Rsquare) * ((nrSamples - 1)/(nrSamples - length(currentRegulators)))

            stats[i, 3] <- meanInModuleCor
            stats[i, 4] <- Rsquare
            stats[i, 5] <- RsquareAjusted
            stats[i, 6] <- homogeneity
            stats[i, 7] <- varEx
        } else {

        }
    }
    dimnames(stats) <- list(rownames(stats, do.NULL = FALSE, prefix = "Module_"), c("nrReg", "nrGen", "MeanInModuleCorr", "Rsquare",
        "RsquareAdjusted", "homogeneity", "condition"))

    return(stats)
}
#' @export
#' @rdname LINKER_runPhase1
#' @param Data Matrix of log-normalized estimated counts of the gene expression data, centered and scaled, containing
#' only the train samples.
#' @param Clusters Clusters generated from the linkerinit function.
LINKER_LearnRegulatoryPrograms <- function(Data, Clusters, RegulatorData, Lambda, alpha, pmax, mode, used_method = "MEAN", NrCores = 1, FDR) {

    RegulatorData_rownames <- rownames(RegulatorData)
    Data_rownames <- rownames(Data)
    NrClusters <- length(unique(Clusters))
    NrGenes <- nrow(Data)
    NrSamples <- ncol(Data)

    y_all <- mat.or.vec(NrClusters, NrSamples)
    ClusterIDs <- unique(Clusters)
    ClusterIDs <- sort(ClusterIDs, decreasing = FALSE)
    cnt <- seq_len(NrClusters)

    # This will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.

    parallClass <- BiocParallel::bpparam()
    parallClass$workers <- NrCores

    # Calculate Bettas
    CalcBettas <- function(i, Data, Clusters, RegulatorData, Lambda, alpha, mode, used_method) {

        CurrentClusterPositions <- which(Clusters %in% ClusterIDs[i])
        nrGenesInClusters <- length(CurrentClusterPositions)

        if (nrGenesInClusters > 1) {
            if (used_method == "LINKER") {
                gene_svd <- svd(Data[CurrentClusterPositions, ])
                y <- gene_svd$v[, 1]
                if (stats::cor(matrixStats::colMeans2(Data[CurrentClusterPositions, ]), y) < 0) {
                  y <- -y
                }
            } else {
                y <- matrixStats::colMeans2(Data[CurrentClusterPositions, ])
            }
        } else {
            y <- Data[CurrentClusterPositions, ]
        }
        X <- RegulatorData

        if (mode == "LASSOmin") {

            fit <- glmnet::cv.glmnet(t(X), y, alpha = alpha)

            nonZeroLambdas <- fit$lambda[which(fit$nzero > 0)]
            nonZeroCVMs <- fit$cvm[which(fit$nzero > 0)]

            if (length(which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))) == 0) {

                # for now: just print a warning, *although* this error WILL cause LINKER to crash in a few steps.
                warnMessage <- paste0("\nOn cluster ", i, " there were no cv.glm results that gave non-zero coefficients.")
                warning(warnMessage)
                bestNonZeroLambda <- fit$lambda.min

            } else {
                bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))]
            }
            b_o <- stats::coef(fit, s = fit$lambda.min)
            b_opt <- c(b_o[2:length(b_o)])  # removing the intercept.

        } else if (mode == "LASSO1se") {

            fit <- glmnet::cv.glmnet(t(X), y, alpha = alpha)

            nonZeroLambdas <- fit$lambda[which(fit$nzero > 0)]
            nonZeroCVMs <- fit$cvm[which(fit$nzero > 0)]

            if (length(which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))) == 0) {

                # for now: just print a warning, *although* this error WILL cause LINKER to crash in a few steps.
                warnMessage <- paste0("\nOn cluster ", i, " there were no cv.glm results that gave non-zero coefficients.")
                warning(warnMessage)
                bestNonZeroLambda <- fit$lambda.min

            } else {
                bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))]
            }
            b_o <- stats::coef(fit, s = fit$lambda.1se)
            b_opt <- c(b_o[2:length(b_o)])  # removing the intercept.

        } else if (mode == 'LASSOparam'){

            fit <- glmnet::cv.glmnet(t(X), y, alpha = alpha)

            fitquants <- stats::quantile(fit$lambda,probs=seq(0,1,1/9))

            nonZeroLambdas <- fit$lambda[which(fit$nzero > 0)]
            nonZeroCVMs <- fit$cvm[which(fit$nzero > 0)]

            if (length(which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))) == 0) {

                # for now: just print a warning, *although* this error WILL cause LINKER to crash in a few steps.
                warnMessage <- paste0("\nOn cluster ", i, " there were no cv.glm results that gave non-zero coefficients.")
                warning(warnMessage)
                bestNonZeroLambda <- fit$lambda.min

            } else {
                bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))]
            }

            b_o <- stats::coef(fit, s = fitquants[Lambda])
            b_opt <- c(b_o[2:length(b_o)])  # removing the intercept.

        } else if (mode == "VBSR") {

            res <- vbsr::vbsr(y, t(X), n_orderings = 15, family = "normal")
            betas <- res$beta
            max_beta <- max(abs(betas))
            betas[res$pval > FDR/(nrow(RegulatorData) * nrow(Data))] <- 0  #Bonferroni
            b_opt <- betas
            b_o <- 0

        } else if (mode == "LM") {
            b_opt <- numeric(length = nrow(RegulatorData))

            for (i in seq_len(nrow(RegulatorData))) {
                x <- X[i, ]
                fit <- stats::lm(y ~ x)
                s <- summary(fit)
                if (s$coefficients[2, "Pr(>|t|)"] < FDR/(nrow(RegulatorData) * nrow(Data))) {
                  b_opt[i] <- s$coefficients[2, 1]
                } else {
                  b_opt[i] <- 0
                }
            }
            b_o <- 0

        } else {
            message("MODE NOT RECOGNIZED")
        }

        return(list(b_opt, y, b_o[1]))

    }

    # Initialize
    BetaY_all <- BiocParallel::bplapply(seq_len(NrClusters), CalcBettas, Data, Clusters, RegulatorData, Lambda, alpha, mode, used_method,
        BPPARAM = parallClass)

    # Transform list of bettas into matrix, as .combine=cbind in foreach
    BetaY_all <- do.call(cbind, BetaY_all)

    Beta <- do.call(cbind, BetaY_all[1, ])
    Beta <- t(Beta)
    colnames(Beta) <- RegulatorData_rownames
    rownames(Beta) <- gsub("result.", "Module_", rownames(Beta))

    y_all <- do.call(cbind, BetaY_all[2, ])
    y_all <- t(y_all)
    rownames(y_all) <- gsub("result.", "Module_", rownames(y_all))

    intercept <- do.call(cbind, BetaY_all[3, ])
    intercept <- t(intercept)
    intercept <- as.numeric(intercept)

    # calculating some statistics

    prediction <- (Beta %*% RegulatorData + intercept)
    error <- y_all - prediction

    result <- list(Beta = Beta, error = error)
    return(result)
}



