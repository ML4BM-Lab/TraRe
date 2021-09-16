#' Cliques generation and plot displaying.
#'
#' By the way the LINKER method works, some highly-correlated driver genes (TFs) may be dropped from the
#' resultant model, as the role they play at the GRN inference process is very similar. Due to this,
#' we propose a method based on cliques (Fully Connected Networks) to recover those dropped drivers.
#'
#' First, `preparedata()` prepares the correlation matrix and the hash table
#' for future uses. Then, `generategraph()` generates an igraph object from genes correlation
#' (threshold dependant) adjacency matrix. After that, `selectmethod()` chooses method
#' in order to remove duplicities generated from the igraph::max_cliques method. Finally,
#' `generatecliques()` generates all the cliques, containing fully connected networks of genes
#' per clique.
#'
#' @param dataset input expression file with genes as rows and samples as columns.
#' @param nassay if SummarizedExperiment object is passed as input to dataset, name of the assay containing
#' the desired matrix. Default: 0
#' @param method method to use in the correlation matrix generation (see stats:cor). Default: 'pearson'
#' @param correlationth threshold to consider edge exists. Default: 0.6
#' @param sparsecorrmatrix boolean variable specifying whether to set to 0 values below threshold or not. Default: TRUE
#' @param numcliques number of cliques to be generated. Default: 'All'
#' @param mandatorygenes array of gene names which is mandatory to include in the returned cliques. Default: c() (none)
#' @param selection integer selecting method. The available options are: 1 - Maximize Genes/Clique,
#' 2 - Maximize Median Correlation Value/Clique, 3- Maximize Avg Variance Correlation Value/Clique or
#' 4 - Maximize Sum(option two, option three).
#' @param table_cliques boolean indicating to return all possible correlated genes instead of groups of cliques.
#' @return List containing the plot generated and a list with all the generated cliques.
#'
#' @examples
#'
#'    ## Suppose we want to recover the drivers LINKER may have dropped out.
#'    ## This method allows to group the highly correlated (above `correlationth`)
#'    ## driver genes so after the GRN generation, we know that if a particular
#'    ## driver gene is taking part of a relationship inside that GRN, all the
#'    ## genes inside this group may be also taking part of the same relationship
#'    ## due to the high correlation. Note that the large this threshold is, the
#'    ## surer we are about this affirmation.
#'
#'
#'    ## For this example, we will only work with the driver genes of a example dataset.
#'
#'    dataset <- readRDS(paste0(system.file('extdata',package='TraRe'),
#'                       '/tfs_linker_example_eg.rds'))
#'
#'    ## Lets select the generated dataset, as the rest of parameters are set by default.
#'
#'    clioutput <- generatecliques(dataset = dataset)
#'
#' @export generatecliques
generatecliques <- function(dataset = NULL, nassay = 1, method = "pearson", correlationth = 0.6,
                            sparsecorrmatrix = TRUE, numcliques = "All",
                            mandatorygenes = c(), selection = 1, table_cliques = TRUE) {

    # All but the last one to print it in a different way.
    # Check for SummarizedExperiment Object

    if (inherits(dataset, "SummarizedExperiment")) {
        dataset <- SummarizedExperiment::assays(dataset)[[nassay]]
    }

    # Unit Tests

    if (is.null(dataset)) {
        stop("dataset field empty")
    }

    if (!(is.matrix(dataset) | is.data.frame(dataset))) {
        stop("matrix or dataframe class is required")
    }

    if (!is.numeric(dataset[1, 1])) {
        stop("non-numeric values inside dataset variable")
    }

    message("Preparing data")
    pdobject <- preparedata(dataset, method, table_cliques, correlationth )

    #check if table_cliques
    if (table_cliques){return(pdobject)}

    message("Generating graph")
    ggobject <- generategraph(correlationth, sparsecorrmatrix, pdobject)

    message("Selecting method")
    smobject <- selectmethod(selection, ggobject, pdobject)

    message("Generate Datasets")
    ml <- smobject$cliques

    # check numcliques is not greater than the amount of modules we have.
    if (numcliques > length(ml) | is.character(numcliques))
        numcliques <- length(ml)

    # ml-> cliques. We assure max variance genes are getting in.
    sortparameter <- vapply(ml[], getMaxVariance, dat = pdobject, FUN.VALUE = 1)
    sortparameter_ix <- sort(sortparameter, decreasing = TRUE, index.return = TRUE)$ix
    sortparameter_ix_numcliques <- sortparameter_ix[seq_len(numcliques)]

    # sort it and select the numcliques
    ml <- ml[sortparameter_ix_numcliques]

    # Get the genes that we need and aren't in our selected ones.
    mandatorygenes_toinclude <- setdiff(mandatorygenes, unlist(ml))


    # Substitute the bottom cliques with the mandatory ones. We can substitute a clique with more than 1 gene for a single gene...

    for (i in seq(length(ml), 1)) {
        if (!length(mandatorygenes_toinclude))
            break

        # check we are not substituting a mandatory gene.
        if (!ml[[i]] %in% mandatorygenes) {

            ml[[i]] <- mandatorygenes_toinclude[1]  #character.pop(0)
            mandatorygenes_toinclude <- mandatorygenes_toinclude[-1]  #character.pop(0)

        }
    }

    # We build a dictionary, where Representatives are the keys and Cliques are the values.

    RC_list <- list()
    RC_list <- ml
    names(RC_list) <- vapply(RC_list, getMaxVarName, dat = pdobject, FUN.VALUE = "String")  #Representatives

    # Plot the selected method.

    SelectedCliquesPlot <- plotcliques(ml, pdobject, sortparameter_ix_numcliques, smobject, numcliques)

    return(list(plot = SelectedCliquesPlot, cliques = RC_list, representatives = names(RC_list)))
}
#' @export
#' @rdname generatecliques
preparedata <- function(dataset, method, table_cliques, correlationth) {

    # Generate variance gene dictionary.
    dataset_variance_dict <- matrixStats::rowVars(as.matrix(dataset))
    names(dataset_variance_dict) <- rownames(dataset)
    dataset_variance_dict <- hash::hash(dataset_variance_dict)

    # Generate the correlation matrix to work with.
    CorrMatrix <- abs(stats::cor(t(dataset), method = method))

    #
    if(table_cliques){

        #retrieve the rownames
        gene_n <- rownames(CorrMatrix)

        #for each gene, retrieve the genes correlated > TH
        gene_corr <- lapply(gene_n,function(gn){
            gene_m <- CorrMatrix[gn,]
            return(names(gene_m[gene_m>correlationth]))
        })
        names(gene_corr) <- gene_n
        #
        return (gene_corr)
    }


    diag(CorrMatrix) <- 0  #eliminate diagonal

    return(list(hash = dataset_variance_dict, mat = CorrMatrix))
}
#' @export
#' @rdname generatecliques
#' @param pdoutput output from preparedata().
generategraph <- function(correlationth, sparsecorrmatrix, pdoutput) {

    message("Generating groups of highly correlated genes and singleton communities")

    # Get the names of high correlated and non- high correlated

    highcorrg_n <- rownames(pdoutput$mat)[apply(pdoutput$mat, 1, function(x) Reduce("|", x > correlationth))]

    nothighcorrg_n <- setdiff(rownames(pdoutput$mat), highcorrg_n)

    # Generate the correlation matrix

    message("Creating the matrix")
    highcorg_g <- pdoutput$mat[highcorrg_n, highcorrg_n]

    # generate sparse matrix
    if (sparsecorrmatrix) {
        highcorg_g[highcorg_g < correlationth] <- 0
    }

    diag(highcorg_g) <- 0  #eliminate diagonal
    # generate graph
    highcorg_g <- igraph::graph_from_adjacency_matrix(highcorg_g, mode = "undirected", weighted = TRUE)

    return(list(graph = highcorg_g, singleton = nothighcorrg_n))

}
#' @export
#' @rdname generatecliques
#' @param ggoutput output from generategraph()
#' @param pdoutput output from preparedata()
selectmethod <- function(selection, ggoutput, pdoutput) {

    # Run igraph::max_cliques() method

    # Cliques<-vapply(igraph::max_cliques(ggoutput$graph)[],names,FUN.VALUE = c(''))
    Cliques <- sapply(igraph::max_cliques(ggoutput$graph)[], names)

    SortingMethod <- switch(selection, length, function(x) getMedianCorrelation(x, dat = pdoutput), function(x) getAvgVariance(x, dat = pdoutput),
        function(x) getMedian_AvgVar(x, dat = pdoutput))

    SortingMethod_tittle <- switch(selection, "Maximizing Genes/Clique", "Maximize Median Correlation Value/Clique", "Maximize Avg Variance Correlation Value/Clique",
        "Maximize Sum(2,3)/Clique")

    SortingMethod_chosen <- vapply(Cliques, SortingMethod, FUN.VALUE = 1)

    # We sort the cliques
    SortedCliques <- sort(SortingMethod_chosen, decreasing = TRUE, index.return = TRUE)$ix

    noduplcliques <- RemoveDuplicities(SortedCliques, Cliques)

    L <- length(noduplcliques)  #We only count the ones with high correlated genes.

    noduplcliques_r <- append(noduplcliques, ggoutput$singleton)  #+SingletonCommunities

    return(list(cliques = noduplcliques_r, tittle = SortingMethod_tittle, L = L))

}
#' @export
#' @rdname generatecliques
#' @param ml cliques from smobject
#' @param pdobject output from preparedata()
#' @param sortparameter_ix_numcliques parameter for sorting cliques.
#' @param smobject output from selectmethod()
plotcliques <- function(ml, pdobject, sortparameter_ix_numcliques, smobject, numcliques) {

    x_axis <- vapply(ml[], getAvgVariance, dat = pdobject, FUN.VALUE = 1)
    y_axis <- vapply(ml[], getMedianCorrelation, dat = pdobject, FUN.VALUE = 1)  #Median Correlation

    Palette <- c("#29bd00", "#ff04d1", "#1714b0", "#b30009", "#030303")

    col_fact <- c()

    mll <- vapply(ml[], length, FUN.VALUE = 1)

    for (i in seq_along(ml)) {

        clique_l <- length(ml[[i]])

        ch_cat <- which(c((clique_l == 1) & (sortparameter_ix_numcliques[i] < smobject$L), (clique_l > 1) & (clique_l <= 4), (clique_l >
            4) & (clique_l <= 7), (clique_l > 7), (sortparameter_ix_numcliques[i] >= smobject$L)))

        col_fact <- c(col_fact, switch(ch_cat, "x == 1", "1 < x <= 4", "4 < x <= 7", "x > 7", "Singleton"))

    }
    col_fact <- as.factor(col_fact)

    SelectedCliquesPlot <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(y_axis, x_axis, color = col_fact), size = 3) + ggplot2::labs(x = "Median Correlation Value",
        y = "Avg Variance Per Clique", title = smobject$tittle, subtitle = paste("TotalNum: ", toString(sum(mll)), ", MedianOfMedian: ",
            toString(round(stats::median(y_axis[y_axis != 1]), 4)), ", Total Singleton Cliques: ", toString(as.double(table(vapply(ml,
                length, FUN.VALUE = 1) == 1)[2])), ", Number of Cliques: ", toString(numcliques), sep = "")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        size = 15, face = "bold"), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12)) + ggplot2::theme(legend.text = ggplot2::element_text(size = 15),
        axis.text = ggplot2::element_text(size = 15), axis.title = ggplot2::element_text(size = 15), legend.title = ggplot2::element_text(size = 15)) +
        ggplot2::scale_color_manual(name = "Genes/Cliques", values = c(`x == 1` = Palette[1], `1 < x <= 4` = Palette[2], `4 < x <= 7` = Palette[3],
            `x > 7` = Palette[4], Singleton = Palette[5]))

    return(SelectedCliquesPlot)
}

#---- Helpers ----

getMaxVariance <- function(clique, dat) {
    return(max(vapply(unlist(clique), function(x) dat$hash[[x]], FUN.VALUE = 1)))
}

getMaxVarName <- function(clique, dat) {
    return(names(which.max(vapply(unlist(clique), function(x) dat$hash[[x]], FUN.VALUE = 1))))
}

getAvgVariance <- function(clique, dat) {
    return(mean(vapply(clique, function(x) dat$hash[[x]], FUN.VALUE = 1)))
}

getMedian_AvgVar <- function(clique, dat) {
    return(getMedianCorrelation(clique) + getAvgVariance(clique, dat))
}

getMedianCorrelation <- function(clique, dat) {

    if (!length(clique))
        return(0) else if (length(clique) == 1)
        return(1)

    # CorrMatrix will change if drivers or targets are run.
    return(stats::median(apply(utils::combn(which(rownames(dat$mat) %in% clique), 2), 2, function(x) dat$mat[x[1], x[2]])))

}

RemoveDuplicities <- function(patron, cliques) {

    mcliques <- cliques[patron]  #modified cliques

    unmcliques <- unlist(mcliques)
    mcliques <- Map("[", mcliques, utils::relist(!duplicated(unmcliques), skeleton = mcliques))

    mcliques <- Filter(Negate(function(X) {
        length(X) == 0
    }), mcliques)

    return(mcliques)

}

