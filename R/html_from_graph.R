#' Html generation
#'
#' From input Gene Regulatory Network, an html is generated containing a table with
#' driver to target phenotype dependent relationships. It is a brief summary containing
#' drivers normalized xor sum, which is the ratio between drivers-targets connections
#' present only in one of both phenotypes over all possible connections, and cliques, which are driver genes
#' that are highly correlated (over a user-decision threshold), and may have been lost during
#' the fitting process of the LINKER method.
#'
#' This functions takes place inside the runrewiring function, so it is not recommended to run it on its own.
#'
#'
#' @param gpath path to the graph object ('refinedsumm.rds'). (RDS format)
#' @param wpath writing path, where the html and txts file will be saved. (Default: temp directory)
#' @param user_mode boolean indicating if this function is called from user or internaly. (Default: TRUE)
#' @param cliquesbool indicating if cliques method should be added to the summary table. (Default: TRUE)
#' @param ... every argument you should pass to generatecliques() in case cliquesbool is TRUE.
#'
#' @return Html and txts containing the mentioned files.
#'
#' @examples
#'
#' ## For this example, we are going to use a generated 'refinedsumm.rds' from the su2c dataset
#' ## (see vignette of TraRe for more information), which is at the external data folder
#' ## of this package.
#'
#' gpath <- paste0(system.file('extdata',package='TraRe'),'/refinedsumm.rds')
#' wpath <- paste0(system.file('extdata',package='TraRe'))
#'
#' ## We are going to use the drivers dataset from the external data folder as well.
#' ## For more information about generatecliques() please check the help page of it.
#'
#' dataset<-readRDS(paste0(system.file('extdata',package='TraRe')
#'                  ,'/tfs_linker_example_eg.rds'))
#' html_from_graph(gpath=gpath,wpath=wpath,dataset=dataset)
#'
#'
#' @export

html_from_graph <- function(gpath = NULL, wpath = paste0(tempdir()), user_mode = TRUE, cliquesbool = TRUE, ...) {

    if (is.null(gpath)) {
        stop("Path to the graph object must be specified")
    }

    if (inherits(try(summary(gpath)$class, TRUE), "try-error")) {
        # check for url object check for weburl check if local url exists

        if (inherits(try(url(gpath), TRUE), "try-error")) {
            if (!file.exists(gpath)) {
                stop("graph object in the specified folder must exist")
            }
        }
    }

    if (!is.logical(cliquesbool)) {
        stop("non-logical variable pass to this cliquesbool argument")
    }

    # List of nodes - XOR-------------------------------------------------------

    # Import graph objects
    graphs <- readRDS(gpath)
    graph_R <- graphs$respond_graph
    graph_NR <- graphs$nonresp_graph

    # Get the edges
    edges_R <- igraph::as_edgelist(graph_R)
    edges_NR <- igraph::as_edgelist(graph_NR)

    # Get the list of nodes in both groups

    # Join driver+target (Responder) and sorting alphabeticaly
    list_edges_NR <- sort(paste(edges_NR[, 2], edges_NR[, 1]))
    # Join driver+target (Non Responder) and sorting alphabeticaly
    list_edges_R <- sort(paste(edges_R[, 2], edges_R[, 1]))
    list_final <- union(list_edges_NR, list_edges_R)

    # Preallocation
    L <- length(list_final)
    R <- rep(0, L)
    NR <- rep(0, L)

    # Fill the table
    NR[list_final %in% list_edges_NR] <- 1
    R[list_final %in% list_edges_R] <- 1

    signature <- as.integer(xor(NR, R))

    # Create Data Frame
    df <- data.frame(do.call(rbind, strsplit(list_final, " ")), NR, R, signature)
    names(df) <- c("Drivers", "Targets", "NR", "R", "XOR")

    # Summary table -----------------------------------------------------------

    # Sum of XOR
    drivers_xor <- split(df$XOR, df$Drivers)
    drivers_xor_sum <- vapply(drivers_xor, sum, FUN.VALUE = c(1))

    # Normalization
    normalized_sum <- round(drivers_xor_sum/lengths(drivers_xor))

    if (cliquesbool) {

        CliquesObject <- generatecliques(...)

        # Add cliques
        Driv <- names(drivers_xor)

        Cliques <- vector("list", length(Driv))
        names(Cliques) <- Driv

        for (x in CliquesObject$cliques) {
            for (y in Driv[Driv %in% x]) {

                Cliques[[y]] <- setdiff(x, y)

            }
        }

        Cliques <- vapply(Cliques, function(x) paste(x, collapse = ", "), FUN.VALUE = c("C"))

        # Reorder (in normalized_sum they are in alphabetical order, in Cliques as they come out from the cliques object)
        Cliques <- Cliques[order(names(Cliques))]


        # Create Data frame with cliques
        sumXOR <- data.frame(normalized_sum, Cliques, check.names = TRUE, check.rows = TRUE)
        colnames(sumXOR) <- c("Normalized XOR Sum", "Cliques")

    } else {

        # Create Data frame without cliques
        sumXOR <- data.frame(normalized_sum, check.names = TRUE, check.rows = TRUE)
        colnames(sumXOR) <- c("Normalized XOR Sum")

    }

    ## Export to txt and html

    # check if this function is called from runrewiring or from user

    if (!user_mode) {
        txts_path <- paste0(wpath, "/txts")
        htmls_path <- paste0(wpath, "/htmls")
    } else {
        # Create directory
        dir.create(wpath)
        # Assign paths
        txts_path <- wpath
        htmls_path <- wpath
    }

    ## Save as .txt

    # Df
    download_df <- paste0(txts_path, "/refined_supermodule_edges.txt")
    utils::write.table(df, file = download_df, sep = "\t", quote = FALSE)

    # SumXOR
    download_sumXOR <- paste0(txts_path, "/refined_supermodule_summary.txt")
    utils::write.table(sumXOR, file = download_sumXOR, sep = "\t", quote = FALSE)

    # Generate path
    path_html <- paste0(htmls_path, "/refined_edges_table.html")

    # Generate html from df and sumXOR
    df_to_html(df, sumXOR, path_html = path_html, download_df = download_df, download_sumXOR = download_sumXOR)

}


# Helper function

df_to_html <- function(df, sumXOR, path_html = "C:/Users/Jesus/Desktop/refined.html", download_df = "", download_sumXOR = "") {


    # Define containers
    write(paste0("<style>\n.floatLeft { width: 25%; float: left; margin-left: 200px;}", "\n.floatRight { width: 25%; float: right; margin-right: 250px;}",
        "\n.container { overflow: hidden; }\n</style>"), file = path_html)


    # First the df file
    write(paste0("<table border = 1 width = '100%'><tr bgcolor = ", "'#AAAAAA'><th><a href = '", download_df, "' target=", "'_blank'>Download Rewiring Edges</a></th><th><a href = '",
        download_sumXOR, "' target='_blank'>Download Summary</a></th></tr></table><br><br>"), file = path_html, append = TRUE)

    write(paste0("<div class='container'>\n", "<div class='floatLeft'>\n", "<table border = 1>"), file = path_html, append = TRUE)

    htmlstr <- table2html_from_graph(df)

    write(htmlstr, file = path_html, append = TRUE)
    write("</div>", file = path_html, append = TRUE)

    # Now the summary
    write(paste0("<div class='floatRight'>\n", "<table border = 1>"), file = path_html, append = TRUE)

    htmlstr <- table2html_from_graph(sumXOR, FALSE)

    write(htmlstr, file = path_html, append = TRUE)
    write("</div>", file = path_html, append = TRUE)

}
