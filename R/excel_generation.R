#' Excel generation
#'
#' From input Gene Regulatory Network, an excel is generated containing a table with
#' driver to target phenotype dependent relationships. Brief summary containing
#' drivers normalized xor sum, which is the ratio between drivers-targets connections
#' present only in one of both phenotypes over all possible connections, and cliques, which are driver genes
#' that are highly correlated (over a user-decision threshold), and may have been lost during
#' the fitting process of the LINKER method.
#'
#'
#' @param gpath path to the graph object ('refinedsumm.rds'). (RDS format)
#' @param wpath writing path, where the excel file ('grnsumm.xlx') will be saved. (Default: working directory)
#' @param cliquesbool indicating if cliques method should be added to the summary table. (Default: TRUE)
#' @param ... every argument you should pass to generatecliques() in case cliquesbool is TRUE.
#'
#' @return Excel containing the mentioned parameters.
#'
#' @examples
#'
#' ## For this example, we are going to use a generated 'refinedsumm.rds' from the su2c dataset
#' ## (see vignette of TraRe for more information), which is at the external data folder
#' ## of this package.
#'
#' gpath <- paste0(system.file('extdata',package='TraRe'),'/refinedsumm.rds')
#' wpath <- system.file('extdata',package='TraRe')
#'
#' ## We are going to use the drivers dataset from the external data folder as well.
#' ## For more information about generatecliques() please check the help page of it.
#'
#' dataset<-readRDS(paste0(wpath,'/tfs_linker_example_eg.rds'))
#' excel_generation(gpath=gpath,wpath=wpath,dataset=dataset)
#'
#'
#' @export

excel_generation <- function(gpath = NULL, wpath = getwd(), cliquesbool=TRUE, ...){

  if (is.null(gpath)){
    stop('Path to the graph object must be specified')
  }
  if (!any(class(gpath)==c('url','connection'))){
    if (!file.exists(gpath)){
      stop('graph object in the specified folder must exist')
    }
  }
  if (!is.logical(cliquesbool)){
    stop('non-logical variable pass to this cliquesbool argument')
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
  list_edges_NR <- sort(paste(edges_NR[,2], edges_NR[,1]))
  # Join driver+target (Non Responder) and sorting alphabeticaly
  list_edges_R <- sort(paste(edges_R[,2], edges_R[,1]))
  list_final <- union(list_edges_NR,list_edges_R)

  # Preallocation
  L <- length(list_final)
  R <- rep(0,L)
  NR <- rep(0,L)

  # Fill the table
  for(i in seq_along(list_final)){
    if (list_final[i]%in%list_edges_NR){
      NR[i] <- 1
    }
    if(list_final[i]%in%list_edges_R){
      R[i] <- 1
    }
  }
  signature <- as.integer(xor(NR,R))

  # Create Data Frame
  excel <- data.frame(do.call(rbind,strsplit(list_final,' ')),NR,R, signature)
  names(excel) <- c("Drivers", "Targets", "NR", "R", "XOR")

  # Summary table -----------------------------------------------------------

  # Sum of XOR
  drivers_xor <- split(excel$XOR,excel$Drivers)
  drivers_xor_sum <- vapply(drivers_xor,sum,FUN.VALUE = c(1))

  # Normalization
  normalized_sum <- vapply(names(drivers_xor),
                           function(x) round(drivers_xor_sum[x]/length(drivers_xor[[x]]),2),
                           FUN.VALUE = c(1.0))

  if (cliquesbool){

    CliquesObject <- generatecliques(...)

    # Add cliques
    Driv <- names(drivers_xor)

    Cliques <- vector('list',length(Driv))
    names(Cliques) <- Driv

    for (x in CliquesObject$cliques){
      for (y in Driv[Driv%in%x]){
        Cliques[[y]]<-setdiff(x,y)
      }
    }

    Cliques<- vapply(Cliques,
                     function(x) paste(x,collapse=','),
                     FUN.VALUE = c('C'))

    # Reorder (in normalized_sum they are in alphabetical order,
    # in Cliques as they come out from the cliques object)
    Cliques <- Cliques[order(names(Cliques))]


    # Create Data frame with cliques
    sumXOR <- data.frame(normalized_sum, Cliques, check.names = TRUE,check.rows = TRUE)
    colnames(sumXOR) <- c("Normalized_XOR_Sum", "Cliques")

  } else{

    # Create Data frame without cliques
    sumXOR <- data.frame(normalized_sum, check.names = TRUE,check.rows = TRUE)
    colnames(sumXOR) <- c("Normalized_XOR_Sum")

  }


  # Export to xlsx ------------------------------------------------------------

  # Create the object
  excelwb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(excelwb, "Table")
  openxlsx::addWorksheet(excelwb, "Summary")

  # Write the data
  openxlsx::writeData(excelwb, "Table", excel, startRow = 1, startCol = 1)
  openxlsx::writeData(excelwb, "Summary", sumXOR, startRow = 1, startCol = 1,rowNames = TRUE)

  openxlsx::saveWorkbook(excelwb, file = paste0(wpath,'/grnsumm.xlsx'),
               overwrite = TRUE)
}
