#' Graph generation from genes correlation adjacency matrix.
#'
#' Graph generation from genes correlation, threshold dependant, adjacency matrix.
#'
#' @param correlationth threshold to consider edge exists.
#' @param sparsecorrmatrix boolean variable specifying whether to fix values below threshold to 0 or not.
#' @param pdoutput output from preparedata().
#'
#' @return List containing an igraph object, containing all the genes which have at least one edge as vertex, and
#' an array of all the genes which do not have correlation values above the threshold.
#'
#' @examples
#'    \dontrun{
#'
#'    #we load the desired dataset.
#'    dataset <- readRDS(".../dataset.rds")
#'
#'    #we select pearson correlation.
#'    method <- "pearson"
#'
#'    pdoutput <- preparedata(dataset,method)
#'    foo <- generategraph(correlationth = 0.6, sparsecorrmatrix = T,pdoutput)
#'    }
#' @export

generategraph<-function(correlationth,sparsecorrmatrix,pdoutput){

  methods::show("Generating groups of highly correlated genes and singleton communities")

  #Get the names of high correlated and non- high correlated

  GenesWhichHaveHighCorrelatedGenes_names<-rownames(pdoutput$mat)[apply(pdoutput$mat,1,function(x) Reduce('|',x>correlationth))]

  GenesWhichDoNotHaveHighCorrelatedGenes_names<-setdiff(rownames(pdoutput$mat),GenesWhichHaveHighCorrelatedGenes_names)

  #Generate the correlation matrix

  methods::show('Creating the matrix')
  GenesWhichHaveHighCorrelatedGenes_graph<-pdoutput$mat[GenesWhichHaveHighCorrelatedGenes_names,GenesWhichHaveHighCorrelatedGenes_names]

  if (sparsecorrmatrix){GenesWhichHaveHighCorrelatedGenes_graph[GenesWhichHaveHighCorrelatedGenes_graph<correlationth]<-0} #generate sparse matrix

  diag(GenesWhichHaveHighCorrelatedGenes_graph)<-0 #eliminate diagonal
  GenesWhichHaveHighCorrelatedGenes_graph<-igraph::graph_from_adjacency_matrix(GenesWhichHaveHighCorrelatedGenes_graph,mode="undirected",weighted=TRUE) #generate graph


  return(list(graph=GenesWhichHaveHighCorrelatedGenes_graph,singleton=GenesWhichDoNotHaveHighCorrelatedGenes_names))

}
