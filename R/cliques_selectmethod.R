#' Method selection for removing duplicities from maximal_cliques method.
#'
#' Method selection in order to remove duplicities generated from the igraph::max_cliques method.
#' Four possibilities are available.
#'
#' @param selection integer selecting method.
#' @param ggoutput output from generategraph()
#' @param pdoutput output from preparedata().
#'
#' @return List containing the cliques generated, without duplicities, filtered by the chosen method, the selected method tittle and
#' the index from which the non-correlated genes start (singleton genes).
#'
#' @examples
#'    \dontrun{
#'    chosenmethod <- 1
#'    g <- igraph::graph object
#'    singleton <- c("Gene1","Gene2","Gene3")
#'
#'    foo <- selectmethod(chosenmethod,g,singleton)
#'    }
#' @export

selectmethod <- function(selection,ggoutput,pdoutput){

  #Run igraph::max_cliques() method
  Cliques<-sapply(igraph::max_cliques(ggoutput$graph)[],names)

  #1 - Maximizing Genes/Clique #sapply(Cliques,length)
  #2 - Maximize Median Correlation Value/Clique #sapply(Cliques,getAverageVarianceVector,VarVector=VarVector)
  #3 - Maximize Avg Variance Correlation Value/Clique #sapply(Cliques,getMedianCorrelation,method=method,DataSet=CorrVector)
  #4 - Maximize Sum(2,3)/Clique #2+3

  SortingMethod<-switch(selection,length,
                        function(x) getMedianCorrelation(x,dat=pdoutput),
                        function(x) getAvgVariance(x,dat=pdoutput),
                        function(x) getMedian_AvgVar(x,dat=pdoutput))

  SortingMethod_tittle<-switch(selection,'Maximizing Genes/Clique','Maximize Median Correlation Value/Clique','Maximize Avg Variance Correlation Value/Clique','Maximize Sum(2,3)/Clique')

  SortingMethod_chosen<-sapply(Cliques,SortingMethod)

  #We sort the cliques
  SortedCliques<-sort(SortingMethod_chosen,decreasing = TRUE,index.return=TRUE)$ix

  CliquesWithoutDuplicities<-RemoveDuplicities(SortedCliques,Cliques)
  L<-length(CliquesWithoutDuplicities) #We only count the ones with high correlated genes.

  CliquesWithoutDuplicities_plusRest<-append(CliquesWithoutDuplicities,ggoutput$singleton) #+SingletonCommunities

  return(list(cliques=CliquesWithoutDuplicities_plusRest,tittle=SortingMethod_tittle,L=L))

}
