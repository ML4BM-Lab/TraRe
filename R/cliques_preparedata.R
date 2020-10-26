#' Prepare data from input dataset.
#'
#' Prepare the correlation matrix and the hash table for future uses.
#'
#' @param dataset input expression file with genes as rows and samples as columns.
#' @param method method to use in the correlation matrix (see stats:cor).
#'
#' @return List containing hash table with genes variance and correlation matrix in the form of adjacency matrix (diag=0)
#'
#' @examples
#'   \dontrun{
#'   #this normally is loaded from outside.
#'   dataset <- readRDS(file="../dataset.rds")
#'
#'   foo <- preparedate (dataset)
#'   }
#' @export

preparedata <- function(dataset,method){

  #Generate variance gene dictionary.
  dataset_variance_dict<-hash::hash() #we use a hash table.
  dataset_variance_dict<-apply(dataset,1,stats::var)

  #Generate the correlation matrix to work with.
  CorrMatrix <- abs(stats::cor(t(dataset),method=method))
  diag(CorrMatrix) <- 0 #eliminate diagonal

  return (list(hash=dataset_variance_dict,mat=CorrMatrix))
}
