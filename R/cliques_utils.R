getMaxVariance <- function(clique,dat){return (max(sapply(unlist(clique),function(x) dat$hash[[x]])))}

getMaxVarName <- function(clique,dat){return (names(which.max(sapply(unlist(clique),function(x) dat$hash[[x]]))))}

getAvgVariance<-function(clique,dat){ return (mean(sapply(clique,function(x) dat$hash[[x]])))}

getMedian_AvgVar<-function(clique,dat){return (getMedianCorrelation(clique)+getAvgVariance(clique,dat))}

getMedianCorrelation<-function(clique,dat){

  if (!length(clique)) return(0)
  else if (length(clique)==1) return(1)

  #CorrMatrix will change if drivers or targets are run.
  return (stats::median(apply(utils::combn(which(rownames(dat$mat)%in%clique),2),2,function(x) dat$mat[x[1],x[2]])))

}

RemoveDuplicities<-function(patron,cliques){

  mcliques<-cliques[patron] #modified cliques

  unmcliques <- unlist(mcliques)
  mcliques<- Map('[', mcliques, utils::relist(!duplicated(unmcliques), skeleton = mcliques))

  mcliques<-Filter(Negate(function(X){length(X)==0}),mcliques)

  return(mcliques)

}


