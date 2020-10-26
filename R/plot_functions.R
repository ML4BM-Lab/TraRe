#' Plotting functions for Gene Regulatory Network.
#'
#' Collection of functions for generating graphs layouts to plot GRN obtained from `NET_run()` method.
#' `return_layout()` generates a layout from the graph object returned by `NET_run()` and `return_layout_phenotype()`
#' plots targets according to the t-statistic from the differential expression analysis of the desired phenotype.
#' `plot_igraph()` takes in the igraph::igraph object and generated layout and generates plot.
#'
#' @param mygraph igraph::igraph object returned from `NET_run()`.
#' @param mytitle Desired tittle.
#' @param titlecol Color for the tittle.
#' @param mylayout desired layout.
#'
#' @return None.
#'
#' @examples
#'    \dontrun{
#'    testNet <- NET_run(lognorm_est_counts, target_filtered_idx,
#'                       regulator_filtered_idx, Gene_set_Collections,
#'    graph_mode=c("VBSR"),FDR=0.05,NrCores=3)
#'
#'    mygraph <- testNet$graphs$VBSR
#'    mylayout <- return_layout(mygraph,c("Reg1","Reg2","Reg3"),c("Target1","Target2","Target3"))
#'
#'    plot_igraph(mygraph,"MyTittle","black",mylayout)
#'    }
#' @export plot_igraph
plot_igraph <- function(mygraph, mytitle, titlecol, mylayout){

  nodecol <- c("darkblue", "darkorange")
  framecol <- c("black", "darkorange")
  shape <- c("circle", "square")
  edge_cscale <- grDevices::colorRamp(c("darkred", "lightgrey", "darkgreen"))

  igraph::E(g)$weight <- rep(1,length(igraph::E(g))) #assign weight 1.
  maxw=max(abs(igraph::E(mygraph)$weight))
  tweight = (igraph::E(mygraph)$weight+maxw)/(2*maxw)
  #tweight <- rep(1, length(E(mygraph)))
  igraph::E(mygraph)$color <- apply(edge_cscale(tweight), 1,
                            function(x) grDevices::rgb(x[1] / 255, x[2] / 255,
                                            x[3] / 255, 0.8))

  degrees = igraph::degree(mygraph, igraph::V(mygraph)$name)
  #show(degrees)
  nodenames = mylayout$genesnames[igraph::V(mygraph)$name]
  #show(nodenames)
  regdegrees = degrees[nodenames]
  regdegrees[which(is.na(regdegrees))]=""
  #show(regdegrees)
  finalnames = apply(cbind(nodenames,regdegrees),1,paste, collapse=" - ")

  plot(mygraph,
       vertex.color = nodecol[as.numeric(igraph::V(mygraph)$type) + 1],
       vertex.shape = shape[as.numeric(igraph::V(mygraph)$type) + 1],
       vertex.label = finalnames,
       vertex.label.cex = 1.5,
       #vertex.label.cex = 3.5,
       vertex.frame.color = framecol[as.numeric(igraph::V(mygraph)$type) + 1],
       #vertex.size = as.numeric(V(mygraph)$type)*10 + 10,
       vertex.size = as.numeric(igraph::V(mygraph)$type)*5 + 5,
       layout = cbind(mylayout$genesx[igraph::V(mygraph)$name],
                      mylayout$genesy[igraph::V(mygraph)$name]
       )
  )
  graphics::title(paste0(mytitle, " ", sum(igraph::V(mygraph)$type==1), "&", sum(igraph::V(mygraph)$type==0)), cex.main = 5, col.main = titlecol)
  graphics::abline(h=0, col=grDevices::rgb(0,0,0,alpha=0.3))
}
#' @export
#' @rdname plot_igraph
#' @param regs regulators name list
#' @param targets targets name list
#' @param namehash dictionary with genes as keys and transcripts as values.
#' If there is no transcripts, build the dictionary with genes as keys and values. ({"g1":"g1","g2":"g2"})
return_layout <- function(regs, targets, namehash){
  nregs <- length(regs)
  myratio <- length(targets) / nregs
  genesx <- c(1:nregs * myratio - myratio / 2, 1:length(targets))
  names(genesx) <- c(regs, targets)
  genesy <- c(rep(c(1, -1), nregs)[1:nregs] * (1 + stats::runif(nregs)),
              rep(0, length(targets)))
  names(genesy) <- c(regs, targets)
  genesnames <- c(namehash[regs], rep("", length(targets)))
  names(genesnames) <- c(regs, targets)
  return(list(genesx = genesx, genesy = genesy,
              genesnames = genesnames))
}
#' @export
#' @rdname plot_igraph
#' @param varfile two column file containing, gene names as rows,
#' t-statistic from the differential expression analysis of the desired phenotype column and
#' a boolean variable for regulator (1) - no regulator (0) column.
return_layout_phenotype <- function(regs, targets, varfile){

  vals = as.numeric(varfile[,"t-stat"])
  genesnames = rownames(varfile)[order(vals)]
  names(genesnames) = genesnames

  genesx = 1:length(vals)
  names(genesx) = genesnames

  orderedregs = names(genesx)[which(varfile[names(genesx),"is-regulator"]==1)]
  absval = max(abs(vals))
  genesy = signif(vals[order(vals)]/absval,3)
  names(genesy) = genesnames
  genesy[orderedregs] = genesy[orderedregs] + rep(c(2,-2),length(regs))[1:length(regs)]

  #my part
  nregs <- length(regs)
  myratio <- length(targets) / nregs
  genesx[orderedregs]= 1:nregs * myratio - myratio / 2

  genesnames[targets] = ""
  #genesnames[regs] = namehash[regs]
  genesnames[regs]= regs

  return(list(genesx=genesx, genesy=genesy, genesnames=genesnames))
}

