#' Plotting functions for Gene Regulatory Network.
#'
#' Collection of functions for generating graphs layouts to plot GRN obtained from `NET_run()` method.
#' `return_layout()` generates a layout from the graph object returned by `NET_run()` and `return_layout_phenotype()`
#' plots targets according to the t-statistic from the differential expression analysis of the desired phenotype.
#' `plot_igraph()` takes in the igraph object and generated layout and generates plot.
#'
#' @param mygraph igraph object returned from `NET_run()`.
#' @param mytitle Desired tittle.
#' @param titlecol Color for the tittle.
#' @param mylayout desired layout.
#'
#' @return plot of the desired single GRN using a specific layout.
#'
#' @examples
#'
#'    ## Assume we have run the rewiring method and the `NET_run()` method to generate the
#'    ## igraph object. We are going to generate and plot both layouts for the example.
#'    ## We are going to generate all the files we need except for the igraph object, which
#'    ## is included as an example file in this package.
#'
#'
#'    ## We load the igraph object that we generated from the `NET_run()` example.
#'    ## Note: the igraph object is inside the list `NET_run()` generates.
#'
#'    graph <- readRDS(paste0(system.file("extdata",package="TraRe"),'/graph_netrun_example.rds'))
#'
#'
#'    ## We first generate the normal layout for the plot.
#'    ## We need the drivers and target names.
#'
#'    drivers <- readRDS(paste0(system.file("extdata",package="TraRe"),'/tfs_cliques_example.rds'))
#'    drivers_n <- rownames(drivers)[1:5]
#'
#'    targets <- readRDS(paste0(system.file("extdata",package="TraRe"),'/targets_linker_example.rds'))
#'    targets_n <- rownames(targets)[1:30]
#'
#'    ## As for this example we are working at gene level (we dont have transcripts inside genes),
#'    ## we will generate a dictionary with genes as keys and values (see param `namehash`)
#'
#'    namehash <- drivers_n
#'    normal_layout <- return_layout(drivers_n,targets_n,namehash)
#'
#'    ## We now generate the phenotype layout and the `varfile` we ned for this layout.
#'    ## (I leave here a way to generate) We need to separate our expression matrix by
#'    ## a binary phenotype, for this case, i will consider the first 40 samples are
#'    ## responding to a treatment (R) and the rest not (NR).
#'
#'    gnames <- c(drivers_n,targets_n)
#'    expmat <-rbind(drivers,targets)
#'    expmat_R <- expmat[,1:40]
#'    expmat_NR <- expmat[,40+(1:28)]
#'
#'
#'    varfile <- t(as.matrix(sapply(gnames,
#'               function(x) c(stats::t.test(expmat_R[x,],expmat_NR[x,])$statistic,
#'               if(x%in%drivers_n) 1 else 0))))
#'
#'    colnames(varfile)<-c("t-stat","is-regulator")
#'
#'    phenotype_layout <- return_layout_phenotype(drivers_n,targets_n,varfile)
#'
#'    plot_igraph(graph,"Normal Layout","black",normal_layout)
#'    plot_igraph(graph,"Phenotype Layout","black",phenotype_layout)
#'
#' @export plot_igraph
plot_igraph <- function(mygraph, mytitle, titlecol, mylayout){

  nodecol <- c("darkblue", "darkorange")
  framecol <- c("black", "darkorange")
  shape <- c("circle", "square")
  edge_cscale <- grDevices::colorRamp(c("darkred", "lightgrey", "darkgreen"))

  igraph::E(mygraph)$weight <- rep(1,length(igraph::E(mygraph))) #assign weight 1.
  maxw <- max(abs(igraph::E(mygraph)$weight))
  tweight = (igraph::E(mygraph)$weight+maxw)/(2*maxw)
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
  genesx <- c(seq_len(nregs) * myratio - myratio / 2, seq_along(targets))
  names(genesx) <- c(regs, targets)
  genesy <- c(rep(c(1, -1), nregs)[seq_len(nregs)] * (1 + stats::runif(nregs)),
              rep(0, length(targets)))
  names(genesy) <- c(regs, targets)
  #genesnames <- c(namehash[regs], rep("", length(targets)))
  genesnames <- c(regs, rep("", length(targets)))
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

  genesx = seq_along(vals)
  names(genesx) = genesnames

  orderedregs = names(genesx)[which(varfile[names(genesx),"is-regulator"]==1)]
  absval = max(abs(vals))
  genesy = signif(vals[order(vals)]/absval,3)
  names(genesy) = genesnames
  genesy[orderedregs] = genesy[orderedregs] + rep(c(2,-2),length(regs))[seq_along(regs)]

  #my part
  nregs <- length(regs)
  myratio <- length(targets) / nregs
  genesx[orderedregs]= seq_len(nregs) * myratio - myratio / 2

  genesnames[targets] = ""
  #genesnames[regs] = namehash[regs]
  genesnames[regs]= regs

  return(list(genesx=genesx, genesy=genesy, genesnames=genesnames))
}

