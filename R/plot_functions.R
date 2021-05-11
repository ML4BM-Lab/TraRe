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
#'    ## We load the igraph object that we generated from the `NET_run()` example.
#'    ## Note: the igraph object is inside the list `NET_run()` generates.
#'
#'    graph <- readRDS(paste0(system.file('extdata',package='TraRe'),
#'                     '/graph_netrun_example.rds'))$graphs$VBSR
#'
#'    ## We first generate the normal layout for the plot.
#'    ## We need the drivers and target names.
#'
#'    drivers <- readRDS(paste0(system.file('extdata',package='TraRe'),'/tfs_linker_example.rds'))
#'    drivers_n <- rownames(drivers)
#'
#'    targets <- readRDS(paste0(system.file('extdata',package='TraRe'),'/targets_linker_example.rds'))
#'    targets_n <- rownames(targets)
#'
#'    ## As for this example we are working at gene level (we dont have transcripts inside genes),
#'    ## we will generate a dictionary with genes as keys and values (see param `namehash`)
#'
#'    normal_layout <- return_layout(drivers_n,targets_n)
#'
#'    ## We now generate the phenotype layout and the `varfile` we ned for this layout.
#'    ## (I leave here a way to generate) We need to separate our expression matrix by
#'    ## a binary phenotype, for this case, i will consider the first 40 samples are
#'    ## responding to a treatment (R) and the rest not (NR).
#'
#'    gnames <- c(drivers_n,targets_n)
#'    expmat <-rbind(drivers,targets)
#'
#'    phenotype <- utils::read.delim(paste0(system.file('extdata',package='TraRe'),
#'                                   '/phenotype_rewiring_example.txt'))
#'
#'    expmat_R <- expmat[,phenotype$Class=='R']
#'    expmat_NR <- expmat[,phenotype$Class=='NR']
#'
#'
#'    varfile <- t(as.matrix(sapply(gnames,
#'               function(x) c(stats::t.test(expmat_R[x,],expmat_NR[x,])$statistic,
#'               if(x%in%drivers_n) 1 else 0))))
#'
#'    colnames(varfile)<-c('t-stat','is-regulator')
#'
#'    phenotype_layout <- return_layout_phenotype(drivers_n,targets_n,varfile)
#'
#'    plot_igraph(graph,mytitle='Normal Layout',titlecol='black',mylayout=normal_layout)
#'    plot_igraph(graph,mytitle='Phenotype Layout',titlecol='black',mylayout=phenotype_layout)
#'
#' @export plot_igraph

plot_igraph <- function(mygraph = NULL, mytitle = "", titlecol = "black", mylayout = NULL) {

    if (is.null(mygraph)) {
        stop("graph object field empty")
    }
    if (is.null(mylayout)) {
        stop("layout field empty")
    }

    if (is.null(igraph::E(mygraph)$weight)) {
        igraph::E(mygraph)$weight <- rep(1, length(igraph::E(mygraph)))
    }

    nodecol <- c("darkblue", "darkorange")
    framecol <- c("black", "darkorange")
    shape <- c("circle", "square")
    edge_cscale <- grDevices::colorRamp(c("darkred", "lightgrey", "darkgreen"))

    maxw <- max(abs(igraph::E(mygraph)$weight))
    tweight <- (igraph::E(mygraph)$weight + maxw)/(2 * maxw)
    igraph::E(mygraph)$color <- apply(edge_cscale(tweight), 1, function(x) grDevices::rgb(x[1]/255, x[2]/255, x[3]/255, 0.8))

    degrees <- igraph::degree(mygraph, igraph::V(mygraph)$name)
    nodenames <- mylayout$genesnames[igraph::V(mygraph)$name]
    regdegrees <- degrees[nodenames]
    regdegrees[which(is.na(regdegrees))] <- ""
    finalnames <- apply(cbind(nodenames, regdegrees), 1, paste, collapse = " - ")

    plot(mygraph, vertex.color = nodecol[as.numeric(igraph::V(mygraph)$type) + 1], vertex.shape = shape[as.numeric(igraph::V(mygraph)$type) +
        1], vertex.label = finalnames, vertex.label.cex = 1.5, vertex.frame.color = framecol[as.numeric(igraph::V(mygraph)$type) +
        1], vertex.size = as.numeric(igraph::V(mygraph)$type) * 5 + 5, layout = cbind(mylayout$genesx[igraph::V(mygraph)$name], mylayout$genesy[igraph::V(mygraph)$name]))
    graphics::title(paste0(mytitle, " ", sum(igraph::V(mygraph)$type == 1), "&", sum(igraph::V(mygraph)$type == 0)), cex.main = 2,
        col.main = titlecol)
    graphics::abline(h = 0, col = grDevices::rgb(0, 0, 0, alpha = 0.3))
}
#' @export
#' @rdname plot_igraph
#' @param regs regulators name list
#' @param targets targets name list
#' @param namehash list containing the drivers genes as names and transcripts as values.
#' If only genes are required, leave it empty.
return_layout <- function(regs = NULL, targets = NULL, namehash = NULL) {

    if (is.null(regs)) {
        stop("regulators field empty")
    }
    if (is.null(targets)) {
        stop("targets field empty")
    }
    if (is.null(namehash)) {
        namehash <- regs
    }


    nregs <- length(regs)
    myratio <- length(targets)/nregs
    genesx <- c(seq_len(nregs) * myratio - myratio/2, seq_along(targets))
    names(genesx) <- c(regs, targets)
    genesy <- c(rep(c(1, -1), nregs)[seq_len(nregs)] * (1 + stats::runif(nregs)), rep(0, length(targets)))
    names(genesy) <- c(regs, targets)
    if (length(names(namehash)) == 0) {
        names(namehash) <- namehash
    }
    genesnames <- c(namehash[regs], rep("", length(targets)))
    names(genesnames) <- c(regs, targets)
    return(list(genesx = genesx, genesy = genesy, genesnames = genesnames))
}
#' @export
#' @rdname plot_igraph
#' @param varfile two column file containing, gene names as rows,
#' t-statistic from the differential expression analysis of the desired phenotype column and
#' a boolean variable for regulator (1) - no regulator (0) column.
return_layout_phenotype <- function(regs = NULL, targets = NULL, varfile = NULL, namehash = NULL) {

    if (is.null(regs)) {
        stop("regulators field empty")
    }
    if (is.null(targets)) {
        stop("targets field empty")
    }
    if (is.null(varfile)) {
        stop("varfile field empty")
    }
    if (is.null(namehash)) {
        namehash <- regs
    }

    # check varfile structure
    if (is.null(rownames(varfile))) {
        stop("genes names must be specified at varfile as rownames")
    }
    if (is.null(colnames(varfile))) {
        stop("colnames must be specified, in particular 'is-regulator' and 't-stat'")
    }
    if (!("is-regulator" %in% colnames(varfile))) {
        stop("varfile must contain the column is-regulator")
    }
    if (!("t-stat" %in% colnames(varfile))) {
        stop("varfile must contain the column t-stat")
    }


    vals <- as.numeric(varfile[, "t-stat"])
    genesnames <- rownames(varfile)[order(vals)]
    names(genesnames) <- genesnames

    genesx <- seq_along(vals)
    names(genesx) <- genesnames

    orderedregs <- names(genesx)[which(varfile[names(genesx), "is-regulator"] == 1)]
    absval <- max(abs(vals))
    genesy <- signif(vals[order(vals)]/absval, 3)
    names(genesy) <- genesnames
    genesy[orderedregs] <- genesy[orderedregs] + rep(c(2, -2), length(regs))[seq_along(regs)]

    # my part
    nregs <- length(regs)
    myratio <- length(targets)/nregs
    genesx[orderedregs] <- seq_len(nregs) * myratio - myratio/2

    genesnames[targets] <- ""
    if (length(names(namehash)) == 0) {
        names(namehash) <- namehash
    }
    genesnames[regs] <- namehash[regs]

    return(list(genesx = genesx, genesy = genesy, genesnames = genesnames))

}
#' @export
#' @rdname plot_igraph
#' @param graph igraph object
#' @param edgelist list containing the edges of the igraph object.
orderGraphWeights <- function(graph, edgelist) {

    weights <- igraph::get.data.frame(igraph::graph.adjacency(as.matrix(igraph::get.adjacency(graph, attr = "weight", type = "upper")),
        weighted = TRUE))
    rownames(weights) <- apply(weights[, seq_len(2)], 1, paste, collapse = "||")
    commonedges <- intersect(edgelist, rownames(weights))
    return(list(commonedges = commonedges, weights = weights[commonedges, "weight"]))

}
#' @export
#' @rdname plot_igraph
#' @param heatm input matrix for plot.
#' @param plotname name of the plot.
#' @param myzlim the range of z values for which colors should be plotted.
#' @param cvec vector of colors for the palette of the plot.
#' @param showRows boolean specifying the option of showing row names.
heatmapplot <- function(heatm, plotname = "", myzlim = c(min(heatm), max(heatm)), cvec = c("red", "white", "blue"), showRows = TRUE) {
    colramp <- (grDevices::colorRampPalette(cvec))(21)
    heatm[heatm < myzlim[1]] <- myzlim[1]
    heatm[heatm > myzlim[2]] <- myzlim[2]
    graphics::image(seq(ncol(heatm)), seq(nrow(heatm)), t(heatm), col = colramp, zlim = myzlim, axes = FALSE, xlab = "", ylab = "")

    graphics::title(main = plotname)
    if (showRows) {
        if (ncol(heatm) > 10) {
            #idxs <- which(1:dim(heatm)[1]%%round(dim(heatm)[1]/10, 0) == 0)
            idxs <- which(seq(nrow(heatm))%%round(seq(nrow(heatm))/10, 0) == 0)
            graphics::axis(2, at = idxs, labels = rownames(heatm)[idxs], las = 2, cex.axis = 0.8, tick = FALSE, col.axis = "black")
        } else {
            #graphics::axis(2, at = 1:dim(heatm)[1], labels = rownames(heatm), las = 2, cex.axis = 0.8, tick = FALSE, col.axis = "black")
            graphics::axis(2, at = seq(nrow(heatm)), labels = rownames(heatm), las = 2, cex.axis = 0.8, tick = FALSE, col.axis = "black")
        }
    }
}
#' @export
#' @rdname plot_igraph
#' @param mymat input matrix for plot.
#' @param rowdesc name uses to specify regulator.
#' @param plotheight height of the plot.
#' @param myshowrows boolean specifying the option of showing row names.
#' @param samps2pheno matrix of sample to phenotype.
#' @param phenostrs strings uses distinguish for phenotypes.
#' @param htmlfile directory to html files.
#' @param imgdir directory for image.
#' @param modnum the number of supermodule.
#' @param plotwidth width of the plot.
#' @param mycvec vector of colors for the palette of the plot.
#' @param plotzlim the range of z values for which colors should be plotted.
plot_expression_row <- function(mymat = NULL, rowdesc = "Regulators", plotheight = 200, myshowrows = TRUE, samps2pheno = NULL, phenostrs = c("nonrespond",
    "responder"), htmlfile = "./", imgdir = "imgs/", modnum = 1, plotwidth = 800, mycvec = c("darkorange", "gray100", "darkblue"),
    plotzlim = c(-10, 10)) {
    write(paste0("<tr>"), htmlfile, append = TRUE)
    for (respstr in phenostrs) {
        myplotname <- paste0("expr", ".mod", modnum, ".", respstr, ".", rowdesc)
        titlename <- paste0(respstr, ".", rowdesc)
        grDevices::png(paste0(imgdir, myplotname, ".png"), width = plotwidth, height = plotheight)
        heatmapplot(mymat[, names(which(samps2pheno == respstr)), drop = FALSE], plotname = titlename, myzlim = plotzlim, cvec = mycvec,
            showRows = myshowrows)
        grDevices::dev.off()
        write(paste0("<td> <img src='", "imgs/", myplotname, ".png", "' alt='", myplotname, "' height='", plotheight, "' width='",
            plotwidth, "'> </td>\n"), htmlfile, append = TRUE)
    }
    write(paste0("</tr>\n"), htmlfile, append = TRUE)
}
#' @export
#' @rdname plot_igraph
#' @param cormats input correlation matrix for plot.
#' @param rowdesc name uses to specify regulator.
#' @param xnames names for the x axis.
#' @param ynames names for the x axis.
#' @param plotheight height of the plot.
#' @param myshowrows boolean specifying the option of showing row names.
#' @param htmlfile directory to html files.
#' @param imgdir directory for image.
#' @param modnum the number of supermodule.
#' @param plotwidth width of the plot.
#' @param mycvec vector of colors for the palette of the plot.
#' @param plotzlim the range of z values for which colors should be plotted.
#' @param plottitle the title of the plot.
plot_correlation_row <- function(cormats = NULL, rowdesc = "regulators", xnames = NULL, ynames = NULL, plotheight = 200, myshowrows = TRUE,
    htmlfile = "./", imgdir = "imgs/", modnum = 1, plotwidth = 200, mycvec = c("darkred", "gray100", "darkgreen"), plotzlim = c(-1,
        1), plottitle = NULL) {

    write(paste0("<tr>"), htmlfile, append = TRUE)
    for (mattype in ls(cormats)) {
        myplotname <- paste0("corr", ".mod.", modnum, ".", mattype, ".", rowdesc)

        titlename <- mattype
        # if(mattype=='cor1'){titlename='nonrespond_only'} if(mattype=='cor2'){titlename='responder_only'}
        if (mattype == "cor1") {
            titlename <- plottitle[1]
        }
        if (mattype == "cor2") {
            titlename <- plottitle[2]
        }
        if (mattype == "corall") {
            titlename <- "all"
        }
        if (mattype == "cordiff") {
            titlename <- "corr_difference"
        }
        titlename <- paste0(titlename, ".", rowdesc)

        # show(myplotname)
        grDevices::png(paste0(imgdir, myplotname, ".png"), width = plotwidth, height = plotheight)
        heatmapplot(cormats[[mattype]][xnames, ynames, drop = FALSE], plotname = titlename, myzlim = plotzlim, cvec = mycvec, showRows = myshowrows)
        grDevices::dev.off()
        write(paste0("<td> <img src='", "imgs/", myplotname, ".png", "' alt='", myplotname, "' height='", plotheight, "' width='",
            plotwidth, "'> </td>\n"), htmlfile, append = TRUE)
    }
}
#' @export
#' @rdname plot_igraph
#' @param pname name for the plot.
#' @param myx the coordinates of points of the first gene.
#' @param myy the coordinates of points of the second gene.
#' @param xgenename the names of the first gene.
#' @param ygenename the names of the second gene.
#' @param mylabels integer class labels.
#' @param alltext text label for the plot.
#' @param plotdir directory for the plot.
plot_gene_pair_scatter <- function(pname, myx, myy, xgenename, ygenename, mylabels, alltext = NULL, plotdir = "") {
    mymax <- max(abs(c(myx, myy)))
    if (is.null(mylabels)) {
        mylabels <- rep(2, length(myx))
    }
    corall <- stats::cor.test(myx, myy)

    grDevices::png(paste0(plotdir, pname, ".png"), 500, 500)

    plot(x = myx, y = myy, main = paste(sep = " ", xgenename, ygenename), xlab = xgenename, ylab = ygenename, xlim = c(mymax * -1,
        mymax), ylim = c(mymax * -1, mymax), type = "p", pch = 16, cex = 1.5, col = scales::alpha(mylabels, 0.5))
    graphics::abline(0, 1, col = "black")
    graphics::abline(0.25, 1, col = "gray80")
    graphics::abline(0.5, 1, col = "gray60")
    graphics::abline(0.75, 1, col = "gray40")
    graphics::abline(-0.25, 1, col = "gray80")
    graphics::abline(-0.5, 1, col = "gray60")
    graphics::abline(-0.75, 1, col = "gray40")

    graphics::abline(h = 0)
    graphics::abline(v = 0)

    if (!is.null(alltext)) {
        laball <- paste0(alltext, "\n", round(corall$estimate, 3), " (", signif(corall$p.value, 3), ")")
        graphics::text(-1 * mymax, 0, labels = laball, col = 1, adj = c(0, 0))
    }

    grDevices::dev.off()
}
#' @export
#' @rdname plot_igraph
#' @param pname name for the plot.
#' @param myx the coordinates of points of the first gene.
#' @param myy the coordinates of points of the second gene.
#' @param xgenename the names of the first gene.
#' @param ygenename the names of the second gene.
#' @param mylabels integer class labels.
#' @param lab1text text label of a class.
#' @param lab2text text label of the other class.
#' @param plotdir directory for the plot.
plot_gene_pair_scatter_by_class <- function(plotdir, pname, myx, myy, xgenename, ygenename, mylabels, lab1text, lab2text, alltext) {

    mymax <- max(abs(c(myx, myy)))
    mycols <- c("darkviolet", "darkgoldenrod")

    corall <- stats::cor.test(myx, myy)
    cor1 <- stats::cor.test(myx[which(mylabels == 1)], myy[which(mylabels == 1)])
    cor2 <- stats::cor.test(myx[which(mylabels == 2)], myy[which(mylabels == 2)])

    treg <- stats::t.test(myx[which(mylabels == 1)], myx[which(mylabels == 2)])
    ttar <- stats::t.test(myy[which(mylabels == 1)], myy[which(mylabels == 2)])

    reglineall <- stats::lm(myy ~ myx)
    regline1 <- stats::lm(myy[which(mylabels == 1)] ~ myx[which(mylabels == 1)])
    regline2 <- stats::lm(myy[which(mylabels == 2)] ~ myx[which(mylabels == 2)])

    grDevices::png(paste0(plotdir, pname, ".png"), 400, 400)
    plot(x = myx, y = myy, main = paste(sep = " ", xgenename, "and", ygenename), xlab = paste(sep = " ", "Regulator:", xgenename),
        ylab = paste(sep = " ", "Target:", ygenename), xlim = c(mymax * -1, mymax), ylim = c(mymax * -1, mymax), type = "p", pch = 16,
        cex = 1.5, col = scales::alpha(mycols[mylabels], 0.5))
    graphics::abline(h = 0, col = grDevices::rgb(0, 0, 0, 0.9))
    graphics::abline(v = 0, col = grDevices::rgb(0, 0, 0, 0.9))

    graphics::abline(h = ttar$estimate[1], lty = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[1])[1]/255, grDevices::col2rgb(mycols[1])[2]/255,
        grDevices::col2rgb(mycols[1])[3]/255, 0.6))
    graphics::abline(h = ttar$estimate[2], lty = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[2])[1]/255, grDevices::col2rgb(mycols[2])[2]/255,
        grDevices::col2rgb(mycols[2])[3]/255, 0.6))
    graphics::points(x = mymax, y = ttar$estimate[1], pch = "+", lty = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[1])[1]/255,
        grDevices::col2rgb(mycols[1])[2]/255, grDevices::col2rgb(mycols[1])[3]/255, 0.6))
    graphics::points(x = mymax, y = ttar$estimate[2], pch = "+", lty = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[2])[1]/255,
        grDevices::col2rgb(mycols[2])[2]/255, grDevices::col2rgb(mycols[2])[3]/255, 0.6))

    graphics::abline(v = treg$estimate[1], lty = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[1])[1]/255, grDevices::col2rgb(mycols[1])[2]/255,
        grDevices::col2rgb(mycols[1])[3]/255, 0.6))
    graphics::abline(v = treg$estimate[2], lty = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[2])[1]/255, grDevices::col2rgb(mycols[2])[2]/255,
        grDevices::col2rgb(mycols[2])[3]/255, 0.6))
    graphics::points(y = -1 * mymax, x = treg$estimate[1], pch = "+", lty = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[1])[1]/255,
        grDevices::col2rgb(mycols[1])[2]/255, grDevices::col2rgb(mycols[1])[3]/255, 0.6))
    graphics::points(y = -1 * mymax, x = treg$estimate[2], pch = "+", lty = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[2])[1]/255,
        grDevices::col2rgb(mycols[2])[2]/255, grDevices::col2rgb(mycols[2])[3]/255, 0.6))

    graphics::abline(reglineall, lty = 5, col = grDevices::rgb(0, 0, 0, 0.6))
    graphics::abline(regline1, lty = 5, lwd = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[1])[1]/255, grDevices::col2rgb(mycols[1])[2]/255,
        grDevices::col2rgb(mycols[1])[3]/255, 0.6))
    graphics::abline(regline2, lty = 5, lwd = 3, col = grDevices::rgb(grDevices::col2rgb(mycols[2])[1]/255, grDevices::col2rgb(mycols[2])[2]/255,
        grDevices::col2rgb(mycols[2])[3]/255, 0.6))

    lab1 <- paste(sep = "", lab1text, "\n", round(cor1$estimate, 3), " (", signif(cor1$p.value, 2), ")")
    lab2 <- paste(sep = "", lab2text, "\n", round(cor2$estimate, 3), " (", signif(cor2$p.value, 2), ")")
    laball <- paste(sep = "", alltext, "\n", round(corall$estimate, 3), " (", signif(corall$p.value, 2), ")")
    labtde <- paste(sep = "", ygenename, "\n", "DE (", signif(ttar$p.value, 2), ")")
    labrde <- paste(sep = "", xgenename, "\n", "DE (", signif(treg$p.value, 2), ")")
    graphics::text(-1 * mymax, mymax, cex = 0.9, labels = lab1, col = mycols[1], adj = c(0, 1))
    graphics::text(-1 * mymax, -1 * mymax, cex = 0.9, labels = lab2, col = mycols[2], adj = c(0, 0))
    graphics::text(-1 * mymax, 0, cex = 0.9, labels = laball, col = 1, adj = c(0, 0))
    graphics::text(mymax, mymax, cex = 0.9, labels = labtde, col = 1, adj = c(1, 1))
    graphics::text(mymax, -1 * mymax, cex = 0.9, labels = labrde, col = 1, adj = c(1, 0))

    grDevices::dev.off()

    return(list(ttar0 = signif(as.numeric(ttar$estimate[1]), 3), ttar1 = signif(as.numeric(ttar$estimate[2]), 3), ttarp = signif(as.numeric(ttar$p.value),
        2), treg0 = signif(as.numeric(treg$estimate[1]), 3), treg1 = signif(as.numeric(treg$estimate[2]), 3), tregp = signif(as.numeric(treg$p.value),
        2)))
}


