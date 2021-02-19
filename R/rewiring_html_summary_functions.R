#' Generate html report from the `runrewiring()` method.
#'
#' @description
#'
#' Contain functions to generate an html report from the rewired linker modules.
#' `create_index_page()` initializes folders and files for the html report. `write_tables_all()` generates
#' a table given stats from the rewiring method. `write_html_table_page()` creates the header for the
#' html file when the table will be placed. `table2html()` gets the `write_html_table_page()` info and
#' return the generated html file. This functions run inside `runrewiring()`
#'
#' @param outdir path for the resultant folder.
#' @param runtag name that is `paste()` to the name folder.
#' @param codedir folder from which `sorttable.js` and `glossary.txt`
#' files are `paste()` into the output folder.
#' @param indexpath name for the index file.
#' @param glossarypath name for the glossary file.
#' @param imgstr path for the image folder.
#' @param txtstr path for the txts folder.
#'
#' @param mytab table to transform into html file.
#' @param tabletype type of table.
#' @param html_idxs index for table rows.
#' @param html_cols index for table columns.
#' @param filestr name of the html file to be generated.
#' @param htmlinfo list containing `indexpath`,`txtstr`, and html dir.
#'
#'
#' @param resultstable modified table from `write_tables_all()`.
#' @param htmlpagefile path of the html file.
#' @param resultsrelpath path of the txt files.
#'
#'
#'
#' @noRd

table2html <- function(resultstable) {

    colheader <- paste(collapse = "</th><th>", c("idx", colnames(resultstable)))
    theader <- paste0("<thead>\n<tr bgcolor='#AAAAAA';><th>", colheader, "</th></tr>\n</thead>\n<tbody>\n")
    head_str = paste0("<script src='sorttable.js'></script>\n<TABLE class=", "'sortable' border =1 >\n", theader)

    numtable = cbind(seq_len(nrow(resultstable)), resultstable)
    rowinnerstrs = apply(numtable, 1, paste, collapse = "</td><td>")
    rowstrs = paste0("<tr><td>", rowinnerstrs, "</td></tr>")
    rows_str = paste(collapse = "\n", rowstrs)

    return(paste0(collapse = "\n", head_str, rows_str, "\n</tbody></table><br>"))

}


write_html_table_page <- function(resultstable, htmlpagefile, resultsrelpath, indexpath = "index.html", glossarypath = "glossary.html") {

    write(paste0("<table border = 1 width = '100%'><tr bgcolor = ", "'#AAAAAA'><th><a href = '", indexpath,
        "' target='_blank'>", "Index</a></th><th><a href = '", glossarypath, "' target=", "'_blank'>Glossary</a></th><th><a href = '",
        resultsrelpath, "' target='_blank'>Download</a></th></tr></table><br><br>"), file = htmlpagefile)
    htmlstr <- table2html(resultstable)
    write(htmlstr, file = htmlpagefile, append = TRUE)
}


write_tables_all <- function(mytab, tabletype = "table", html_idxs = seq_len(nrow(mytab)), html_cols = colnames(mytab),
    filestr = "data", htmlinfo = list(htmldir = "html/", indexpath = "index.html", txtstr = "txts/"),extradir='') {

    htmlpath <- paste0(extradir,filestr, "_", tabletype, ".html")
    resultspath <- paste0(extradir,htmlinfo$txtstr, filestr, "_", tabletype, ".txt")
    message("Writing table: ", resultspath)
    utils::write.table(mytab, paste0(htmlinfo$htmldir, resultspath), sep = "\t", row.names = FALSE, col.names = TRUE,
        quote = FALSE)
    write(paste0("<a href = \"", htmlpath, "\" target=\"_blank\">", tabletype, "</a><br>"), file = paste0(htmlinfo$htmldir,
        htmlinfo$indexpath), append = TRUE)
    write_html_table_page(resultstable = mytab[html_idxs, html_cols], htmlpagefile = paste0(htmlinfo$htmldir,
        htmlpath), resultsrelpath = resultspath, indexpath = htmlinfo$indexpath)
}

create_index_page <- function(outdir = "./", runtag = "run", codedir = "code/", indexpath = "index.html", glossarypath = "glossary.html",
    imgstr = "imgs/", txtstr = "txts/") {

    htmldir <- paste0(outdir, runtag, "/")
    dir.create(file.path(htmldir))

    file.copy(from = paste0(codedir, "sorttable.js"), to = htmldir)
    dir.create(file.path(paste0(htmldir, imgstr)))
    dir.create(file.path(paste0(htmldir, txtstr)))

    glossary <- as.matrix(utils::read.table(paste0(codedir, "glossary.txt"), header = TRUE, sep = "\t", quote = ""))
    file.copy(from = paste0(codedir, "glossary.txt"), to = htmldir)
    abspath <- paste0(htmldir, indexpath)

    write(paste0("<br>"), file = abspath)

    return(list(htmldir = htmldir, indexpath = indexpath, imgstr = imgstr, txtstr = txtstr, glossarypath = glossarypath,
        abspath = abspath))
}

violinPlots <- function(norm_expr_mat_keep, keepsamps, keeplabels, nodesumm, edgesumm, appendmat, htmlinfo) {
    modhtmlfile = paste0(htmlinfo$htmldir, htmlinfo$indexpath)
    imgdir = paste0(htmlinfo$htmldir, htmlinfo$imgstr)

    normexpmat = norm_expr_mat_keep[, keepsamps]
    nodesumm2 = type.convert(as.data.frame(nodesumm, stringsAsFactors=F))
    colnames(nodesumm2) = make.names(colnames(nodesumm2))
    # edgesumm2 = type.convert(as.data.frame(cbind(edgesumm, appendmat), stringsAsFactors=F))
    # colnames(edgesumm2) = make.names(colnames(edgesumm))

    reg_info = nodesumm2[nodesumm2$is.regulator==1,c("t.stat","t.pval")]
    reg_info = reg_info[order(reg_info$t.stat, decreasing=T),]
    reg_vec = rownames(reg_info)

    # --- regulator expression distributions ---
    write(paste0("<table style='width:100%' bgcolor='gray'><tr><td><h1>",
            "Regulator Expression Levels",
            "</h1></td></tr></table><br>\n"),
            modhtmlfile, append = T)

    respond_plotlist = list()
    nonresp_plotlist = list()
    for(reg in reg_vec){
      respond_plotlist[[reg]] = as.numeric(normexpmat[reg,keeplabels==1])
      nonresp_plotlist[[reg]] = as.numeric(normexpmat[reg,keeplabels==0])
      show(c(reg, signif(mean(respond_plotlist[[reg]])), signif(mean(nonresp_plotlist[[reg]]))))
    }

    plotwidth = 500
    plotheight = 800
    myplotname = paste0("violin.regexp")

    grDevices::png(paste0(imgdir, myplotname, ".png"),
        width = plotwidth, height = plotheight)
    vioplot::vioplot(x=respond_plotlist, col = "palevioletred", side = "right",
            names = paste(sep="||", rownames(reg_info), signif(reg_info$t.pval,2)),
            plotCentre = "line", horizontal = T, cex.axis=0.8)
    vioplot::vioplot(x=nonresp_plotlist, col = "lightblue", side = "left", add = T,
            plotCentre = "line", horizontal = T)
    abline(v=-5:5*2, col="lightgrey")
    title(xlab = "Expression Value (CQN)", ylab = "Regulator")
    legend("topleft", fill = c("palevioletred", "lightblue"), legend = c("R", "NR"), title = "Outcome", cex=.5)
    grDevices::dev.off()

    rankdf = data.frame(matrix(vector(), length(reg_vec), 0,
                      dimnames=list(sample(reg_vec), c())),
                      stringsAsFactors=F)
    rankdfcols = c()
    rankdfcols = c(rankdfcols, "regulator", "reg.expr.tpval", "reg.expr.rank")
    rankdf = cbind(rankdf[rownames(reg_info),], rownames(reg_info), signif(reg_info$t.pval,2), rank(reg_info$t.pval))
    colnames(rankdf) = rankdfcols

    write(paste0("<img src='", htmlinfo$imgstr, myplotname, ".png",
               "' alt='", myplotname,
               "' height='", plotheight, "' width='", plotwidth, "'><br>\n"),
          modhtmlfile, append = T)

    # --- target expression diff distributions ---
    write(paste0("<table style='width:100%' bgcolor='gray'><tr><td><h1>",
             "Target Expression Diff Distributions",
             "</h1></td></tr></table><br>\n"),
            modhtmlfile, append = T)
    write("<table style='width:100%'> ", modhtmlfile, append = T)
    write(paste0("<tr>"), modhtmlfile, append = T)

    for (cname in c("all", "nonresp", "respond")){
      nonzeroedges = edgesumm[edgesumm[,paste0(cname,".weights")]!=0,]
      plotlist = list()
      for(reg in reg_vec){
        targets = as.character(nonzeroedges[nonzeroedges$reg==reg,"target"])
        plotlist[[reg]] = as.numeric(nodesumm2[targets,"t.stat"])
        show(c(reg, length(targets), signif(mean(plotlist[[reg]]))))
      }
      plotlist = plotlist[names(sort(unlist(lapply(plotlist, mean)), decreasing=T))]
      plotlist[["supmod"]] = as.numeric(nodesumm2[,"t.stat"])

      # added png
      plotwidth = 500
      plotheight = 800
      myplotname0 = paste0("violin.tarexpdiff.", cname)

      grDevices::png(paste0(imgdir, myplotname0, ".png"),
          width = plotwidth, height = plotheight)

      vioplot::vioplot(x=plotlist, col = "goldenrod", side = "right", main = cname,
              names = paste(sep="||", names(plotlist),
                                      unlist(lapply(plotlist, length)),
                                      signif(unlist(lapply(plotlist, mean)),3)),
              plotCentre = "line", horizontal = T, cex.axis=0.6)
      abline(v=-5:5*2, col="lightgrey")
      title(xlab = "T Statistic: Higher in Responders <-> Not Signficant <-> Higher in Non-Responders", ylab = "Regulator Targets")
      grDevices::dev.off()

      if(cname=="all"){
        rankdfcols = c(rankdfcols, "all.targets", "tar.diffexpr.avg", "tar.diffexpr.rank")
        tmpvals = signif(unlist(lapply(plotlist, mean)),3)[reg_vec]
        rankdf = cbind(rankdf[reg_vec,], unlist(lapply(plotlist, length))[reg_vec], tmpvals, rank(tmpvals))
        colnames(rankdf) = rankdfcols
      }

      write(paste0(" <img src='", htmlinfo$imgstr, myplotname0, ".png",
                   "' alt='", myplotname0,
                   "' height='", plotheight, "' width='", plotwidth, "'> &emsp; \n"),
            modhtmlfile, append = T)
    }
    write(paste0("</tr>\n"), modhtmlfile, append = T)
    write(paste0("</table>"), modhtmlfile, append = T)

    # --- target regulator correlation distributions ---
    write(paste0("<table style='width:100%' bgcolor='gray'><tr><td><h1>",
                    "Target Regulator Correlation Distributions",
                    "</h1></td></tr></table><br>\n"),
              modhtmlfile, append = T)
    write("<table style='width:100%'> ", modhtmlfile, append = T)
    write(paste0("<tr>"), modhtmlfile, append = T)

    # calculate module correlation
    modmat = t(normexpmat[rownames(nodesumm2),])
    corall = cor(modmat)
    corR = cor(modmat[keeplabels==1, ])
    corNR = cor(modmat[keeplabels==0, ])
    corR[is.na(corR)]=0
    corNR[is.na(corNR)]=0
    cordiff = corNR - corR

    for (loopmode in c("all", "specific")){
      # plot regulator target set correlation
      respond_plotlist = list()
      nonresp_plotlist = list()
      mypvals = rep(2, length(reg_vec))
      names(mypvals) = reg_vec
      for(reg in reg_vec){
        respondtars = as.character(edgesumm$target[edgesumm$respond.weights!=0 & edgesumm$reg==reg])
        nonresptars = as.character(edgesumm$target[edgesumm$nonresp.weights!=0 & edgesumm$reg==reg])
        if(loopmode=="all"){
          respondtars = as.character(edgesumm$target[edgesumm$all.weights!=0 & edgesumm$reg==reg])
          nonresptars = as.character(edgesumm$target[edgesumm$all.weights!=0 & edgesumm$reg==reg])
        }
        respond_plotlist[[reg]] = as.numeric(corR[reg,respondtars])
        nonresp_plotlist[[reg]] = as.numeric(corNR[reg,nonresptars])
        if(length(respondtars) < 2 | length(nonresptars) < 2){
          mypvals[reg] = 1
          respond_plotlist[[reg]] = 0
          nonresp_plotlist[[reg]] = 0
        } else if (var(respond_plotlist[[reg]]) == 0 | var(nonresp_plotlist[[reg]]) == 0) {
          mypvals[reg] = 1
        } else {
          tres = t.test(respond_plotlist[[reg]], nonresp_plotlist[[reg]])
          mypvals[reg] = signif(tres$p.value, 1)
        }
        show(c(reg, length(respondtars), length(nonresptars), signif(mean(respond_plotlist[[reg]])), signif(mean(nonresp_plotlist[[reg]])), as.numeric(mypvals[reg])))
      }
      respond_plotlist= respond_plotlist[names(sort(mypvals,decreasing = T))]
      nonresp_plotlist= nonresp_plotlist[names(sort(mypvals,decreasing = T))]
    }

    for (loopmode in c("all", "specific")){
      # plot histogram of target/regulator correlation
      plotwidth = 500
      plotheight = 800
      myplotname = paste0("violin.tarregcorr.", loopmode)

      grDevices::png(paste0(imgdir, myplotname, ".png"),
          width = plotwidth, height = plotheight)
      vioplot::vioplot(x=respond_plotlist, col = "palevioletred", side = "right", main = loopmode,
              names = paste(sep="||", names(respond_plotlist),
                            unlist(lapply(respond_plotlist, length)),
                            unlist(lapply(nonresp_plotlist, length)),
                            mypvals[names(respond_plotlist)]),
              plotCentre = "line", horizontal = T, cex.axis=0.6)
      vioplot::vioplot(x=nonresp_plotlist, col = "lightblue", side = "left", add = T,
              plotCentre = "line", horizontal = T)
      abline(v=-5:5/5, col="lightgrey")
      title(xlab = "Pearson Correlation", ylab = "Regulator Targets")
      legend("topleft", fill = c("palevioletred", "lightblue"), legend = c("R", "NR"), title = "Outcome", cex=.5)
      grDevices::dev.off()

      if(loopmode=="all"){
        rankdfcols = c(rankdfcols, "regtar.corr.tpval", "regtar.corr.rank")
        rankdf = cbind(rankdf[names(mypvals),], mypvals, rank(mypvals))
        colnames(rankdf) = rankdfcols
      }

      write(paste0("<img src='", htmlinfo$imgstr, myplotname, ".png",
                   "' alt='", myplotname,
                   "' height='", plotheight, "' width='", plotwidth, "'> &emsp; \n"),
            modhtmlfile, append = T )
    }
    write(paste0("</tr>\n"), modhtmlfile, append = T)
    write(paste0("</table>"), modhtmlfile, append = T)

    return(rankdf)
}

bipartiteGraphsSumm <- function(cluster_num, modsumm, modmeth, htmlinfo, outdir){
    modhtmlfile = paste0(htmlinfo$htmldir, htmlinfo$indexpath)
    imgdir = paste0(htmlinfo$htmldir, htmlinfo$imgstr)

    write(
      paste0(
        "<table style='width:100%' bgcolor='gray'><tr><td><h1>",
        "Modules Summary",
        "</h1></td></tr></table><br>\n"
      ),
      modhtmlfile, append = T
    )

    # write plot to index page
    write(paste0("<img src='", outdir,'/supermodule_',cluster_num,'/imgs/',pname, ".png",
                   "' alt='", pname,
                   "' height='", 750, "' width='", 1500, "'> &emsp; <br>\n"),
            modhtmlfile,
            append = TRUE)

    sortidxs <- sort(as.numeric(modsumm$nodesumm[, "t-pval"]),
                     decreasing = F, index.return = T)$ix
    write_tables_all(modsumm$nodesumm[sortidxs, ],
                     tabletype = paste0(modmeth, "_nodesumm"),
                     filestr = "data", html_idxs = 1:dim(modsumm$nodesumm)[1],
                     htmlinfo = htmlinfo)

    sortidxs <- sort(as.numeric(modsumm$fulledgesumm[, "all.weights"]),
                     decreasing = F, index.return = T)$ix
    write_tables_all(modsumm$fulledgesumm[sortidxs, ],
                     tabletype = paste0(modmeth, "_edgesumm"),
                     filestr = "data", html_idxs = 1:dim(modsumm$fulledgesumm)[1],
                     htmlinfo = htmlinfo)

}

geneOrder <- function(modsumm, keepsamps, keeplabels, norm_mat_keep) {
    modregs = levels(unique(modsumm$fulledgesumm$reg))
    modtargs = levels(unique(modsumm$fulledgesumm$target))
    keepfeat = unique(c(modregs, modtargs))
    mat = t(norm_mat_keep[keepfeat, keepsamps])
    
    # calculate correlation matrices
    corall = cor(mat)
    cor1 = cor(mat[(keeplabels + 1) == 1, 1:ncol(mat)])
    cor2 = cor(mat[(keeplabels + 1) == 2, 1:ncol(mat)])
    cor1[is.na(cor1)]=0
    cor2[is.na(cor2)]=0
    cordiff = cor1 - cor2
    cormats = list(
      corall = corall,
      cordiff = cordiff,
      cor1 = cor1,
      cor2 = cor2
    )

    # order genes by clustering of corrdiff
    targetorder = labels(as.dendrogram(hclust(dist(cordiff[modtargs, modtargs]))))
    regorder = modregs
    if (length(modregs) > 1) {
      regorder = labels(as.dendrogram(hclust(dist(cordiff[modregs, modregs]))))
    }

    return(list("modregs"=modregs, "modtargs"=modtargs, "mat"=mat, "targetorder"=targetorder, "regorder"=regorder, "cormats"=cormats))
}


superModuleStatistics <- function(modregs, modtargs, mat, keeplabels, htmlinfo) {
    modhtmlfile = paste0(htmlinfo$htmldir, htmlinfo$indexpath)

    regnames = paste(collapse = ", ", modregs)
    split = as.numeric(table(keeplabels))
    stats = c(
      length(modtargs),
      length(modregs),
      regnames,
      dim(mat),
      split
    )

    statsnames = c(
      "num-targets",
      "num-regulators",
      "regulator-names",
      "num-samples",
      "num-genes",
      "num-nonresp",
      "num-respond"
    )

    stattab = cbind(statsnames, stats)

    write(
      paste0(
        "<table style='width:100%' bgcolor='gray'><tr><td><h1>",
        "SuperModule Statistics",
        "</h1></td></tr></table><br>\n"
      ),
      paste0(htmlinfo$htmldir, htmlinfo$indexpath),
      append = T
    )
    write(table2html(stattab), modhtmlfile, append = T)
}

createLegendPlot <- function(htmlinfo){
    imgdir = paste0(htmlinfo$htmldir, htmlinfo$imgstr)

    plotwidth=450
    plotheight=150
    myplotname = paste0("correlation_colorscale")
    grDevices::png(paste0(imgdir, myplotname,".png"), width=plotwidth, height=plotheight)
    ngroups = 9
    vals = seq(-1,1,length.out=ngroups)
    colramp = colorRampPalette(c("darkred","gray100","darkgreen"))(ngroups)
    image(vals, 1, as.matrix(vals,ngroups,1), col=colramp, axes=F, main="", xlab="", ylab="")
    axis(1, vals, cex.axis=2)
    grDevices::dev.off()

    myplotname = paste0("expression_colorscale")
    grDevices::png(paste0(imgdir, myplotname,".png"), width=plotwidth, height=plotheight)
    ngroups = 9
    vals = seq(-10,10,length.out=ngroups)
    colramp = colorRampPalette(c("darkorange","gray100","darkblue"))(ngroups)
    image(vals, 1, as.matrix(vals,ngroups,1), col=colramp, axes=F, main="", xlab="", ylab="")
    axis(1, vals, cex.axis=2)
    grDevices::dev.off()
}

# plot target correlation
correlationOfModuleGene <- function(mymod, regorder, targetorder, mat, cormats, keeplabels, htmlinfo, phenotype_class_vals){
    modhtmlfile = paste0(htmlinfo$htmldir, htmlinfo$indexpath)
    imgdir = paste0(htmlinfo$htmldir, htmlinfo$imgstr)

    write(
      paste0(
        "<table style='width:100%' bgcolor='gray'><tr><td><h1>",
        "Correlation Of Module Genes ",
        "</h1></td></tr></table><br>\n"
      ),
      modhtmlfile,
      append = T
    )

    write(
      paste0(
        "<img src='",
        htmlinfo$imgstr,
        "correlation_colorscale",
        ".png",
        "' alt='",
        "correlation_colorscale",
        "' height='",
        150,
        "' width='",
        450,
        "'><br>\n"
      ),
      modhtmlfile,
      append = T
    )

    write("<table style='width:100%'> ", modhtmlfile, append = T)
    plot_correlation_row(
      cormats = cormats,
      rowdesc = "targets",
      xnames = targetorder,
      ynames = targetorder,
      plotheight = 450,
      myshowrows = TRUE,
      htmlfile = modhtmlfile,
      imgdir = imgdir,
      modnum = mymod,
      plotwidth = 450,
      mycvec = c("darkred", "gray100", "darkgreen"),
      plotzlim = c(-1, 1),
      plottitle = phenotype_class_vals
    )


    # plot regulator correlation
    plot_correlation_row(
      cormats = cormats,
      rowdesc = "regulators",
      xnames = regorder,
      ynames = regorder,
      plotheight = 200,
      myshowrows = TRUE,
      htmlfile = modhtmlfile,
      imgdir = imgdir,
      modnum = mymod,
      plotwidth = 200,
      mycvec = c("darkred", "gray100", "darkgreen"),
      plotzlim = c(-1, 1),
      plottitle = phenotype_class_vals
    )

    # plot regulator by target correlation
    plot_correlation_row(
      cormats = cormats,
      rowdesc = "all",
      xnames = targetorder,
      ynames = regorder,
      plotheight = 450,
      myshowrows = TRUE,
      htmlfile = modhtmlfile,
      imgdir = imgdir,
      modnum = mymod,
      plotwidth = 200,
      mycvec = c("darkred", "gray100", "darkgreen"),
      plotzlim = c(-1, 1),
      plottitle = phenotype_class_vals
    )

    write(paste0("</table>"), modhtmlfile, append = T)

}

expressionTableOfModuleGenes <- function(mymod, nodesumm, htmlinfo) {
    modhtmlfile = paste0(htmlinfo$htmldir, htmlinfo$indexpath)

    write(
      paste0(
        "<table style='width:100%' bgcolor='gray'><tr><td><h1>",
        "Expression Table of Module Genes",
        "</h1></td></tr></table><br>\n"
      ),
      modhtmlfile,
      append = T
    )

    sortidxs = sort(as.numeric(nodesumm[, "t-pval"]),
                    decreasing = F,
                    index.return = T)$ix
    htmlidxs = union(sortidxs[1:25], which(nodesumm[, "is-regulator"] == "1"))
    write(table2html(nodesumm[htmlidxs, ]), modhtmlfile, append = T)

    # TODO : fix link
    sortidxs <- sort(as.numeric(nodesumm[, "t-stat"]),
                     decreasing = F, index.return = T)$ix
    write_tables_all(nodesumm[sortidxs, ],
                     tabletype = paste0(mymod, "_mod_node_summ"),
                     filestr = "data", html_idxs = 1:dim(nodesumm)[1],
                     htmlinfo = htmlinfo)

    # resultspath = paste(sep = ".", "mod_node_summ", mymod, "txt")
    # show(paste0("Writing table: ", resultspath))
    # write(
    #   paste0(
    #     "<a href = '",
    #     htmlinfo$txtstr,
    #     resultspath,
    #     "' target='_blank'>Download Full Module Node Table</a> <br>\n"
    #   ),
    #   modhtmlfile,
    #   append = T
    # )
}

expressionPlotsOfModuleGenes <- function(mymod, regorder, targetorder, mat, samps2pheno, phenostrs, htmlinfo) {
    modhtmlfile = paste0(htmlinfo$htmldir, htmlinfo$indexpath)
    imgdir = paste0(htmlinfo$htmldir, htmlinfo$imgstr)

    write(
        paste0(
          "<table style='width:100%' bgcolor='gray'><tr><td><h1>",
          "Expression Plots of Module Genes ",
          "</h1></td></tr></table><br>\n"
        ),
        modhtmlfile,
        append = T
      )

    write(
      paste0(
        "<img src='",
        htmlinfo$imgstr,
        "expression_colorscale",
        ".png",
        "' alt='",
        "correlation_colorscale",
        "' height='",
        150,
        "' width='",
        450,
        "'><br>\n"
      ),
      modhtmlfile,
      append = T
    )

    write("<table style='width:100%'> ", modhtmlfile, append = T)

    # plot regulator expression
    plot_expression_row(
      mymat = t(mat)[regorder, , drop = F],
      rowdesc = "regulators",
      plotheight = 200,
      myshowrows = TRUE,
      samps2pheno = samps2pheno,
      phenostrs = phenostrs,
      htmlfile = modhtmlfile,
      imgdir = imgdir,
      modnum = mymod,
      plotwidth = 800,
      mycvec = c("darkorange", "gray100", "darkblue"),
      plotzlim = c(-10, 10)
    )

    # plot target expression
    plot_expression_row(
      mymat = t(mat)[targetorder, ],
      rowdesc = "targets",
      plotheight = 600,
      myshowrows = TRUE,
      samps2pheno = samps2pheno,
      phenostrs = phenostrs,
      htmlfile = modhtmlfile,
      imgdir = imgdir,
      modnum = mymod,
      plotwidth = 800,
      mycvec = c("darkorange", "gray100", "darkblue"),
      plotzlim = c(-10, 10)
    )
    write(paste0("</table>"), modhtmlfile, append = T)
}

nullDistributionOfRewiringStatistic <- function(mat, keeplabels, modmeth, mymod, htmlinfo){
    modhtmlfile = paste0(htmlinfo$htmldir, htmlinfo$indexpath)
    imgdir = paste0(htmlinfo$htmldir, htmlinfo$imgstr)

    result = NULL
    result = rewiring_test_pair_detail(mat, keeplabels + 1, perm = 1000)

    write(
      paste0(
        "<table style='width:100%' bgcolor='gray'><tr><td><h1>",
        "Null Distribution of Rewiring Statistic",
        "</h1></td></tr></table><br>\n"
      ),
      modhtmlfile,
      append = T
    )

    # plot histogram of dave stats
    plotwidth = 400
    plotheight = 400
    myplotname = paste0("hist.", modmeth, ".mod", mymod)
    grDevices::png(paste0(imgdir, myplotname, ".png"),
        width = plotwidth,
        height = plotheight)
    rewiring_score = result$T_star
    hist(main = paste0("True Val = ", signif(result$T, 3)), rewiring_score)
    abline(v = result$T, col = "red")
    grDevices::dev.off()

    write(
      paste0(
        "<img src='", htmlinfo$imgstr, myplotname, ".png",
        "' alt='", myplotname,
        "' height='", plotheight, "' width='", plotwidth, "'><br>\n"
      ),
      modhtmlfile, append = T)
}

regulatorSummaryAndRank <- function(rankdf, htmlinfo){
    modhtmlfile = paste0(htmlinfo$htmldir, htmlinfo$indexpath)
  
    write(paste0("<table style='width:100%' bgcolor='gray'><tr><td><h1>",
                       "Regulator Summary And Rank",
                       "</h1></td></tr></table><br>\n"),
                modhtmlfile, append = T )
    write(table2html(rankdf), modhtmlfile, append = T)
}

