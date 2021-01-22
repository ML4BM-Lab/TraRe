#' Create HTML report
#'
#' Contains all the necessary functions to run `create_html_summary()`.
#'
#' @noRd

linker_summarize_rungraphs <- function(rungraphs = NULL, iso_table = NULL, weighted_chip_evidence = NULL, graphstr = "data", 
    htmlinfo = list(htmldir = "html/", indexpath = "index.html", txtstr = "txts/"), runinfo = list(nregs = "1", 
        ntargets = "5", nboots = "10")) {
    
    message("Processing ", graphstr, " graph method ...")
    write(paste0("<br><br><b>", graphstr, "</b><br><br>"), file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), 
        append = TRUE)
    
    # summarize edge results
    edgesinfo <- edgeinfo_from_graphs(rungraphs = rungraphs, iso_table = iso_table, weighted_chip_evidence = weighted_chip_evidence)
    # edgesinfo <<- edgesinfo
    
    # write edge results summary
    write_tables_all(edgesinfo, tabletype = "edges", html_cols = c("edgekey", "weight", "reg-origid", "target-origid", 
        "reg-geneid", "target-geneid", "chip-evidence", "num-chip-peaks"), html_idxs = seq_len(min(1000, dim(edgesinfo)[1])), 
        filestr = graphstr, htmlinfo = htmlinfo)
    
    # regulators summary
    regsinfo <- genetable_summary(filter = "reg-origid", filterNeigh = "target-origid", tabletype = "regulators", 
        nboots = as.integer(runinfo$nboots), minsupport = 4, edgesinfo = edgesinfo, graphstr = graphstr, htmlinfo = htmlinfo)
    
    # targets summary
    targetsinfo <- genetable_summary(filter = "target-origid", filterNeigh = "reg-origid", tabletype = "targets", 
        nboots = as.integer(runinfo$nboots), minsupport = 4, edgesinfo = edgesinfo, graphstr = graphstr, htmlinfo = htmlinfo)
    
    # summarize edge results
    message("Summarizing ", graphstr, " graph method ...")
    ngraphmods <- length(rungraphs)
    myweights <- as.numeric(edgesinfo[, "weight"])
    myevidence <- as.numeric(edgesinfo[, "chip-evidence"])
    
    write(paste0("<br>", ngraphmods, " graphs with an average of ", signif(sum(myweights)/ngraphmods, 3), " edges each.<br>"), 
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
    write(paste0(sum(myweights), " (", signif(sum(myweights)/(runinfo$nregs * runinfo$ntargets * runinfo$nboots) * 
        100, 3), "%) of ", runinfo$nregs, "*", runinfo$ntargets, "*", runinfo$nboots, " possible edges across bootstraps.<br>"), 
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
    write(paste0(length(myweights), " (", signif(length(myweights)/(runinfo$nregs * runinfo$ntargets) * 100, 
        3), "%) of ", runinfo$nregs, "*", runinfo$ntargets, " unique edges found.<br><br>"), file = paste0(htmlinfo$htmldir, 
        htmlinfo$indexpath), append = TRUE)
    
    evid_weight_raw <- cbind(c(myweights), c(myevidence))
    if (length(unique(myweights)) == 1 || length(unique(myweights)) == 1) {
        evid_weight_raw <- rbind(evid_weight_raw, c(0, -1))
    }
    evidence_by_weight_final <- cumSumTable(evid_weight_raw)
    message(evidence_by_weight_final)
    
    summary_final <- cbind(rownames(evidence_by_weight_final), evidence_by_weight_final)
    if (dim(summary_final)[2] == 6) {
        colnames(summary_final) <- c("Support", "nEdges", "cumEdges", "%NA", "%NoPeak", "%Peaks")
    } else {
        colnames(summary_final) <- c("Support", "nEdges", "cumEdges")
    }
    
    htmlstr <- table2html(summary_final)
    write(htmlstr, file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
    return(summary_final)
}

linker_create_index_page <- function(outdir = "./", runtag = "run", codedir = "./", indexpath = "index.html", 
    glossarypath = "glossary.html", imgstr = "imgs/", txtstr = "txts/") {
    
    dir.create(file.path(outdir), showWarnings = FALSE)
    htmldir <- paste0(outdir, runtag, "/")
    dir.create(file.path(htmldir), showWarnings = FALSE)
    file.copy(from = paste0(codedir, "sorttable.js"), to = htmldir)
    dir.create(file.path(paste0(htmldir, imgstr)), showWarnings = FALSE)
    dir.create(file.path(paste0(htmldir, txtstr)), showWarnings = FALSE)
    
    glossary <- as.matrix(utils::read.table(paste0(codedir, "glossary.txt"), header = TRUE, sep = "\t", quote = ""))
    file.copy(from = paste0(codedir, "glossary.txt"), to = htmldir)
    abspath <- paste0(htmldir, indexpath)
    
    write(paste0("<br>"), file = abspath)
    write_html_table_page(glossary, paste0(htmldir, glossarypath), glossarypath)
    
    return(list(htmldir = htmldir, indexpath = indexpath, imgstr = imgstr, txtstr = txtstr, glossarypath = glossarypath, 
        abspath = abspath))
}

genetable_summary <- function(filter, filterNeigh, tabletype, nboots = 10, minsupport = 4, edgesinfo, graphstr, 
    htmlinfo) {
    message(tabletype, " processing...")
    mygids <- unique(edgesinfo[, filter])
    # only extract top 1000 targets
    if (tabletype == "targets") {
        tmptable <- table(as.data.frame(edgesinfo[, c("target-origid", "weight")]))
        topgenes <- rownames(tmptable)
        names(topgenes) <- rownames(tmptable)
        if (nboots > 1) {
            scaledtable <- t(t(as.matrix(tmptable)) * as.numeric(colnames(tmptable)))
            keepsupport <- nboots:minsupport
            colstrs <- intersect(colnames(scaledtable), as.character(keepsupport))
            topgenes <- rowSums(scaledtable[, as.character(colstrs)])
        }
        sortidxs <- sort(topgenes, decreasing = TRUE, index.return = TRUE)$ix
        mygids <- rownames(tmptable)[sortidxs[seq_len(min(1000, length(topgenes)))]]
    }
    ### 
    genetable <- t(sapply(mygids, extract_gene_row, filter, filterNeigh, edgesinfo, nboots))
    colnames(genetable) <- c("gid", "nEdges", "nNeigh", "nMax-Conf-Neigh", "topNeigh", "suppTable")
    sortval <- (as.numeric(genetable[, "nMax-Conf-Neigh"]) * 1e+09 + as.numeric(genetable[, "nEdges"]))
    sortidxs <- sort(sortval, decreasing = TRUE, index.return = TRUE)$ix[seq_len(min(1000, length(sortval)))]
    write_tables_all(genetable, tabletype = tabletype, html_idxs = sortidxs, filestr = graphstr, htmlinfo = htmlinfo)
    return(genetable)
}

edgeinfo_from_graphs <- function(rungraphs, iso_table, weighted_chip_evidence) {
    
    message("Extracting edges from all graphs...")
    edgelist_list <- lapply(rungraphs, extract_edge_strings_from_graph)
    edgestable <- table(unlist(edgelist_list))
    
    reg_origids <- unlist(lapply(strsplit(names(edgestable), "\\|\\|"), "[[", 2))
    tar_origids <- unlist(lapply(strsplit(names(edgestable), "\\|\\|"), "[[", 1))
    message("reg_origids: ")
    showfirstlast(reg_origids)
    message("tar_origids: ")
    showfirstlast(tar_origids)
    
    # get gene ids
    reg_gids <- as.character(iso_table[reg_origids, "iso_ensgs"])
    tar_gids <- as.character(iso_table[tar_origids, "iso_ensgs"])
    message("reg_gids: ")
    showfirstlast(reg_gids)
    message("tar_gids: ")
    showfirstlast(tar_gids)
    
    # get chip evidence
    raw_evidence <- rep(-1, length(edgestable))
    if (!is.null(weighted_chip_evidence)) {
        message("...Extracting chip evidence")
        raw_evidence <- weighted_chip_evidence[cbind(gsub("-", ".", reg_gids), gsub("-", ".", tar_gids))]
    }
    binary_evidence <- raw_evidence
    binary_evidence[binary_evidence > 1] <- 1
    methods::show(table(binary_evidence))
    sortval <- edgestable * 1e+06 + raw_evidence
    
    # combine edge summary table
    edgesinfo <- cbind(names(edgestable), edgestable, reg_origids, tar_origids, reg_gids, tar_gids, binary_evidence, 
        raw_evidence, sortval)
    colnames(edgesinfo) <- c("edgekey", "weight", "reg-origid", "target-origid", "reg-geneid", "target-geneid", 
        "chip-evidence", "num-chip-peaks", "sort-val")
    
    return(edgesinfo[sort(as.numeric(edgesinfo[, "sort-val"]), decreasing = TRUE, index.return = TRUE)$ix, ])
}

# Helpers -----------------------------------------------------------------

showfirstlast <- function(mynames) {
    methods::show(c(length(mynames), length(unique(mynames)), sort(mynames)[c(1, length(mynames))]))
}


cumSumTable <- function(mydataframe) {
    dftable <- table(as.data.frame(as.matrix(mydataframe)))
    idxnums <- sort(as.numeric(rownames(dftable)), decreasing = TRUE)
    dftable <- dftable[as.character(idxnums), ]
    if (is.null(dim(dftable))) {
        # if only one column or row
        nCount <- dftable
        cumCount <- cumsum(dftable)
        # show(tmptab)
        return(cbind(nCount, cumCount))
    } else {
        cumul_table <- apply(dftable, 2, cumsum)
        nCount <- rowSums(dftable)
        cumCount <- rowSums(cumul_table)
        dftable_final <- cbind(nCount, cumCount, signif(cumul_table/rowSums(cumul_table) * 100, 3))
        # show(dftable_final)
        return(dftable_final)
    }
}


extract_edge_strings_from_graph <- function(mygraph) {
    elist <- igraph::get.edgelist(mygraph)
    # small graphs lose node names when extracted
    if (is.numeric(elist)) {
        return()
    } else {
        return(apply(elist, 1, paste, collapse = "||"))
    }
}

write_tables_all <- function(mytab, tabletype = "table", html_idxs = seq_len(dim(mytab)[1]), html_cols = colnames(mytab), 
    filestr = "data", htmlinfo = list(htmldir = "html/", indexpath = "index.html", txtstr = "txts/")) {
    htmlpath <- paste0(filestr, "_", tabletype, ".html")
    resultspath <- paste0(htmlinfo$txtstr, filestr, "_", tabletype, ".txt")
    message("Writing table: ", resultspath)
    utils::write.table(mytab, paste0(htmlinfo$htmldir, resultspath), sep = "\t", row.names = F, col.names = T, 
        quote = F)
    write(paste0("<a href = \"", htmlpath, "\" target=\"_blank\">", tabletype, "</a><br>"), file = paste0(htmlinfo$htmldir, 
        htmlinfo$indexpath), append = T)
    write_html_table_page(resultstable = mytab[html_idxs, html_cols], htmlpagefile = paste0(htmlinfo$htmldir, 
        htmlpath), resultsrelpath = resultspath, indexpath = htmlinfo$indexpath)
}

extract_gene_row <- function(mygid, myfield, myfield2, edgesinfo, nboots = 10, maxneigh = 3) {
    # mygid = 'ENSG00000196092' myfield = 'reg-geneid' myfield2 = 'target-origid' mygid = 'ENSG00000181804'
    # myfield2 = 'reg-origid' myfield = 'target-geneid'
    
    gididxs <- which(edgesinfo[, myfield] == mygid)
    mysubtab <- edgesinfo[gididxs, c(myfield2, "weight", "num-chip-peaks", "chip-evidence")]
    if (length(gididxs) == 1) {
        mysubtab <- t(as.matrix(mysubtab))
    }
    
    evid_weight_tab <- ""
    if (length(gididxs) > 1) {
        evid_weight_tab <- cumSumTable(mysubtab[, c("weight", "chip-evidence")])
        evid_weight_tab <- evid_weight_tab[c(1, dim(evid_weight_tab)[1]), ]
    }
    
    return(c(mygid, sum(as.numeric(mysubtab[, "weight"])), length(gididxs), sum(as.numeric(mysubtab[, "weight"]) == 
        nboots), tab2cell(mysubtab[seq_len(min(maxneigh, length(gididxs))), ], keeprnames = FALSE), tab2cell(evid_weight_tab)))
}


tab2cell <- function(mytab, keepheader = TRUE, keeprnames = TRUE) {
    if (keepheader) {
        mytab <- rbind(colnames(mytab), mytab)
    }
    if (keeprnames) {
        mytab <- cbind(rownames(mytab), mytab)
    }
    return(paste(collapse = "<br>", apply(mytab, 1, paste, collapse = " | ")))
}

