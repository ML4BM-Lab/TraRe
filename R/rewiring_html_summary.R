#' Contain neccesary functions to generate rewiring html report
#'
#'
#' @noRd

maxrepeat <- function(invector) {
    sort(table(invector), decreasing = TRUE)[1]
}

module_node_summary <- function(norm_mat_keep, runmoddata, name2idx, nonrespond_idxs, responder_idxs) {
    # get graph method free edge data
    
    genenames = c(runmoddata$regulators, runmoddata$target_genes)
    regidxs = name2idx[runmoddata$regulators]
    targetidxs = name2idx[runmoddata$target_genes]
    geneidxs = c(regidxs, targetidxs)
    isreg = c(rep(1, length(regidxs)), rep(0, length(targetidxs)))
    
    allmat = norm_mat_keep[geneidxs, , drop = FALSE]
    allmat0 = norm_mat_keep[geneidxs, nonrespond_idxs, drop = FALSE]
    allmat1 = norm_mat_keep[geneidxs, responder_idxs, drop = FALSE]
    
    # count of most frequent value in gene
    repAll = apply(allmat, 1, maxrepeat)
    rep0 = apply(allmat0, 1, maxrepeat)
    rep1 = apply(allmat1, 1, maxrepeat)
    
    # find average value in gene
    avgAll = rowMeans(allmat)
    avg0 = rowMeans(allmat0)
    avg1 = rowMeans(allmat1)
    
    avgdiff = avg0 - avg1
    
    test_vals = NULL
    for (gene in genenames) {
        tres = stats::t.test(allmat0[gene, ], allmat1[gene, ])
        wres = stats::wilcox.test(allmat0[gene, ], allmat1[gene, ])
        myrow = c(signif(tres$p.value, 3), signif(tres$statistic, 3), signif(wres$p.value, 3), signif(wres$statistic, 
            3))
        test_vals = rbind(test_vals, myrow)
    }
    nodesumm = cbind(genenames, isreg, signif(avgAll, 3), signif(avg0, 3), signif(avg1, 3), signif(avgdiff, 
        3), test_vals, repAll, rep0, rep1)
    colnames(nodesumm) = c("gene-name", "is-regulator", "avg-all", "avg-nonresp", "avg-respond", "avg-diff", 
        "t-pval", "t-stat", "wilcox-pval", "wilcox-stat", "num-repeats-all", "num-repeats-nonresp", "num-repeats-respond")
    rownames(nodesumm) = genenames
    return(nodesumm)
}


module_edge_summary <- function(norm_mat_keep, runmoddata, name2idx, nonrespond_idxs, responder_idxs) {
    # get graph method free edge data
    
    regidxs = name2idx[runmoddata$regulators]
    targetidxs = name2idx[runmoddata$target_genes]
    
    regmat = norm_mat_keep[runmoddata$regulators, , drop = FALSE]
    regmat0 = norm_mat_keep[runmoddata$regulators, nonrespond_idxs, drop = FALSE]
    regmat1 = norm_mat_keep[runmoddata$regulators, responder_idxs, drop = FALSE]
    
    tarmat = norm_mat_keep[runmoddata$target_genes, ]
    tarmat0 = norm_mat_keep[runmoddata$target_genes, nonrespond_idxs]
    tarmat1 = norm_mat_keep[runmoddata$target_genes, responder_idxs]
    
    edgesumm = NULL
    for (myreg in runmoddata$regulators) {
        for (mytar in runmoddata$target_genes) {
            pcorall = stats::cor.test(regmat[myreg, ], tarmat[mytar, ])
            pcor0 = stats::cor.test(regmat0[myreg, ], tarmat0[mytar, ])
            pcor1 = stats::cor.test(regmat1[myreg, ], tarmat1[mytar, ])
            spcorall = stats::cor.test(regmat[myreg, ], tarmat[mytar, ], method = "spearman")
            spcor0 = stats::cor.test(regmat0[myreg, ], tarmat0[mytar, ], method = "spearman")
            spcor1 = stats::cor.test(regmat1[myreg, ], tarmat1[mytar, ], method = "spearman")
            myrow = c(paste0(mytar, "||", myreg), myreg, mytar, pcorall$estimate, pcorall$p.value, pcor0$estimate, 
                pcor0$p.value, pcor1$estimate, pcor1$p.value, spcorall$estimate, spcorall$p.value, spcor0$estimate, 
                spcor0$p.value, spcor1$estimate, spcor1$p.value)
            edgesumm = rbind(edgesumm, myrow)
        }
    }
    colnames(edgesumm) = c("key", "reg", "target", "all-pearson", "all-pearson-pvalue", "nonresp-pearson", "nonresp-pearson-pvalue", 
        "respond-pearson", "respond-pearson-pvalue", "all-spearman", "all-spearman-pvalue", "nonresp-spearman", 
        "nonresp-spearman-pvalue", "respond-spearman", "respond-spearman-pvalue")
    rownames(edgesumm) = edgesumm[, "key"]
    return(edgesumm)
}


# summarize module
summarize_module <- function(norm_expr_mat_keep, runmoddata, name2idx, nonrespond_idxs, responder_idxs) {
    summarized_list <- NULL
    
    # node summary
    nodesumm = module_node_summary(norm_expr_mat_keep, runmoddata, name2idx, nonrespond_idxs, responder_idxs)
    
    # edge summary
    edgesumm = module_edge_summary(norm_expr_mat_keep, runmoddata, name2idx)
    
    modregs = runmoddata$regulators
    modtargs = runmoddata$target_genes
    
    regidxs = name2idx[modregs]
    targetidxs = name2idx[modtargs]
    
    
    # get chip evidence and support data
    appendmat = matrix(0, dim(edgesumm)[1], 6)
    rownames(appendmat) = rownames(edgesumm)
    colnames(appendmat) = c("support", "chip-evidence", "num-chip-peaks", "all-weights", "nonresp-weights", 
        "respond-weights")
    
    # compute graphs and extract VBSR edge weights
    
    cut <- FALSE  #define variable to generate grah plots after trycatch graph generations.
    
    full_graph <- tryCatch({
        
        fg <- NET_compute_graph_all_VBSR(norm_expr_mat_keep, regidxs, targetidxs)
        
        # extract VBSR edge weights
        orderweights <- orderGraphWeights(fg, rownames(appendmat))
        appendmat[orderweights$commonedges, "all-weights"] <- orderweights$weights
        
        fg
        
    }, error = function(x) {
        methods::show(x)
        methods::show("cant generate full graph")
        cut <- TRUE
        return(NULL)
    })
    
    respond_graph <- tryCatch({
        rg <- NET_compute_graph_all_VBSR(norm_expr_mat_keep[, responder_idxs], regidxs, targetidxs)
        
        # extract VBSR edge weights
        orderweights <- orderGraphWeights(rg, rownames(appendmat))
        appendmat[orderweights$commonedges, "respond-weights"] <- orderweights$weights
        
        rg
        
    }, error = function(x) {
        methods::show(x)
        methods::show("cant generate responder graph")
        return(NULL)
    })
    
    nonresp_graph <- tryCatch({
        
        ng <- NET_compute_graph_all_VBSR(norm_expr_mat_keep[, nonrespond_idxs], regidxs, targetidxs)
        
        # extract VBSR edge weights
        orderweights <- orderGraphWeights(ng, rownames(appendmat))
        appendmat[orderweights$commonedges, "nonresp-weights"] <- orderweights$weights
        
        ng
        
    }, error = function(x) {
        
        methods::show(x)
        methods::show("cant generate nonresponder graph")
        return(NULL)
        
    })
    
    
    fulledgesumm = utils::type.convert(as.data.frame(cbind(edgesumm, appendmat), stringsAsFactors = FALSE))
    colnames(fulledgesumm) = make.names(colnames(fulledgesumm))
    
    return(list(nodesumm = nodesumm, fulledgesumm = fulledgesumm, full_graph = full_graph, respond_graph = respond_graph, 
        nonresp_graph = nonresp_graph, cut = cut))
}
