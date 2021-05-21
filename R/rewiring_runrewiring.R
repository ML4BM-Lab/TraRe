#' GRN Modules Rewiring Method.
#'
#' Gene Regulatory Network modules Rewiring method. It performs a permutation test,
#' (what we call rewiring test) and generates an html report containing a correlation matrix
#' with the higher scores obtained from the rewiring test. This matrix is shown in the way of a
#' heatmap, and its sorted by a hierarchical clustering for better interpretation.
#'
#' @param ObjectList Output from `preparerewiring()`containing some required parameters.
#'
#' @return It creates a folder (in tempdir() by default) containing the files explained above.
#'
#' @examples
#'
#' ## Lets assume that we have already generated the ObjectList, we will load it from
#' ## the folder containing the examples files. After this, it is all straight forward.
#'
#' objectlist <- readRDS(file=paste0(system.file('extdata',package='TraRe'),
#'                       '/prepared_rewiring_example.rds'))
#'
#'
#' ## We are going to create the folder containing
#' ## the graphs, reports, etc, and then we are deleting it.
#' ## If you want to keep it, do not run the last line.
#'
#' ## We are modifying output directory for this example.
#' objectlist$outdir <- paste(getwd(),'examplefolder',sep='/')
#'
#' runrewiring(ObjectList = objectlist)
#' unlink(objectlist$outdir,recursive = TRUE)
#'
#'
#' @export
runrewiring <- function(ObjectList) {

    # helper directory

    codedir <- paste0(system.file("extdata", package = "TraRe"), "/RewiringReport/")

    # initialize common parameters
    regulator_info_col_name <- ObjectList$regulator_info_col_name
    phenotype_class_vals <- ObjectList$phenotype_class_vals
    phenotype_class_vals_label <- ObjectList$phenotype_class_vals_label
    outdir <- paste0(ObjectList$outdir, "/supermod_rewiring")
    orig_test_perms <- ObjectList$orig_test_perms
    retest_thresh <- ObjectList$retest_thresh
    retest_perms <- ObjectList$retest_perms
    logfile <- ObjectList$logfile

    # set up output html page, we use the first argv.
    dir.create(ObjectList$outdir)
    indexpageinfo <- create_index_page(outdir = outdir, runtag = "", codedir = codedir)
    imgdir <- paste0(indexpageinfo$htmldir, indexpageinfo$imgstr)


    # Initialize log file
    logfile_p <- paste0(indexpageinfo$htmldir, "logfile.txt")

    write("----Prepare Rewiring summary----", logfile_p)
    write(logfile, logfile_p, append = TRUE)
    write("----End of Prepare Rewiring summary----", logfile_p, append = TRUE)

    # set.seed(1)
    dqrng::dqset.seed(1)

    # we create rundata and combine the modules of both parsers. Lets generate the dupla (dataset's method - dataset's number)

    duplas <- unlist(lapply(seq_along(ObjectList$datasets), function(i) paste(names(ObjectList$datasets[[i]]$rundata$modules), i)))

    # Initialize c_allstats and statsnames variables

    c_allstats <- c()  #This is for the combined heatmap
    c_module_membership_list <- hash::hash()  #This is for the combined heatmap

    statsnames <- c("module-method", "module-index", "orig-pval", "revised-pvalue", "num-targets", "num-regulators", "regulator-names",
        "target-names", "num-samples", "num-genes", "num-class1", "num-class2")

    for (dupla in duplas) {

        # Initialize rewired module's hash table
        module_membership_list <- hash::hash()

        # Initialize allstats array
        allstats <- NULL

        # For instance: 'VBSR X' to 'VBSR' and 'X'
        modmeth_i <- unlist(strsplit(dupla, " "))

        modmeth <- modmeth_i[1]
        i <- as.numeric(modmeth_i[2])

        # Output to the user which dupla we are working with
        message(modmeth, " ", i)

        rundata <- ObjectList$datasets[[i]]$rundata
        norm_expr_mat_keep <- ObjectList$datasets[[i]]$norm_expr_mat_keep
        keepsamps <- ObjectList$datasets[[i]]$keepsamps
        keeplabels <- ObjectList$datasets[[i]]$keeplabels
        class_counts <- ObjectList$datasets[[i]]$class_counts
        final_signif_thresh <- ObjectList$datasets[[i]]$final_signif_thresh
        responder <- ObjectList$datasets[[i]]$responder
        gene_info_df_keep <- ObjectList$datasets[[i]]$gene_info_df_keep
        name2idx <- ObjectList$datasets[[i]]$name2idx
        allregs <- ObjectList$datasets[[i]]$allregs
        alltargs <- ObjectList$datasets[[i]]$alltargs


        # This will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.

        parallClass <- BiocParallel::bpparam()
        parallClass$workers <- ObjectList$NrCores

        GenerateStats <- function(mymod) {

            signify <- NULL
            modregs <- unique(rundata$modules[[modmeth]][[mymod]]$regulators)
            modtargs <- unique(rundata$modules[[modmeth]][[mymod]]$target_genes)
            regnames <- paste(collapse = ", ", modregs)
            targnames <- paste(collapse = ", ", modtargs)
            keepfeats <- unique(c(modregs, modtargs))
            modmat <- t(norm_expr_mat_keep[keepfeats, keepsamps])
            modmeth_i_c <- paste(modmeth_i, collapse = " ")

            orig_pval <- rewiring_test(modmat, keeplabels + 1, perm = orig_test_perms)
            new_pval <- orig_pval
            stats <- c(modmeth_i_c, mymod, signif(orig_pval, 3), signif(new_pval, 3), length(modtargs), length(modregs), regnames,
                targnames, dim(modmat), class_counts)
            if (orig_pval < retest_thresh | orig_pval == 1 | mymod%%300 == 0) {
                # methods::show(paste(c('ModNum and NumGenes', mymod, length(keepfeats))))
                result <- rewiring_test_pair_detail(modmat, keeplabels + 1, perm = retest_perms)
                new_pval <- result$pval
                stats <- c(modmeth_i_c, mymod, signif(orig_pval, 3), signif(new_pval, 3), length(modtargs), length(modregs), regnames,
                  targnames, dim(modmat), class_counts)

                if (new_pval <= final_signif_thresh | new_pval == 1) {
                  # save as list
                  modname <- paste0(modmeth, ".", i, ".mod.", mymod)
                  # module_membership_list[[modname]] <- keepfeats
                  signify <- list(modname, keepfeats)
                }

            }
            return(list(stats, signify))

        }
        foreach_allstats <- BiocParallel::bplapply(seq_along(rundata$modules[[modmeth]]), GenerateStats, BPPARAM = parallClass)

        for (elements in foreach_allstats) {

            # now we recover first allstats matrix
            foreach_stats <- elements[[1]]
            allstats <- rbind(allstats, foreach_stats)

            # and then update the module_membership dictionary
            hashtable <- elements[[2]]

            if (!is.null(hashtable)) {

                module_membership_list[[hashtable[[1]]]] <- hashtable[[2]]

            }
        }

        # generate txt file with pvals
        colnames(allstats) <- statsnames
        rownames(allstats) <- allstats[, "module-index"]
        allstats[, "orig-pval"] <- signif(as.numeric(allstats[, "orig-pval"]), 3)
        allstats[, "revised-pvalue"] <- signif(as.numeric(allstats[, "revised-pvalue"]), 3)

        message("Generating Heatmap")

        clusters <- gen_heatmap(ObjectList, module_membership_list, allstats,
                                imgdir=imgdir,outdir=outdir,
                                indexpageinfo=indexpageinfo, modmeth = modmeth, i = i)


        # concatenate allstats and module_membership_list to generate the combinated heatmap
        # check if module_membership_list is empty
        if (length(module_membership_list)!=0){

        c_allstats <- rbind(c_allstats, allstats)
        c_module_membership_list <- hash::hash(c(hash::as.list.hash(c_module_membership_list),
                                               hash::as.list.hash(module_membership_list)))


        # Create multiplicity table
        supermod_regs_list <- NULL
        supermod_targs_list <- NULL


        # Detection of clusters (all but the last one)
        message(length(clusters$clusters) - 1, " clusters detected!")

        # select every cluster we have found except the last one.

        for (numclus in seq_along(clusters$clusters)) {

            message("Cluster number: ", numclus)

            # Name the folder of supermodule numclus within method modmeth
            foldername_p <- paste0('supermodule',i,'.',modmeth,'.', numclus)

            # Create dir for cluster numclus
            dir.create(paste0(outdir, "/", foldername_p))

            # Create txts folder in numclus's folder
            dir.create(paste0(outdir, "/", foldername_p, "/txts"))
            # Create imgs folder in numclus's folder
            dir.create(paste0(outdir, "/", foldername_p, "/imgs"))
            # Create html folder in numclu's folder
            dir.create(paste0(outdir, "/", foldername_p, "/htmls"))

            # Write to the logfile
            mods <- clusters$clusters[[numclus]]
            mods <- split(mods, cut(seq_along(mods), max(length(mods)/5, 2), labels = FALSE))
            mods <- paste0(vapply(mods, FUN = paste0, collapse = "|", FUN.VALUE = ""), collapse = "\n")

            # Add to logfile
            write(paste0("\nSupermodule ", numclus, ", mode ", modmeth, " ", i, " : ", mods, "\n"), logfile_p, append = TRUE)

            # add reference index to main index page
            write(paste0("<a href = '",foldername_p, "/index.html'>Cluster",numclus,"_Summary</a><br>"),
                paste0(indexpageinfo$htmldir, indexpageinfo$indexpath), append = TRUE)

            # add reference index to main index page
            write(paste0("<a href = '../index.html'>Return to Main Rewiring Summary</a><br>"),
                paste0(indexpageinfo$htmldir, foldername_p, "/", indexpageinfo$indexpath), append = TRUE)

            for (clusmod in clusters$clusters[[numclus]]) {

                clusmod_vec <- unlist(strsplit(clusmod, "\\."))
                # methods::show(clusmod_vec)
                mregs <- unique(rundata$modules[[clusmod_vec[1]]][[as.numeric(clusmod_vec[3])]]$regulators)
                mtargs <- unique(rundata$modules[[clusmod_vec[1]]][[as.numeric(clusmod_vec[3])]]$target_genes)
                supermod_regs_list <- c(supermod_regs_list, mregs)
                supermod_targs_list <- c(supermod_targs_list, mtargs)

            }

            reg_multiplicity <- sort(table(supermod_regs_list), decreasing = TRUE)
            targ_multiplicity <- sort(table(supermod_targs_list), decreasing = TRUE)
            multitab <- NULL
            for (i in sort(unique(c(reg_multiplicity, targ_multiplicity)), decreasing = TRUE)) {
                multitab <- rbind(multitab, c(i, sum(reg_multiplicity == i), paste(collapse = ", ", names(reg_multiplicity)[reg_multiplicity ==
                  i]), sum(targ_multiplicity == i), paste(collapse = ", ", names(targ_multiplicity)[targ_multiplicity == i])))
            }
            colnames(multitab) <- c("multiplicity", "number-of-regs", "reg-name-list", "number-of-targs", "targ-name-list")

            # Multipliticy Table tittle
            write(paste0("<table style='width:100%' bgcolor='gray'><tr><td><h1>", "Multiplicity Table", "</h1></td></tr></table><br>\n"),
                paste0(indexpageinfo$htmldir, foldername_p, "/", indexpageinfo$indexpath), append = TRUE)
            # write multitab table to index.html
            write(table2html(multitab), paste0(indexpageinfo$htmldir, "/", foldername_p, "/", indexpageinfo$indexpath), append = TRUE)

            alllabels <- responder[keepsamps]
            samps2pheno <- alllabels
            samps2pheno[which(alllabels == phenotype_class_vals_label[2])] <- phenotype_class_vals[2]
            samps2pheno[which(alllabels == phenotype_class_vals_label[1])] <- phenotype_class_vals[1]

            nonrespond_idxs <- names(samps2pheno)[which(samps2pheno == phenotype_class_vals[1])]
            responder_idxs <- names(samps2pheno)[which(samps2pheno == phenotype_class_vals[2])]

            rawrunmoddata <- list(regulators = names(reg_multiplicity), target_genes = names(targ_multiplicity))

            message("Generating raw graph")
            rawsumm <- summarize_module(norm_expr_mat_keep, rawrunmoddata, name2idx, nonrespond_idxs, responder_idxs)

            # If we have failed generating the raw full graph,
            # we cant generate refined graphs so we generate the summaries without the graphs
            # and we finish.

            if (rawsumm$cut){
                rawsummary(indexpageinfo,rawrunmoddata,rawsumm, norm_expr_mat_keep, outdir,foldername_p,modmeth,cut=TRUE)
                next
            }

            rawsummary(indexpageinfo,rawrunmoddata,rawsumm, norm_expr_mat_keep, outdir,foldername_p,modmeth)
            refinedmod <- unique(names(igraph::V(rawsumm$full_graph)))
            refinedrunmoddata <- list(regulators = refinedmod[(gene_info_df_keep[refinedmod, regulator_info_col_name] == "1")], target_genes = refinedmod[(gene_info_df_keep[refinedmod,
                regulator_info_col_name] == "0")])
            message("Generating refined graph")
            refinedsumm <- summarize_module(norm_expr_mat_keep, refinedrunmoddata, name2idx, nonrespond_idxs, responder_idxs)

            # summary of refined
            write(paste0("<table style='width:100%' bgcolor='gray'><tr><td><h1>", "Refined Modules Summary", "</h1></td></tr></table><br>\n"),
                paste0(indexpageinfo$htmldir, foldername_p, "/", indexpageinfo$indexpath), append = TRUE)

            pname <- paste(sep = ".", "igraphs.refined.graphs")
            grDevices::png(paste0(outdir, "/", foldername_p, "/imgs/", pname, ".png"), 1500, 750)
            graphics::par(mfrow = c(1, 3))

            mylayout <- return_layout_phenotype(refinedrunmoddata$regulators, refinedrunmoddata$target_genes, refinedsumm$nodesumm,
                rownames(norm_expr_mat_keep))


            try(plot_igraph(refinedsumm$full_graph, paste0(ncol(norm_expr_mat_keep), " Samples"), "black", mylayout))
            try(plot_igraph(refinedsumm$nonresp_graph, paste0(length(nonrespond_idxs), " Phenotype1"), "darkviolet", mylayout))
            try(plot_igraph(refinedsumm$respond_graph, paste0(length(responder_idxs), " Phenotype2"), "darkgoldenrod", mylayout))
            grDevices::dev.off()

            # write plot to index page
            write(paste0("<img src='", "imgs/", pname, ".png", "' alt='", pname, "' height='", 750, "' width='", 1500, "'> &emsp; <br>\n"),
                paste0(indexpageinfo$htmldir, foldername_p, "/", indexpageinfo$indexpath), append = TRUE)

            # Write tables for refinedsumm
            sortidxs <- sort(as.numeric(refinedsumm$nodesumm[, "t-pval"]), decreasing = FALSE, index.return = TRUE)$ix
            write_tables_all(refinedsumm$nodesumm[sortidxs, ], tabletype = paste0(modmeth, "_refined_nodesumm"), filestr = "data",
                html_idxs = seq_len(nrow(refinedsumm$nodesumm)), htmlinfo = indexpageinfo, extradir = paste0(foldername_p, "/"))

            sortidxs <- sort(as.numeric(refinedsumm$fulledgesumm[, "all.weights"]), decreasing = FALSE, index.return = TRUE)$ix
            write_tables_all(refinedsumm$fulledgesumm[sortidxs, ], tabletype = paste0(modmeth, "_refined_edgesumm"), filestr = "data",
                html_idxs = seq_len(nrow(refinedsumm$fulledgesumm)), htmlinfo = indexpageinfo, extradir = paste0(foldername_p, "/"))

            # Write raw and refined r object nodesumm, fulledgesumm, full_graph, respond_graph, nonresp_graph
            saveRDS(refinedsumm, file = paste0(outdir, "/", foldername_p, "/refinedsumm.rds"))

            tryCatch({

                # Generate the htmls for every supermodule from the refinedsumm.rds
                html_from_graph(gpath = paste0(outdir, "/", foldername_p, "/refinedsumm.rds"), wpath = paste0(outdir, "/", foldername_p),
                  user_mode = FALSE, dataset = norm_expr_mat_keep[allregs, ])

            }, error = function(x) {
                message("Html summary could not be generated for this supermodule")
            })

        }

        }

    }

    if (length(duplas) > 1) {

        message("Generating joint heatmap")

        # Generate the combined heatmap

        joint_clusters <- gen_heatmap(ObjectList, c_module_membership_list, c_allstats,
                                      imgdir=imgdir,outdir=outdir,
                                      indexpageinfo=indexpageinfo, cmp=TRUE, modmeth="combined")


        message('Joint cluster generated')

        #Read the combined sigmodules txt
        cmb_txt <- utils::read.delim(paste0(outdir,'/txts/sigmodules_combined.txt'),sep='\t',header=FALSE)
        cmb_txt_n <- gsub('[a-zA-Z]+.[0-9].mod.','',cmb_txt[,1])

        #Create hash table to store modname-cluster pairs
        clust_hash <- hash::hash()

        #Assign each modnumber to the cluster it is assigned to.
        for (i in seq_along(joint_clusters$clusters)){
            clustnames <- gsub('[a-zA-Z]+.[0-9].','',joint_clusters$clusters[[i]])
            for (name in clustnames){
                clust_hash[[name]] <- i
            }
        }

        #Add that info into the sigmodules_combined txt
        cmb_txt <- cbind(cmb_txt,sapply(cmb_txt_n,function(x){
            if (!identical(clust_hash[[x]],NULL)) clust_hash[[x]] else -1
            }))

        colnames(cmb_txt) <- c('Modname','Cluster')
        #Write back
        utils::write.table(cmb_txt,paste0(outdir,'/txts/sigmodules_combined.txt'),
                    sep='\t',quote=FALSE,row.names=FALSE)
    }
}


#' @export
#' @rdname runrewiring
#' @param module_membership_list hash list containing the pair module-significance
#' from the rewiring test.
#' @param allstats array containing information about module's rewiring.
#' @param imgdir path to image directory.
#' @param outdir path to the output directory.
#' @param indexpageinfo list containing information about index.html page.
#' @param cmp boolean indicating the joint heatmap case.
#' @param modmeth method of evaluation from LINKER_run.
#' @param i integer containing the dupla being evaluated.
gen_heatmap <- function(ObjectList, module_membership_list, allstats, imgdir, outdir, indexpageinfo, cmp = FALSE, modmeth = "", i=0) {

    # calculate overlap between signif mod sets
    fisher_tbl <- NULL
    fisher_cols <- c("user_gene_set", "property_gene_set", "universe_count", "user_count", "property_count", "overlap_count",
                     "pval")

    if (cmp) {
        # If compare mode enable, select the max of dataset's genes as universe_size
        universe_size <- max(vapply(ObjectList$datasets, function(x) nrow(x$norm_expr_mat_keep), FUN.VALUE = 1))
    } else {
        universe_size <- nrow(ObjectList$datasets[[i]]$norm_expr_mat_keep)
    }

    all_modules <- names(module_membership_list)

    methods::show(paste(c("Significant Modules: ", all_modules)))
    if (!length(all_modules))
        return(message("No significant modules found."))

    # save significant modules as .txt
    if (cmp) {
        sigmodules_p <- paste0("txts/sigmodules_", modmeth, ".txt")
    }else{
        sigmodules_p <- paste0("txts/sigmodules_", paste(modmeth, i, sep = "_"), ".txt")
    }
    utils::write.table(all_modules, file = paste(outdir, sigmodules_p, sep = "/"), quote = FALSE, sep = "\n", row.names = FALSE,
                       col.names = FALSE)

    module_pairs <- gtools::combinations(length(all_modules), 2, all_modules, repeats = TRUE)

    # fisher test (hipergeometric distribution)

    pvals <- NULL
    for (pair_idx in seq_len(nrow(module_pairs))) {

        mod1genes <- module_membership_list[[module_pairs[pair_idx, ][1]]]
        mod2genes <- module_membership_list[[module_pairs[pair_idx, ][2]]]
        stats <- c(module_pairs[pair_idx, ][1], module_pairs[pair_idx, ][2], universe_size, length(mod1genes), length(mod2genes),
                   length(intersect(mod1genes, mod2genes)))
        # show(stats)

        # building contigency table.
        contig_tbl <- as.table(matrix(c(length(intersect(mod1genes, mod2genes)),
                                        length(setdiff(mod1genes, mod2genes)),
                                        length(setdiff(mod2genes, mod1genes)),
                                        universe_size - length(mod2genes) - length(mod1genes) + length(intersect(mod1genes, mod2genes))),
                                        ncol = 2, byrow = TRUE))
        # show(contig_tbl)
        res <- stats::fisher.test(contig_tbl, alternative = "g")
        fisher_tbl <- rbind(fisher_tbl, c(stats, res$p.value))
        pvals <- c(pvals, res$p.value)

    }

    # generate txt file with pvals
    colnames(fisher_tbl) <- fisher_cols
    fisher_tbl <- as.data.frame(fisher_tbl)
    fisher_tbl$pval <- signif(log10(as.numeric(fisher_tbl$pval)) * -1, 3)
    # Maximum pval is 300
    fisher_tbl$pval[fisher_tbl$pval > 300] <- 300

    # Bonferroni correction (to implement)

    # get the max(fisher_tbl$pval) for scaling heatmap colors
    pval_max <- ceiling(max(fisher_tbl$pval))

    simmat <- matrix(-0.1, length(all_modules), length(all_modules), dimnames = list(all_modules, all_modules))
    simmat[cbind(fisher_tbl$user_gene_set, fisher_tbl$property_gene_set)] <- fisher_tbl$pval
    simmat[cbind(fisher_tbl$property_gene_set, fisher_tbl$user_gene_set)] <- fisher_tbl$pval

    rownames(simmat) <- gsub(paste0("mod."), "", rownames(simmat))
    colnames(simmat) <- gsub(paste0("mod."), "", colnames(simmat))


    pvc_result <- pvclust::pvclust(simmat, method.dist = "cor", method.hclust = "average", nboot = 1000)

    clusters <- pvclust::pvpick(pvc_result)

    if (cmp){
        myplotname <- paste0("mod_sim.", modmeth)
    }else{
        myplotname <- paste0("mod_sim.", modmeth, ".", i)
    }

    grDevices::png(paste0(imgdir, myplotname, ".dendro.png"), width = 8 * 300, height = 4 * 300, res = 300, pointsize = 8)
    plot(pvc_result, hang = -1, cex = 1)
    pvclust::pvrect(pvc_result, alpha = 0.9)
    grDevices::dev.off()

    # create heatmap plot
    if (length(all_modules) > 1) {

        # Generate the color palette
        my_palette <- (grDevices::colorRampPalette(c("darkgrey", "yellow", "green")))(n = 299)

        # Generate a const to compensate pvals
        compensate_const <- pval_max/300

        # Fix a lower bound to 0.01 for darkgrey and apply compensate_const if that lower bound is not reached
        lower_bound <- max(2, 4.9 * compensate_const)

        # Generate new_compensate_const based on the lower_bound
        new_compensate_const <- lower_bound/4.9

        # Correct higherbound based on the new compensate const
        higher_bound <- 199 * new_compensate_const

        # for darkgrey, yellow and green respectively
        my_col_breaks <- c(seq(0, lower_bound, length = 100), seq(lower_bound + 0.1, higher_bound, length = 100), seq(higher_bound +
                                                                                                                          1, 300 * new_compensate_const, length = 100))

        grDevices::png(paste0(imgdir, myplotname, ".heatm.png"), width = 8 * 300, height = 8 * 300, res = 300, pointsize = 8)
        row_order <- labels(stats::as.dendrogram(pvc_result$hclust))
        # show(row_order)
        heatm <- simmat[row_order, row_order]
        gplots::heatmap.2(heatm, Rowv = FALSE, Colv = FALSE, scale = "none", col = my_palette, breaks = my_col_breaks, dendrogram = "none",
                          cellnote = round(heatm, 0), trace = "none", density.info = "none", notecol = "black", margins = c(8, 8), cexRow = 1.5,
                          cexCol = 1.5, lmat = rbind(c(4, 4), c(2, 1), c(0, 3)), lwid = c(1, 4), lhei = c(1, 8, 0.1))
        grDevices::dev.off()

        # write plots to index page
        write(paste0("<table style='width:100%' bgcolor='gray'><tr><td><h1>", paste0("Rewiring Summary for Dataset", i, " using ", modmeth), "</h1></td></tr></table><br>\n"),
            paste0(indexpageinfo$htmldir, indexpageinfo$indexpath), append = TRUE)
        write(paste0("<img src='", indexpageinfo$imgstr, myplotname, ".dendro.png", "' alt='", myplotname, "' height='", 300,
                     "' width='", 600, "'> &emsp; <br>\n"), paste0(indexpageinfo$htmldir, indexpageinfo$indexpath), append = TRUE)
        write(paste0("<img src='", indexpageinfo$imgstr, myplotname, ".heatm.png", "' alt='", myplotname, "' height='", 600,
                     "' width='", 600, "'> &emsp; <br>\n"), paste0(indexpageinfo$htmldir, indexpageinfo$indexpath), append = TRUE)

    } else warning("Unique module found, heatmap can not be generated")
    # end module simularity if

    # write all stats table
    sortidxs <- sort(as.numeric(allstats[, "revised-pvalue"]), decreasing = FALSE, index.return = TRUE)$ix

    if (cmp){
        write_tables_all(allstats[sortidxs, ], tabletype = paste0(modmeth, "_mod_rewiring_scores"), filestr = "data", html_idxs = seq_len(nrow(allstats)),
                         htmlinfo = indexpageinfo)
    }else{
        write_tables_all(allstats[sortidxs, ], tabletype = paste0(paste(modmeth, i, sep = "_"), "_mod_rewiring_scores"), filestr = "data", html_idxs = seq_len(nrow(allstats)),
                         htmlinfo = indexpageinfo)
    }


    # write module similarity table
    sortidxs <- sort(as.numeric(fisher_tbl[, "pval"]), decreasing = FALSE, index.return = TRUE)$ix

    if (cmp){
        write_tables_all(fisher_tbl[sortidxs, ], tabletype = paste0(modmeth, "_sigmod_overlap"), filestr = "data", html_idxs = seq_len(nrow(fisher_tbl)),
                         htmlinfo = indexpageinfo)
    }else{
        write_tables_all(fisher_tbl[sortidxs, ], tabletype = paste0(paste(modmeth, i, sep = "_"), "_sigmod_overlap"), filestr = "data", html_idxs = seq_len(nrow(fisher_tbl)),
                         htmlinfo = indexpageinfo)
    }


    return(clusters)
}
#' @export
#' @rdname runrewiring
#' @param rawrunmoddata list containing regulators and target genes from multiplicity table.
#' @param rawsumm list containing generated graph related information and cut variable.
#' @param foldername_p path to the current foldername.
#' @param norm_expr_mat_keep gene expression matrix from ObjectList.
#' @param cut boolean to avoid errors from previous generated graphs.
rawsummary <- function(indexpageinfo, rawrunmoddata,rawsumm, norm_expr_mat_keep, outdir,foldername_p,modmeth,cut = FALSE) {

    if (!cut) {
        # Variable to check if previously graphs were generated, to avoid errors.

        # summary of raw tittle
        write(paste0("<table style='width:100%' bgcolor='gray'><tr><td><h1>", "Raw Modules Summary", "</h1></td></tr></table><br>\n"),
              paste0(indexpageinfo$htmldir, foldername_p, "/", indexpageinfo$indexpath), append = TRUE)

        pname <- paste(sep = ".", "igraphs.raw.full_graph")
        grDevices::png(paste0(outdir, "/", foldername_p, "/imgs/", pname, ".png"), 1500, 750)
        mylayout <- return_layout_phenotype(rawrunmoddata$regulators, rawrunmoddata$target_genes, rawsumm$nodesumm, rownames(norm_expr_mat_keep))

        try(plot_igraph(rawsumm$full_graph, paste0(ncol(norm_expr_mat_keep), " Samples"), "black", mylayout))
        grDevices::dev.off()

        # write plot to index page
        write(paste0("<img src='", "imgs/", pname, ".png", "' alt='", pname, "' height='", 750, "' width='", 1500, "'> &emsp; <br>\n"),
              paste0(indexpageinfo$htmldir, foldername_p, "/", indexpageinfo$indexpath), append = TRUE)
    }

    # Write tables for rawsumm
    sortidxs <- sort(as.numeric(rawsumm$nodesumm[, "t-pval"]), decreasing = FALSE, index.return = TRUE)$ix
    write_tables_all(rawsumm$nodesumm[sortidxs, ], tabletype = paste0(modmeth, "_raw_nodesumm"), filestr = "data", html_idxs = seq_len(nrow(rawsumm$nodesumm)),
                     htmlinfo = indexpageinfo, extradir = paste0(foldername_p, "/"), glossarypath = "../glossary.html")

    sortidxs <- sort(as.numeric(rawsumm$fulledgesumm[, "all.weights"]), decreasing = FALSE, index.return = TRUE)$ix
    write_tables_all(rawsumm$fulledgesumm[sortidxs, ], tabletype = paste0(modmeth, "_raw_edgesumm"), filestr = "data",
                     html_idxs = seq_len(nrow(rawsumm$fulledgesumm)), htmlinfo = indexpageinfo, extradir = paste0(foldername_p, "/"),
                     glossarypath = "../glossary.html")

    # Write raw and refined r object nodesumm, fulledgesumm, full_graph, respond_graph, nonresp_graph
    if (!cut)
        saveRDS(rawsumm, file = paste0(outdir, "/", foldername_p, "/rawsumm.rds"))

}
