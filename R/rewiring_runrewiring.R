#' GRN modules Rewiring method.
#'
#' Gene Regulatory Network modules Rewiring method.
#'
#' @param ObjectList Output from preparerewiring function.
#' @param orig_test_perms Initial permutations for first test (default: 100) .
#' @param retest_thresh Threshold if a second test is performed (default: 1000) .
#' @param retest_perms Permutations if a second test is performed (default: 0.08) .
#'
#' @return None. Generate an html report containing a hierarchical clustering of the rewired modules,
#' a correlation matrix containing the rewiring scores of the modules as aheatmap and all the rewired
#' modules drivers and targets statistical information about the rewiring process.
#'
#' @examples
#' \dontrun{
#' rewiringinput<-preparerewiring(name="example",linker_rds_v,expr_mat_v,gene_info_v,
#'                                phenotype_file_v,final_signif_thresh)
#' runrewiring(rewiringinput)
#' }
#'
#' @export
runrewiring<- function(ObjectList,orig_test_perms=100,retest_thresh=0.08,retest_perms=1000){

# helper directory

  codedir <- paste0(system.file("extdata",package="TraRe"),"/")
  methods::show(codedir)

  #Here we can start the new method.----

  # set up output html page, we use the first argv.
  indexpageinfo <- create_index_page(outdir = ObjectList[[1]]$outdir, runtag = "",
                                     codedir = codedir)
  imgdir <- paste0(indexpageinfo$htmldir, indexpageinfo$imgstr)

  set.seed(1)
  module_membership_list <- hash::hash()

  #we create rundata and combine the modules of both parsers.

  for (modmeth in names(ObjectList[[1]]$rundata$modules)) { #VBSR
       methods::show(paste(c("ModuleMethod", modmeth)))
       allstats <- NULL
       statsnames <- c("module-method", "module-index", "orig-pval",
                       "revised-pvalue", "num-targets", "num-regulators",
                       "regulator-names","target-names", "num-samples", "num-genes",
                       "num-class1", "num-class2")
      for (i in 1:length(ObjectList)){ #1:3 if extramodule stablished.

        modmeth_i<-paste(modmeth,i)
        methods::show(paste('VBSR:',i))
        rundata<-ObjectList[[i]]$rundata
        norm_expr_mat_keep<-ObjectList[[i]]$norm_expr_mat_keep
        keepsamps<-ObjectList[[i]]$keepsamps
        keeplabels<-ObjectList[[i]]$keeplabels
        class_counts<-ObjectList[[i]]$class_counts
        final_signif_thresh<-ObjectList[[i]]$final_signif_thresh

        for (mymod in 1:length(rundata$modules[[modmeth]])) {

            modregs <- unique(rundata$modules[[modmeth]][[mymod]]$regulators)
            modtargs <- unique(rundata$modules[[modmeth]][[mymod]]$target_genes)
            regnames <- paste(collapse = ", ", modregs)
            targnames <- paste(collapse = ", ", modtargs)
            keepfeats <- unique(c(modregs, modtargs))
            modmat <- t(norm_expr_mat_keep[keepfeats, keepsamps]) #this has to be modified.

            orig_pval <- dave_test(modmat, keeplabels + 1, perm = orig_test_perms) #this too
            new_pval <- orig_pval
            stats <- c(modmeth_i, mymod, signif(orig_pval, 3), signif(new_pval, 3),
                       length(modtargs), length(modregs), regnames,targnames, dim(modmat),
                       class_counts)
            if (orig_pval < retest_thresh | orig_pval == 1 | mymod %% 300 == 0) {
                methods::show(paste(c("ModNum and NumGenes", mymod, length(keepfeats))))
                result <- dave_test_pair_detail(modmat, keeplabels + 1, #this too
                                                perm = retest_perms)
                new_pval <- result$pval
                stats <- c(modmeth_i, mymod, signif(orig_pval, 3),
                           signif(new_pval, 3), length(modtargs),
                           length(modregs), regnames,targnames, dim(modmat), class_counts) #this too

                if (new_pval <= final_signif_thresh | new_pval == 1) {
                    # keep as significant
                    modname <- paste0(modmeth,'.',i, ".mod.", mymod) #we include the 'i' for the parse.
                    module_membership_list[[modname]] <- keepfeats
                }
            }
            allstats <- rbind(stats, allstats)
        } # end mymod
      }

      # generate txt file with pvals
      colnames(allstats) <- statsnames
      rownames(allstats) <- allstats[, "module-index"]
      allstats[, "orig-pval"] <- signif(as.numeric(allstats[, "orig-pval"]), 3)
      allstats[, "revised-pvalue"] <- signif(as.numeric(allstats[, "revised-pvalue"]), 3)

      # calculate overlap between signif mod sets
      fisher_tbl <- NULL
      fisher_cols <- c("user_gene_set", "property_gene_set", "universe_count",
                       "user_count", "property_count", "overlap_count", "pval")

      universe_size <- length(rownames(ObjectList[[1]]$norm_expr_mat_keep)) #we select the Promote one. (it contains all of the S+some of the supermodule)
      all_modules <- names(module_membership_list)
      methods::show(paste(c("Significant Modules: ", all_modules)))

      methods::show(length(all_modules))
      methods::show(all_modules)

      #save significant modules as .txt
      utils::write.table(all_modules,file=paste(ObjectList[[1]]$outdir,'sigmodules.txt',sep="/"),quote=F,sep="\n",row.names = F,col.names = F)

      module_pairs <- gtools::combinations(length(all_modules), 2, all_modules,
                                   repeats = TRUE)
      pvals <- NULL
      for (pair_idx in 1:nrow(module_pairs)) {
          mod1genes <- module_membership_list[[module_pairs[pair_idx, ][1]]]
          mod2genes <- module_membership_list[[module_pairs[pair_idx, ][2]]]
          stats <- c(module_pairs[pair_idx, ][1], module_pairs[pair_idx, ][2],
                       universe_size, length(mod1genes), length(mod2genes),
                       length(intersect(mod1genes, mod2genes)))
          #show(stats)

          contig_tbl <- as.table(matrix(c(length(intersect(mod1genes, mod2genes)),
                                          length(setdiff(mod1genes, mod2genes)),
                                          length(setdiff(mod2genes, mod1genes)),
                                          universe_size - length(mod2genes) -
                                          length(mod1genes) +
                                          length(intersect(mod1genes, mod2genes))),
                                        ncol = 2, byrow = TRUE))
          #show(contig_tbl)
          res <- stats::fisher.test(contig_tbl, alternative = "g")
          fisher_tbl <- rbind(fisher_tbl, c(stats, res$p.value))
          pvals <- c(pvals, res$p.value)
      }

      # generate txt file with pvals
      colnames(fisher_tbl) <- fisher_cols
      fisher_tbl <- as.data.frame(fisher_tbl)
      fisher_tbl$pval <- signif(log10(as.numeric(fisher_tbl$pval)) * -1, 3)
      fisher_tbl$pval[fisher_tbl$pval > 300] <- 300

      simmat <- matrix(-0.1, length(all_modules), length(all_modules),
                       dimnames = list(all_modules, all_modules))
      simmat[cbind(fisher_tbl$user_gene_set, fisher_tbl$property_gene_set)] <- fisher_tbl$pval
      simmat[cbind(fisher_tbl$property_gene_set, fisher_tbl$user_gene_set)] <- fisher_tbl$pval

      rownames(simmat) <- gsub(paste0("mod."), "", rownames(simmat))
      colnames(simmat) <- gsub(paste0("mod."), "", colnames(simmat))

      pvc_result <- pvclust::pvclust(simmat, method.dist = "cor", method.hclust = "average",
                            nboot = 1000)
      clusters <- pvclust::pvpick(pvc_result)
      myplotname <- paste0("mod_sim.", modmeth)
      grDevices::png(paste0(imgdir, myplotname, ".dendro.png"),
          width = 8 * 300,
          height = 4 * 300,
          res = 300, pointsize = 8)
      plot(pvc_result, hang = -1, cex = 1)
      pvclust::pvrect(pvc_result, alpha = 0.9)
      grDevices::dev.off()

      # create heatmap plot
      if (length(all_modules) > 1) {
          myplotname <- paste0("mod_sim.", modmeth)
          my_palette <- grDevices::colorRampPalette(c("darkgrey",
                                           "yellow", "green"))(n = 299)
          my_col_breaks <- c(seq(0, 4.9, length = 100), # for darkgrey
                             seq(5, 199, length = 100), # for yellow
                             seq(200, 300, length = 100)) # for green

          grDevices::png(paste0(imgdir, myplotname, ".heatm.png"),
              width = 8 * 300,
              height = 8 * 300,
              res = 300, pointsize = 8)
          row_order <- labels(stats::as.dendrogram(pvc_result$hclust))
          #show(row_order)
          heatm <- simmat[row_order, row_order]
          gplots::heatmap.2(heatm,
                    Rowv = FALSE,
                    Colv = FALSE,
                    scale = "none",
                    col = my_palette,
                    breaks = my_col_breaks,
                    dendrogram = "none",
                    cellnote = round(heatm, 0),
                    trace = "none",
                    density.info = "none",
                    notecol = "black",
                    margins = c(8, 8),
                    cexRow = 1.5,
                    cexCol = 1.5,
                    lmat = rbind(c(4, 4), c(2, 1), c(0, 3)),
                    lwid = c(1, 4),
                    lhei = c(1, 8, .1)
          )
          grDevices::dev.off()

          # write plot to index page
          write(paste0("<img src='", indexpageinfo$imgstr, myplotname, ".dendro.png",
                       "' alt='", myplotname,
                       "' height='", 300, "' width='", 600, "'> &emsp; <br>\n"),
                paste0(indexpageinfo$htmldir, indexpageinfo$indexpath),
                append = T)
          write(paste0("<img src='", indexpageinfo$imgstr, myplotname, ".heatm.png",
                       "' alt='", myplotname,
                       "' height='", 600, "' width='", 600, "'> &emsp; <br>\n"),
                paste0(indexpageinfo$htmldir, indexpageinfo$indexpath),
                append = T)
      } # end module simularity if

      # write all stats table
      sortidxs <- sort(as.numeric(allstats[, "revised-pvalue"]),
                       decreasing = F, index.return = T)$ix
      write_tables_all(allstats[sortidxs, ],
                       tabletype = paste0(modmeth, "_mod_rewiring_scores"),
                       filestr = "data", html_idxs = 1:dim(allstats)[1],
                       htmlinfo = indexpageinfo)

      # write module similarity table
      sortidxs <- sort(as.numeric(fisher_tbl[, "pval"]),
                       decreasing = F, index.return = T)$ix
      write_tables_all(fisher_tbl[sortidxs, ],
                       tabletype = paste0(modmeth, "_sigmod_overlap"),
                       filestr = "data", html_idxs = 1:dim(fisher_tbl)[1],
                       htmlinfo = indexpageinfo)
  }
}
