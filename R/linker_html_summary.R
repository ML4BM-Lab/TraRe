#' Create HTML summary
#'
#' Given results of Linker runs and other annotations, builds an interconnected
#' html website that summarizes graph edges by support and ChIP evidence (if
#' provided).
#'
#' @param rfiles, character vector of the string file names where results of linker
#'      runs have been saved using saveRDS, e.g., linkerResult <- LINKER_run(...),
#'      saveRDS(linkerResult, file = "StringFileName").
#' @param tagstr, a string name to tag related html summary results
#' @param mapfile, a string name of file that contains identifier and regulator status information, must contain the columns
#'      "uniq_isos" - contains name of rows in LINKER_run `lognorm_est_counts` input matrix
#'      "iso_ensts" - transcript ids, but can be left blank
#'      "iso_ensgs" - contains gene identifiers that match identifiers available in the evidfile
#'      "iso_ensgvs" - versioned gene identifiers, can be left blank
#'      "iso_gnames" - gene name for display
#'      "iso_descs" - extra information about row, can be left blank
#'      "filtered" - whether to include row, can be left blank
#'      "regulator" - whether row is regulator or target (has values 1 or 0 respectively)
#' @param outdir, a string name of directory location to save html
#' @param evidfile, a string name of file that contain ChIP evidence information, must be a matrix where
#'      the column names are the "iso_ensgs" names of the targets
#'      the row names are the "iso_ensgs" names of the regulators
#'      the matrix values are
#'          -1: meaning that information is missing, the regulator was not chipped, or the gene names were not mapped
#'          0: meaning that the ChIP'ed regulator did not have a peak in the targets genomic region (gene body +/- 20KB)
#'          positive integer: number of peaks of the ChIP'ed regulator in the targets genomic region (gene body +/- 20KB)
#' @return allsummaries, a data frame that contains for each module and graph
#'      the number of graph edges at each possible level of support and the
#'      percentage of cumulative edges with each type of ChIP evidence
#'
#' @examples
#'
#'    ## We are going to use files from the example folder.
#'    ## `create_html_summary()` only needs the paths, so thats what
#'    ## we are giving it.
#'
#'    evidpath <- paste0(system.file("extdata",package="TraRe"),'/ChIP',
#'                       '/Tumor_OV50_intersectBed.weighted_evidence.txt')
#'
#'    rfiles <- c(paste0(system.file("extdata",package="TraRe"),'/ChIP',
#'                       '/Tumor_OV50.tar8855_reg638.VBSR.m100_b10.rds'))
#'
#'    tagstr <- "Tumor_OV50.tar8855_reg638"
#'
#'    mapfile <- paste0(system.file("extdata",package="TraRe"),'/ChIP',
#'                       '/Tumor_OV50.gene_info.txt')
#'
#'    ## We are going to create a folder for this example
#'
#'    \dontrun{
#'    dir.create(paste0(getwd(),'/summaryfolder'),showWarnings=FALSE)
#'    outdir <- paste0(getwd(),'/summaryfolder/')
#'
#'
#'    create_html_summary(rfiles,tagstr,mapfile,outdir,evidpath)
#'    }
#' @export
create_html_summary <- function(rfiles,tagstr,mapfile,outdir = getwd(),evidfile){
  
  runinfo <- list()
  runinfo$rfiles <- rfiles
  runinfo$tagstr <- tagstr
  runinfo$mapfile <- mapfile
  runinfo$outdir <- outdir
  runinfo$evidfile <- evidfile
  
  #################### load information on gene identifiers ####################
  
  methods::show(paste0("Loading id conversion table: ", runinfo$mapfile, "..."))
  iso_table <- as.matrix(utils::read.table(runinfo$mapfile, header = TRUE, sep = "\t",
                                           quote = ""))
  rownames(iso_table) <- iso_table[, "uniq_isos"]
  all_reg_genes <- iso_table[which(iso_table[, "regulator"] == 1),
                             "iso_gnames"]
  all_tar_genes <- iso_table[which(iso_table[, "regulator"] == 0),
                             "iso_gnames"]
  runinfo$nregs <- length(all_reg_genes)
  runinfo$ntargets <- length(all_tar_genes)
  
  ################ create overall summary for html index page #################
  indexpath <- paste(sep = ".", "index",  runinfo$tagstr, "html")
  htmlinfo <- linker_create_index_page(outdir = runinfo$outdir,
                                       runtag = runinfo$tagstr,
                                       indexpath = indexpath,
                                       codedir = paste0(system.file("extdata",package="TraRe"),'/'))
  
  # write gene info stats
  write(paste0("<br>", length(all_tar_genes), " target isoforms covering ",
               length(unique(all_tar_genes)), " genes.<br>"),
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
  write(paste0(length(all_reg_genes), " regulator isoforms covering ",
               length(unique(all_reg_genes)), " genes.<br><br>"),
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
  
  ################# load information on chip evidence, if exists ###############
  weighted_chip_evidence <- NULL
  if (file.exists(runinfo$evidfile)){
    methods::show(paste0("Loading chip evidence table: ", runinfo$evidfile, "..."))
    weighted_chip_evidence <- as.matrix(utils::read.table(runinfo$evidfile, header = TRUE,
                                                          row.names = 1, sep = "\t",
                                                          quote = ""))
    methods::show("chip evidence regulators: ")
    showfirstlast(rownames(weighted_chip_evidence))
    methods::show("chip evidence regulators: ")
    showfirstlast(colnames(weighted_chip_evidence))
    binary_chip_evidence <- weighted_chip_evidence
    binary_chip_evidence[binary_chip_evidence > 1] <- 1
    binary_chip_summary <- table(as.numeric(binary_chip_evidence))
    methods::show(signif(binary_chip_summary / ( dim(weighted_chip_evidence)[1] *
                                                   dim(weighted_chip_evidence)[2]) * 100, 3))
    
    nonzeroregs <- sum(rowSums(binary_chip_evidence) >= 0)
    nonzerotargs <- sum(colSums(binary_chip_evidence) !=
                          dim(binary_chip_evidence)[1] * -1)
    
    write(paste0(nonzeroregs, " regulators and ", nonzerotargs,
                 " targets with possible chip evidence."),
          file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
    write(paste0(binary_chip_summary["1"], " (",
                 signif(binary_chip_summary["1"] / (nonzeroregs *
                                                      nonzerotargs) * 100, 3),
                 "%) of ", nonzeroregs, "*", nonzerotargs,
                 " possible chip edges have at least one peak.<br>"),
          file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
    bin_summ <- cumSumTable(cbind(as.numeric(binary_chip_evidence),
                                  1 - as.numeric(binary_chip_evidence)))
    myrnames <- c("Peaks", "noPeak", "NA")
    bin_summ <- cbind(myrnames, bin_summ)
    colnames(bin_summ) <- c("RowType", "nCount", "cumCount", "%RowPeaks",
                            "%RowNoPeak", "%RowNAs")
    rownames(bin_summ) <- myrnames
    htmlstr <- table2html(bin_summ)
    write(htmlstr, file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = TRUE)
  }
  
  ###################### for each saved linker run result ########################
  allsummaries <- NULL
  
  for (rfile in runinfo$rfiles){
    
    methods::show(paste0("Processing ", rfile, "..."))
    rundata <- readRDS(rfile)
    
    runinfo$nboots <- 1
    runtypestr <- "single_gene"
    mymodmeths <- ls(rundata$graphs)
    if (!is.null(rundata$raw_results)){
      runtypestr <- "multi_module"
      mymodmeths <- ls(rundata$modules)
    }
    
    for ( mymodmeth in mymodmeths){
      # mymodmeth <- mymodmeths[1]
      
      graphslist <- ls(rundata$graphs)
      if (runtypestr != "single_gene"){
        graphslist <- ls(rundata$graphs[[mymodmeth]])
        runinfo$nboots <- length(rundata$raw_results[[mymodmeth]]$bootstrapResults)
      }
      
      # for each graph
      for ( graphmeth in graphslist){
        # graphmeth="VBSR"
        
        graphstr <- paste(sep = ".", runinfo$tagstr, mymodmeth, runtypestr,
                          graphmeth)
        
        if (runtypestr == "single_gene"){
          methods::show(paste0("****************Single Gene Network ", graphmeth))
          rungraphs <- list(rundata$graphs[[graphmeth]])
        } else {
          methods::show(paste0("****************", mymodmeth, " Method: ", graphmeth))
          rungraphs <- rundata$graphs[[mymodmeth]][[graphmeth]]
        }
        final_summary <- linker_summarize_rungraphs(rungraphs = rungraphs,
                                                    iso_table = iso_table,
                                                    weighted_chip_evidence = weighted_chip_evidence,
                                                    graphstr = graphstr,
                                                    runinfo = runinfo,
                                                    htmlinfo = htmlinfo)
        labeled_table <- cbind(runtypestr, mymodmeth, graphmeth, final_summary)
        allsummaries <- rbind(allsummaries, labeled_table)
        resultspath <- paste0(htmlinfo$txtstr,
                              runinfo$tagstr, ".all_summaries.txt")
        methods::show(paste0("Writing table: ", resultspath))
        utils::write.table(allsummaries, paste0(htmlinfo$htmldir, resultspath),
                           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      } # end graphmeth
    } # end modmeth
  } # end rfile
  return(allsummaries)
}

