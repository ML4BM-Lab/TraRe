#' Prepare rewiring data for running the method.
#'
#' Prepare neccessary files for running `runrewiring()`
#'
#' @param name Desired name of the folder which is generated. The chosen
#' threshold will be `paste()` to the folder's name.
#' @param TraReObj TrareObj generated during preprocessing step.
#' @param linker_output linker output file.
#' @param final_signif_thresh Significance threshold for the rewiring method. The lower the threshold, the restrictive the method.
#' @param orig_test_perms Initial permutations for first test (default: 100) .
#' @param retest_thresh Threshold if a second test is performed (default: 0.08) .
#' @param retest_perms Permutations if a second test is performed (default: 1000) .
#' @param use_graphs Boolean specifying the use of graphs as linkeroutput modules. (TRUE by default).
#' @param outdir Directory for the output folder to be located (default: tempdir())
#' @param nrcores Number of cores to run the parallelization within the rewiring test (default: 3).
#' @param last_cluster Boolean specifying whether to include the last_cluster in the rewiring or not. (default: FALSE)
#'
#'
#' @return Return a list containing: LINKER's output, expression matrix, boolean array from phenotype file,
#' array containing number of c(R,NR) samples, significance threshold and output directory.
#'
#' @examples
#'
#' ## Load the linker output
#' linker_output <- readRDS(paste0(system.file('extdata',package='TraRe'),
#'                                  '/linker_rewiring_example.rds'))
#'
#' ## Load the phenotype file
#' phenotype_p <- paste0(system.file('extdata',package='TraRe'),
#'                      '/phenotype_rewiring_example.txt')
#' phenotype <- read.delim(phenotype_p, row.names = 1)
#' colnames(phenotype) <- 'phenotype'
#' 
#' ## We will be using the same file we generated for the LINKER_run phase.
#' ## TraReObj <- readRDS(paste0(system.file('extdata',package='TraRe'),'/TraReObj.rds'))
#' 
#' ## Add the phenotype if it was not in the SummarizedExperiment when we did the preprocessing.
#' ## TraReObj <- rewiring_add_phenotype(TraReObj, phenotype)
#'
#' ## outdir <- system.file('extdata',package='TraRe')
#' ## prepared <- preparerewiring(name='example', linker_output= linker_output,
#' ##                            TraReObj = TraReObj, 
#' ##                            final_signif_thresh=0.05,
#' ##                            nrcores=1,outdir=outdir)
#' 
#' @export
preparerewiring <- function(name = "defaultname", linker_output = NULL, TraReObj = NULL, 
                            final_signif_thresh = 0.001,
                            orig_test_perms = 100, retest_thresh = 0.08, retest_perms = 1000,
                            use_graphs = TRUE, outdir = tempdir(), nrcores = 3, last_cluster = FALSE) {

    # checks
    if (is.null(linker_output)) {
        stop("linker_output file required.")
    }

    # check if TraReObj has been provided.
    if (is.null(TraReObj)) {
        stop("Trare object is required")

    }

    #check if phenotype is configured
    if (length(TraReObj@pheno) == 1){
        stop('phenotype has not been added, please refer to TraRe::rewiring_add_phenotype')
    }

    # Check for comparison mode
    linker_files <- 1
    if (linker_files > 1) {
        warning("Data comparison mode selected, only heatmap will be generated.")
    }

    # Create folder name
    foldername <- paste(name, paste(final_signif_thresh, collapse = "_"), sep = "_")

    # Concatenate with outdir path
    outdir <- paste(outdir, foldername, sep = "/")

    rewobjects <- list()
    rewobjects$datasets <- list()

    for (i in seq_along(linker_files)) {

        rewobject <- list()

        # retrieve data from TraRe object
        lognorm_est_counts <- TraReObj@lognorm_counts
        geneinfo <- rownames(lognorm_est_counts)[TraReObj@regulator_idx]
        # phenotype <- TraReObj@pheno
        regs <- geneinfo
        targs <- rownames(lognorm_est_counts)[TraReObj@target_idx]
        phenosamples <- colnames(lognorm_est_counts)[!is.na(TraReObj@pheno)]
	      phenotype <- as.logical(TraReObj@pheno[!is.na(TraReObj@pheno)])
        
        #generate list with name index
        name2idx <- seq_len(nrow(lognorm_est_counts))
        names(name2idx) <- rownames(lognorm_est_counts)

        # check if NR/R proportions are similar to ensure property functioning of the method.
        klzero <- sum(phenotype == 0, na.rm = TRUE)
        klone <- sum(phenotype == 1, na.rm = TRUE)

        if (use_graphs){
            #Format linkeroutput modules to graphs
            message("\nFrom here on, graphs will be taken into account for Rewiring and some of them may be dropped out.\nHence, it is recommended to take preparerewiring object's run data to proceed with further analysis\n")
            rundata <- graph_to_modules(linker_output, geneinfo)
        }

        if (min(klone, klzero)/max(klone, klzero) < 0.8) {
            warning(paste0("phenotype samples proportions imbalance ", toString(c(klzero, klone)), " (<80%)."))
        }

        #generate line for 
        NumRegsTargs <- c("NumRegs and NumTargs: [", length(regs), ",", length(targs), "]")
        message(NumRegsTargs)

        SampleNames <- c("\nSample Names: [", phenosamples, "]\n\nNumber of samples: ", length(phenosamples))
        message(SampleNames)

        class_counts <- as.numeric(table(phenotype))
        ClassPerCounts <- c("Class Per Counts: ", "(", paste(class_counts, collapse = ","), ")")
        message(ClassPerCounts)


        rewobject$rundata <- rundata
        rewobject$lognorm_est_counts <- lognorm_est_counts
        rewobject$final_signif_thresh <- final_signif_thresh
        rewobject$name2idx <- name2idx
        rewobject$regs <- regs
        rewobject$targs <- targs
        rewobject$pheno <- phenotype
        rewobject$phenosamples <- phenosamples
        rewobject$class_counts <- class_counts
        rewobjects$datasets[[i]] <- rewobject
    }

    rewobjects$outdir <- outdir
    rewobjects$orig_test_perms <- orig_test_perms
    rewobjects$retest_thresh <- retest_thresh
    rewobjects$retest_perms <- retest_perms
    rewobjects$NrCores <- nrcores
    rewobjects$last_cluster <- last_cluster

    # Create logfile

    logfile <- list(NumRegsTargs, SampleNames, ClassPerCounts)
    logfile <- vapply(logfile, FUN = paste0, collapse = "", FUN.VALUE = "")

    rewobjects$logfile <- logfile

    return(rewobjects)

}

# Helper function
graph_to_modules <- function(linkeroutput, geneinfo){

    ## The structure we want to get is
    ## linkeroutput$modules[[link_mode]][[graph_mode]][[num_module]]$(target or regs)

    #names of linkeroutput
    linkeroutput_names <- names(linkeroutput$modules)

    linkeroutput <- lapply(linkeroutput_names,function(x){

        graph_modes <- names(linkeroutput$graphs[[x]])

        selected <- graph_modes[1]

        ## Select VBSR if more than 1 is available

        if (length(graph_modes)>1 & 'VBSR'%in%graph_modes){
            selected <- 'VBSR'
        }

        module_list <- lapply(seq_along(linkeroutput$graphs[[x]][[selected]]),function(y){

            graph <- linkeroutput$graphs[[x]][[selected]][[y]]

            totgenes <- unique(names(igraph::V(graph)))
            regulators <- intersect(totgenes, geneinfo)

            if (identical(regulators,character(0))){

                message('Module number ',y,' has been deleted')
                return(NULL)

            }

            list(regulators = regulators,
                 target_genes = setdiff(totgenes, regulators),
                 bootstrap_idx = linkeroutput$modules[[x]][[y]]$bootstrap_idx)

        })

        #generate old index
        orig_index <- unlist(sapply(seq_along(module_list),function(x){

            if (!is.null(module_list[[x]])){return(x)}

        }))


        module_list <- Filter(Negate(function(X) {

            length(X) == 0

        }),module_list)

        return(list(modules=module_list,graphs=linkeroutput$graphs[[x]][[selected]][orig_index]))

    })

    names(linkeroutput) <- linkeroutput_names

    #initialize the final linkeroutut
    new_linkeroutput <- list(modules=list(),graphs=list())

    for (method in names(linkeroutput)){

        #assign the modules
        new_linkeroutput[['modules']][[method]] <- linkeroutput[[method]][['modules']]
        #assign the graphs
        new_linkeroutput[['graphs']][[method]] <- linkeroutput[[method]][['graphs']]
    }

    return(new_linkeroutput)

}
