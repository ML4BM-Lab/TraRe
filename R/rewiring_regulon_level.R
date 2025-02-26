#' Perform the rewiring at the regulon level
#'
#' Perform the rewiring test to regulons using a hypergeometric test.
#' @param linker_output Output from LINKER_run function (recommended with 50 bootstraps).
#' @param TraReObj the TrareObj generated during preprocessing step before GRN inference.
#' @param fpath Desired path for the rewiring file to be generated. If fpath not provided, it will create a rewiring module list file in the current directory with the name "rewiring_gene_level_fs_<final_signif_thresh>.txt". If file already exists it will skip the fast rewiring step.
#' @param rewired Bool indicating whether the analysis is to be done regulons from rewired modules (default= TRUE)
#' @param final_signif_thresh Significance threshold for the rewiring method. The lower the threshold, the more restrictive the method. Default set to 0.05.
#' @param fpath Path of the output file from rewired_gene_level function with rewired modules list.
#' @param outdir Directory for the output folder to be located (default: tempdir()).
#'
#' @return Return an object list with a data.frame for each regulon with filtered targets by pvalue, multiplicity and the igraph object.
#' 
#' @examples 
#' ## We will be using an output file we generated with LINKER_run and the 
#' ## same expression datainput  object.
#' 
#' ## linker_output <- readRDS(paste0(system.file('extdata',package='TraRe'),
#' ##                                     'linkeroutput_rewiring_example.rds'))
#'
#' ## TraReObj <- readRDS(paste0(system.file('extdata',package='TraRe'),
#' ##                                                        '/TraReObj.rds'))
#' 
#' ## Add the phenotype to TraReObj if it was not before.
#' ## TraReObj <- rewiring_add_phenotype(TraReObj, phenotype)
#' 
#' ## Select directory for output file
#' ## outdir <- system.file('extdata',package='TraRe')
#
#' ## regulons <- rewiring_regulon_level(linker_output = linker_output, 
#' ##                                    TraReObj = TraReObj,
#' ##                                    rewired = TRUE,
#' ##                                    final_signif_thresh = 0.05,
#' ##                                    outdir = outdir)
#' @export

rewiring_regulon_level <- function(linker_output,
                                   TraReObj,
                                   fpath='',
                                   rewired = TRUE,
                                   final_signif_thresh = 0.05,
                                   outdir=tempdir()){

  ### Let's work
  #Prepare rewiring creating object
  preparedrewiring <- TraRe::preparerewiring(TraReObj = TraReObj,
                                             linker_output= linker_output,
                                             final_signif_thresh = final_signif_thresh,
                                             outdir=outdir)
  
  ## Get regulons
  regulons_all_modules <- get_regulons(preparedrewiring,
                                       fpath = fpath,
                                       rewired = TRUE)

  ## Rewire regulons
  finals_adj <- rewiring_regulon(regulons_all_modules,
                                 preparedrewiring,
                                 final_signif_thresh)

  ## Merge regulons
  regulons_uniq <- merge_regulons(finals_adj)

  ## Filter regulons
  order_filtr_regulons <- filter_regulons(regulons_uniq)

  return(order_filtr_regulons)
}

### Helper functions
## function get_regulons:
# Takes as input prepareobj object, the number of rewired modules' path (if wanted to get the regulons from only those modules).
# Retrieves an object with the regulon subgraphs of the modules of linkeroutput and the boostrap idx
get_regulons <- function(prepareobj,
                         fpath = '',
                         rewired = rewired){

  if (rewired){
    
    if (!file.exists(fpath)){
      
      #define fpath (check default value)
      if (fpath==''){
        final_signif_thresh <- prepareobj$datasets[[1]]$final_signif_thresh
        fpath <- paste0(getwd(),'/rewiring_gene_level','_fs_',final_signif_thresh,'.txt')
        
      }else{
        
        message('Adding significant threshold to output file name')
        
        #(check extension to add final sig th)
        if (substr(fpath,nchar(fpath)-3,nchar(fpath)-3)=='.'){
          fpath <- paste0(substr(fpath,1,nchar(fpath)-4),'_fs_',final_signif_thresh,'.txt')
        }else{
          fpath <-paste0(fpath,'_fs_',final_signif_thresh,'.txt')
        }
        
      }
      
      message('Performing fast rewiring')
      #Do the fast rewiring
      rew_list_names <- fast_rew(prepareobj)
      
      #write rew_list_names
      utils::write.table(rew_list_names, file = fpath, quote = FALSE, sep = '\t',
                         row.names = FALSE, col.names = FALSE)
      message(paste0('Rewired modules list generated in ',fpath))
    }else{
      message('Fast rewiring file found to use reiwired module list')
    }
    #Read the rewired modules from the 50 bootstraps
    sigmod <- utils::read.delim(fpath, header = FALSE)[,1]
    # Rewired graphs
    message('Extracting regulons from rewired modules')
    graph_objs <- lapply(sigmod,function(x) prepareobj$datasets[[1]]$rundata$graphs$VBSR[[x]])
    #Get the bootstraps of the rewired graphs
    bootst_rewired_mods <- sapply(sigmod,function(x) prepareobj$datasets[[1]]$rundata$modules$VBSR[[x]]$bootstrap_idx)
  }else{
    message('Extracting all regulons')
    graph_objs <- prepareobj$datasets[[1]]$rundata$graphs$VBSR
  }

  #Get the regulons from all graphs
  regulons_all_modules <- lapply(seq_along(graph_objs),function(idx){
    g <- graph_objs[[idx]]

    # Regulators in module
    regs_mod <- igraph::V(g)[igraph::get.vertex.attribute(g,'type')]$name

    out <- lapply(regs_mod,function(reg){
      v <- igraph::V(g)[igraph::neighbors(g,reg)]
      regulons <- igraph::induced_subgraph(g,c(names(v),reg))
    })

    return(list(regulons=out,bootsts_idx=bootst_rewired_mods[idx]))
  })

  #Change regulons names
  names(regulons_all_modules) <- paste0("module_",sigmod)
  total_regulons <- sum(unlist(lapply(regulons_all_modules, function(x) length(x$regulons))))

  if (rewired) {
    message(paste(total_regulons,"regulons extracted from",length(sigmod), "rewired modules"))
  }else{message(paste(length(regulons_all_modules),"regulons extracted from all modules"))}

  return(regulons_all_modules)
}

## function rewiring_single_module
rewiring_single_module <- function(ObjectList,module,mymod=1){

  if (inherits(module,'igraph')){
    genes <- names(igraph::V(module))
  }
  #We are working with just 1
  i<-1
  #Needed info
  orig_test_perms<-ObjectList$orig_test_perms
  retest_thresh<-ObjectList$retest_thresh
  retest_perms<-ObjectList$retest_perms
  lognorm_est_counts <- ObjectList$datasets[[i]]$lognorm_est_counts
  pheno <- ObjectList$datasets[[i]]$pheno
  phenosamples <- ObjectList$datasets[[i]]$phenosamples

  #Now rewiring test!

  modregs <- intersect(ObjectList$datasets[[i]]$allregs,genes)
  modtargs <- intersect(ObjectList$datasets[[i]]$alltargs,genes)
  keepfeats <- unique(c(modregs, modtargs))
  modmat <- t(lognorm_est_counts[keepfeats, phenosamples])
  orig_pval <- TraRe::rewiring_test(modmat, pheno + 1, perm = orig_test_perms)
  new_pval <- orig_pval

  if (orig_pval < retest_thresh | orig_pval == 1 | mymod %% 300 == 0) {
    result <- TraRe::rewiring_test_pair_detail(modmat, pheno + 1,perm = retest_perms)
    new_pval <- result$pval}
  return(new_pval)

}

## function rewiring_regulon:
# Takes as input the output of get_regulons, preparedrewiring object,
# and the significant threshold of the rewiring method.
# Retreives an object that includes the regulon information and the pvalues obtained in the rewiring method


### THIS CAN BE PARALLELIZED
rewiring_regulon <- function(regulons_all_modules,
                             prepareobj,
                             final_signif_thresh){

  #Unlist the regulons
  graphs_regulons <- lapply(regulons_all_modules, function(x) x$regulons)
  unlisted_regulons <- unlist(graphs_regulons, recursive = FALSE, use.names = FALSE)

  # Rewiring of regulon (slow)
  message("Calculating the rewiring of the regulons")
  #Compute pvals for regulons on every rewired module
  pvals <- lapply(unlisted_regulons, function(r){
    rewiring_single_module(prepareobj,r)
  })

  # Calculate padj
  message('Calculating padjust values')
  padj <- stats::p.adjust(unlist(pvals),"BH")

  # Track of module number and bootstrap
  # Compute n? regulons/module:
  nr_regulons_per_module<- sapply(regulons_all_modules,function(x)length(x$regulons))

  # Create a vector with the number of regulons/modules
  module_nr_vector <- rep(names(nr_regulons_per_module),nr_regulons_per_module)

  # Extract the bootstrap of the regulon
  boost_rewired_mods <- as.numeric(sapply(regulons_all_modules,function(x) x$bootsts_idx))

  # Create a vector with the number of bootstrap of the regulon
  boots_idx_vector <- rep(boost_rewired_mods,nr_regulons_per_module)


  # Compile everything in  list
  finals_adj <- lapply(seq_along(unlisted_regulons),function(r){

    list(driver=igraph::V(unlisted_regulons[[r]])[igraph::get.vertex.attribute(unlisted_regulons[[r]],'type')]$name,

         targets=igraph::V(unlisted_regulons[[r]])[!igraph::get.vertex.attribute(unlisted_regulons[[r]],'type')]$name,

         p_value=unlist(pvals)[r],

         padj=padj[r],

         module_nr=module_nr_vector[r],

         boots_idx=boots_idx_vector[r],

         G=unlisted_regulons[[r]])
  })
  return(finals_adj)
}

## function fishermeth:
# Calculates the pvalue using Fisher's method
fishermeth <- function(pvals_v,showmess=FALSE,method='NA'){
  if (showmess){
    print(pvals_v)
  }
  if (method!='NA'){
    pvals_v <- stats::p.adjust(pvals_v,method=method)
  }
  #degrees of freedom
  df <- 2*length(pvals_v)
  #approximation by fisher
  newvar <- (-2) * sum(log(pvals_v))
  #Why upper tail??
  stats::pchisq(newvar,df,lower.tail = FALSE)
}


## function merge_regulons:
# Takes as input the output of rewirin_regulon and retrieves an object containing the info for
# the unified regulon (merged regulons withe the same driver).

merge_regulons <- function(finals_adj){
  message("Calculating unique regulons")
  #((repeat))
  targs <- sapply(finals_adj,function(x) x$targets)
  # should be equal to targslist: identical(targslist,targs)
  padjs <- sapply(finals_adj,function(x) x$padj)
  # should be equal to padj: identical(padj,padjs)

  #Calculate drivers within regulons
  uni_drivers <- unique(sapply(finals_adj,function(x) x$driver))

  #Evaluate each driver that appear in regulons within rewired drivers
  regulons_uniq <- lapply(uni_drivers, function(x){

    # Idx where there is a regulon for that driver
    bool_driv<- sapply(finals_adj, function(y) x%in%y$driver)

    iter_targs <- targs[bool_driv]
    iter_padjs <- padjs[bool_driv]

    # List of targets of this driver
    targs_vec <- unique(unlist(iter_targs))

    # For each target get the pvals if there is an edge
    info_edge <- t(sapply(targs_vec,function(j){
      bool <- sapply(iter_targs, function(k) j%in%k )
      m <- sum(bool) #multiplicity
      pval_targ <- iter_padjs[bool]
      pfish <- fishermeth(unlist(pval_targ))
      list(pfish_targ=pfish, multiplicity=m)
    }))

    p_value <- as.numeric(info_edge[,1])
    multiplicity <- as.numeric(info_edge[,2])

    list(driver=x,target=targs_vec, p_value=p_value,multiplicity=multiplicity)
  })

  #Give name to the drivers

  names(regulons_uniq) <- uni_drivers
  return(regulons_uniq)
}

## function filter_regulons: takes as input the output of merge_regulons and outputs an object
# with a data.frame for each regulon with filtered targets by pvalue and its corresponding info and the graph

filter_regulons <- function(regulons_uniq){
  message('Filtering results')
  # Check which regulons have no sig edges and get a bool
  bool_pval <- sapply(regulons_uniq, function(x) x$p_value<=0.05)
  bool_reg <- sapply(bool_pval,sum)!=0
  pval_th <- 0.05
  if (length(bool_pval)==length(regulons_uniq)){
    message("No regulon was found to be significantly rewired")
    pval_th=1
    bool_reg <- seq(regulons_uniq)
  }
  
  # Filter edges whose pval<0.05 and output ordered by multiplicity data frame
  dfs <- lapply(regulons_uniq[bool_reg], function(x){
    y <- as.data.frame(x)
    y <- y[y$p_value<=pval_th,]
    y <- y[order(y$multiplicity,decreasing = TRUE),]
    rownames( y ) <- seq_len( nrow( y ) )
    return(y)
  })

  graphs <- lapply(dfs,function(x){
    g <- igraph::graph_from_data_frame(x[,1:2], directed = FALSE)
    igraph::set_edge_attr(g, "weight",index = igraph::E(g), x$multiplicity)
  })
  return(list(regulons=dfs, graph=graphs))
}
