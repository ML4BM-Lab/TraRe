#' Perform the rewiring at the regulon level
#'
#' Perform the rewiring test to regulons using a hypergeometric test.
#'
#' @param linker_output_p Path of output file from linker function. RDS format is required.
#' @param lognorm_est_counts_p Path of lognorm counts of the gene expression matrix.
#' @param gene_info_p Path of a two-column file containing genes and 'regulator' boolean variable.
#' @param phenotype_p Path of a two-column file containing used samples and Responder or No Responder 'Class' (NR,R).
#' @param final_signif_thresh Significance threshold for the rewiring method. The lower the threshold, the restrictive the method. Default set to 0.05.
#' @param rewired Bool indicating whether the analysis is to be done regulons from rewired modules (default= TRUE)
#' @param sigmodules_p Path of the output file from rewired_gene_level function with rewired modules list.
#'
#'
#' @return Return an object with a data.frame for each regulon with filtered targets by pvalue, multiplicity and the igraph object.
#'
#' @export

rewiring_regulon_level <- function(linker_output_p,
                                   lognorm_est_counts_p,
                                   gene_info_p,
                                   phenotype_p,
                                   sigmodules_p,
                                   rewired = TRUE,
                                   final_signif_thresh = 0.05){

  ### Let's work
  prepareobj <- TraRe::preparerewiring(linker_output_p = linker_output_p,
                                       lognorm_est_counts_p = lognorm_est_counts_p,
                                       gene_info_p = gene_info_p,
                                       phenotype_p = phenotype_p,
                                       final_signif_thresh = final_signif_thresh,
                                       use_graphs = T)

  ## Get regulons
  regulons_all_modules <- get_regulons(prepareobj,
                                       sigmodules_p = sigmodules_p,
                                       rewired = TRUE)

  ## Rewire regulons
  finals_adj <- rewiring_regulon(regulons_all_modules, prepareobj)

  ## Merge regulons
  regulons_uniq <- merge_regulons(finals_adj)

  ## Filter regulons
  order_filtr_regulons <- filter_regulons(regulons_uniq)

  return(order_filtr_regulons)
}

## function get_regulons:
# Takes as input prepareobj object, the number of rewired modules' path (if wanted to get the regulons from only those modules).
# Retrieves an object with the regulon subgraphs of the modules of linkeroutput and the boostrap idx
get_regulons <- function(prepareobj,
                         sigmodules_p = NULL,
                         rewired = rewired){

  if (rewired){
    message('Extracting regulons from rewired modules')
    #Read the rewired modules from the 50 bootstraps
    sigmod <- utils::read.delim(sigmodules_p,header = F)[,1]
    # Rewired graphs
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
  norm_expr_mat_keep<-ObjectList$'datasets'[[i]]$norm_expr_mat_keep
  keepsamps<-ObjectList$'datasets'[[i]]$keepsamps
  keeplabels<-ObjectList$'datasets'[[i]]$keeplabels

  #Now rewiring test!

  modregs <- intersect(ObjectList$datasets[[i]]$allregs,genes)
  modtargs <- intersect(ObjectList$datasets[[i]]$alltargs,genes)
  keepfeats <- unique(c(modregs, modtargs))
  modmat <- t(norm_expr_mat_keep[keepfeats, keepsamps])
  orig_pval <- TraRe::rewiring_test(modmat, keeplabels + 1, perm = orig_test_perms)
  new_pval <- orig_pval

  if (orig_pval < retest_thresh | orig_pval == 1 | mymod %% 300 == 0) {
    result <- TraRe::rewiring_test_pair_detail(modmat, keeplabels + 1,perm = retest_perms)
    new_pval <- result$pval}
  return(new_pval)

}

## function rewiring_regulon:
# Takes as input the output of get_regulons, preparedrewiring object,
# and the significant threshold of the rewiring method.
# Retreives an object that includes the regulon information and the pvalues obtained in the rewiring method


### THIS CAN BE PARALLELIZED
rewiring_regulon <- function(regulons_all_modules,
                             prepareobj=prepareobj,
                             final_signif_thresh=final_signif_thresh){

  #Unlist the regulons
  graphs_regulons <- lapply(regulons_all_modules, function(x) x$regulons)
  unlisted_regulons <- unlist(graphs_regulons,recursive = F,use.names = F)

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

  # Filter edges whose pval<0.05 and output ordered by multiplicity data frame
  dfs <- lapply(regulons_uniq[bool_reg], function(x){
    y <- as.data.frame(x)
    y <- y[y$p_value<=0.05,]
    y[order(y$multiplicity,decreasing = TRUE),]
    rownames( y ) <- seq_len( nrow( y ) )
    return(y)
  })

  graphs <- lapply(dfs,function(x){
    g <- igraph::graph_from_data_frame(x[,1:2],directed = F)
    igraph::set_edge_attr(g, "weight",index = igraph::E(g), x$multiplicity)
  })
  return(list(regulons=dfs, graph=graphs))
}
