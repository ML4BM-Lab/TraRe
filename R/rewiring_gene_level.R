#' Perform the rewiring at the gene level
#'
#' Perform the rewiring test to check if regulators are enriched in rewiring modules using a hypergeometric test.
#'
#' @param filepath Desired path for the rewiring file to be generated.
#' @param linker_output_p Output file from linker function path. RDS format is required.
#' @param lognorm_est_counts_p Lognorm counts of the gene expression matrix path.
#' @param gene_info_p Path of a two-column file containing genes and 'regulator' boolean variable.
#' @param phenotype_p Path of a two-column file containing used samples and Responder or No Responder 'Class' (NR,R).
#' @param final_sig_th Significance threshold for the rewiring method. The lower the threshold, the restrictive the method. Default set to 0.05.
#' @param include_cliques Boolean specifying to include cliques in the returned score matrix. Default set to FALSE.
#' @param ImpTH Threshold for the refinement of the returned list. Default set to 0.05.
#' @param cliquesTH Correlation threshold if include_cliques is set to TRUE. Default set to 0.8.
#' @param nrcores Number of cores to run the parallelization within the rewiring test (default: 3).
#'
#'
#' @return Return a matrix containing, for each gene (or genes if include_cliques is set to TRUE) the pvalue and odds ratio
#' from the hypergeometric test within three categories.
#'
#' @export
rewiring_gene_level <- function(linker_output_p,
                               lognorm_est_counts_p,
                               gene_info_p,
                               phenotype_p,
                               filepath='',
                               final_sig_th=0.05,
                               include_cliques=FALSE,
                               ImpTH=0.05,
                               cliquesTH=0.8,
                               nrcores=3){

  #Do fast rewiring
  fast_rew <- function(ObjectList){

    #initialize common parameters
    regulator_info_col_name<-ObjectList$regulator_info_col_name
    phenotype_class_vals<-ObjectList$phenotype_class_vals
    phenotype_class_vals_label<-ObjectList$phenotype_class_vals_label
    outdir<-ObjectList$outdir
    orig_test_perms<-ObjectList$orig_test_perms
    retest_thresh<-ObjectList$retest_thresh
    retest_perms<-ObjectList$retest_perms
    logfile <- ObjectList$logfile

    #set.seed(1)
    dqrng::dqset.seed(1)

    #we create rundata and combine the modules of both parsers.
    #Lets generate the dupla (dataset's method - dataset's number)

    duplas <- unlist(lapply(seq_along(ObjectList$datasets),
                            function(i) paste(names(ObjectList$datasets[[i]]$rundata$modules),i)))

    #Initialize c_allstats and statsnames variables
    c_allstats <- c() #This is for the combined heatmap
    c_module_membership_list <- c() #This is for the combined heatmap

    statsnames <- c("module-method", "module-index", "orig-pval",
                    "revised-pvalue", "num-targets", "num-regulators",
                    "regulator-names","target-names", "num-samples", "num-genes",
                    "num-class1", "num-class2")

    for (dupla in duplas) {


      #Initialize rewired module's hash table
      module_membership_list <- hash::hash()

      #Initialize allstats array
      allstats <- NULL

      #For instance: 'VBSR X' to 'VBSR' and 'X'
      modmeth_i <- unlist(strsplit(dupla,' '))

      modmeth <- modmeth_i[1]
      i <- as.numeric(modmeth_i[2])

      #Output to the user which dupla we are working with
      message(modmeth,' ',i)

      rundata<-ObjectList$'datasets'[[i]]$rundata
      norm_expr_mat_keep<-ObjectList$'datasets'[[i]]$norm_expr_mat_keep
      keepsamps<-ObjectList$'datasets'[[i]]$keepsamps
      keeplabels<-ObjectList$'datasets'[[i]]$keeplabels
      class_counts<-ObjectList$'datasets'[[i]]$class_counts
      final_signif_thresh<-ObjectList$'datasets'[[i]]$final_signif_thresh
      responder<-ObjectList$'datasets'[[i]]$responder
      gene_info_df_keep<-ObjectList$'datasets'[[i]]$gene_info_df_keep
      name2idx<-ObjectList$'datasets'[[i]]$name2idx


      # This will register nr of cores/threads, keep this here
      # so the user can decide how many cores based on
      # their hardware.

      parallClass <- BiocParallel::bpparam()
      parallClass$workers <- ObjectList$NrCores

      GenerateStats <- function(mymod){

        signify <- NULL
        modregs <- unique(rundata$modules[[modmeth]][[mymod]]$regulators)
        modtargs <- unique(rundata$modules[[modmeth]][[mymod]]$target_genes)
        regnames <- paste(collapse = ", ", modregs)
        targnames <- paste(collapse = ", ", modtargs)
        keepfeats <- unique(c(modregs, modtargs))
        modmat <- t(norm_expr_mat_keep[keepfeats, keepsamps])
        modmeth_i_c <- paste(modmeth_i,collapse=' ')

        orig_pval <- TraRe::rewiring_test(modmat, keeplabels + 1, perm = orig_test_perms)
        new_pval <- orig_pval
        stats <- c(modmeth_i_c, mymod, signif(orig_pval, 3), signif(new_pval, 3),
                   length(modtargs), length(modregs), regnames,targnames, dim(modmat),
                   class_counts)
        if (orig_pval < retest_thresh | orig_pval == 1 | mymod %% 300 == 0) {
          #methods::show(paste(c("ModNum and NumGenes", mymod, length(keepfeats))))
          result <- TraRe::rewiring_test_pair_detail(modmat, keeplabels + 1,perm = retest_perms)
          new_pval <- result$pval
          stats <- c(modmeth_i_c, mymod, signif(orig_pval, 3),
                     signif(new_pval, 3), length(modtargs),
                     length(modregs), regnames,targnames, dim(modmat), class_counts)

          if (new_pval <= final_signif_thresh | new_pval == 1) {
            # save as list
            modname <- paste0(modmeth,'.',i, ".mod.", mymod)
            #module_membership_list[[modname]] <- keepfeats
            signify <- list(modname,keepfeats)
          }

        }
        return(list(stats,signify))

      }
      foreach_allstats <-BiocParallel::bplapply(seq_along(rundata$modules[[modmeth]]),
                                                GenerateStats, BPPARAM = parallClass)

      for (elements in foreach_allstats){

        #now we recover first allstats matrix
        foreach_stats <- elements[[1]]
        allstats<-rbind(allstats,foreach_stats)

        #and then update the module_membership dictionary
        hashtable <- elements[[2]]

        if (!is.null(hashtable)){

          module_membership_list[[hashtable[[1]]]] <- hashtable[[2]]

        }
      }

    }

    return(module_membership_list)

  }

  #generate dataframe with driver scores
  score_driv <- function(driver,linker_output,rew_list_names,showmess=FALSE){

    bootdist <- sapply(linker_output,function(x) x$bootstrap_idx)
    nbootstraps <- max(bootdist)

    oddsratio <- function(contig_tbl){
      #+0.5 using the Haldane Anscombe correction
      #(C/D) / (A/B)
      oddsr <- ((contig_tbl[2,1]+0.5) / (contig_tbl[2,2]+0.5)) / ((contig_tbl[1,1]+0.5) / (contig_tbl[1,2]+0.5))
      if (showmess){
        message(oddsr)
      }
      return(oddsr)
    }

    #go through every bootstrap
    pvals <- sapply(seq(nbootstraps),function(x) {

      #message('Bootstrap: ', x)

      #All the modules within a bootstrap
      pos <- which(bootdist==x)

      #A (rewired modules)
      A_modnames <- intersect(pos,rew_list_names)

      #Ac (non rewired modules)
      Ac_modnames <- setdiff(pos,rew_list_names)

      rew_mods <- sapply(A_modnames,function(y) sum(driver%in%linker_output[[y]]$regulators))
      non_rew_mods <- sapply(Ac_modnames,function(y) sum(driver%in%linker_output[[y]]$regulators))


      A_0 <- sum(rew_mods==0)
      A_1 <- sum(rew_mods!=0)
      Ac_0 <- sum(non_rew_mods==0)
      Ac_1 <- sum(non_rew_mods!=0)

      contig_tbl <- as.table(matrix(c(A_0,Ac_0,A_1,Ac_1),
                                    ncol = 2, byrow = FALSE))
      if (showmess){
        message('Every bootstrap')
        print(contig_tbl)
      }

      #show(contig_tbl)
      #Null hypothesis tail!
      res <- stats::fisher.test(contig_tbl,alternative = 'l')

      #Pvalue!
      c(res$p.value, oddsratio(contig_tbl))

    })

    #evaluate all bootstraps at once
    pvals_all <- function(){

      #A (rewired modules)
      A_modnames <- rew_list_names

      #Ac (non rewired modules)
      Ac_modnames <- setdiff(seq_along(bootdist),rew_list_names)


      rew_mods <- sapply(A_modnames,function(y) sum(driver%in%linker_output[[y]]$regulators))
      non_rew_mods <- sapply(Ac_modnames,function(y) sum(driver%in%linker_output[[y]]$regulators))


      A_0 <- sum(rew_mods==0)
      A_1 <- sum(rew_mods!=0)
      Ac_0 <- sum(non_rew_mods==0)
      Ac_1 <- sum(non_rew_mods!=0)

      contig_tbl <- as.table(matrix(c(A_0,Ac_0,A_1,Ac_1),
                                    ncol = 2, byrow = FALSE))

      if (showmess){
        message('All bootstraps')
        print(contig_tbl)
      }

      #show(contig_tbl)
      #Null hypothesis tail!
      res <- stats::fisher.test(contig_tbl,alternative='l')

      #Pvalue!
      c(res$p.value, oddsratio(contig_tbl))


    }

    #evaluate every 5 bootstraps
    pvals_b_5 <- sapply(seq(nbootstraps/5),function(x) {

      #message('Bootstrap: ', x)

      #All the modules within a bootstrap
      pos <- which(bootdist%in%seq(x*5-4,x*5))

      #A (rewired modules)
      A_modnames <- intersect(pos,rew_list_names)

      #Ac (non rewired modules)
      Ac_modnames <- setdiff(pos,rew_list_names)


      rew_mods <- sapply(A_modnames,function(y) sum(driver%in%linker_output[[y]]$regulators))
      non_rew_mods <- sapply(Ac_modnames,function(y) sum(driver%in%linker_output[[y]]$regulators))


      A_0 <- sum(rew_mods==0)
      A_1 <- sum(rew_mods!=0)
      Ac_0 <- sum(non_rew_mods==0)
      Ac_1 <- sum(non_rew_mods!=0)

      contig_tbl <- as.table(matrix(c(A_0,Ac_0,A_1,Ac_1),
                                    ncol = 2, byrow = FALSE))

      if (showmess){
        message('Every 5 bootstraps')
        print(contig_tbl)
      }

      #Null hypothesis tail!
      res <- stats::fisher.test(contig_tbl,alternative='l')

      #Pvalue and odds ratio #(c/d) / (a/b)
      c(res$p.value, oddsratio(contig_tbl))

    })

    fishermeth <- function(pvals_v,showmess=FALSE){

      if (showmess){
        print(pvals_v)
      }

      #degrees of freedom
      df <- 2*length(pvals_v)

      #approximation by fisher
      newvar <- (-2) * sum(log(pvals_v))

      #Why upper tail??
      stats::pchisq(newvar,df,lower.tail = FALSE)

    }

    #Reformat
    pvals <- c(fishermeth(pvals[1,]),
               mean(pvals[2,]))
    pvals_b_5 <- c(fishermeth(pvals_b_5[1,]),
                   mean(pvals_b_5[2,]))

    results_l <- list(single_b = pvals,
                      all_b = pvals_all(),
                      every5_b = pvals_b_5)
    #return as a dataframe
    return(t(data.frame(do.call(c,results_l))))

  }


  #Include cliques (Without duplicities)
  add_cliques <- function(regs_mrm,score_matrix,cliquesTH){

    #Generate cliques
    gcliques <- generatecliques(utils::read.delim(lognorm_est_counts_p)[regs_mrm,],
                                correlationth=cliquesTH)

    #retrieve gene-level drivers
    rew_gl_drivers <- rownames(score_matrix)

    cliqued_rew_gl <- sapply(rew_gl_drivers,function(drive){paste0(gcliques[[drive]],collapse='||')},
                             USE.NAMES = FALSE)

    #change the rownames
    rownames(score_matrix) <- cliqued_rew_gl

    return(score_matrix)

    }

  message('Preparing rewiring object')

  #Prepare rewiring creating object
  preparedrewiring <- TraRe::preparerewiring(linker_output_p = linker_output_p,
                                             lognorm_est_counts_p = lognorm_est_counts_p,
                                             gene_info_p = gene_info_p,
                                             phenotype_p = phenotype_p,
                                           final_signif_thresh = final_sig_th,
                                           nrcores=nrcores)
  if (!file.exists(fpath)){

    #define fpath (check default value)
    if (filepath==''){

      fpath <- paste0(getwd(),'/rewiring_gene_level','_fs_',final_sig_th,'.txt')

    }else{

      message('Adding significant threshold to output file name')

      #(check extension to add final sig th)
      if (substr(filepath,nchar(filepath)-3,nchar(filepath)-3)=='.'){
        fpath <- paste0(substr(filepath,1,nchar(filepath)-4),'_fs_',final_sig_th,'.txt')
      }else{
        fpath <-paste0(filepath,'_fs_',final_sig_th,'.txt')
      }

    }

    message('Performing fast rewiring')
    #Do the fast rewiring
    rew_list <- fast_rew(preparedrewiring)
    rew_list_names <- sub('[a-zA-Z]+.[0-9].mod.','',ls(rew_list))

    #write rew_list_names
    utils::write.table(rew_list_names,file=fpath,quote=F,sep='\t',
                row.names = FALSE,col.names=FALSE)
  }else{
    message('File found, reading from there')
  }

  #read file
  rew_list_names <- utils::read.delim(fpath,header=FALSE)[,1]

  #Intersect genes to obtain TFs
  gene_info_tfs <- utils::read.delim(gene_info_p)[,c('uniq_isos','regulator')]
  drivers_from_geneinfo <- gene_info_tfs[gene_info_tfs$regulator==1,'uniq_isos']
  regs_mrm <- intersect(rownames(utils::read.delim(lognorm_est_counts_p)),drivers_from_geneinfo)

  #confirm regulator genes
  #rew_driv <- unique(unlist(sapply(rew_list_names,function(x) preparedrewiring$datasets[[1]]$rundata$modules$VBSR[[x]]$regulators)))

  #define the linker output
  linker_output <- preparedrewiring$datasets[[1]]$rundata$modules$VBSR

  #Initialize the parallel process
  parallClass <- BiocParallel::bpparam()
  parallClass$workers <- preparedrewiring$NrCores

  message('Building score matrix')
  #Build score matrix
  score_matrix <- BiocParallel::bplapply(regs_mrm,score_driv,linker_output=linker_output,
                                        rew_list_names=rew_list_names,BPPARAM = parallClass)
  #Reformat
  score_matrix <- as.data.frame(do.call(rbind,score_matrix),row.names = regs_mrm)

  #Add colnames
  colnames(score_matrix) <- c('1B - FS','1B - OR',
                    'AllB - FS','AllB - OR',
                    '5B - FS','5B - OR')

  #List refinement
  imp <- as.numeric(which(apply(score_matrix,1,function(x) any(x[c(1,3,5)]<ImpTH))))

  #filter by refinement
  score_matrix_refined <- score_matrix[imp[order(score_matrix[imp,3])],]

  if (include_cliques){

    old_genes <- length(rownames(score_matrix_refined))
    score_matrix_refined <- add_cliques(regs_mrm = regs_mrm,
                                        score_matrix = score_matrix_refined,
                                        cliquesTH = cliquesTH)

    new_genes <- length(unique(unlist(sapply(rownames(score_matrix_refined),
                                             function(i) unlist(strsplit(i,split='\\|\\|'))))))

    message('Drivers db augmented from ',old_genes, ' to ',new_genes)

  }

  return(score_matrix_refined)

}

