#' Transcriptional Rewiring package for R.
#'
#' @section Description of TraRe:
#'
#' *TraRe* (Transcriptional Rewiring) is an R package which contains the necessary tools
#' to carry out: Identification of module-based gene regulatory networks (GRN); score-based
#' classification of these modules via a rewiring test; visualization of rewired modules
#' to analyze condition-based GRN deregulation and drop out genes recovering via cliques
#' methodology. For each tool, html report can be generated containing useful information
#' about the generated GRN and statistical data about the performed tests. This tools have
#' been developed considering sequenced data (RNA-Seq).
#'
#' *TraRe* offers the following tools
#'
#' \itemize{
#'
#'   \item *LINKER* generates, from an initial RNA-Seq dataset where drivers
#'   (Transcription Factors) and targets genes are provided, GRN modules in three different
#'   forms: as raw results, from the phase I output modules; as modules from the phase II
#'   output modules and in the form of bipartite graphs, where drivers and targets relationships
#'   are defined. This is done using the `LINKER_run()` function. `?LINKER_run` for more information.
#'
#'   \item *LINKER* also provides a way of generating a single GRN from specified
#'   list of driver and target genes. This ease the task of analyzing relationships between
#'   drivers and targets by constraining all the provided genes to be only in a single GRN.
#'   This is done by `NET_run()` function. Type `?NET_run` for more information
#'
#'   \item The *Rewiring* test performs a permutation test over a certain condition to infer
#'   if that condition is producing any deregulation on our generated GRN. Bootstrapping plays
#'    an important role, as the non-convex nature of this biological events makes necessary to
#'    ensure that a certain behavior is repeated across bootstraps, and to confirm this event
#'    does not come from a particular realization. As bootstrapping has been performed in LINKER,
#'    this step will take advantage of them and will try to group highly scored modules, to infer
#'    modules similar behavior GRN across bootstraps. It will outputs a correlation matrix in the
#'    form of a heatmap (sorted by hierarchical clustering to ease interpretation), containing
#'    similar highly scored modules. `preparerewiring()` will return an object containing the
#'    necessary information for calling `runrewiring()` and generate graph objects, reports and graphs.
#'    Type `?preparerewiring` and `?runrewiring` for more information.
#'
#'  \item The *Visualization* module contains a graphical way of detecting condition-dependent
#'  deregulation on the selected rewired modules containing gene regulatory networks. Once we
#'  have selected a cluster of modules that across bootstraps have similar behavior as GRN,
#'  we can generate single GRNs of the genes that belong to those modules, filtering by samples
#'  that belong to the condition we want to evaluate. Check `plot_igraph`,`return_layout` and
#'  `return_layout_phenotype`.
#'
#'  \item The *Results* module generate an excel file is containing drivers-targets relationships
#'  and cliques. The way LINKER method works can make some highly-correlated driver genes (TFs)
#'  may be dropped from the resultant model, as the role they play at the GRN inference process is
#'  very similar. Due to this, we propose a method based on cliques (Fully Connected Networks) to
#'  recover those dropped drivers. Check `excel_generation` and `generatecliques` for more information.
#'
#'
#' }
#' @docType package
#' @keywords internal
#' @name TraRe
"_PACKAGE"



