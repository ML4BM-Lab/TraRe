% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rewiring_regulon_level.R
\name{rewiring_regulon_level}
\alias{rewiring_regulon_level}
\title{Perform the rewiring at the regulon level}
\usage{
rewiring_regulon_level(
  linker_output,
  TraReObj,
  fpath = "",
  rewired = TRUE,
  final_signif_thresh = 0.05,
  outdir = tempdir()
)
}
\arguments{
\item{linker_output}{Output from LINKER_run function (recommended with 50 bootstraps).}

\item{TraReObj}{the TrareObj generated during preprocessing step before GRN inference.}

\item{fpath}{Path of the output file from rewired_gene_level function with rewired modules list.}

\item{rewired}{Bool indicating whether the analysis is to be done regulons from rewired modules (default= TRUE)}

\item{final_signif_thresh}{Significance threshold for the rewiring method. The lower the threshold, the more restrictive the method. Default set to 0.05.}

\item{outdir}{Directory for the output folder to be located (default: tempdir()).}
}
\value{
Return an object list with a data.frame for each regulon with filtered targets by pvalue, multiplicity and the igraph object.
}
\description{
Perform the rewiring test to regulons using a hypergeometric test.
}
\examples{
## We will be using an output file we generated with LINKER_run and the 
## same expression datainput  object.

## linker_output <- readRDS(paste0(system.file('extdata',package='TraRe'),
##                                     'linkeroutput_rewiring_example.rds'))

## TraReObj <- readRDS(paste0(system.file('extdata',package='TraRe'),
##                                                        '/TraReObj.rds'))

## Add the phenotype to TraReObj if it was not before.
## TraReObj <- rewiring_add_phenotype(TraReObj, phenotype)

## Select directory for output file
## outdir <- system.file('extdata',package='TraRe')
## regulons <- rewiring_regulon_level(linker_output = linker_output, 
##                                    TraReObj = TraReObj,
##                                    rewired = TRUE,
##                                    final_signif_thresh = 0.05,
##                                    outdir = outdir)
}
