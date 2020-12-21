# TraRe

*TraRe* (Transcriptional Rewiring) is an R package which contains the necessary tools to carry out: 
Identification of module-based gene regulatory networks (GRN); score-based classification of these modules 
via a rewiring test; visualization of rewired modules to analyze condition-based GRN deregulation and drop 
out genes recovering via cliques methodology. For each tool, html report can be generated containing useful 
information about the generated GRN and statistical data about the performed tests. These tools have been 
developed considering sequenced data (RNA-Seq).

![TraRe Package](https://github.com/ubioinformat/TraRe/blob/Version_1_3_0/vignettes/Trare.png)

# Installation 
```{r, eval=FALSE}

 if (!requireNamespace("BiocManager"))
 install.packages("BiocManager")
 BiocManager::install("TraRe")

```