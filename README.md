# TraRe

*TraRe* (Transcriptional Rewiring) is an R package which contains the necessary tools to carry out: 
Identification of module-based gene regulatory networks (GRN); score-based classification of these modules 
via a Rewiring test; visualization of rewired modules to analyze condition-based GRN deregulation and drop 
out genes recovering via cliques methodology. For each tool, html report can be generated containing useful 
information about the generated GRN and statistical data about the performed tests. These tools have been 
developed considering for RNA-Seq data.

![TraRe Package](https://github.com/ubioinformat/TraRe/blob/master/vignettes/Trare.png)

For more information see:
*Bayesian Machine Learning Enables Identification of Transcriptional Network Disruptions Associated with Drug-Resistant Prostate Cancer* 
Charles Blatti, Jesús de la Fuente, Huanyao Gao, Irene Marín-Goñi, Zikun Chen, Sihai D. Zhao, Winston Tan, Richard Weinshilboum, Krishna R. Kalari, Liewei Wang, Mikel Hernaez; 
*Cancer Research* 15 April 2023; 83 (8): 1361–1380. https://doi.org/10.1158/0008-5472.CAN-22-1910

# Installation 
```{r, eval=FALSE}

 if (!requireNamespace("devtools"))
 install.packages("devtools")
devtools::install_github('ubioinformat/TraRe')

```