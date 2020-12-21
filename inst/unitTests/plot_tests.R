plot_tests<- function(){

  graph <- readRDS(paste0(system.file("extdata",package="TraRe"),'/graph_netrun_example.rds'))

  #plot_function

  obs <-tryCatch(plot_igraph(mygraph=NULL,mylayout=NULL),error=conditionMessage)

  RUnit::checkIdentical("graph object field empty", obs)

  obs <-tryCatch(plot_igraph(mygraph=graph,mylayout=NULL),error=conditionMessage)

  RUnit::checkIdentical("layout field empty", obs)

  #return_layout

  obs <-tryCatch(return_layout(regs=NULL,
                               targets=c("gene4","gene5","gene6")),error=conditionMessage)

  RUnit::checkIdentical("regulators field empty", obs)

  obs <-tryCatch(return_layout(regs=c("gene1","gene2","gene3"),
                               targets=NULL),error=conditionMessage)

  RUnit::checkIdentical("targets field empty", obs)

  #return_layout_phenotype

  obs <-tryCatch(return_layout_phenotype(regs=NULL,
                               targets=c("gene4","gene5","gene6"),
                               varfile=NULL),error=conditionMessage)

  RUnit::checkIdentical("regulators field empty", obs)

  obs <-tryCatch(return_layout_phenotype(regs=c("gene1","gene2","gene3"),
                               targets=NULL,
                               varfile=NULL),error=conditionMessage)

  RUnit::checkIdentical("targets field empty", obs)

  obs <-tryCatch(return_layout_phenotype(regs=c("gene1","gene2","gene3"),
                               targets=c("gene4","gene5","gene6"),
                               varfile=NULL),error=conditionMessage)

  RUnit::checkIdentical("varfile field empty", obs)


  #return_layout_phenotype varfile data

  varfile <- matrix(seq_len(12),6,2)

  obs <-tryCatch(return_layout_phenotype(regs=c("gene1","gene2","gene3"),
                                         targets=c("gene4","gene5","gene6"),
                                         varfile=varfile),error=conditionMessage)

  RUnit::checkIdentical("genes names must be specified at varfile as rownames", obs)

  rownames(varfile)<-c("gene1","gene2","gene3","gene4","gene5","gene6")

  obs <-tryCatch(return_layout_phenotype(regs=c("gene1","gene2","gene3"),
                                         targets=c("gene4","gene5","gene6"),
                                         varfile=varfile),error=conditionMessage)

  RUnit::checkIdentical("colnames must be specified, in particular 'is-regulator' and 't-stat'", obs)

  colnames(varfile)<-c("category1","category2")

  obs <-tryCatch(return_layout_phenotype(regs=c("gene1","gene2","gene3"),
                                         targets=c("gene4","gene5","gene6"),
                                         varfile=varfile),error=conditionMessage)

  RUnit::checkIdentical("varfile must contain the column is-regulator", obs)

  colnames(varfile)<-c("is-regulator","category2")

  obs <-tryCatch(return_layout_phenotype(regs=c("gene1","gene2","gene3"),
                                         targets=c("gene4","gene5","gene6"),
                                         varfile=varfile),error=conditionMessage)

  RUnit::checkIdentical("varfile must contain the column t-stat", obs)


}
