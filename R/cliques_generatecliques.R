#' Cliques generation and plot displaying.
#'
#' Cliques generation and plot displaying of the selected method.
#'
#' @param dataset input expression file with genes as rows and samples as columns.
#' @param method method to use in the correlation matrix (see stats:cor).
#' @param correlationth threshold to consider edge exists.
#' @param sparsecorrmatrix boolean variable specifying whether to fix values below threshold to 0 or not.
#' @param numcliques number of cliques to be generated.
#' @param mandatorygenes array of gene names which is mandatory to include in the returned cliques.
#' @param selection integer selecting method.
#'
#' @return List containing the plot generated and a list with all the generated cliques.
#'
#' @examples
#'    \dontrun{
#'    mat <- matrix(c(0,0,0,0,0,1,0,1,0),3,3)
#'    g <- graph_from_adjacency_matrix(mat)
#'    cli <- max_cliques(g)
#'    tittle <- "Maximizing genes/cliques"
#'    L <-  350
#'    numcliques <- 'All'
#'    mandatorygenes <- c()
#'
#'    generatecliques(cli,tittle,L,numcliques,mandatorygenes)
#'    }
#' @export
generatecliques<-function(dataset,method,correlationth,sparsecorrmatrix,numcliques,mandatorygenes,selection=1){ #All but the last one to print it in a different way.

  methods::show("Preparing data")
  pdobject <- preparedata(dataset,method)

  methods::show("Generating graph")
  ggobject <- generategraph(correlationth,sparsecorrmatrix,pdobject)

  methods::show("Selecting method")
  smobject <- selectmethod(selection,ggobject,pdobject)

  methods::show("Generate Datasets")

  ml<-smobject$cliques
  #check numcliques is not greater than the amount of modules we have.
  if (numcliques>length(ml)|is.character(numcliques)) numcliques<-length(ml)

  #ml-> cliques.
  sortparameter<-sapply(ml[],getMaxVariance,dat=pdobject) #We assure max variance genes are getting in.
  sortparameter_ix<-sort(sortparameter,decreasing = TRUE,index.return=TRUE)$ix
  sortparameter_ix_numcliques<-sortparameter_ix[1:numcliques]

  ml<-ml[sortparameter_ix_numcliques] #sort it and select the numcliques

  mandatorygenes_toinclude<-setdiff(mandatorygenes,unlist(ml)) #Get the genes that we need and aren't in our selected ones.

  #Substitute the bottom cliques with the mandatory ones.
  #We can substitute a clique with more than 1 gene for a single gene...

  for (i in length(ml):1){
    if (!length(mandatorygenes_toinclude)) break

    if (!ml[[i]]%in%mandatorygenes){ #check we are not substituting a mandatory gene.

      ml[[i]]<-mandatorygenes_toinclude[1]                   #character.pop(0)
      mandatorygenes_toinclude<-mandatorygenes_toinclude[-1] #character.pop(0)

    }
  }

  #We build a dictionary, where Representatives are the keys and Cliques are the values.

  RC_list<-list()
  RC_list<-ml
  names(RC_list)<- sapply(RC_list,getMaxVarName,dat=pdobject) #Representatives

  #Plot the selected method.

  x_axis<-sapply(ml[],getAvgVariance,dat=pdobject)
  y_axis<-sapply(ml[],getMedianCorrelation,dat=pdobject) #Median Correlation


  Palette<-c("#29bd00","#ff04d1","#1714b0","#b30009","#030303")

  col_fact <- c()

  mll<-sapply(ml[],length)

  for (i in 1:length(ml)){

    clique_l<-length(ml[[i]])

    ch_cat<-which(c((clique_l==1)&(sortparameter_ix_numcliques[i]<smobject$L),(clique_l>1)&(clique_l<=4),(clique_l>4)&(clique_l<=7),(clique_l>7),(sortparameter_ix_numcliques[i]>=smobject$L)))

    col_fact<-c(col_fact,switch(ch_cat,"x == 1","1 < x <= 4", "4 < x <= 7","x > 7","Singleton"))

  }
  col_fact<-as.factor(col_fact)

  SelectedCliquesPlot<-ggplot2::ggplot()+

    ggplot2::geom_point(ggplot2::aes(y_axis,x_axis,color=col_fact),size=3)+

    ggplot2::labs(
      x="Median Correlation Value",
      y="Avg Variance Per Clique",
      title = smobject$tittle,
      subtitle = paste("TotalNum: ",toString(sum(mll)),
                       ", MedianOfMedian: ",
                       toString(round(stats::median(y_axis[y_axis!=1]),4)),
                       ", Total Singleton Cliques: ",
                       toString(as.double(table(sapply(ml,length)==1)[2])),
                       ", Number of Cliques: ",
                       toString(numcliques),
                       sep="")
    )+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size=15,face = "bold"),plot.subtitle = ggplot2::element_text(hjust = 0.5,size=12))+
    ggplot2::theme(legend.text = ggplot2::element_text(size = 15),axis.text=ggplot2::element_text(size=15),axis.title = ggplot2::element_text(size=15),legend.title = ggplot2::element_text(size = 15))+
    ggplot2::scale_color_manual(name="Genes/Cliques",values = c("x == 1"=Palette[1], "1 < x <= 4"=Palette[2],"4 < x <= 7"=Palette[3],"x > 7"=Palette[4],"Singleton"=Palette[5]))

  return(list(plot=SelectedCliquesPlot,cliques=RC_list,representatives=names(RC_list)))
}
