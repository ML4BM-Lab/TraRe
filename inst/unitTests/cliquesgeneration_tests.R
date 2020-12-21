cliquesgeneration_tests<- function(){

  obs <-tryCatch(generatecliques(dataset=matrix(c("1","2","3","4"),2,2),
                          method="spearman",
                          correlationth=0.7,
                          sparsecorrmatrix=FALSE),error=conditionMessage)

  RUnit::checkIdentical("non-numeric values inside dataset variable", obs)


  obs <-tryCatch(generatecliques(dataset=NULL,
                                 method="spearman",
                                 correlationth=0.7,
                                 sparsecorrmatrix=FALSE),error=conditionMessage)

  RUnit::checkIdentical("dataset field empty", obs)

  obs <-tryCatch(generatecliques(dataset="text",
                                 method="spearman",
                                 correlationth=0.7,
                                 sparsecorrmatrix=FALSE),error=conditionMessage)

  RUnit::checkIdentical("matrix or dataframe class is required", obs)

}
