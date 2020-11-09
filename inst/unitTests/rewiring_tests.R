rewiring_tests<- function(){


  obs <-tryCatch(preparerewiring(expr_matrix_file="path2",
                                 gene_info_file="path3",
                                 phenotype_file="path4"),error=conditionMessage)

  RUnit::checkIdentical("linker_saved_file field required", obs)

  obs <-tryCatch(preparerewiring(linker_saved_file="path1",
                                 gene_info_file="path3",
                                 phenotype_file="path4"),error=conditionMessage)

  RUnit::checkIdentical("expr_matrix_file field required", obs)

  obs <-tryCatch(preparerewiring(linker_saved_file="path1",
                                 expr_matrix_file="path2",
                                 phenotype_file="path4"),error=conditionMessage)

  RUnit::checkIdentical("gene_info_file field required", obs)

  obs <-tryCatch(preparerewiring(linker_saved_file="path1",
                                 expr_matrix_file="path2",
                                 gene_info_file="path3"),error=conditionMessage)

  RUnit::checkIdentical("phenotype_file field required", obs)

}
