runlinker_testmatrix<- function(){

  obs <-tryCatch(LINKER_run(lognorm_est_counts=NULL,
                            target_filtered_idx=seq_len(15),
                            regulator_filtered_idx=15+seq_len(3)),error=conditionMessage)

  RUnit::checkIdentical("lognorm_est_counts field empty", obs)


  obs <-tryCatch(LINKER_run(lognorm_est_counts="text",
                            target_filtered_idx=seq_len(15),
                            regulator_filtered_idx=15+seq_len(3)),error=conditionMessage)

  RUnit::checkIdentical("matrix class is required for input dataset", obs)

  obs <-tryCatch(LINKER_run(lognorm_est_counts=matrix(vapply(seq_len(18),toString,FUN.VALUE = c("1")),6,3),
                            target_filtered_idx=seq_len(15),
                            regulator_filtered_idx=15+seq_len(3)),error=conditionMessage)

  RUnit::checkIdentical("non-numeric values inside lognorm_est_counts variable", obs)

  obs <-tryCatch(LINKER_run(lognorm_est_counts=matrix(seq_len(18),6,3),
                            target_filtered_idx=seq_len(15),
                            regulator_filtered_idx=15+seq_len(3)),error=conditionMessage)

  RUnit::checkIdentical("null field detected in row names or column names, check lognorm_est_counts matrix", obs)

}
