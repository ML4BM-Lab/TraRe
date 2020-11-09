runlinkernet_testindex<- function(){


  lognorm_est_counts<-utils::read.table(file=paste0(system.file("extdata",package="TraRe"),
                                                    '/expression_rewiring_example.txt'),
                                        header=TRUE,sep="\t",stringsAsFactors = FALSE)

  obs <-tryCatch(NET_run(lognorm_est_counts=lognorm_est_counts,
                            target_filtered_idx=NULL,
                            regulator_filtered_idx=150+seq_len(1000)),error=conditionMessage)

  RUnit::checkIdentical("target_filtered_idx field empty", obs)


  obs <-tryCatch(NET_run(lognorm_est_counts=lognorm_est_counts,
                            target_filtered_idx=seq_len(150),
                            regulator_filtered_idx=NULL),error=conditionMessage)

  RUnit::checkIdentical("regulator_filtered_idx field empty", obs)

  obs <-tryCatch(NET_run(lognorm_est_counts=lognorm_est_counts,
                            target_filtered_idx=seq_len(150),
                            regulator_filtered_idx=150+seq_len(1000-1)),error=conditionMessage)

  RUnit::checkIdentical("the total number of genes is not equal to the sum of target_filtered_idx and regulatory_filtered_idx lengths", obs)

  obs <-tryCatch(NET_run(lognorm_est_counts=lognorm_est_counts,
                            target_filtered_idx=vapply(seq_len(150),toString,FUN.VALUE = c("1")),
                            regulator_filtered_idx=150+seq_len(1000)),error=conditionMessage)

  RUnit::checkIdentical("targets and regulators index arrays must be numeric", obs)

}
