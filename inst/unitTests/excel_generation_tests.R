excel_generation_tests <- function(){

  obs <-tryCatch(excel_generation('Not working path'),error=conditionMessage)

  RUnit::checkIdentical('refinedsumm.rds in the specified folder must exist', obs)

  workingpath <- system.file('extdata',package='TraRe')

  obs <-tryCatch(excel_generation(workingpath,cliquesbool='non-logical'),error=conditionMessage)

  RUnit::checkIdentical('non-logical variable pass to this cliquesbool argument', obs)
}
