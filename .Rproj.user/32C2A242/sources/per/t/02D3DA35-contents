#' Generate html report.
#'
#' @description
#'
#' Contain functions to generate an html report from the output modules. `table2html()`
#' creates html summary of a set of linker runs on the same data. `write_html_table_page()` place all
#' pngs on same html page, for run directory. `write_tables_all()` create html summary of a set of
#' linker runs on the same data. `create_index_page()` generates the rest of pages of the html report.
#'
#'
#' @param resultstable example of param.
#' @param htmlpagefile example of param.
#' @param resultsrelpath example of param.
#' @param indexpath example of param.
#' @param glossarypath example of param.
#' @param mytab example of param.
#' @param tabletype example of param.
#' @param html_idxs example of param.
#' @param html_cols example of param.
#' @param filestr example of param.
#' @param htmlinfo example of param.
#' @param outdir example of param.
#' @param runtag example of param.
#' @param codedir example of param.
#' @param imgstr example of param.
#'
#'
#' @examples
#' \dontrun{
#' example of example.
#' }
#'
#' @export

table2html <- function(resultstable){

  colheader <- paste(collapse = "</th><th>", c("idx", colnames(resultstable)))
  theader <- paste0("<thead>\n<tr bgcolor='#AAAAAA';><th>", colheader,
                    "</th></tr>\n</thead>\n<tbody>\n")
  head_str = paste0("<script src='sorttable.js'></script>\n<TABLE class=",
                    "'sortable' border =1 >\n", theader)

  numtable = cbind(1:dim(resultstable)[1], resultstable)
  rowinnerstrs = apply(numtable, 1, paste, collapse = "</td><td>")
  rowstrs = paste0("<tr><td>", rowinnerstrs, "</td></tr>")
  rows_str = paste(collapse="\n", rowstrs)

  return(paste0(collapse="\n", head_str, rows_str, "\n</tbody></table><br>"))
}


write_html_table_page <- function(resultstable, htmlpagefile, resultsrelpath,
                                  indexpath="index.html",
                                  glossarypath="glossary.html"){

  write(paste0("<table border = 1 width = '100%'><tr bgcolor = ",
               "'#AAAAAA'><th><a href = '", indexpath, "' target='_blank'>",
               "Index</a></th><th><a href = '", glossarypath, "' target=",
               "'_blank'>Glossary</a></th><th><a href = '", resultsrelpath,
               "' target='_blank'>Download</a></th></tr></table><br><br>"),
        file = htmlpagefile)
  htmlstr = table2html(resultstable)
  write(htmlstr, file = htmlpagefile, append = T)
}


write_tables_all <- function(mytab, tabletype="table",
                             html_idxs=1:dim(mytab)[1],
                             html_cols=colnames(mytab),
                             filestr="data",
                             htmlinfo=list(htmldir="html/",
                                           indexpath="index.html",
                                           txtstr="txts/")){
  htmlpath = paste0(filestr, "_", tabletype, ".html")
  resultspath = paste0(htmlinfo$txtstr, filestr, "_", tabletype, ".txt")
  methods::show(paste0("Writing table: ", resultspath))
  utils::write.table(mytab, paste0(htmlinfo$htmldir, resultspath), sep='\t',
              row.names=F, col.names=T, quote=F)
  write(paste0('<a href = "',htmlpath,'" target="_blank">',
               tabletype,'</a><br>'),
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append=T)
  write_html_table_page(resultstable=mytab[html_idxs,html_cols],
                        htmlpagefile=paste0(htmlinfo$htmldir, htmlpath),
                        resultsrelpath=resultspath,
                        indexpath=htmlinfo$indexpath)
}

create_index_page <- function(outdir="./", runtag="run", codedir="code/",
                              indexpath="index.html",
                              glossarypath="glossary.html",
                              imgstr="imgs/", txtstr="txts/"){

  htmldir <- paste0(outdir, runtag, "/")
  dir.create(file.path(htmldir))

  file.copy(from = paste0(codedir, "sorttable.js"), to = htmldir)
  dir.create(file.path(paste0(htmldir, imgstr)))
  dir.create(file.path(paste0(htmldir, txtstr)))

  glossary <- as.matrix(utils::read.table(paste0(codedir, "glossary.txt"),
                                   header = TRUE, sep = "\t", quote=""))
  file.copy(from = paste0(codedir, "glossary.txt"), to = htmldir)
  abspath <- paste0(htmldir, indexpath)

  write(paste0("<br>"), file = abspath)
  #  write(paste0("<a href = '",indexpath,"' target='_blank'>Index</a><br><br>"), file = abspath)
  #  write_html_table(paste0(htmldir, glossarypath), glossary, glossarypath)

  return(list(htmldir = htmldir, indexpath = indexpath, imgstr = imgstr,
              txtstr = txtstr, glossarypath = glossarypath, abspath = abspath))
}
