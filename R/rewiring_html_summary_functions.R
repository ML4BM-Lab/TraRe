#' Generate html report from the `runrewiring()` method.
#'
#' @description
#'
#' Contain functions to generate an html report from the rewired linker modules.
#' `create_index_page()` initializes folders and files for the html report. `write_tables_all()` generates
#' a table given stats from the rewiring method. `write_html_table_page()` creates the header for the
#' html file when the table will be placed. `table2html()` gets the `write_html_table_page()` info and
#' return the generated html file. This functions run inside `runrewiring()`
#'
#' @param outdir path for the resultant folder.
#' @param runtag name that is `paste()` to the name folder.
#' @param codedir folder from which `sorttable.js` and `glossary.txt`
#' files are `paste()` into the output folder.
#' @param indexpath name for the index file.
#' @param glossarypath name for the glossary file.
#' @param imgstr path for the image folder.
#' @param txtstr path for the txts folder.
#'
#' @param mytab table to transform into html file.
#' @param tabletype type of table.
#' @param html_idxs index for table rows.
#' @param html_cols index for table columns.
#' @param filestr name of the html file to be generated.
#' @param htmlinfo list containing `indexpath`,`txtstr`, and html dir.
#'
#'
#' @param resultstable modified table from `write_tables_all()`.
#' @param htmlpagefile path of the html file.
#' @param resultsrelpath path of the txt files.
#'
#'
#'
#' @noRd

table2html <- function(resultstable){

  colheader <- paste(collapse = "</th><th>", c("idx", colnames(resultstable)))
  theader <- paste0("<thead>\n<tr bgcolor='#AAAAAA';><th>", colheader,
                    "</th></tr>\n</thead>\n<tbody>\n")
  head_str = paste0("<script src='sorttable.js'></script>\n<TABLE class=",
                    "'sortable' border =1 >\n", theader)

  numtable = cbind(seq_len(nrow(resultstable)), resultstable)
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
  write(htmlstr, file = htmlpagefile, append = TRUE)
}


write_tables_all <- function(mytab, tabletype="table",
                             html_idxs=seq_len(nrow(mytab)),
                             html_cols=colnames(mytab),
                             filestr="data",
                             htmlinfo=list(htmldir="html/",
                                           indexpath="index.html",
                                           txtstr="txts/")){
  htmlpath = paste0(filestr, "_", tabletype, ".html")
  resultspath = paste0(htmlinfo$txtstr, filestr, "_", tabletype, ".txt")
  message("Writing table: ", resultspath)
  utils::write.table(mytab, paste0(htmlinfo$htmldir, resultspath), sep='\t',
              row.names=FALSE, col.names=TRUE, quote=FALSE)
  write(paste0('<a href = "',htmlpath,'" target="_blank">',
               tabletype,'</a><br>'),
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append=TRUE)
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

  return(list(htmldir = htmldir, indexpath = indexpath, imgstr = imgstr,
              txtstr = txtstr, glossarypath = glossarypath, abspath = abspath))
}
