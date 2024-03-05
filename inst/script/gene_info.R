## code to prepare `gene_info` dataset goes here

# Download Human Protein Atlas (current version v.22)
temp <- tempfile()
download.file("https://v22.proteinatlas.org/download/proteinatlas.tsv.zip", temp)
human_protein_atlas_v22 <- read.delim(unz(temp, "proteinatlas.tsv"))
unlink(temp)

# Filter data columns and get Transcription factor information
tf_idx <- c()
i <- 0
for (protein in strsplit(human_protein_atlas_v22$Protein.class,split = ", ")){
  i <- i+1
  if ("Transcription factors"%in% protein){
    tf_idx <- c(tf_idx,i)
  }
}
regulator <- integer(length=nrow(human_protein_atlas_v22))
regulator[tf_idx] <- 1

gene_info <- cbind(human_protein_atlas_v22[,1:4],regulator)
gene_info <- gene_info[,c(1,3,5,2,4)]
# write.table(gene_info, file = paste0(system.file("data-raw", package="TraRe"), "/gene_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
usethis::use_data(gene_info, overwrite = TRUE)
