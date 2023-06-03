## code to prepare `PRAD_SU2C` dataset goes here

# Download count matrix
fpkm_polya_counts <- readr::read_tsv("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/prad_su2c_2019/data_mrna_seq_fpkm_polya.txt",)
# # Download clinical data (for phenotype)
# clinical_data <- readr::read_tsv("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/prad_su2c_2019/data_clinical_patient.txt",skip = 4,trim_ws = TRUE)
# sample_data <- readr::read_tsv("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/prad_su2c_2019/data_clinical_sample.txt",skip = 4, trim_ws = TRUE)
# sample_data$PATIENT_ID <- as.character(sample_data$PATIENT_ID)
# # Merge clinical data
# full_clinical_data <- dplyr::left_join(sample_data,clinical_data) %>% 
#   dplyr::filter(SAMPLE_ID%in%colnames(fpkm_polya_counts[,2:ncol(fpkm_polya_counts)])) %>% 
#   dplyr::filter(ABI_ENZA_EXPOSURE_STATUS == "Naive")


# Download clinical data (for phenotype) 2
temp <- tempfile()
download.file("https://www.pnas.org/doi/suppl/10.1073/pnas.1902651116/suppl_file/pnas.1902651116.sd01.xlsx",temp,mode="wb")
clinical1 <- readxl::read_xlsx(temp,sheet = 1)
clinical2 <- readxl::read_xlsx(temp,sheet = 2)
unlink(temp)
# Merge clinical data
clinical_data <- dplyr::full_join(clinical1,clinical2) %>% 
                 dplyr::filter(`Sample ID`%in%colnames(fpkm_polya_counts[,-1])) %>% # get only samples for RNAseq
                 dplyr::filter(`Abiraterone (ABI) and Enzalutamide (ENZA) Exposure Status` == "Naive") %>% # get only ARSI naive samples
                 dplyr::filter(!is.na(`OS from first-line ARSI start`))  # get only those samples that have response phenotype
                      

# Categorize phenotype
OS_quartiles <- quantile(clinical_data$`OS from first-line ARSI start`)
subset_clinical_data <- clinical_data %>% 
                      dplyr::mutate(phenotype = `OS from first-line ARSI start` < OS_quartiles[4] & `OS from first-line ARSI start` > OS_quartiles[2] ) %>% 
                      dplyr::mutate(phenotype = na_if(phenotype,T)) %>% # NA for quartiles 2 and 3
                      dplyr::mutate(phenotype = phenotype +`OS from first-line ARSI start` > OS_quartiles[4]) # Responder = TRUE (quartile 4); Non respodner FALSE (quartile 1)
# Save clinical data
write.table(subset_clinical_data, file = "data-raw/phenotype.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset_clinical_data[,c(2,16)], "inst/extdata/new_phenotype_rewiring_example.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Load gene_info
gene_info <- read.delim("data-raw/gene_info.tsv", sep = "\t")
gene_order <- match(as.data.frame(fpkm_polya_counts)[,1],gene_info_new$Gene)

# Subset matrix
ass <- as.matrix(fpkm_polya_counts[,full_clinical_data$`Sample ID`])
rownames(ass) <- as.data.frame(fpkm_polya_counts)[,1]
max(ass)
min(ass)
# Exploratory plots
# plot(rowMeans(ass),rowSds(ass), xlab="mean",ylab="sds",main="meanSdplot")
# boxplot(ass)

vst_out <- DESeq2::varianceStabilizingTransformation(round(ass,0))
max(vst_out)
min(vst_out)
# Exploratory plots
# plot(rowMeans(vst_out),rowSds(vst_out), xlab="mean",ylab="sds",main="meanSdplot")
# boxplot(vst_out)


# Subset genes
# set.seed(123)
# gene_idx <- sample(nrow(gene_info), 2500)
# write.table(fpkm_polya_counts,"inst/extdata/new_expression_example.txt",sep = "\t", quote = FALSE, row.names = FALSE)


PRAD_SU2C <- SummarizedExperiment::SummarizedExperiment(assays = vst_out),
                                                        rowData = gene_info[gene_idx,],
                                                        colData = subset_clinical_data)

gene_order <- match(as.data.frame(fpkm_polya_counts)[,1],gene_info_new$Gene)
PRAD_SU2C <- SummarizedExperiment::SummarizedExperiment(assays = as.matrix(fpkm_polya_counts[,full_clinical_data$`Sample ID`]),
                                                        rowData = gene_info_new[gene_order,],
                                                        colData = full_clinical_data)

ass <- assays(PRAD_SU2C)[[1]]
plot(rowMeans(ass),rowSds(ass), xlab="mean",ylab="sds",main="meanSdplot")
max(ass)
min(ass)
boxplot(ass)
vst_out <- DESeq2::varianceStabilizingTransformation(round(ass,0))
plot(rowMeans(vst_out),rowSds(vst_out), xlab="mean",ylab="sds",main="meanSdplot")
max(vst_out)
min(vst_out)
boxplot(vst_out)


read_vars <- rowVars(as.matrix(fpkm_polya_counts[,-1]))
raw_counts <- raw_counts[which(read_vars > 1e-10),]

# Save dataset
saveRDS(PRAD_SU2C, "inst/extdata/PRAD_SU2C.rds")
usethis::use_data(PRAD_SU2C, overwrite = TRUE)
