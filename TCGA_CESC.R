setwd("E:/Qian/R Studio/TCGA_CESC_NEW")
getwd()

library(rjson)
library(tidyverse)

#read meta.data
json <- jsonlite::fromJSON("metadata.cart.2025-02-12.json")

#get sample names and file names
sample_id <- sapply(json$associated_entities, function(x){x[,1]})
sample_id[1:10]
file_sample <- data.frame(sample_id, file_name=json$file_name)

#get counts location
count_file <-list.files("gdc_download_20250212_234100.651941/",
                        pattern = "*.tsv", recursive = T)
#count_file[1:10]

#get every file name
count_file_name <- strsplit(count_file, split = '/')
#count_file_name[1:10]
count_file_name <- sapply(count_file_name, function(x){x[2]})
#head(count_file_name)

#construct a empty data matrix
matrix <- data.frame(matrix(nrow=60666, ncol=0))

#read one by one and combine
for (i in 1:length(count_file)){
  path = paste0("gdc_download_20250212_234100.651941//", count_file[i]) # the folder name of counts 
  data <- read.delim(path, fill = T, header = F, row.names =1)
  colnames(data) <- data[2,]
  
  #3:unstranded, counts; 4:stranded_first; 5:stranded_second; 6:tpm_unstranded; 7:fpkm_unstranded; 8:fpkm_yq_unstranded
  data <- data[3]
  #data <- data[6]
  #data <- data[7]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix, data)
}

#transfer to gene_symbol
matrix0 <- matrix[-c(1:6),]
path =paste0("gdc_download_20250212_234100.651941//", count_file[1])
data <- as.matrix(read.delim(path, fill = T, header=F, row.names=1))
gene_name <- data[-c(1:6),1]
matrix1 <-cbind(gene_name, matrix0)
gene_type <- data[-c(1:6),2]
matrix2 <- cbind(gene_type, matrix1)

#remove repeat genes and keep max, min or mean
matrix2 <- aggregate(. ~gene_name, data=matrix2, max)
table(gene_name)

#keep mRNA
matrix3 <- subset(x=matrix2, gene_type =="protein_coding")
table(gene_type)

#make gene_name as col and transfer to derive format
rownames(matrix3) <- matrix3[,1]
matrix3 <- matrix3[,-c(1,2)]

write.table(matrix3, "TCGA_CESC_count.txt", sep="\t", quote=F, row.names=T)




