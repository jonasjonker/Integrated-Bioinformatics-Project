install.packages("stringr")
library("stringr")
data=fread("C:\\Users\\user\\Desktop\\pbmc_unsorted_3k_atac_peak_annotation.tsv")
df <- data.frame(peaks=character(),
                 gene=character(), 
                 distance=integer(), 
                 peak_type=character()) 
for (i in 1:length(data$peak)){
  if (str_detect(data$peak_type[i],";")==TRUE){
    print("entered")
    gene_sep=str_split(data$gene[i],";")[[1]]
    dis_sep=str_split(data$distance[i],";")[[1]]
    peak_type_sep=str_split(data$peak_type[i],";")[[1]]
    for (j in 1:length(gene_sep)){
      row_extra <- data.frame(data$peak[i], gene_sep[j], dis_sep[j],    peak_type_sep[j]) 
      names(row_extra) <- c("peak", "gene", "distance", "peak_type") 
      df <- rbind(df, row_extra)
    }
  }
  else{
    rbind(df,data[i])
  }
}

write.table(df, file='pre_processed_pbmc_unsorted_3k_atac_peak_annotation.tsv', quote=FALSE, sep='\t',row.names = FALSE)

