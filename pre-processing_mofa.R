install.packages("stringr")
library("stringr")
data=fread("C:\\Users\\user\\Desktop\\pbmc_unsorted_3k_atac_peak_annotation.tsv")
# df <- data.frame(peaks=character(),
#                  gene=character(), 
#                  distance=integer(), 
#                  peak_type=character()) 
# for (i in 1:length(data$peak)){
#   if (str_detect(data$peak_type[i],";")==TRUE){
#     print("entered")
#     gene_sep=str_split(data$gene[i],";")[[1]]
#     dis_sep=str_split(data$distance[i],";")[[1]]
#     peak_type_sep=str_split(data$peak_type[i],";")[[1]]
#     for (j in 1:length(gene_sep)){
#       row_extra <- data.frame(data$peak[i], gene_sep[j], dis_sep[j],    peak_type_sep[j]) 
#       names(row_extra) <- c("peak", "gene", "distance", "peak_type") 
#       df <- rbind(df, row_extra)
#     }
#   }
#   else{
#     rbind(df,data[i])
#   }
# }
# # also maybe remove duplictate lines
# df_unique <- df[!duplicated(df), ] 


for ( i in 1: length(data$peak)){
  row_values=str_split(data$peak[i],"_")[[1]]
  temp_row_values=paste(row_values[1],":",row_values[2],"-",row_values[3],sep='')
  data$peak[i]=temp_row_values
}


write.table(data, file='pre_processed_pbmc_unsorted_3k_atac_peak_annotation_temp.tsv', quote=FALSE, sep='\t',row.names = FALSE)
