# Loading Libraries
library(ArchR)
library(dplyr)
library(Seurat)

# _____________________________________________
# ---------------- Gene Scores ----------------
# _____________________________________________

# navigate to PBMC project with filtered doublets
PATH_TO_DATA <- readline(prompt="Enter path to data: ")  # /Users/nrank/Desktop/BioInf/IntegratedProject/TrainedMOFA_ArchR/
FILENAME_ARROWFILE <- readline(prompt="Enter file name of arrow file: ")   # 10k_sorted.arrow
FILENAME_SEURAT <- readline(prompt="Enter file name of seurat gene expression file: ") # seurat_10k_sorted.rds

setwd(PATH_TO_DATA)

# get GeneScoreMatrix summarised experiment
gsm_se <- getMatrixFromArrow(FILENAME_ARROWFILE, 'GeneScoreMatrix')

# get gsm from assays
gsm = assays(gsm_se)[['GeneScoreMatrix']]

# getting gene names
row_data = rowData(gsm_se)
genes = row_data["name"]

# convert to dataframe
gsm_df <- as.data.frame(as.matrix(gsm))
row.names(gsm_df) <- genes$name

# remove heading from colnames
# colnames(gsm_df) <- substring(colnames(gsm_df), 11, 30)
colnames(gsm_df) <- substring(colnames(gsm_df), nchar(colnames(gsm_df))-18+1, nchar(colnames(gsm_df)))    # NEW

# _____________________________________________
# ---------- GEX ------------
# _____________________________________________

# Load the PBMC dataset
seurat_obj <- readRDS(file = FILENAME_SEURAT)

data_matrix = GetAssayData(object=seurat_obj, assay="RNA", slot="data")
gex_df = as.data.frame(data_matrix)

# _____________________________________________
# --- Ordering and Filtering ---
# _____________________________________________
# check which genes are in both datasets
gsm_df_rn <- rownames(gsm_df)
gex_df_rn <- rownames(gex_df)

#intersect(gex_df_rn,gsm_df_rn) %>% length # number of genes that are in both datasets
#setdiff(gsm_df_rn, intersect(gex_df_rn,gsm_df_rn)) %>% length
#setdiff(gex_df_rn, intersect(gex_df_rn,gsm_df_rn))

# keep only rows and cells that are in the intersect 
rows.to.keep_gex<-which(rownames(gex_df) %in% intersect(gex_df_rn,gsm_df_rn))
rows.to.keep_gsm<-which(rownames(gsm_df) %in% intersect(gex_df_rn,gsm_df_rn))

# all genes that left over after the first intersect (isoforms, etc.)
# left.rownames.gex <- gex_df[-rows.to.keep_gex, ] %>% rownames() # NEW
# left.rownames.gsm <- gsm_df[-rows.to.keep_gsm, ] %>% rownames() # NEW
# gex_df_left <- gex_df[-rows.to.keep_gex, ]
# gsm_df_left <- gsm_df[-rows.to.keep_gsm, ]
# gex_df_left <-  gex_df_left[order(row.names(gex_df_left)), ]
# gsm_df_left <-  gsm_df_left[order(row.names(gsm_df_left)), ]


# genes that are orignially in the intersect of GEX, ATAQ
gex_df <- gex_df[rows.to.keep_gex, intersect(colnames(gex_df), colnames(gsm_df))]
gsm_df <- gsm_df[rows.to.keep_gsm, intersect(colnames(gex_df), colnames(gsm_df))]

# order rows and columns in increasing order
gex_df <-  gex_df[, order(names(gex_df))]
gex_df <-  gex_df[order(row.names(gex_df)), ]
gsm_df <-  gsm_df[, order(names(gsm_df))]
gsm_df <-  gsm_df[order(row.names(gsm_df)), ]



# check correlation between ranks of genes - all genes
# taus <- c()
# 
# for (cn in c(1:ncol(gsm_df))){
#   print(cn)
#   a <- cor.test(gsm_df[, cn], gex_df[, cn], method='kendall')
#   tau <- a$estimate[[1]]
#   taus <- c(taus, tau)
# }
# 
# hist(taus)
# saveRDS(taus, 'taus_3k_sorted.rds')

#a <- cor.test(gsm_df[, 1], gex_df[, 1], method='kendall')
#a$estimate[[1]]
#cor.test(x, y, method = "spearman") in its "stats" package (also cor(x, y, method = "spearman") will work


# check correlation between ranks of genes -  upper 15% of the genes
# comparison of upper 15% ranks
taus_allCells <- c()
common_prop_allCells <- c()

for (cell in c(1:ncol(gsm_df))){
# order genes by count/score per cell
#cell <-  1
gsm_allCells <- gsm_df
gex_allCells <- gex_df

gsm_allCells$rown <- rownames(gsm_allCells)
gex_allCells$rown <- rownames(gex_allCells)

gsm_cell <- gsm_allCells %>% select(c(cell, 'rown'))
gex_cell <- gex_allCells %>% select(c(cell, 'rown'))

cellId <- colnames(gsm_cell)[1]
colnames(gsm_cell) <- c('cellID', 'rown')
colnames(gex_cell) <- c('cellID', 'rown')
gsm_cell <- gsm_cell %>% arrange(desc(cellID))
gex_cell <- gex_cell %>% arrange(desc(cellID))

# cut off the first 15 %
gsm_cell_cut <- gsm_cell[1:round(nrow(gsm_cell)*0.15), ]
gex_cell_cut <- gex_cell[1:round(nrow(gex_cell)*0.15), ]


intersect(gsm_cell_cut$rown, gex_cell_cut$rown)

gsm_cell_int <- gsm_cell_cut %>% filter(rown %in% intersect(gsm_cell_cut$rown, gex_cell_cut$rown))
gex_cell_int <- gex_cell_cut %>% filter(rown %in% intersect(gsm_cell_cut$rown, gex_cell_cut$rown))
gsm_cell_int <-  gsm_cell_int[order(row.names(gsm_cell_int)), ]
gex_cell_int <-  gex_cell_int[order(row.names(gex_cell_int)), ]

common_prop <- length(intersect(gsm_cell_cut$rown, gex_cell_cut$rown))/nrow(gsm_cell_cut)
tau_cell <- cor.test(gsm_cell_int[, 'cellID'], gex_cell_int[, 'cellID'], method='kendall')
tau_cell <- tau_cell$estimate[[1]]
taus_allCells <- c(taus_allCells, tau_cell)
common_prop_allCells <- c(common_prop_allCells, common_prop)
}

par(mfrow=c(1,1))
hist(taus_allCells, main = "", xlab = expression(tau), cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)

par(mfrow=c(1,1))
hist(common_prop_allCells, main = " ", xlab = 'proportion of intersecting genes', cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)

sd_taus <- sd(taus_allCells)
sd_common_prop <- sd(common_prop_allCells)
mean_taus <- mean(taus_allCells)
mean_common_prop <- mean(common_prop_allCells)
cat('sd_taus: ', sd_taus)
cat('mean_taus: ', mean_taus)
cat('sd_common_prop: ', sd_common_prop)
cat('mean_common_prop: ', mean_common_prop)


## NOTES
# ### NEW
# # cut off the .1, .2, ... from each gene name/ rowname
# gex_df_left$rown <- rownames(gex_df_left)
# gex_df_left <- gex_df_left %>% rowwise %>% dplyr::mutate(rn_cut = strsplit(rown, '[.]')[[1]][1]) %>% ungroup
# gex_df_left <- gex_df_left %>% select(-rown) 
# gex_df_left <- setDT(gex_df_left)[, lapply(.SD, sum) , by = .(rn_cut)]
# rownames(gex_df_left) <- gex_df_left$rn_cut
# gex_df_left <- gex_df_left %>% select(-rn_cut)
# 
# 
# # getting all genes that are in ATAQ and isoforms of GEX
# # keep only rows and cells that are in the intersect 
# gsm_df_rn2 <- rownames(gsm_df_left)
# gex_df_rn2 <- rownames(gex_df_left)
# rows.to.keep_gex2<-which(rownames(gex_df_left) %in% intersect(gex_df_rn2,gsm_df_rn2))
# rows.to.keep_gsm2<-which(rownames(gsm_df_left) %in% intersect(gex_df_rn2,gsm_df_rn2))
# 
# gex_left <- gex_df[rows.to.keep_gex2, intersect(colnames(gex_df_left), colnames(gsm_df_left))]
# gsm_left <- gsm_df[rows.to.keep_gsm2, intersect(colnames(gex_df_left), colnames(gsm_df_left))]
# 
# ncol(gex_left)
# ncol(gsm_left)
# 
# nrow(gex_left)
# nrow(gsm_left)
# 
# 
# intersect(gsm_df_rn2, gex_df_rn2)
# 
# 
# 
# 
# #rownames(gex_df_left) <- gex_df_left$rn_cut
# 
# 
# 
# gex_df_left2 <- gex_df_left[c(1:100), ] %>% group_by(rn_cut) %>% summarise_each(funs(sum))
# gex_df_left2 <- gex_df_left %>% group_by(rn_cut) %>% summarise_each(funs(sum))
# 
# gex_df_left2 <- gex_df_left[c(1:1000), c(1:2000, ncol(gex_df_left))] %>% group_by(rn_cut) %>% summarise_each(funs(sum))
# 
# 
# gex_df_left$rn_cut %>% unique %>% length
# 
# library(data.table)
# b <- setDT(gex_df_left)[, lapply(.SD, sum) , by = .(rn_cut)]
# 
# 
# c <- data.frame(a =c(1,2,2,3,4,5,6,5,4), b=c(1,1,1,2,2,2,3,3,7), c=c(1,1,0,0,2,2,3,3,7))
# setDT(c)[, lapply(.SD, sum) , by = .(a)]
# 
# 
# library(dplyr)
# df1 %>%
#   group_by(elevation, distance) %>% 
#   summarise_each(funs(sum))
# 
# 
# gex_df_left %>% select(rown, rn_cut, c(1:10)) %>% head(20) %>% View
# 
# a <- gex_df_left %>% head(20)
# 
# rownames(gex_df_left) <- strsplit(rownames(gex_df_left), '[.]')[[1]][1]
# 
# 
# 
# # getting all genes that are in ATAQ and isoforms of GEX
# # keep only rows and cells that are in the intersect 
# gsm_df_rn2 <- rownames(gsm_df_left)
# gex_df_rn2 <- rownames(gex_df_left)
# rows.to.keep_gex2<-which(rownames(gex_df_left) %in% intersect(gex_df_rn2,gsm_df_rn2))
# rows.to.keep_gsm2<-which(rownames(gsm_df_left) %in% intersect(gex_df_rn2,gsm_df_rn2))
# 
# gex_left <- gex_df[rows.to.keep_gex2, intersect(colnames(gex_df_left), colnames(gsm_df_left))]
# gsm_left <- gsm_df[rows.to.keep_gsm2, intersect(colnames(gex_df_left), colnames(gsm_df_left))]
# 
# 
# 
# 
# CGATTCCTCAGCACCA-1
# 
# 
# left.rownames.gsm <- gsm_df[-rows.to.keep_gsm, ] %>% rownames()
# rows.not.keep <- which(!(rownames(gsm_df) %in% intersect(gex_df_rn,gsm_df_rn)))
# gsm_df[-rows.to.keep_gsm, ] %>% rownames() %>% length
# length(rows.to.keep_gsm)
# length(rows.not.keep)
# nrow(gsm_df)
# order(left.rownames.gex)
# 
# # aggregate by rowname
# df %>% group_by(row.names) %>% summarise_each(sum)
# 
# strsplit("AC000032.1", '[.]')[[1]][1]
# 
