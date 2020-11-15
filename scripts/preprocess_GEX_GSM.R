# Loading Libraries
library(ArchR)
library(dplyr)
library(Seurat)
library(ggplot2)

# ++++++++++++++++++++++++++++++++++ 
# ++++++++++ Gene Scores +++++++++++
# ++++++++++++++++++++++++++++++++++

# navigate to PBMC project with filtered doublets
PATH_TO_DATA <- readline(prompt="Enter path to data: ")
FILENAME_ARROWFILE <- readline(prompt="Enter file name of arrow file: ")
FILENAME_SEURAT <- readline(prompt="Enter file name of seurat expression file: ")

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
colnames(gsm_df) <- substring(colnames(gsm_df), 11, 30)

# ++++++++++++++++++++++++++++++++++ 
# ++++++++ Gene Expression +++++++++
# ++++++++++++++++++++++++++++++++++ 

# Load the PBMC dataset
seurat_obj <- readRDS(file = FILENAME_SEURAT)

data_matrix = GetAssayData(object=seurat_obj, assay="RNA", slot="data")
gex_df = as.data.frame(data_matrix)

# ++++++++++++++++++++++++++++++++++
# +++++ Ordering and Filtering +++++
# ++++++++++++++++++++++++++++++++++

# check which genes are in both datasets
gsm_df_rn <- rownames(gsm_df)
gex_df_rn <- rownames(gex_df)

intersect(gex_df_rn, gsm_df_rn) %>% length # number of genes that are in both datasets
setdiff(gsm_df_rn, intersect(gex_df_rn, gsm_df_rn))
setdiff(gex_df_rn, intersect(gex_df_rn, gsm_df_rn))

# keep only rows and cells that are in the intersect
rows.to.keep_gex<-which(rownames(gex_df) %in% intersect(gex_df_rn, gsm_df_rn))
rows.to.keep_gsm<-which(rownames(gsm_df) %in% intersect(gex_df_rn, gsm_df_rn))

gex_df <- gex_df[rows.to.keep_gex, intersect(colnames(gex_df), colnames(gsm_df))]
gsm_df <- gsm_df[rows.to.keep_gsm, intersect(colnames(gex_df), colnames(gsm_df))]

# order rows and columns in increasing order
gex_df <-  gex_df[, order(names(gex_df))]
gex_df <-  gex_df[order(row.names(gex_df)), ]
gsm_df <-  gsm_df[, order(names(gsm_df))]
gsm_df <-  gsm_df[order(row.names(gsm_df)), ]

# some basic exploratory plots
gsm_df_flipped = t(gsm_df)
gex_df_flipped = t(gex_df)
d1 = density(gsm_df_flipped[ ,2000])
plot(d1)
max(gex_df_flipped[,210])
d2 = hist(gex_df_flipped[ ,2000])
plot(d2)

# dataframe with the average expected, observed, and difference for each gene
exp_obv_diff = data.frame(exp=double(), obv=double(), abs_diff=double(), rel_diff=double())

for(i in 1:ncol(gsm_df_flipped)) {
  exp_i = sum(gsm_df_flipped[,i])/nrow(gsm_df_flipped)
  obv_i = sum(gex_df_flipped[,i])/nrow(gex_df_flipped)
  abs_diff_i = sum(c(abs(exp_i - obv_i)))
  rel_diff_i = abs_diff_i/obv_i
  
  exp_obv_diff[i, 1] = exp_i
  exp_obv_diff[i, 2] = obv_i
  exp_obv_diff[i, 3] = abs_diff_i
  exp_obv_diff[i, 4] = rel_diff_i
}

row.names(exp_obv_diff) <- colnames(gsm_df_flipped)

# function to plot expected and observed and relative difference for any number of genes
plot_diffs <- function(dataframe = exp_obv_diff, n = 20) {
  ggplot(exp_obv_diff[1:n,]) +
    geom_point(aes(x=row.names(exp_obv_diff[1:n,]), y = exp), colour='red', alpha=0.7) +
    geom_point(aes(x=row.names(exp_obv_diff[1:n,]), y = obv), colour='blue', alpha=0.7) +
    geom_point(aes(x=row.names(exp_obv_diff[1:n,]), y = abs_diff), colour='green', alpha=0.7) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

plot_diffs(n=10)

# plot density estimate
ggplot(exp_obv_diff, aes(diffs)) + geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + xlim(0,5)
