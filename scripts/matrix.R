# load ArchR
library(ArchR)

# navigate to PBMC project with filtered doublets
setwd("/home/james/Documents/leuven/second-year/IBP/ArchR/ProjPBMC/ArrowFiles")

# get GeneScoreMatrix summarised experiment for 3k unsorted
gsm_se <- getMatrixFromArrow('3k_unsorted.arrow', 'GeneScoreMatrix')

# view summarised experiment data
rowData(gsm_se)
colData(gsm_se)

# get gsm from assays
gsm = assays(gsm_se)[['GeneScoreMatrix']]
fix(gsm)
  