library(Seurat)
library(data.table)
library(purrr)

#basedir <- "/Users/ricard/data/10x_rna_atac/original"
#outfile <- "/Users/ricard/data/10x_rna_atac/seurat.rds"


###############
## Load data ##
###############

# counts <- Read10X_h5(paste0(basedir,"/filtered_feature_bc_matrix.h5"))
counts <- Read10X(paste0("C:\\Users\\user\\Desktop\\pbmc_unsorted_3k_filtered_feature_bc_matrix\\","filtered_feature_bc_matrix"))


###################
## Create Seurat ##
###################

seurat <- CreateSeuratObject(
  counts = counts["Gene Expression"][[1]],
  project = "scRNA+scATAC PBMC",
  min.cells = 1
)

# Add ATAC modality
seurat[["ATAC"]] <- CreateAssayObject(counts = counts["Peaks"][[1]])

seurat

##################
## Add metadata ##
##################

#metadata <- fread("C:\\Users\\user\\Desktop\\sample_metadata.csv") %>% .[,barcode:=gsub("-1","",barcode)]
metadata <- fread("C:\\Users\\user\\Desktop\\metadata_file_3k_unsorted.csv")
rename_celltypes <- c(
  "B cell progenitor"="Lymphoid",
  "CD4 Memory"="Lymphoid",
  "CD8 Naive"="Lymphoid",
  "NK cell"="Lymphoid",
  "CD14 Monocytes"= "Myeloid",
  "CD4 Naive"="Lymphoid",
  "Dendritic cell"="Lymphoid",
  "pDC"="Lymphoid",
  "CD16 Monocytes"="Myeloid",
  "CD8 effector"="Lymphoid",
  "Double negative T cell"="Lymphoid",
  "pre-B cell"="Lymphoid"
)
metadata$broad_celltype <- stringr::str_replace_all(metadata$celltype,rename_celltypes)

dt <- data.table(barcode=colnames(seurat)) %>%
  merge(metadata,by="barcode",all.x= TRUE) %>%
  .[,c("pass_rnaQC","pass_accQC"):=FALSE] %>%
  .[!is.na(celltype),c("pass_rnaQC","pass_accQC"):=TRUE] %>%
  tibble::column_to_rownames("barcode")


seurat <- AddMetaData(seurat,dt)

head(seurat@meta.data)

