library(Seurat)
library(data.table)
library(purrr)

basedir <- "/Users/ricard/data/10x_rna_atac/original"
outfile <- "/Users/ricard/data/10x_rna_atac/seurat.rds"


###############
## Load data ##
###############

# counts <- Read10X_h5(paste0(basedir,"/filtered_feature_bc_matrix.h5"))
counts <- Read10X(paste0(basedir,"/filtered_feature_bc_matrix"))

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

metadata <- fread("/Users/ricard/data/10x_rna_atac/sample_metadata.csv") %>%
  .[,barcode:=gsub("-1","",barcode)]

rename_celltypes <- c(
  "naive CD4 T cells" = "Lymphoid",
  "memory CD4 T cells" = "Lymphoid",
  "naive CD8 T cells" = "Lymphoid",
  "CD56 \\(bright\\) NK cells" = "Lymphoid",
  "CD56 \\(dim\\) NK cells" = "Lymphoid",
  "memory B cells" = "Lymphoid",
  "naive B cells" = "Lymphoid",
  "effector CD8 T cells" = "Lymphoid",
  "MAIT T cells" = "Lymphoid",
  "non-classical monocytes" = "Myeloid",
  "intermediate monocytes" = "Myeloid",
  "classical monocytes" = "Myeloid",
  "myeloid DC" = "Myeloid",
  "plasmacytoid DC" = "Lymphoid"
)
metadata$broad_celltype <- stringr::str_replace_all(metadata$celltype,rename_celltypes)

dt <- data.table(barcode=colnames(seurat)) %>%
  merge(metadata,by="barcode", all.x=TRUE) %>%
  .[,c("pass_rnaQC","pass_accQC"):=FALSE] %>%
  .[!is.na(celltype),c("pass_rnaQC","pass_accQC"):=TRUE] %>%
  tibble::column_to_rownames("barcode")

seurat <- AddMetaData(seurat, dt)

head(seurat@meta.data)

##########
## Save ##
##########

saveRDS(seurat, outfile, compress = FALSE)
