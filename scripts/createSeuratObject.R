'Create Seurat Object

Usage:
  createSeuratObject.R <feature_matrix> <cell_annotation> <peak_data> <outfile> 

Options:
  -h --help     Show this screen.
  --version     Show version.
' -> doc
library(docopt)
args <- docopt(doc, version = '1.0')

# preconditions
if (!dir.exists(args$feature_matrix))   stop(args$feature_matrix, " doesn't exist.")
if (!file.exists(args$cell_annotation)) stop(args$cell_annotation, " doesn't exist.")
if (!file.exists(args$peak_data))       stop(args$peak_data, " doesn't exist.")

# load libraries
library(Seurat)
library(data.table)
library(purrr)
library(reticulate)
library(data.table)
library(ggplot2)
library(Seurat)
library(Signac)
library(msigdbr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MOFA2)

###############################################################################
message("1) read data")########################################################
###############################################################################
counts <- Seurat::Read10X(args$feature_matrix)


###############################################################################
message("2) create seurat object (with RNAseq and ATAC assay)")################
###############################################################################
seurat <- CreateSeuratObject(counts    = counts["Gene Expression"][[1]],
                             project   = "scRNA+scATAC PBMC",
                             min.cells = 1)
seurat[["ATAC"]] <- CreateAssayObject(counts = counts["Peaks"][[1]])


###############################################################################
message("3) add broad celltypes [skipped]")####################################
###############################################################################
metadata <- data.table::fread(args$cell_annotation) 
metadata$barcode <- metadata$V1
metadata$V1 <- NULL

rename_celltypes <- c(
    "B cell progenitor"      = "Lymphoid",
    "CD4 Memory"             = "Lymphoid",
    "CD8 Naive"              = "Lymphoid",
    "NK cell"                = "Lymphoid",
    "CD4 Naive"              = "Lymphoid",
    "Dendritic cell"         = "Lymphoid",
    "pDC"                    = "Lymphoid",
    "CD8 effector"           = "Lymphoid",
    "Double negative T cell" = "Lymphoid",
    "pre-B cell"             = "Lymphoid",
    "CD14+ Monocytes"         = "Myeloid",
    "CD16+ Monocytes"         = "Myeloid"
)
metadata$celltype       <- metadata$predicted.id
metadata$predicted.id   <- NULL
metadata$broad_celltype <- stringr::str_replace_all(metadata$celltype,
                                                    rename_celltypes)

###############################################################################
message("4) Quality Control")##################################################
###############################################################################
# DT[, col := val]   # update (or add) a column called "col" with value "val".
# DT[i, col := val]  # same as above, but only for those rows specified in i.
###############################################################################
dt <- data.table(barcode = colnames(seurat$RNA)) %>%
      merge(metadata,
            by = "barcode",
            all.x = TRUE) %>%
      .[,c("pass_rnaQC", "pass_accQC"):=FALSE] %>%
      .[!is.na(broad_celltype), c("pass_rnaQC", "pass_accQC"):=TRUE] %>%
      tibble::column_to_rownames("barcode")
seurat <- AddMetaData(seurat, dt)
# seurat <- seurat %>% . [,seurat@meta.data$pass_accQC==TRUE & 
#                          seurat@meta.data$pass_rnaQC==TRUE]


###############################################################################
message("5) add feature and peak metadata")####################################
###############################################################################
feature_metadata      <- fread(file.path(args$feature_matrix, 
                                         "features.tsv.gz")) %>% 
                         setnames(c("ens_id",
                                    "gene",
                                    "view",
                                    "chr",
                                    "start",
                                    "end"))
feature_metadata.rna  <- feature_metadata[view == "Gene Expression"]
feature_metadata.atac <- feature_metadata[view == "Peaks"] %>% 
                         .[,ens_id:=NULL] %>%
                         setnames("gene","peak")
feature_peakdata.atac <- fread(args$peak_data) %>%
                         .[,c("peak","peak_type")] %>% 
                         .[peak_type %in% c("distal", "promoter")]
feature_peakdata.atac <- feature_peakdata.atac[,peak:=sub("_",":",peak)]
feature_peakdata.atac <- feature_peakdata.atac[,peak:=sub("_","-",peak)]
feature_metadata.atac <- feature_metadata.atac %>%
                         merge(feature_peakdata.atac,
                               by = "peak",
                               all.x = TRUE)

###############################################################################
message("6) add chromatin assay")##############################################
###############################################################################
pfm <- getMatrixSet(JASPAR2020, opts = list(species="Homo sapiens"))

for (i in c("distal","promoter")) {
    # Create GRanges
    peaks.granges <- feature_metadata.atac %>%
                   .[peak_type == i] %>%
                   .[,c("chr", "start", "end", "peak")] %>%
                   makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                                            ignore.strand = TRUE) 
    peaks.granges <- keepStandardChromosomes(peaks.granges,
                                             pruning.mode="coarse")

    # create motif matrix
    motif.matrix <- CreateMotifMatrix(
                      features   = peaks.granges,
                      pwm        = pfm,
                      genome     = BSgenome.Hsapiens.UCSC.hg38,
                      use.counts = FALSE) %>% 
                    as.matrix

    # AddChromatinAssay to the Seurat object
    seurat@assays[[paste0("ATAC_",i)]] <- CreateChromatinAssay(
    seurat@assays$ATAC@counts[peaks.granges$peak,], 
    ranges = peaks.granges,
    motifs = CreateMotifObject(motif.matrix, pfm)
    ) 
}

###############################################################################
message("7) normalize data")###################################################
###############################################################################
seurat <- NormalizeData(seurat,
                        normalization.method = "LogNormalize",
                        assay = "RNA")
seurat <- ScaleData(seurat,
                    do.center = TRUE,
                    do.scale = FALSE)
for (i in c("ATAC_distal","ATAC_promoter")) {
    seurat <- RunTFIDF(seurat, assay = i)
}

###############################################################################
message("8) feature analysis")#################################################
###############################################################################
seurat <- FindVariableFeatures(seurat, 
                               selection.method = "vst", 
                               nfeatures        = 5000,
                               assay            = "RNA",
                               verbose          = FALSE)
for (i in c("ATAC_distal","ATAC_promoter")) {
    seurat <- FindTopFeatures(seurat, 
                              assay      = i,
                              min.cutoff = 2000)
}

###############################################################################
message("9) save seurat object")###############################################
###############################################################################
saveRDS(seurat, paste0(args$outfile, "_seurat.rds"))