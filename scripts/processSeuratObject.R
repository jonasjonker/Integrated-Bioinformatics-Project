'Create Seurat Object

Usage:
  processSeuratObject.R <feature_matrix> <meta_data> <outfile> 

Options:
  -h --help     Show this screen.
  --version     Show version.
' -> doc
library(docopt)
args <- docopt(doc, version = '1.0')

# preconditions
if (!file.exists(args$seurat_obj))      stop(args$meta_data, " doesn't exist.")

# load libraries
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

use_python("C:\\Users\\user\\anaconda3", required=TRUE)

seurat <- readRDS(paste0(args$outfile, "_seurat.rds"))

seurat <- seurat %>% . [,seurat@meta.data$pass_accQC==TRUE & 
                         seurat@meta.data$pass_rnaQC==TRUE]


pfm <- getMatrixSet(JASPAR2020,opts=list(species="Homo sapiens"))

feature_metadata      <- fread(args$meta_data) %>% 
                         setnames(c("ens_id","gene","view","chr","start","end"))
feature_metadata.rna  <- feature_metadata[view=="Gene Expression"]

# atac metadata
feature_metadata.atac <- feature_metadata[view=="Peaks"] %>% 
                         .[,ens_id:=NULL] %>%
                         setnames("gene","peak")
feature_peakdata.atac <- fread() %>% 
                         .[,c("peak","peak_type")] %>% 
                         .[peak_type %in% c("distal", "promoter")]
feature_metadata.atac <- feature_metadata.atac %>%
                         merge(feature_peakdata.atac,
                               by = "peak",
                               all.x = TRUE)

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

seurat <- NormalizeData(seurat,
                        normalization.method = "LogNormalize",
                        assay = "RNA")
seurat <- ScaleData(seurat,
                    do.center = TRUE,
                    do.scale = FALSE)


for (i in c("ATAC_distal","ATAC_promoter")) {
    seurat <- RunTFIDF(seurat, assay = i)
}

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

saveRDS(seurat, paste0(args$outfile, "_seurat.rds"))
