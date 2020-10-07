install.packages('devtools', clean=TRUE, quiet=TRUE)
install.packages('remotes', clean=TRUE, quiet=TRUE)

install.packages("ggpubr")
# installing bioconductor and dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# automatically install Bioconductor dependencies
setRepositories(ind=1:2)

# core packages needed
install.packages('data.table', clean=TRUE, quiet=TRUE)
install.packages('ggplot2', clean=TRUE, quiet=TRUE)

# automatically install Bioconductor dependencies
setRepositories(ind=1:2)

# installing packages from Bioconductor
BiocManager::install(c(
  'AnnotationFilter',
  'BiocGenerics',
  'GenomeInfoDb',
  'GenomicFeatures',
  'GenomicRanges',
  'IRanges',
  'Rsamtools',
  'S4Vectors',
  'JASPAR2020',
  'TFBSTools',
  'ggbio',
  'motifmatchr',
  'AnnotationDbi',
  'Seurat',
  'Signac',
  'msigdbr',
  'BSgenome.Hsapiens.UCSC.hg38')
)

BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
# installing MOFA2


devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))

install.packages("reticulate")
library(reticulate)
use_python("C:\\Users\\user\\anaconda3", required=TRUE)

#use_condaenv(condaenv = "C:\\Users\\user\\Anaconda3\\Lib\\site-packages", conda = "auto", required = FALSE)

library(data.table)
library(ggplot2)
library(Seurat)
library(Signac)

# for GSEA analysis
library(msigdbr)

library(ggpubr)

# For motif enrichment analysis
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

# MOFA
library(MOFA2)


# load the data
# for now we are going to ftp the counts from ricard
s <- readRDS(url("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/seurat.rds"))
seurat=s
seurat



head(seurat@meta.data[,c("celltype","broad_celltype","pass_rnaQC","pass_accQC")])
# this will be the count of each cell type
table(seurat@meta.data$celltype)
# if we want the broad cell type count
table(seurat@meta.data$broad_celltype)

# now we have to do some pre-processing to check for the quality of the reads
seurat <- seurat %>% . [,seurat@meta.data$pass_accQC==TRUE & seurat@meta.data$pass_rnaQC==TRUE]

seurat
# therefore this expression of RNA is 29732 genes and 10032 cells
# we can also list the features -- these will be gene names 
# they can we entered into ensembl hg to see what they are
seurat@assays[["RNA"]]

# https://www.youtube.com/watch?v=OjFFPIENrWQ 
# here he discusses how in the secondary analysis, they create a matrix of cells by peaks where the values will represent how open the region of the genome was per cell
# to see how many peaks we have:
seurat@assays[["ATAC"]] 

## Load additional information
# we will need a position specific matirix from the jaspar database
#this will get a list of motif position freq from jaspar
pfm=getMatrixSet(JASPAR2020,opts=list(species="Homo sapiens"))


# load metadata
# metadata will have information on both RNA and ATAC peaks
feature_metadata =fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/filtered_feature_bc_matrix/features.tsv.gz") %>% setnames(c("ens_id","gene","view","chr","start","end"))
# fetch the rna metadata
feature_metadata.rna =feature_metadata[view=="Gene Expression"]
head(feature_metadata.rna,n=3)
# fetch the atac metadata
feature_metadata.atac=feature_metadata[view=="Peaks"] %>% .[,ens_id:=NULL] %>% setnames("gene","peak")
head(feature_metadata.atac,n=3)

# now we will classify the peaks into promoter overlapping and distal

foo <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/atac_peak_annotation.tsv") %>% .[,c("peak","peak_type")] %>% .[peak_type%in%c("distal", "promoter")]

feature_metadata.atac <- feature_metadata.atac %>% merge(foo,by="peak",all.x=TRUE)

table(feature_metadata.atac$peak_type)



# Parse the Seurat Object
# here we are splitting based on peaktype and creating a chromatin assay.
# so you are making 3 assays from the ATAC assay: ATAC, ATAC_distal and ATAC_promoter
for (i in c("distal","promoter")) {
  
  # Create GRanges
  peaks.granges <- feature_metadata.atac %>% .[peak_type==i] %>% .[,c("chr","start","end","peak")] %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)
  
  peaks.granges <- keepStandardChromosomes(peaks.granges, pruning.mode="coarse")
  
  # Scan motifs throughout the DNA sequence of each peak and create a binary matrix of motif-peak presence.
  # scan the DNA sequence of each peak for the presnce of each motif
  motif.matrix <- CreateMotifMatrix(
    features = peaks.granges,
    pwm = pfm,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    use.counts = FALSE) %>% as.matrix
  
  # AddChromatinAssay to the Seurat object
  seurat@assays[[paste0("ATAC_",i)]] <- CreateChromatinAssay(
    seurat@assays$ATAC@counts[peaks.granges$peak,], 
    ranges = peaks.granges,
    motifs = CreateMotifObject(motif.matrix, pfm)
  )
  
}
seurat

## SECTION 6: NORMALIZATION

# now we need to normalize

## now for RNA
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", assay = "RNA")
seurat <- ScaleData(seurat, do.center = TRUE, do.scale = FALSE)

## now for ATAC
for (i in c("ATAC_distal","ATAC_promoter")) {
  seurat <- RunTFIDF(seurat, assay = i)
}

## SECTION 7: FEATURE SELECTION!!!
# 7.1 RNA finding most variable genes using VST 
seurat <- FindVariableFeatures(seurat, 
                               selection.method = "vst", 
                               nfeatures = 5000,
                               assay = "RNA",
                               verbose = FALSE
)
# 7.2 ATAC finding variable peaks (ATAC seq dta gives peaks and not genes like RNA does)
for (i in c("ATAC_distal","ATAC_promoter")) {
  seurat <- FindTopFeatures(seurat, assay=i, min.cutoff = 2000)
  print(length(seurat[[i]]@var.features))
}

# so there are 8456 top features for distal and 11198 for promoter

## SECTION 8 : TRAIN THE MODEL
# 8.1 create the mofa object ( this mofa fuction can take a seurat object)
mofa <- create_mofa(seurat, assays = c("RNA","ATAC_distal","ATAC_promoter"))# this will throw an error its okay
mofa
# this function is to check the groups(columns) and their dimensions
# if there are grey bars, they are the missing information
# remember: MOFA , once trained, can impute misssing values
plot_data_overview(mofa)
data_opt=get_default_data_options(mofa)
data_opt
# we can see that the groups(columns) arent scaled to have same total varience and same with the views(rows)
# peaks.granges@seqnames@values

slotNames(mofa)
names(mofa@data)

# 8.2: define mofa options
# the most important parameter is the number of factors
# this NEEDS to be specified and usually needs to be large enough

model_opts=get_default_model_options(mofa)
model_opts$num_factors=15

mofa=prepare_mofa(mofa,model_options = model_opts)

# 8.3: run mofa and possibly fail check reticulate
# mofa=run_mofa(mofa)
mofa <- readRDS(url("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/mofa.rds"))

# SECTION 9: MOFA downstream analysis

#9.1 Add the metadata to the model
samples_metadata(mofa) <- seurat@meta.data %>% 
  tibble::rownames_to_column("sample") %>%
  as.data.table
#9.2 Assess correlation between factors
# this is very imp we dont want correlation between our factors this is a primary assumption of principal component analysis
plot_factor_cor(mofa)
# 9.3: variance decomposition
plot_variance_explained(mofa,max_r2=4)
plot_variance_explained(mofa,plot_total = TRUE)[[2]]

# 9.4 Characterisation of factors
# 9.4.1 using covariates
correlate_factors_with_covariates(mofa,covariates=c("nFeatures_RNA","nFeature_ATAC"))
# we see that factor 1 is the strongest source of variation becuase it is associated to the number of expressed genes per cell and the number of accessible peaks per cell
# now we want to visualize
plot_factor(mofa, factors=1, group_by = "celltype", color_by="broad_celltype") +
  theme(
    axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
  )
# 9.4.2 visualize usign feauture weights
plot_weights(mofa, 
             view = "RNA", 
             factors = 1, 
             nfeatures = 20, 
             text_size = 4
)
# visualisation of covar pattersn
# using RNA data
plot_data_scatter(mofa, 
                  view = "RNA", 
                  factor = 1, 
                  features = 6,
                  color_by = "broad_celltype",
                  add_lm = T,
                  dot_size = 1
)
# using ATAC data

plot_top_weights(mofa, 
                 view = "ATAC_distal", 
                 factors = 1, 
                 sign = "positive",
                 nfeatures = 15,
)
plot_top_weights(mofa, 
                 view = "ATAC_promoter", 
                 factors = 1, 
                 sign = "positive",
                 nfeatures = 15,
)
# we can also see this with a heatmap
plot_data_heatmap(mofa, 
                  view = "ATAC_promoter", 
                  factor = 1, 
                  features = 50,
                  show_rownames = F, show_colnames = F, 
                  cluster_rows = T, cluster_cols = F,
                  annotation_samples = "broad_celltype"
)
# we can do denoising 
plot_data_heatmap(mofa, 
                  view = "ATAC_promoter", 
                  factor = 1, 
                  features = 50,
                  show_rownames = F, show_colnames = F, 
                  cluster_rows = T, cluster_cols = F,
                  annotation_samples = "broad_celltype",
                  denoise = TRUE
)
# 9.4.3 now we can do factor 2
plot_factor(mofa, factors=2, group_by = "celltype", color_by="broad_celltype") +
  theme(
    axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
  )
# plot the assocaietd RNA weights
plot_weights(mofa, 
             view = "RNA", 
             factors = 2, 
             nfeatures = 10, 
             text_size = 4
)
# now color by BANK1 expression
plot_factor(mofa, factors=2, group_by = "celltype", color_by="BANK1") +
  theme(
    axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
  )
# and visualize the ATAC signatures
plot_data_heatmap(mofa, 
                  view = "ATAC_promoter", 
                  factor = 2, 
                  features = 50,
                  show_rownames = F, show_colnames = F, 
                  cluster_rows = T, cluster_cols = F,
                  annotation_samples = "celltype",
                  denoise = TRUE
)

# 9.4.4 characterisation of factor 4
# this could be a source of technical variation
# we could do another type of normalization
plot_factor(mofa, factors=4, group_by = "celltype", color_by="nFeature_ATAC") +
  theme(
    axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
  )
plot_weights(mofa, 
             view = "ATAC_promoter", 
             factors = 4, 
             nfeatures = 0, 
             text_size = 4
)
# 9,5 non linear dimensionality reduction
# 9.5.1 using MOFA factors
factors <- 1:get_dimensions(mofa)[["K"]]
factors <- factors[!factors%in%c(4,7)]

mofa <- run_umap(mofa, 
                 factors = factors, 
                 n_neighbors = 15,  
                 min_dist = 0.30
)

temp1=plot_dimred(mofa, 
            method = "UMAP", 
            color_by = "celltype", 
            label = TRUE, 
            stroke=0.05, 
            dot_size = 1, 
            legend = FALSE
) + scale_fill_manual(values=colors)
print(temp1)

for (i in paste0("Factor",1:3)) {
  p <- plot_dimred(mofa, 
                   method = "UMAP", 
                   color_by = i, 
                   stroke = 0.05, 
                   dot_size = 1
  )
  print(p)
}

K=15

# 9.5.2 use rna alone
DefaultAssay(seurat) <- "RNA"
seurat <- RunPCA(seurat, npcs = K, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = 'pca', dims = 1:K, verbose = FALSE)

# okay so this command:
# The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
## DISCUSS THIS!!!

DimPlot(seurat, label = TRUE, reduction="umap") + 
  NoLegend() + NoAxes() + scale_fill_manual(values=colors)

# 9.5.3 use the atac data alone
DefaultAssay(seurat) <- "ATAC_distal"
seurat <- RunSVD(seurat, n = K, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 1:K, verbose = FALSE)
DimPlot(seurat, label = TRUE, reduction="umap") + 
  NoLegend() + NoAxes() + scale_fill_manual(values=colors)

# 9.6 Gene set enrichment analysis
# 9.6.1 load the gene set annotations
msgidb.matrix <- msigdbr(
  species = "Homo sapiens",
  category = "C5", 
  subcategory = "BP"
) %>% as.data.table %>% .[,id:=1] %>%
  dcast(gs_name~gene_symbol, value.var="id", fill=0) %>% 
  matrix.please
# matrix please:P

# 9.6.2 run GSEA
# GSEA on positive weights
gsea.positive <- run_enrichment(mofa, 
                                feature.sets = msgidb.matrix, 
                                view = "RNA",
                                sign = "positive"
)

# GSEA on negative weights
gsea.negative <- run_enrichment(mofa, 
                                feature.sets = msgidb.matrix, 
                                view = "RNA",
                                sign = "negative"
)

names(gsea.positive)

# 9.6.3 visualise GSEA results
plot_enrichment(gsea.positive, factor = 1, max.pathways = 15)
plot_enrichment(gsea.negative, factor = 1, max.pathways = 15)
plot_enrichment_detailed(gsea.positive,
                         factor = 1,
                         max.genes = 10,
                         max.pathways = 5
)

# 9.7 motif enrichment
# 9.7.1 Run
# define motif matrix
motif.matrix <- t(as.matrix(seurat[["ATAC_distal"]]@motifs@data))

# Run GSEA enrichment analysis using the motif-peak matrix, (+) weights
motif.enrichment.positive <- run_enrichment(mofa,
                                            view = "ATAC_distal", 
                                            factors = 1:2,
                                            feature.sets = motif.matrix,
                                            sign = "positive"
)

# Run GSEA enrichment analysis using the motif-peak matrix, (-) weights
motif.enrichment.negative <- run_enrichment(mofa,
                                            view = "ATAC_distal", 
                                            factors = 1:2,
                                            feature.sets = motif.matrix,
                                            sign = "negative"
)

# 9.7.2 Visualise
plot_enrichment(motif.enrichment.positive, factor = 1, max.pathways = 15)
plot_enrichment(motif.enrichment.negative, factor = 1, max.pathways = 15)

sig.motifs.positive <- motif.enrichment.positive$pval.adj[,"Factor1"] %>%
  sort %>% head(n=6) %>% names
MotifPlot(seurat[["ATAC_distal"]], motifs = sig.motifs.positive)

# 9.7.3 Validation using ChromVar
seurat <- RunChromVAR(
  object = seurat, 
  assay = "ATAC_distal",
  genome = BSgenome.Hsapiens.UCSC.hg38, 
  new.assay.name = "chromvar"
)

DefaultAssay(seurat) <- "chromvar"

motifs.to.plot <- c(sig.motifs.positive[1:2], sig.motifs.negative[1:2])

FeaturePlot(seurat,
            features = motifs.to.plot,
            reduction = "umap",
            combine = TRUE
) & NoLegend() & NoAxes()


# packageVersion()/ sessionInfo()