library(Seurat)

#read data
multiome.data = Read10X("../data/filtered_feature_bc_matrix_pbmc_unsorted_10k")
multiome_RNA_Seurat = CreateSeuratObject(counts = multiome.data[['Gene Expression']], project = "MultiOmeRna", min.cells = 3, min.features = 0)

#EDIT(!) after re-reading this script --> It might be better to not filter cells in this step (otherwise some cells wont have an annotation), so better set min.features at 0

annotated_Seurat = readRDS("../data/pbmc_10k_v3.rds")

#normalize data (can take a long time)
pbmc.list = c(multiome_RNA_Seurat, annotated_Seurat)
for (i in 1:length(pbmc.list)) {
  pbmc.list[[i]] = SCTransform(pbmc.list[[i]], verbose = TRUE)
}

#select top 3000 highly variable features
pbmc.features = SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)

#prepare for integration
pbmc.list = PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features, verbose = TRUE)

#we want to query the multiome data for cell types in the annotated dataset
pbmc.query = pbmc.list[[1]]
pbmc.reference = pbmc.list[[2]]

#find transfer achnors
pbmc.anchors = FindTransferAnchors(reference = pbmc.reference, query = pbmc.query, dims = 1:30, normalization.method = 'SCT', features = pbmc.features)

#predict cell types based on transfer anchors
predictions = TransferData(anchorset = pbmc.anchors, refdata = pbmc.reference$celltype, dims = 1:30)
pbmc.query = AddMetaData(pbmc.query, metadata = predictions)

#check if cell type annotation makes sense
pbmc.query = RunPCA(pbmc.query, features = pbmc.features)
pbmc.query = RunUMAP(pbmc.query, dims = 1:30)
pdf('UMAP_w_Predicted_cell_types_10k_unsorted.pdf')
DimPlot(pbmc.query, reduction = "umap", group.by = 'predicted.id', label = TRUE)
dev.off()

#write meta_data to file
write.csv(predictions, file = '../data/CellType_annotation_pbmc_unsorted_10k.csv', row.names = TRUE)