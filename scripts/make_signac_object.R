#!/usr/bin/env Rscript
'
Usage:
    create_signac_obj.R <dir10x> <fragments> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.

Arguments:
    10xdir     directory containing: barcodes.tsv.gz features.tsv.gz matrix.mtx.gz
    fragments  fragment file
    output     output file
' -> doc

library(docopt)
message("Parse command line arguments.")
args <- docopt(doc, version = 'make_signac_object 0.1')
if (!dir.exists(args$dir10x)) { stop(args$dir10x, " does not exists") }
if (!file(args$fragments)) { stop(args$fragments, " does not exists") }

message("Load libraries.")
library(Signac)
library(Seurat)

message("Read 10x data.")
raw_feature <- Seurat::Read10X(args$dir10x)

message("Read peak file.")
fragments <- args$fragments  

message("create seurat object from RNAseq data.")
snare <- Seurat::CreateSeuratObject(counts = raw_feature$`Gene Expression`)

message("Add ATAC data to seurat object.")
snare[['ATAC']] <- Signac::CreateChromatinAssay(
  counts = raw_feature$Peaks,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = fragments
)

message("Save to .rds file.")
saveRDS(snare, args$output)

message("Done.")
