'Find Doublets.

Usage:
  doublet_inference.R <data_dir>

Options:
  -h --help     Show this screen.
  --version     Show version.
' -> doc
library(docopt)
args <- docopt(doc, version = 'doublet_inference 1.0')

# preconditions
if (!dir.exists(args$data_dir)) stop(args$data_dir, " doesn't exist.")

# load libraries
library(ArchR)

# file locations
ArrowFiles <- file.path(args$data_dirlist.files(file.path(args$data_dir, "ArrowFiles")))

# find doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
