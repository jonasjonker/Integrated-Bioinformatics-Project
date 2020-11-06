library(ArchR)

ArrowFiles <- c("../data/3k_unsorted.arrow",
                "../data/3k_sorted.arrow",
                "../data/10k_unsorted.arrow", 
                "../data/10k_sorted.arrow")

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
