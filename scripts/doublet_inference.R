library(ArchR)

arrow_3k_unsorted <- "/home/james/Documents/leuven/second-year/IBP/ArchR/arrow_files/3k_unsorted.arrow"
arrow_3k_sorted <- "/home/james/Documents/leuven/second-year/IBP/ArchR/arrow_files/3k_sorted.arrow"
arrow_10k_unsorted <- "/home/james/Documents/leuven/second-year/IBP/ArchR/arrow_files/10k_unsorted.arrow"
arrow_10k_sorted <- "/home/james/Documents/leuven/second-year/IBP/ArchR/arrow_files/10k_sorted.arrow"

ArrowFiles <- c(arrow_3k_unsorted, arrow_3k_sorted, arrow_10k_unsorted, arrow_10k_sorted)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
