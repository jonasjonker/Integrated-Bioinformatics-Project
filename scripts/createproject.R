# load ArchR
library(ArchR)

# set number of threads (usually set to half the number of available cores)
addArchRThreads(threads = 1)

# add genome
addArchRGenome("hg38")

# assign paths to arrow files
arrow_3k_unsorted <- "/home/james/Documents/leuven/second-year/IBP/ArchR/3k_unsorted.arrow"
arrow_3k_sorted <- "/home/james/Documents/leuven/second-year/IBP/ArchR/3k_sorted.arrow"
arrow_10k_unsorted <- "/home/james/Documents/leuven/second-year/IBP/ArchR/10k_unsorted.arrow"
arrow_10k_sorted <- "/home/james/Documents/leuven/second-year/IBP/ArchR/10k_sorted.arrow"
ArrowFiles <- c(arrow_3k_unsorted, arrow_3k_sorted, arrow_10k_unsorted, arrow_10k_sorted)

# make an ArchR project
proj_pbmc <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ProjPBMC",
  copyArrows = TRUE # This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# create a backup of the Archr Project before filtering doublets
# saveArchRProject(ArchRProj = proj_pbmc, outputDirectory = "Save-ProjPBMC", load = FALSE)

# Defining a fixed filterDoublets function for R v4.0.3
filterDoublets <- function(ArchRProj = NULL, cutEnrich = 1, cutScore = -Inf, filterRatio = 1){
  
  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for (i in seq_along(fn)) {
    tryCatch({
      eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
  }
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = cutEnrich, name = "cutEnrich", valid = c("numeric"))
  .validInput(input = cutScore, name = "cutScore", valid = c("numeric"))
  .validInput(input = filterRatio, name = "filterRatio", valid = c("numeric"))
  
  if(any(grepl("filterDoublets", names(ArchRProj@projectSummary)))){
    stop("Already ran filterDoublets on ArchRProject! Cannot be re-ran on an ArchRProject!")
  }
  
  df <- getCellColData(ArchRProj, c("Sample", "DoubletEnrichment", "DoubletScore"))
  splitDF <- split(seq_len(nrow(df)), as.character(df$Sample))
  
  cellsFilter <- lapply(splitDF, function(y){
    
    x <- df[y, ,drop = FALSE]
    
    n <- nrow(x)
    
    x <- x[order(x$DoubletEnrichment, decreasing = TRUE), ]
    
    if(!is.null(cutEnrich)){
      x <- x[which(x$DoubletEnrichment >= cutEnrich), ]
    } 
    
    if(!is.null(cutScore)){
      x <- x[which(x$DoubletScore >= cutScore), ]
    } 
    
    if(nrow(x) > 0){
      head(rownames(x), filterRatio * n * (n / 100000))
    }else{
      NULL
    }
    
  }) %>% unlist(use.names=FALSE)
  
  message("Filtering ", length(cellsFilter), " cells from ArchRProject!")
  tabRemove <- table(df[cellsFilter,]$Sample)
  tabAll <- table(df$Sample)
  samples <- unique(df$Sample)
  for(i in seq_along(samples)){
    if(!is.na(tabRemove[samples[i]])){
      message("\t", samples[i], " : ", tabRemove[samples[i]], " of ", tabAll[samples[i]], " (", round(100 * tabRemove[samples[i]] / tabAll[samples[i]], 1),"%)")
    }else{
      message("\t", samples[i], " : ", 0, " of ", tabAll[samples[i]], " (0%)")
    }
  }
  
  if(length(cellsFilter) > 0){
    
    ArchRProj@cellColData <- ArchRProj@cellColData[rownames(ArchRProj@cellColData) %ni% cellsFilter,,drop=FALSE]
    
  }
  
  ArchRProj
  
}

# filtering doublets in the proj_pbmc project
proj_pbmc_filtered <- filterDoublets(proj_pbmc)
