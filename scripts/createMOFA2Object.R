'Create MOFA2 Object

Usage:
  createMOFA2Object.R <factors> <outfile> 

Options:
  -h --help     Show this screen.
  --version     Show version.
' -> doc
library(docopt)
args <- docopt(doc, version = '1.0')

seurat <- readRDS(paste0(args$outfile, "_seurat.rds"))

library(MOFA2)

mofa <- create_mofa(seurat, assays = c("RNA",
                                       "ATAC_distal",
                                       "ATAC_promoter")) # this will throw an error its okay

data_opt   <- get_default_data_options(mofa)
model_opts <- get_default_model_options(mofa)

model_opts$num_factors=as.integer(args$factors)

mofa <- prepare_mofa(mofa, model_options = model_opts)

trainedmofa <- run_mofa(mofa, outfile = paste0(args$outfile, "_trained_mofa.rds"))