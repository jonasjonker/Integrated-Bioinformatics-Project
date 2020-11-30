# load libraries
library(mygene)
library(ggplot2)
library(splines)
library(dataPreparation)
library(biomaRt)
library(dplyr)

# loading in the data
PATH_TO_DATA <- readline(prompt="Enter path to data: ")
GSM_FILE <- readline(prompt="Enter file name of GSM file: ")
GEX_FILE <- readline(prompt="Enter file name of GEX file: ")

setwd(PATH_TO_DATA)
gsm_df = read.table(GSM_FILE,sep=" ")
gex_df = read.table(GEX_FILE,sep=" ")

#++++++++++++++++++++++++++++++++++++++++++
#++++++++++++ UTILITY FUNCTIONS +++++++++++
#++++++++++++++++++++++++++++++++++++++++++

# generate a data frame with gene score and expression data for a given gene
getGeneData <- function(geneName) {
  geneDF <- data.frame(matrix(nrow=ncol(gsm_df), ncol=2))
  colnames(geneDF) <- c("GSM", "GEX")
  geneDF$GSM <- t(gsm_df[geneName,])
  geneDF$GEX <- t(gex_df[geneName,])
  return(geneDF)
}

# split data into training and test sets
splitData <- function(data, prop=0.8) {
  train_index <- sample(1:nrow(data), prop * nrow(data))
  test_index <- setdiff(1:nrow(data), train_index)
  train <- data[train_index, ]
  test <- data[test_index, ]
  splitData = list(train, test)
}

# get gene length and add it to a matrix
getGeneLengths <- function(genelist) {
  ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  lengthData <- getBM(mart = ensembl_human,
                               attributes = c("external_gene_name","start_position","end_position"),
                               filters = "external_gene_name",
                               values = genelist)
  lengthData$length = lengthData$end_position - lengthData$start_position
  return(lengthData)
}

# plot a scatter plot for gene score vs expression for a given gene and fit a local regression on the data
plotSingleGene <- function(geneName) {
  geneDF = getGeneData(geneName)
  gg <- ggplot(geneDF, aes(x=GSM, y=GEX)) + 
    geom_point() + 
    geom_smooth(method="loess", se=T) + 
    labs(subtitle="GSM vs GEX", 
         y="GEX", 
         x="GSM", 
         title=geneName) 
  plot(gg)
}

# Get gene ontology terms for all genes
getOntologies <- function(go_data, numGO){
  ontologies <- data.frame(matrix(nrow=length(go_data$go.BP), ncol=numGO+1))
  colnames(ontologies) <- c("Gene", paste("GO", c(1:numGO), sep=""))
  ontologies[is.na(ontologies)] <- 0
  
  for (i in 1:length(go_data$go.BP)) {
    ontologies$Gene[i] <- go_data$query[i]
    if (is.null(go_data$go.BP[[i]]$term)==FALSE){
      for (j in 1:numGO){
        if (is.null(go_data$go.BP[[i]]$term[j]) == FALSE) {
          ontologies[i,j+1] <- go_data$go.BP[[i]]$term[j]
        }
        else{
          ontologies[i,j+1]=0
        }
      }
    }
    ontologies[is.na(ontologies)] <- 0
  }
  return(ontologies)
}

getGenesWithGO <- function(geneOntology, ontologies) {
  genes = ontologies[apply(ontologies, 1, function(x) any(grepl(geneOntology, x))), ][1]
  filtered_gsm = data.frame(GSM = c(t(gsm_df[genes$Gene,])), stringsAsFactors=FALSE)
  filtered_gex = data.frame(GEX = c(t(gex_df[genes$Gene,])), stringsAsFactors=FALSE)
  ontologyData = as.data.frame(cbind(filtered_gsm, filtered_gex))
  return(ontologyData)
}

# get gene data for a given ontology
getOntologyData <- function(geneOntology, numGOs) {
  genelist <- row.names(gsm_df)
  go_data <- queryMany(genelist, scopes='symbol', fields=c('go'), species='human')
  ontologies <- getOntologies(go_data, numGOs)
  results <-getGenesWithGO(geneOntology, ontologies)
  return(results)
}

# plot a scatter plot for gene score vs expression for a given gene ontology and fit a local regression on the data
plotGeneOntology <- function(geneOntology, results, numcells=5000) {
  gene_GO_DF = getGenesFromGO(geneOntology, results)
  gg <- ggplot(gene_GO_DF[1:numcells,], aes(x=GSM, y=GEX)) + 
    geom_point() + 
    geom_smooth(method="loess", se=T) + 
    labs(subtitle="GSM vs GEX", 
         y="GEX", 
         x="GSM", 
         title=geneOntology) 
  plot(gg)
}

# getting length data
gene_lengths = getGeneLengths(genelist)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++ Single Gene Models +++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fit a smoothing spline to data from a single gene
singleGeneSpline <- function(geneName, cv=FALSE, df=5) {
  data = getGeneData(geneName)
  split_data = splitData(data, prop=0.8)
  train = split_data[[1]]
  test = split_data[[2]]
  
  plot(train$GSM, train$GEX, cex=.5, main=paste(geneName, " Smoothing Spline"))
  if (cv==TRUE) {
    fit = smooth.spline(train$GSM, train$GEX, cv=TRUE)
  }
  else {
    fit = smooth.spline(train$GSM, train$GEX, df=df)
  }
  lines(fit, col="red", lwd=2)
  pred = predict(fit, test$GSM, se=T)
  eval_data = data.frame(expected=test$GEX, predicted=pred[2])
  correlation = cor(eval_data[,1], eval_data[,2])
  plot(eval_data)
  return(correlation)
}

# fitting using local regression
singleGeneLoess <- function(geneName, span=0.5) {
  data = getGeneData(geneName)
  gsmlims = range(data$GSM)
  gsm.grid=seq(from=gsmlims[1], to=gsmlims[2])
  
  split_data = splitData(data, prop=0.8)
  train = split_data[[1]]
  test = split_data[[2]]
  
  plot(train$GSM, train$GEX, cex=.5, main=paste(geneName, " Local Regression"))
  fit = loess(train$GEX~train$GSM, span=span, data=data)
  lines(fit, col="red", lwd=2)
  
  preds = predict(fit, test$GSM, se=T)
  eval_data = data.frame(expected=test$GEX, predicted=pred[2])
  correlation = cor(eval_data[,1], eval_data[,2])
  plot(eval_data)
  return(correlation)
}

#+++++++++++++++++++++++++++++++++++
#++++++++ ONTOLOGY MODELS ++++++++++
#+++++++++++++++++++++++++++++++++++

ontologySpline <- function(geneOntology, cv=FALSE, df=5) {
  data = getGeneData(geneOntolo)
  split_data = splitData(data, prop=0.8)
  train = split_data[[1]]
  test = split_data[[2]]
  
  plot(train$GSM, train$GEX, cex=.5, main=paste(geneName, " Smoothing Spline"))
  if (cv==TRUE) {
    fit = smooth.spline(train$GSM, train$GEX, cv=TRUE)
  }
  else {
    fit = smooth.spline(train$GSM, train$GEX, df=df)
  }
  lines(fit, col="red", lwd=2)
  pred = predict(fit, test$GSM, se=T)
  eval_data = data.frame(expected=test$GEX, predicted=pred[2])
  correlation = cor(eval_data[,1], eval_data[,2])
  plot(eval_data)
  return(correlation)
}