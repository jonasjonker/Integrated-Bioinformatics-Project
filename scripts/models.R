# load libraries
library(mygene)
library(ggplot2)
library(splines)
library(dataPreparation)
library(biomaRt)
library(dplyr)
library(ggthemes)
library(scales)

# setting color scheme
stata_pal(scheme = "s2color")

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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++ Single Gene Models +++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# calculate correlation for a given gene
singleGeneCorrelation <- function(geneName) {
  data = getGeneData(geneName)
  cor(data[,1], data[,2])
}

# fit a linear regression to data from a single gene
singleGeneRegression <- function(geneName) {
  data = getGeneData(geneName)
  split_data = splitData(data, prop=0.8)
  train = split_data[[1]]
  test = split_data[[2]]
  
  model <- lm(GEX~GSM, data=train)
  preds <- predict(model, newdata=test)
  
  conf.int <-predict(model, interval="confidence")
  lmdata <- cbind(train, conf.int)
  
  p <- ggplot(lmdata, aes(GSM, GEX)) +
       ggtitle(paste(geneName, "Linear Regression")) +
       geom_point() +
       stat_smooth(method=lm) +
       geom_line(aes(y=lwr), color="red", linetype="dashed") +
       geom_line(aes(y=upr), color="red", linetype="dashed")
  print(p)
  
  eval_data = data.frame(expected=test$GEX, predicted=preds)
  correlation = cor(eval_data[,1], eval_data[,2])
  return(correlation)
}

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

# fit a linear regression to data from a single gene
ontologyRegression <- function(geneOntology, numGO=2) {
  data = getOntologyData(geneOntology, numGO)
  split_data = splitData(data, prop=0.8)
  train = split_data[[1]]
  test = split_data[[2]]
  
  model <- lm(GEX~GSM, data=train)
  preds <- predict(model, newdata=test)
  
  conf.int <-predict(model, interval="confidence")
  lmdata <- cbind(train, conf.int)
  
  p <- ggplot(lmdata, aes(GSM, GEX)) +
    ggtitle(paste(geneOntology, "Linear Regression")) +
    geom_point() +
    stat_smooth(method=lm) +
    geom_line(aes(y=lwr), color="red", linetype="dashed") +
    geom_line(aes(y=upr), color="red", linetype="dashed")
  print(p)
  
  eval_data = data.frame(expected=test$GEX, predicted=preds)
  correlation = cor(eval_data[,1], eval_data[,2])
  return(correlation)
}

ontologySpline <- function(geneOntology, cv=FALSE, df=5, plot=FALSE, numGO=2) {
  data = getOntologyData(geneOntology, numGO)
  split_data = splitData(data, prop=0.8)
  train = split_data[[1]]
  test = split_data[[2]]
  
  if (plot==TRUE) {
    plot(train$GSM, train$GEX, cex=.5, main=paste(geneOntology, " Smoothing Spline"))
  }
  if (cv==TRUE) {
    fit = smooth.spline(train$GSM, train$GEX, cv=TRUE)
  }
  else {
    fit = smooth.spline(train$GSM, train$GEX, df=df)
  }
  
  if (plot==TRUE) {
    lines(fit, col="red", lwd=2)
  }
  pred = predict(fit, test$GSM, se=T)
  eval_data = data.frame(predicted=pred[2], expected=test$GEX)
  correlation = cor(eval_data[,1], eval_data[,2])
  return(correlation)
}

#######################################
############ ANALYSIS #################
#######################################

###### SINGLE GENE CORRELATION ######
# testing a number of genes for correlation

# LYMPHOID
BCL11B_cor <- singleGeneCorrelation("BCL11B")
CD8A_cor <- singleGeneCorrelation("CD8A")
CCR7_cor <- singleGeneCorrelation("CCR7")
lymphoid_cor <- c(BCL11B_cor, CD8A_cor, CCR7_cor)

# MYELOID
LYN_cor <- singleGeneCorrelation("LYN")
SLC8A1_cor <- singleGeneCorrelation("SLC8A1")
PLXDC2_cor <- singleGeneCorrelation("PLXDC2")
LRMDA_cor <- singleGeneCorrelation("LRMDA")
NAMPT_cor <- singleGeneCorrelation("NAMPT")
myeloid_cor <- c(LYN_cor, SLC8A1_cor, PLXDC2_cor, LRMDA_cor, NAMPT_cor)

# BROAD CELL-TYPE
broad_cell_type_cor <- c(lymphoid_cor, myeloid_cor)

# Non-specific
NEAT1_cor <- singleGeneCorrelation("NEAT1")
ZEB2_cor <- singleGeneCorrelation("ZEB2")
IL12RB1_cor <- singleGeneCorrelation("IL12RB1")
PLCB1_cor <- singleGeneCorrelation("PLCB1")
SERINC5_cor <- singleGeneCorrelation("SERINC5")
non_specific_cor <- c(NEAT1_cor, ZEB2_cor, IL12RB1_cor, PLCB1_cor, SERINC5_cor)

# plotting specific vs non-specific
names_cor = c("broad cell-type", "non-specific")
barplot(c(mean(broad_cell_type_cor), mean(non_specific_cor)), names=names_cor, col='darkred', ylab="correlation", width=c(0.1,0.1), space=c(0.2,0.2))

###### SINGLE GENE REGRESSION ######

# LYMPHOID
BCL11B <- singleGeneRegression("BCL11B")
CD8A <- singleGeneRegression("CD8A")
CCR7 <- singleGeneRegression("CCR7")
lymphoid <- c(BCL11B, CD8A, CCR7)

# MYELOID
LYN <- singleGeneRegression("LYN")
SLC8A1 <- singleGeneRegression("SLC8A1")
PLXDC2 <- singleGeneRegression("PLXDC2")
LRMDA <- singleGeneRegression("LRMDA")
NAMPT <- singleGeneRegression("NAMPT")
myeloid <- c(LYN, SLC8A1, PLXDC2, LRMDA, NAMPT)

# BROAD CELL-TYPE
broad_cell_type <- c(lymphoid, myeloid)

# B-Cell
# BANK1 <- singleGeneRegression("BANK1")
# BCL11B <- singleGeneRegression("BCL11B")
# BLK <- singleGeneRegression("BLK")
# CD79A <- singleGeneRegression("CD79A")
# bcell <- c(BANK1, BCL11B, BLK, CD79A)
# 
# # DENDRITIC
# PLD4 <- singleGeneRegression("PLD4")
# GZMB <- singleGeneRegression("GZMB")
# dendritic <- c(PLD4, GZMB)
# 
# # GRANULOCYTES
# CCR3 <- singleGeneRegression("CCR3")
# MS4A3 <- singleGeneRegression("MS4A3")
# FCGR3B <- singleGeneRegression("FCGR3B")
# granulocytes <- c(CCR3, MS4A3, FCGR3B)
# 
# # NATURAL KILLER
# KLRF1 <- singleGeneRegression("KLRF1")
# XCL2 <- singleGeneRegression("XCL2")
# natural_killer <- c(KLRF1, XCL2)
# 
# # T-CELLS
# CCR7 <- singleGeneRegression("CCR7")
# CD8B <- singleGeneRegression("CD8B")
# CCL5 <- singleGeneRegression("CCL5")
# t_cells <- c(CCR7, CD8B, CCL5)

# Non-specific
NEAT1 <- singleGeneRegression("NEAT1")
ZEB2 <- singleGeneRegression("ZEB2")
IL12RB1 <- singleGeneRegression("IL12RB1")
PLCB1 <- singleGeneRegression("PLCB1")
SERINC5 <- singleGeneRegression("SERINC5")
non_specific <- c(NEAT1, ZEB2, IL12RB1, PLCB1, SERINC5)

# plotting specific vs non-specific
names = c("broad cell-type", "non-specific")
barplot(c(mean(broad_cell_type), mean(non_specific)), names=names, col='darkred', ylab="correlation", width=c(0.1,0.1), space=c(0.2,0.2))

# lymphoid, myeloid, broad cell-type, and non-specific
names = c("lymphoid", "myeloid", "broad cell-type", "non-specific")
barplot(c(mean(lymphoid), mean(myeloid), mean(broad_cell_type), mean(non_specific)), names=names)

###### SINGLE GENE SPLINE ######
# LYMPHOID
BCL11B_spline <- singleGeneSpline("BCL11B")
CD8A_spline <- singleGeneSpline("CD8A")
CCR7_spline <- singleGeneSpline("CCR7")
lymphoid_spline <- c(BCL11B_spline, CD8A_spline, CCR7_spline)

# MYELOID
LYN_spline <- singleGeneSpline("LYN")
SLC8A1_spline <- singleGeneSpline("SLC8A1")
PLXDC2_spline <- singleGeneSpline("PLXDC2")
LRMDA_spline <- singleGeneSpline("LRMDA")
NAMPT_spline <- singleGeneSpline("NAMPT")
myeloid_spline <- c(LYN_spline, SLC8A1_spline, PLXDC2_spline, LRMDA_spline, NAMPT_spline)

# BROAD CELL-TYPE
broad_cell_type_spline <- c(lymphoid_spline, myeloid_spline)

# # B-Cell
# BANK1_spline <- singleGeneSpline("BANK1")
# BCL11B_spline <- singleGeneSpline("BCL11B")
# BLK_spline <- singleGeneSpline("BLK")
# CD79A_spline <- singleGeneSpline("CD79A")
# bcell_spline <- c(BANK1_spline, BCL11B_spline, BLK_spline, CD79A_spline)
# 
# # DENDRITIC
# # PLD4_spline <- singleGeneSpline("PLD4")
# GZMB_spline <- singleGeneSpline("GZMB")
# dendritic_spline <- c(GZMB_spline)
# 
# # GRANULOCYTES
# CCR3_spline <- singleGeneSpline("CCR3")
# MS4A3_spline <- singleGeneSpline("MS4A3")
# FCGR3B_spline <- singleGeneSpline("FCGR3B")
# granulocytes_spline <- c(CCR3_spline, MS4A3_spline, FCGR3B_spline)
# 
# # NATURAL KILLER
# KLRF1_spline <- singleGeneSpline("KLRF1")
# XCL2_spline <- singleGeneSpline("XCL2")
# natural_killer_spline <- c(KLRF1_spline, XCL2_spline)
# 
# # T-CELLS
# CCR7_spline <- singleGeneSpline("CCR7")
# CD8B_spline <- singleGeneSpline("CD8B")
# CCL5_spline <- singleGeneSpline("CCL5")
# t_cells_spline <- c(CCR7_spline, CD8B_spline, CCL5_spline)

# Non-specific
NEAT1_spline <- singleGeneSpline("NEAT1")
ZEB2_spline <- singleGeneSpline("ZEB2")
IL12RB1_spline <- singleGeneSpline("IL12RB1")
PLCB1_spline <- singleGeneSpline("PLCB1")
SERINC5_spline <- singleGeneSpline("SERINC5")
non_specific_spline <- c(NEAT1_spline, ZEB2_spline, IL12RB1_spline, PLCB1_spline, SERINC5_spline)

# plotting specific vs non-specific
names = c("broad cell-type", "non-specific")
barplot(c(mean(broad_cell_type_spline), mean(non_specific_spline)), names=names, col='darkred', ylab="correlation")

# lymphoid, myeloid, broad cell-type, and non-specific
names = c("lymphoid", "myeloid", "broad cell-type", "non-specific")
barplot(c(mean(lymphoid_spline), mean(myeloid_spline), mean(broad_cell_type_spline), mean(non_specific_spline)), names=names)

###### SINGLE GENE MODEL COMPARISONS ######
# plotting broad cell vs non-specific for correlation, linear, splines
names = c("broad cell-type", "non-specific","broad cell-type_spline", "non-specific_spline")
barplot(c(mean(broad_cell_type), mean(non_specific), mean(broad_cell_type_spline), mean(non_specific_spline)), names=names, col='darkred', ylab="correlation")

######## GENE ONTOLOGY #######
# running some sample gene ontology models to show poor performance
angiogenesis = ontologyRegression("angiogenesis", 2)
angiogenesis_spline = ontologySpline("angiogenesis", 2)

mcd = ontologyRegression("myeloid cell differentiation", 3)
mcd_spline = ontologySpline("cell development", 3)
 