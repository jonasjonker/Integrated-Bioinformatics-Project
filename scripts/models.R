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
singleGeneRegression <- function(geneName, plot=FALSE) {
  data = getGeneData(geneName)
  split_data = splitData(data, prop=0.8)
  train = split_data[[1]]
  test = split_data[[2]]
  
  model <- lm(GEX~GSM, data=train)
  preds <- predict(model, newdata=test)
  
  conf.int <-predict(model, interval="confidence")
  lmdata <- cbind(train, conf.int)
  
  if (plot==TRUE) {
    p <- ggplot(lmdata, aes(GSM, GEX)) +
      ggtitle(paste(geneName, "Linear Regression")) +
      geom_point() +
      stat_smooth(method=lm) +
      geom_line(aes(y=lwr), color="red", linetype="dashed") +
      geom_line(aes(y=upr), color="red", linetype="dashed")
    print(p) 
  }
  
  eval_data = data.frame(expected=test$GEX, predicted=preds)
  correlation = cor(eval_data[,1], eval_data[,2])
  return(correlation)
}

# fit a smoothing spline to data from a single gene
singleGeneSpline <- function(geneName, cv=TRUE, df=5, plot=FALSE) {
  data = getGeneData(geneName)
  split_data = splitData(data, prop=0.8)
  train = split_data[[1]]
  test = split_data[[2]]
  
  if (plot==TRUE) {
    plot(train$GSM, train$GEX, cex=.5, main=paste(geneName, " Smoothing Spline")) 
  }
  
  if (cv==TRUE) {
    fit = smooth.spline(train$GSM, train$GEX, cv=TRUE)
  }  else {
    fit = smooth.spline(train$GSM, train$GEX, df=df)
  }
  
  if (plot==TRUE) {
    lines(fit, col="red", lwd=2)
  }

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
cor_bct <- c(lymphoid_cor, myeloid_cor)

# Non-specific
NEAT1_cor <- singleGeneCorrelation("NEAT1")
ZEB2_cor <- singleGeneCorrelation("ZEB2")
IL12RB1_cor <- singleGeneCorrelation("IL12RB1")
PLCB1_cor <- singleGeneCorrelation("PLCB1")
SERINC5_cor <- singleGeneCorrelation("SERINC5")
cor_ns <- c(NEAT1_cor, ZEB2_cor, IL12RB1_cor, PLCB1_cor, SERINC5_cor)

###### SINGLE GENE REGRESSION ######

# initialise empty vectors to store results
reg_bct <- c()
reg_ns <- c()

for (i in 1:10) {
  # Broad cell-type
  BCL11B <- singleGeneRegression("BCL11B")
  CD8A <- singleGeneRegression("CD8A")
  CCR7 <- singleGeneRegression("CCR7")
  LYN <- singleGeneRegression("LYN")
  SLC8A1 <- singleGeneRegression("SLC8A1")
  PLXDC2 <- singleGeneRegression("PLXDC2")
  LRMDA <- singleGeneRegression("LRMDA")
  NAMPT <- singleGeneRegression("NAMPT")
  bct <- mean(c(BCL11B, CD8A, CCR7, LYN, SLC8A1, PLXDC2, LRMDA, NAMPT))
  reg_bct <- c(reg_bct, bct)
  
  # Non-specific
  NEAT1 <- singleGeneRegression("NEAT1")
  ZEB2 <- singleGeneRegression("ZEB2")
  IL12RB1 <- singleGeneRegression("IL12RB1")
  PLCB1 <- singleGeneRegression("PLCB1")
  SERINC5 <- singleGeneRegression("SERINC5")
  ns <- mean(c(NEAT1, ZEB2, IL12RB1, PLCB1, SERINC5))
  reg_ns <- c(reg_ns, ns) 
}

###### SINGLE GENE SPLINE ######

# initialise empty vectors to store results
spl_bct <- c()
spl_ns <- c()

for (i in 1:10) {

  # LYMPHOID
  BCL11B_spline <- singleGeneSpline("BCL11B")
  CD8A_spline <- singleGeneSpline("CD8A")
  CCR7_spline <- singleGeneSpline("CCR7")
  LYN_spline <- singleGeneSpline("LYN")
  SLC8A1_spline <- singleGeneSpline("SLC8A1")
  PLXDC2_spline <- singleGeneSpline("PLXDC2")
  LRMDA_spline <- singleGeneSpline("LRMDA")
  NAMPT_spline <- singleGeneSpline("NAMPT")
  bct <- mean(c(BCL11B_spline, CD8A_spline, CCR7_spline, LYN_spline, SLC8A1_spline, PLXDC2_spline, LRMDA_spline, NAMPT_spline))
  spl_bct <- c(spl_bct, bct)
 
  # Non-specific
  NEAT1_spline <- singleGeneSpline("NEAT1")
  ZEB2_spline <- singleGeneSpline("ZEB2")
  IL12RB1_spline <- singleGeneSpline("IL12RB1")
  PLCB1_spline <- singleGeneSpline("PLCB1")
  SERINC5_spline <- singleGeneSpline("SERINC5")
  ns <- mean(c(NEAT1_spline, ZEB2_spline, IL12RB1_spline, PLCB1_spline, SERINC5_spline))
  spl_ns <- c(spl_ns, ns)
}

###### SINGLE GENE MODEL COMPARISONS ######
barplot_df <- data.frame(model=c("No model", "No model", "Linear", "Linear", "Spline", "Spline"),
                         genes=c("Broad cell-type", "Non-specific", "Broad cell-type", "Non-specific", "Broad cell-type", "Non-specific"),
                         corr=round(c(mean(cor_bct), mean(cor_ns), mean(reg_bct), mean(reg_ns), mean(spl_bct), mean(spl_ns)),digit=2),
                         sd=c(sd(cor_bct), sd(cor_ns), sd(reg_bct), sd(reg_ns), sd(spl_bct), sd(spl_ns)))

ggplot(data=barplot_df, aes(x=model, y=corr, fill=genes)) +
  geom_bar(stat="identity", position=position_dodge()) +
  # geom_errorbar(aes(ymin=corr-sd, ymax=corr+sd), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=corr), vjust=-0.3, position=position_dodge(0.9), size=3.5) +
  theme_classic() +
  theme(legend.position="top", legend.title=element_blank()) +
  scale_fill_calc() +
  xlab(element_blank()) +
  ylab("Pearson Correlation")
  

# plotting broad cell vs non-specific for correlation, linear, splines
names = c("Broad", "NS", "broad linear", "NS linear","broad spline", "NS spline")
barplot(c(mean(cor_bct), mean(cor_ns), mean(reg_bct), mean(reg_ns), mean(spl_bct), mean(spl_ns)), names=names, col='darkred', ylab="correlation")


######## GENE ONTOLOGY #######
# running some sample gene ontology models to show poor performance
angiogenesis = ontologyRegression("angiogenesis", 2)
angiogenesis_spline = ontologySpline("angiogenesis", 2)

mcd = ontologyRegression("myeloid cell differentiation", 3)
mcd_spline = ontologySpline("cell development", 3)
 