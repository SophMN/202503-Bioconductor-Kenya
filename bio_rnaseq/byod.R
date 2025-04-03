#Downloading data directly from GEO NCBI database
BiocManager::install("GEOquery")
BiocManager::install("macrophage")
library(macrophage)
data("gse")
library(GEOquery)
library(SummarizedExperiment)
library(DESeq2)
library(vsn)
library(ggplot2)
library(apeglm)
library(tidyverse)
library(airway)
library(cowplot)
library(ExploreModelMatrix)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

options(timeout = 1e6)
geo_data <- getGEO("GSE234712", GSEMatrix = FALSE, destdir = "output")
class(geo_data)
rm(geo_data)
geo_data <- getGEO("GSE195521", GSEMatrix = TRUE, destdir = "output")

#Examine the class of the geo data object
class(geo_data)
geo_data[[1]]
class(geo_data[[1]])

#Since geo_data is of the expression set class, convert an expression set into a SE object
se <- makeSummarizedExperimentFromExpressionSet(geo_data[[1]])
se
rm(se)

#Examine metadata
colData(se) |> View()
dim(colData(se))

#Check for the unique levels within the variables in your sample metadata
unique(se$tissue.ch1)
unique(se$treatment.ch1)

#Look at the 
table(se$treatment.ch1)

#Examine rowData
dim(rowData(se))

#Examine the counts data 
assay(se)
head(assay(se))
dim(assay(se))

##Bring your own data day: Install the airway package
#Load the data
se <- data("airway")
class(airway)
se <- airway
rm(airway)
se

#Access the assay, sample and features metadata
assay(se)
dim(assay(se))
colData(se)
dim(colData(se))
rowData(se)
dim(rowData(se))

unique(se$dex)
table(se$dex)

#Exploratory data analysis and quality control
#Remove unexpressed genes
se <- se[rowSums(assay(se, "counts")) > 5, ]
se$dex <- factor(se$dex, levels = c("untrt", "trt"))


#Calculate the library sizes in preparation for DESeq2 analysis
se$libSize <- colSums(assay(se))
se
colData(se)

#Create a DESeqDataSet
dds <- DESeqDataSet(se, design = ~0 + dex) #Remove intercept
dds <- DESeq(dds)
resultsNames(dds)

#Estimate size factors
dds <- estimateSizeFactors(dds)
dds

#Variance stabilisation
#Relationship between mean and variance of the counts
meanSdPlot(assay(dds), ranks = FALSE)

#Perform variance stabilisation transformation
vsd <- vst(dds, blind = TRUE)
meanSdPlot(assay(vsd), ranks = FALSE)

#Quality control with principal component analysis
pcaData <- plotPCA(vsd, intgroup = "dex", returnData = TRUE)
View(pcaData)

#Plot the variance in the experimental design
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData |>
  ggplot(aes(x = PC1, y = PC2, color = dex)) +
  geom_point() +
  xlab(paste0(x = "PC1:", percentVar[1], "%variance")) +
  ylab(paste0(y = "PC2:", percentVar[2], "%variance"))

#Differential expression analysis
#Estimate gene-wise dispersions
dds <- estimateDispersions(dds)

#Visualise the dispersion estimates
plotDispEsts(dds)
dds
rowData(dds)

#Statistical testing and modelling
dds <- nbinomWaldTest(dds)
dds
rowData(dds)
results <- results(dds)
resultsNames(dds)
summary(results)

#Extract the results for the dex treatment
res
res_dex <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
summary(res_dex)
View(res_dex)

BiocManager::install("ashr")
library(ashr)
#Shrink the LFC of genes with low mean normalised counts
res_dexLFC <- lfcShrink(dds, contrast = c("dex", "trt", "untrt"), 
                        res = res_dex, type = "ashr")
resultsNames(dds)

#Filter the results by LFC and adjusted p-values
fold_filtered <- subset(res_dex, abs(log2FoldChange) > 1 & padj < 0.05)
summary(fold_filtered)
fold_filtered
fold_filtered_top <- fold_filtered[order(fold_filtered$padj),]
fold_filtered_top

#Visualise the output
plotMA(res_dex)
plotMA(res_dexLFC)
plotMA(fold_filtered)

#Bind output results
head(as.data.frame(res_dex))
head(as.data.frame(rowData(dds)))
dds_dex <- cbind(as.data.frame(rowData(dds)), as.data.frame(res_dex))

#Create a vector of the upregulated genes
res_dex
up <- res_dex |>
  as.data.frame() |>
  filter(log2FoldChange > 1, !is.na(padj), padj < 0.05)
nrow(up)  
View(up)

#Filter the gene names symbols from the SE
rowData(se)
up_symbols <- rowData(se) |>
  as.data.frame()
rm(up_symbols)
  rownames_to_column("gene") |>
  filter(gene_name %in% up) 
  
#Import macrophage package summarized experiment object
data("gse")

#Examine the assay, sample and feature metadata
assay(gse)
dim(assay(gse))
colData(gse)
colnames(colData(gse))
dim(colData(gse))
rowData(gse)
colnames(rowData(gse))
dim(rowData(gse))

##Exploratory data analysis and quality control
#Remove the unexpresed genes
gse <- gse[rowSums(assay(gse, "counts") > 5, ]



