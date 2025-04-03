#Load packages
library(SummarizedExperiment)
library(DESeq2)
library(vsn)
library(ggplot2)
library(ComplexHeatmap)
library(hexbin)
library(iSEE)
library(ExploreModelMatrix)
library(apeglm)
library(cowplot)

#Load SE
se <- readRDS("output/GSE96870_se.rds")
se

#Load the normalised DESeqDataSet
dds <- readRDS("output/GSE96870_dds.rds")
dds

#Consider the variation in gene counts (dispersions in our dataset)
dds <- estimateDispersions(dds)
dds
colData(dds)
rowData(dds)

#Visualise the dispersion estimates
plotDispEsts(dds)

#Testing
dds <- nbinomWaldTest(dds)
rowData(dds)

#Explore the results for specific contrasts in our experimental design (sex + time)
resTime <- results(dds, contrast = c("time", "Day8", "Day0"))
summary(resTime)

#Filter our data by p-value <0.05
resTime2 <- results(dds, contrast = c("time", "Day8", "Day0"), alpha = 0.05)
summary(resTime2)

#Filter our differential expression output by log2 fold change and adjusted p-values
resTime2
fold_filtered <- subset(resTime2, abs(log2FoldChange) > 1 & padj < 0.05)
summary(fold_filtered)

#Visualise the output
plotMA(resTime2)
plotMA(fold_filtered)

#Reduce or shrink the noise from the genes with low mean count
fold_filtered_fc <- lfcShrink(dds, coef = "time_Day8_vs_Day0", res = resTime2)
plotMA(fold_filtered_fc)

#Load the variance-stabilised data set
vsd
head(as.data.frame(resTime2))
class(resTime2)
head(as.data.frame(rowData(dds)))

#Combine and export results
dds_fc <- cbind(as.data.frame(rowData(dds)), as.data.frame(resTime2))
write.csv(dds_fc, file = "output/GSE96870_dds_fc.csv")


