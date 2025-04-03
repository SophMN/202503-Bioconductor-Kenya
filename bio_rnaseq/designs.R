#Load packages
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

#Load the data
se <- readRDS("output/GSE96870_se.rds")
se

#Differential expression analysis
##Day4 vs Day0
#Create a DESeqDataSet in prep for DEA
dds <- DESeqDataSet(se, design = ~sex + time)

#Run DESeq to estimate size factors, dispersions and calculate test stats with Wald Test
dds <- DESeq(dds)

#Get results for Day4 vs Day0
results_d4 <- results(dds, contrast = c("time", "Day4", "Day0"))

#Day4 vs Day0 - Female
#Create a new DESeq data set that explores interaction btwn sex and time
dds2 <- DESeqDataSet(se, design = ~ sex * time)

#Create a new column in our SE object that takes into account this interaction
se$sex_time <- paste(se$sex, se$time, sep = "_")
unique(se$sex_time)

#Create a new DESeq data set with a new experimental design taking sex and time into account 
dds2 <- DESeqDataSet(se, design = ~sex_time)

#Perform DEA
dds2 <- DESeq(dds2)

#Extract the results for sex_time focusing on female mice
results_d4_female <- results(dds2, contrast = c("sex_time", "Female_Day4", "Female_Day0"))
head(results_d4_female)

#Goal: create a vector of upregulated DEGs
up <- results_d4_female |>
  as.data.frame() |>
  filter(log2FoldChange > 1,
         !is.na(padj),
         padj < 0.05)
#Get the number of DEGs
nrow(up)

#Get the IDs of the up-regulated DEGs
up <- rownames(up)
up

#Overrepresentation analysis
rowData(se)

#Extract the ENTREZ IDs from the SE object
up_entrez <- rowData(se) |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  filter(gene %in% up) |>
  pull(ENTREZID)
up_entrez
head(up_entrez)
nrow(up_entrez)

#Perform functional annotation analysis using the GO database
go_day4 <- enrichGO(gene = up_entrez,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP")
#Convert to a data frame to view the results
go_day4
as.data.frame(go_day4) |>
  View()

#Visualise the results with the enrichplot package
barplot(go_day4, showCategory = 52)
dotplot(go_day4)










