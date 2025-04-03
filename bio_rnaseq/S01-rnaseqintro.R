#Import data
download.file(url = "https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_counts_cerebellum.csv", destfile = "data/GSE96870_counts_cerebellum.csv")
download.file(url = "https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_coldata_cerebellum.csv", destfile = "data/GSE96870_coldata_cerebellum.csv")
download.file(url = "https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_coldata_all.csv", destfile = "data/GSE96870_coldata_all.csv")
download.file(url = "https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_rowranges.tsv", destfile = "data/GSE96870_rowranges.tsv")

#Load the required packages
library(AnnotationDbi)
library(org.Mm.eg.db)
library(hgu95av2.db)
library(SummarizedExperiment)
library(tidyverse)

#Load the files to R
#CSV
coldata <- read.csv("data/GSE96870_coldata_cerebellum.csv", row.names = 1)

#CSV
counts <- read.csv("data/GSE96870_counts_cerebellum.csv", row.names = 1)

#TSV 
rowranges <- read.delim("data/GSE96870_rowranges.tsv", sep = "\t",
                        colClasses = c(ENTREZID = "character"),
                        header = TRUE, quote = "", row.names = 5)

#Evaluate the dimensions of the variables
dim(counts)
dim(coldata)
dim(rowranges)
View(counts)
View(coldata)
View(rowranges)

#Create a SummarizedExperiment that includes all the data sets of interest 
?SummarizedExperiment()
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  rowRanges = as(rowranges, "GRanges"),
  colData = coldata
)

#View the SummarizedExperiment object 
se


head(assay(se))
colData(se)
head(rowData(se))

#Save the SummarizedExperiment object
saveRDS(se, "output/GSE96870_se.RDS")
rm(se)

#Loading your SummarisedExperiment
se <- readRDS("output/GSE96870_se.RDS")








