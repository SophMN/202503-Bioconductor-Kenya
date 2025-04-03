#Load packages
library(SummarizedExperiment)
library(DESeq2)
library(vsn)
library(ggplot2)
library(ComplexHeatmap)
library(hexbin)
library(iSEE)

#Load the expression data
se <- readRDS("output/GSE96870_se.rds")
se

#Extracting data from the SE object
rowData(se)
colData(se)
assay(se)

#Subsetting SE objects based on columns (colData)
se_female <- se[ , se$sex == "Female"]
se_female

se_female_infected_t8 <- se[ , se$sex == "Female" & 
                               se$infection == "InfluenzaA" &
                               se$time == "Day8"]
se_female_infected_t8

#Subsetting SE objects based on rows (rowData)
se2 <- se[rowData(se)$gbkey == "ncRNA", ]
se2
idx <- rowData(se)$gbkey == "ncRNA"
idx
se[idx, ]

#Subset the se object to keep only samples from male mice at time point 4
colData(se)
se
se_male <- se[ , se$sex == "Male" & se$time == "Day4"]
se_male
ncol(se_male)

#Filtering genes out with counts of less than 5
##Example 1: housekeeping genes (high mean expression; low variance)
##Example 2: unexpressed genes
se <- se[rowSums(assay(se, "counts")) > 5, ]
se
nrow(se)

#Check for library size differences
se$libSize <- colSums(assay(se))
se
colData(se)

#Visualise the difference in library size by sex
colData(se) |>
  as.data.frame() |> #convert the SummarizedExperiment object into a data frame for visualisation
  ggplot(aes(x = geo_accession, y = libSize / 1e6, fill = sex)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Total counts in millions") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#Save SummarizedExperiment 
saveRDS(se, "output/GSE96870_se.rds")
rm(se)

#Estimate size factors by creating a DESeq2 object
dds <- DESeqDataSet(se, design = ~ sex + time)
dds
dds <- estimateSizeFactors(dds)
dds
colData(dds)

#Variance-stabilising transformation
meanSdPlot(assay(dds), ranks = FALSE)
vsd <- vst(dds, blind = TRUE) #Transformation doesn't take the experimental design into account 
vsd
meanSdPlot(assay(vsd), ranks = FALSE)

#Principal component analysis (PCA)
plotPCA(vsd, intgroup = c("sex", "time"))

#Convert the variance-stabilised data into a data frame to visualise the difference in sex and time
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("sex", "time"), returnData = TRUE)

#Convert PCA data into 

#Save outputs
saveRDS(se, "output/GSE96870_se.rds")
rm(se)
saveRDS(dds, "output/GSE96870_dds.rds")
rm(dds)
saveRDS(vsd, "output/GSE96870_vsd.rds")
rm(vsd)
write.csv(pcaData, "output/GSE96870_pca.csv")
rm(pcaData)
