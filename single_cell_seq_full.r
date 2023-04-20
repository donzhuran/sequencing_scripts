# working directory ----
#Setting the working directory
setwd("/Users/siddhantkalra/Desktop/Decode_Workshop/Sc_RNAseq")
#pdf(file="Plots.pdf")
# Packages ----
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("AnnotationDbi")
library("AnnotationDbi") 

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Hs.eg.db")

library("org.Hs.eg.db")

#install.packages("dplyr")
library("dplyr")

# install.packages("tidyr")
library(tidyr)

#install.packages("data.table")
library(data.table)

#install.packages("Seurat")
library(Seurat)


'''
### HANDLING THE DATA ###############################################################
'''

# data was taken from the experiment data available online


# Importing data ----
f <- read.csv("EXP0001_PCG_beforeQC.txt", sep="\t", row.names = 1)
View(f[1:10,1:10])

# removing the first row
f <- f[-1,]
View(f[1:10,1:10])

# Annotating the genes
f$Gene <- mapIds(org.Hs.eg.db, # create new column ("Gene"), take the library homo sapiens,
                   keys=row.names(f), # take row names from file as "keys", 
                   column="SYMBOL",  # and map them to human genome, column "symbol"
                   keytype="ENSEMBL", # return me mapped but transformed symbols as ensembl symbols
                   multiVals="first") # if there are more values, give me the first

# Setting Gene column as first
f <- f %>% select(Gene, everything())

View(f[1:10,1:10])
#View(f$Gene)

# removing all rows with distinct gene names, and which are not NA
f <- distinct(f,Gene, .keep_all= TRUE) # keeping uniques
f <- f %>% drop_na(Gene) # removing NAs

# convert this column as row ids
row.names(f)<-f$Gene
f$Gene <- NULL # remove the first column


pb <- CreateSeuratObject(counts = f, 
                         min.cells = 3, # keep only rows where at least in 3 columns the value is >0
                         min.features = 200) # keep only those columns that have values >0 in at least 200 genes
pb

# quality control - expression of mitochondrial genes
pb[["percent.mt"]] <- PercentageFeatureSet(pb, pattern = "^MT-") 
pb[["percent.mt"]] # 5-10% is ok

# Visualize QC metrics as a violin plot
pdf(file="Plot1.pdf") # create pdf file in the background
VlnPlot(pb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off() # close the pdf file
# they should be as high as possible, except percent.mt, which should be as small as possible

# Normalization
pb <- NormalizeData(pb)

'''
### ANALYZING GENE EXPRESSIONS #################################################
'''

# Variable features (genes that have different (variable) expressions among cells)
pb <- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)

list_of_variable_features <- VariableFeatures(pb)
list_of_variable_features <- as.data.frame(list_of_variable_features)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pb), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pb)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# scaling the data
all.genes <- rownames(pb)
pb <- ScaleData(pb, features = all.genes)

# PCA analysis
pb <- RunPCA(pb, features = VariableFeatures(object = pb))
VizDimLoadings(pb, dims = 1:2, reduction = "pca")
DimPlot(pb, reduction = "pca")

# heatmap for first 15 principal components
DimHeatmap(pb, dims = 1:15, cells = 500, balanced = TRUE)

# compute p-values for all the genes in specific principal components
pb <- JackStraw(pb, num.replicate = 100)
pb <- ScoreJackStraw(pb, dims = 1:20)
JackStrawPlot(pb, dims = 1:15) # data must be above the dashed line
ElbowPlot(pb)

# clustering genes
pb <- FindNeighbors(pb, dims = 1:10) # clustering based on gene expression
pb <- FindClusters(pb, resolution = 0.5) 
head(Idents(pb), 5)
pb <- RunUMAP(pb, dims = 1:10)
DimPlot(pb, reduction = "umap")

# finding genes (markers) for differentiating those clusters
pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pb.markers,file="Markers_info.csv",
            sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

# interested only in some genes
VlnPlot(pb, features = c("FTL","RAP1A","ISG15"))
VlnPlot(pb, features = c("ISG15"))

# which cells have specific expression (on UMAP plot)
FeaturePlot(pb, features = "STATH")
FeaturePlot(pb, features = "ISG15")
FeaturePlot(pb, features = "ISG15", pt.size = 2)

# which cells belong to which cluster
check <- Idents(pb)
check <- as.data.frame(check)

check<-setDT(check, keep.rownames = TRUE)[] # take "check" dataset, make the row name as the first column (Cell id)
colnames(check) <- c("Cell_Id", "Cluster")

cluster2 <- subset(check, check$Cluster==2)

# filtering gene lists to custom gene list
merge_data <- merge(f1,f2,by.x="Gene", by.y="Gene")

# metadata
meta_data <- read.csv("SraRunTable.txt", sep=",")
meta_data <- read.csv("filereport_read_run_PRJNA400576_tsv.txt", sep=",")
meta_data <- meta_data[,c(1,5,21,25,27)]
