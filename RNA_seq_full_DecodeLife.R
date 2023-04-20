#---------Bulk RNAseq Analysis

#Set up the working directory----
setwd("Desktop/Decode_Workshop/Bulk_RNAseq/")

#Loading required packages----
#install.packages("R.utils")
library(R.utils)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsubread")
library(Rsubread)

#install.packages("data.table")
library(data.table)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("RUVSeq")
library(RUVSeq)

#(if the above installation of RUVseq didn't work, try this)
#source("http://bioconductor.org/biocLite.R")
#biocLite("RUVSeq")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
library(DESeq2)

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

#install.packages("pheatmap")
library(pheatmap)

#install.packages("RColorBrewer")
library(RColorBrewer)

#install.packages("ggplot2")
library(ggplot2)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rqc")
# library(Rqc)

#Downloading FastQ files----
url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_1.fastq.gz"
destination<-"SRR5924196_1.fastq.gz"
download.file(url,destination)

url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_2.fastq.gz"
destination<-"SRR5924196_2.fastq.gz"
download.file(url,destination)

url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_1.fastq.gz"
destination<-"SRR5924198_1.fastq.gz"
download.file(url,destination)

url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_2.fastq.gz"
destination<-"SRR5924198_2.fastq.gz"
download.file(url,destination)


#Downloading Genome file----
url<-"ftp://ftp.ensembl.org/pub/release-96/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
destination<-"Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
download.file(url,destination)
gunzip(destination)

#Downloading GTF file---- (which part of the gene corresponds to which gene); file must be from the same source as Genome file (Ensembl)
url<-"ftp://ftp.ensembl.org/pub/release-96/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
destination<-"Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
download.file(url,destination)
gunzip(destination)


"""
CHECKING QUALITY OF FASTQ FILES -------------------------------------
- 
"""

#install.packages("fastqcr")
#Your system should also have JAVA installed.
#visit www.java.com for installation
library(fastqcr)
fastqc_install()
fastqc()
# save all results in a single folder
qc <- qc_aggregate("FASTQC/")
qc
view(qc)


"""
MAPPING TO REFERENCE GENOME -------------------------------------
- 
"""

#Building Index---- (hash table from Genome - for alignment)
buildindex("Sc_full_index_rsubread",
           "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
           indexSplit=F)
# choosing the read from both directions (pair-end reads)
reads1 <- list.files(pattern = "_1.fastq.gz$" )
reads2 <- list.files(pattern = "_2.fastq.gz$" )
# checking if we have all the reads from both ends
all.equal(length(reads1),length(reads2))


#Performing Alignment----
align(index="Sc_full_index_rsubread",
      readfile1=reads1,
      readfile2=reads2,
      input_format="gzFASTQ",
      output_format="BAM",
      nthreads=10)

#Checking the BAM files generated----
bam.files <- list.files(pattern = ".BAM$", full.names = TRUE)
bam.files # one BAM file for both FASTQ files

#Checking the mapping quality----
props <- propmapped(files=bam.files)
props # 0.8/0.85 or more is okay


"""
QUANTIFICATION -------------------------------------
- 
"""

#Generating Feature counts----
fcLim <- featureCounts(files = bam.files,
                       annot.ext="Saccharomyces_cerevisiae.R64-1-1.96.gtf",
                       GTF.featureType="exon",
                       GTF.attrType="gene_id",
                       isGTFAnnotationFile=TRUE,
                       isPairedEnd = TRUE)

fc <- data.frame(fcLim[["counts"]])
colnames(fc) <- c("Normal", "Tumor") # adding feature characteristics


# since this is a simple data, we will use another dataset to explain the rest
#Working on actual data----
f <- read.csv("GSE143630_RCC_htseq_counts.txt",
              sep =" ",
              row.names = 1)
View(f[1:10,1:10])

#T1 stage count----
count_of_T1=length(grep(x = colnames(f),
                            pattern = "^T1.")) # count how many patients are T1
count_of_T1

#T2 stage count----
count_of_T2=length(grep(x = colnames(f),
                           pattern = "^T2.")) # count how many patients are T2
count_of_T2


"""
GENE EXPRESSION ANALYSIS -------------------------------------

"""

#Filtering only genes that have values higher than 0 in both rows
filter <- apply(f, 1, function(x) length(x[x>0])>=2)
filtered <- f[filter,]   
genes <- rownames(filtered) 

# creating tables explaining what is what
t1<-rep(c("T1"),each=count_of_T1)
t2<-rep(c("T2"),each=count_of_T2)
x<-c(t1,t2)
x<-as.factor(x)
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x,
                                                  row.names=colnames(filtered)))
set

#Setting Color theme
colors <- brewer.pal(3, "Set2")

#Plotting basic graphs
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

differences <- makeGroups(x)
set3 <- RUVs(set, genes, k=1, differences)
pData(set3)

#Computing Differential expression
dds <- DESeqDataSetFromMatrix(countData = counts(set3),
                              colData = pData(set3),
                              design= ~ W_1 + x)

dds <- DESeq(dds)

res <- results(dds)
head(res)

write.table(res, file = "RUVseq_analysis.csv",
            sep = ",", col.names = NA,
            qmethod = "double")

df <- read.csv("RUVseq_analysis.csv")
df <- subset(df, df$padj <= 0.05) # taking only statistically significant fold changes

#Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')


#Heatmap
f1 <- f[1:4,] # taking only 4 genes
pheatmap(f1, cluster_cols = FALSE)


#Boxplot
f <- t(f)
f <- as.data.frame(f)
View(f[1:10,1:10])

Tumor_stage <- c(rep(c("T1"),each=count_of_T1),
                 rep(c("T2"),each=count_of_T2))
f <- cbind(Tumor_stage,f)
View(f[,1:10])

p<-ggplot(f, aes(x=Tumor_stage, y=A1CF, fill=Tumor_stage)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
p

#Violin Plot
p<-ggplot(f, aes(x=Tumor_stage, y=A1CF, fill=Tumor_stage)) +
  geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2))
p
