#https://github.com/compgenome365/TCGA-Assembler-2
install.packages(c("HGNChelper", "RCurl", "httr", "stringr", "digest", "bitops"), dependencies=T)
getwd()
setwd("C:/Users/shiva/Desktop/cancer-databases/TCGA-Assembler")
source("./Module_A.R")##to download data
source("./Module_B.R")##for data processing

##########Tutorial
vPatientID <- c("TCGA-A7-A13F", "TCGA-AO-A12B", "TCGA-AR-A1AP", "TCGA-
AR-A1AQ", "TCGA-AR-A1AS", "TCGA-AR-A1AV", "TCGA-AR-A1AW", "TCGA-BH-
A0BZ", "TCGA-BH-A0DD", "TCGA-BH-A0DG")
#Downloading of data (RNA Seq and CNA)#we use Module A functions to download data of seven different platforms for the samples and save the data files. There are four input arguments taken by every data downloading
#function, including cancerType, assayPlatform, inputPatientIDs, and saveFolderName (specifying the directory)
path_geneExp <- DownloadRNASeqData(cancerType = "BRCA", assayPlatform
                                   = "gene.normalized_RNAseq", inputPatientIDs = vPatientID,
                                   saveFolderName = ".")
path_copyNumber <- DownloadCNAData(cancerType = "BRCA", assayPlatform
                                   = "cna_cnv.hg19", inputPatientIDs = vPatientID, saveFolderName = ".")

path_miRNAExp <- DownloadmiRNASeqData(cancerType = "BRCA", assayPlatform = "mir_HiSeq.hg19.mirbase20", inputPatientIDs = vPatientID, saveFolderName = ".")

path_somaticMutation <- DownloadSomaticMutationData(cancerType = "BRCA", assayPlatform = "somaticMutation_DNAseq", inputPatientIDs = vPatientID, saveFolderName = ".")

#Data processing (RNA Seq and CNA)

#These functions perform quality control on data (such as validate/correct the gene symbols and draw box plot of data for identifying outliers), remove redundant genomic feature descriptions, and separate different types of measurements into their own data tables (such as separating the raw read counts and the reads per million miRNA mapped (RPM) of miRNA-seq data into two tables). 
#These functions transform data into the matrix format and ensure that each row in the matrix corresponds to a unique genomic feature and each column corresponds to a sample.
list_geneExp <- ProcessRNASeqData(inputFilePath = path_geneExp[1],
                                  outputFileName = paste("BRCA", "geneExp", sep = "__"), dataType =
                                    "geneExp", outputFileFolder = "./ManualExampleData/")
list_copyNumber <- ProcessCNAData(inputFilePath = path_copyNumber[1],
                                  outputFileName = paste("BRCA", "copyNumber", sep = "__"),
                                  refGenomeFile = "./SupportingFiles/Hg19GenePosition.txt",
                                  outputFileFolder = "./ManualExampleData/")


#Advanced Data processing (RNA Seq and CNA)

#First,we organize the data of each platform into a list object that includes three elements, i.e. Des ,Data , and dataType . Des is the description of genomic features. Data is the data. dataType indicates the data platform
#CombineMultiPlatformData: This command identifies the samples that are covered by all seven platforms and merges their data
l_geneExp <- list(Des = list_geneExp$Des, Data = list_geneExp$Data,
                  dataType = "geneExp")
l_copyNumber <- list(Des = list_copyNumber$Des, Data =
                       list_copyNumber$Data, dataType = "copyNumber")
inputDataList <- vector("list", 2)
inputDataList[[1]] <- l_copyNumber
inputDataList[[2]] <- l_geneExp
list_CombinedData <- CombineMultiPlatformData(inputDataList =
                                                inputDataList)
write.table(cbind(list_CombinedData$Des, list_CombinedData$Data), file
            = paste(".", "CombinedMultiPlatformData.txt", sep = "/"), quote =
              FALSE, sep = "\t", na = "", col.names = TRUE, row.names = FALSE)

######Download BRCA RNA seq for all patients and normal
geneExp <- DownloadRNASeqData(cancerType = "BRCA", assayPlatform
                                   = "gene.normalized_RNAseq",
                                   saveFolderName = ".")

######Download BRCA RNA seq for all patients (TP)

geneExp2 <- DownloadRNASeqData(cancerType = "BRCA", assayPlatform
                              = "gene.normalized_RNAseq",tissueType = "TP",
                              saveFolderName = ".")
