library(chromPlot)
library(tidyverse)
library(ggpubr)
library(rstatix)

setwd("C:\\Users\\HP\\Desktop\\Assignment01")
dataFile=read.delim("GRCh38_hg38_variants_2020-02-25.txt",header = T,sep = "\t")

 
# Genomic elements
GenomicDataFrame=data.frame(Chrom=dataFile$chr,Start=dataFile$start,End=dataFile$end)

GenomicDataFrame$Colors <- "red"

vec=c(1:10,"X","Y")

chromPlot(bands=GenomicDataFrame, gaps=hg_gap, chr=c(8),figCols = 8)




#Histogram
GenomicDataFrame$Colors <- NULL
chromPlot(gaps=hg_gap, bands=GenomicDataFrame, annot1=GenomicDataFrame, chr=c(8), figCols=6)

#Histograms for different sub variants-------------------------------------------------------------------
#01.Gains

library(dplyr)
# Subset the dataframe using logical indexing
Gains <- dataFile[dataFile$variantsubtype == "gain", ]
GainsDataFrame=data.frame(Chrom=Gains$chr,Start=Gains$start,End=Gains$end)

GainsDataFrame$Colors <- "red"
chromPlot(bands=GainsDataFrame, gaps=hg_gap, chr=c(19:21), figCols=3)

GainsDataFrame$Colors <- "black"
chromPlot(gaps=hg_gap, bands=GainsDataFrame, annot1=GainsDataFrame, chr=c(1:3,19:21), figCols=6)


#02.Losses
# Subset the dataframe using logical indexing
Losses <- dataFile[dataFile$variantsubtype == "loss", ]
LossDataFrame=data.frame(Chrom=Losses$chr,Start=Losses$start,End=Losses$end)

LossDataFrame$Colors <- "red"
chromPlot(bands=LossDataFrame, gaps=hg_gap, chr=c(1:3,19:21), figCols=3)

LossDataFrame$Colors <- "black"
chromPlot(gaps=hg_gap, bands=LossDataFrame, annot1=LossDataFrame, chr=c(1:3,19:21), figCols=6)

#03.Duplication
# Subset the dataframe using logical indexing
Duplication <- dataFile[dataFile$variantsubtype == "duplication", ]
DuplicationDataFrame=data.frame(Chrom=Duplication$chr,Start=Duplication$start,End=Duplication$end)

DuplicationDataFrame$Colors <- "red"
chromPlot(bands=DuplicationDataFrame, gaps=hg_gap, chr=c(19:21), figCols=3)

DuplicationDataFrame$Colors <- "black"
chromPlot(gaps=hg_gap, bands=DuplicationDataFrame, annot1=DuplicationDataFrame, chr=c(1:3,19:21), figCols=6)

#04.Deletion
# Subset the dataframe using logical indexing
Deletion <- dataFile[dataFile$variantsubtype == "deletion", ]
DeletionDataFrame=data.frame(Chrom=Deletion$chr,Start=Deletion$start,End=Deletion$end)

DeletionDataFrame$Colors <- "red"
chromPlot(bands=DeletionDataFrame, gaps=hg_gap, chr=c(19:21), figCols=3)

DeletionDataFrame$Colors <- "black"
chromPlot(gaps=hg_gap, bands=DeletionDataFrame, annot1=DeletionDataFrame, chr=c(1:3,19:21), figCols=6)

#05.Insertion
# Subset the dataframe using logical indexing
Insertion <- dataFile[dataFile$variantsubtype == "insertion", ]

InsertionDataFrame=data.frame(Chrom=Insertion$chr,Start=Insertion$start,End=Insertion$end)

InsertionDataFrame$Colors <- "red"
chromPlot(bands=InsertionDataFrame, gaps=hg_gap, chr=c(19:21), figCols=3)

InsertionDataFrame$Colors <- "black"
chromPlot(gaps=hg_gap, bands=InsertionDataFrame, annot1=InsertionDataFrame, chr=c(1:3,19:21), figCols=6)




# Large stacked segments
# chromPlot(gaps=hg_gap, segment=GenomicDataFrame, noHist=TRUE, annot1=hg_gap,chrSide=c(-1,1,1,1,1,1,1,1), chr=c(21), stack=TRUE, figCol=3,bands=GenomicDataFrame)          
   


#Trying to do an ANOVA
#Hypothesis

#Alternative:
#Lengths of Gains,Losses,Insertions,Deletions,Duplications are not equal to each other

#Null:
#Lengths of Gains,Losses,Insertions,Deletions,Duplications are equal to each other

#Getting targeted variables into one table
Insertion <- dataFile[dataFile$variantsubtype == "insertion", ]
Deletion<- dataFile[dataFile$variantsubtype == "deletion", ]
Duplication=dataFile[dataFile$variantsubtype == "duplication", ]
Gains=dataFile[dataFile$variantsubtype == "gain", ]
Losses=dataFile[dataFile$variantsubtype == "loss", ]

MainDataFrame=rbind(Insertion,Deletion,Duplication,Gains,Losses)

LenDataFrame=data.frame(subvar=MainDataFrame$variantsubtype,LengthOfCnv=MainDataFrame$end-MainDataFrame$start)
CorrectedLenDataFrame=LenDataFrame[LenDataFrame$LengthOfCnv!=0,]

# install.packages("tidyverse")
# install.packages("ggpubr")
# install.packages("rstatix")



#Summary Statistics
sumStat=CorrectedLenDataFrame%>%
  group_by(subvar)%>%
  get_summary_stats(LengthOfCnv,type="common")

#Visualization
ggboxplot(CorrectedLenDataFrame,x="subvar",y="LengthOfCnv",fill = "subvar")

res.Kruskal=CorrectedLenDataFrame%>%
  kruskal_test(LengthOfCnv~subvar)
res.Kruskal

res.dunn<- CorrectedLenDataFrame %>%
  dunn_test(LengthOfCnv~subvar, p.adjust.method = "bonferroni")
res.dunn






       