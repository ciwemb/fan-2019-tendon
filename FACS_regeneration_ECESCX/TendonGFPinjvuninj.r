#set working directory
setwd("/Volumes/sequence/harvey/git/TendonGFP_injvuninj_DESeq")

#Load metadata table and store it
TendonGFPs_metadata<-read.table("TendonGFPs.metadata.tsv", header=TRUE)
TendonGFPs_metadata

#metadata data summary
summary(TendonGFPs_metadata)

#Load DESeq library
library("DESeq")

#Define the count table??

#Load in the read count data
cds<-newCountDataSetFromHTSeqCount(TendonGFPs_metadata, directory="/Volumes/sequence/harvey/git/TendonGFP_injvuninj_DESeq/")
cds

#Normalizing Data
cds<-estimateSizeFactors(cds)
sizeFactors(cds)

#Estimating Dispersion; parametricDispersion failed so had to do local fitType
cds<-estimateDispersions(cds, fitType = c("local"), sharingMode= c("gene-est-only"))
plotDispEsts(cds)
head(fData(cds))
write.csv(counts(cds, normalized=TRUE), file="/Volumes/sequence/harvey/git/TendonGFP_injvuninj_DESeq/normalized_cds_counts.csv")

#Adjust Variance for PCA
vsd<-varianceStabilizingTransformation(cds)
vsd

#Make PCA plot with DESeq
plotPCA(vsd)

#Calculating Differential Expression; write out file as .csv
GFPinj_vs_GFP <- nbinomTest(cds, "GFPinj", "GFP" )
GFPinj_vs_GFP

write.csv(GFPinj_vs_GFP, file = "/Volumes/sequence/harvey/git/TendonGFP_injvuninj_DESeq/GFPinj_vs_GFP2.csv")

#GOI list for plot subsetting
Canonical <- c("Bgn","Ctgf","Dcn","Fmod","Fn1","Tnc","Tnmd","Egr1","Scx","Mkx")
Canonical
Tenocyte10X <- c("Abi3bp","Cilp2","Clec11a","Clu","CrispId2","Ctsk","Dkk3","Ecm2","Fibin","Kera","Leprel1","Lum","Pdgfrl","Prelp","Prg4","Serpinf1","Sfrp2","Sparc","Thbs1","Thbs2","Thbs4")
Tenocyte10X
Sheath10X <- c("Has1","Ly6a","Npdc1","Phlda3","Plin2","S100a10")
Sheath10X
Collagen <- c("Col1a1","Col1a2","Col3a1","Col5a1","Col5a2","Col6a1","Col6a2","Col6a3","Col8a1","Col12a1")
Collagen
#Plot Differential Expression MA plots
plotMA(GFPinj_vs_GFP[order(GFPinj_vs_GFP$padj),])

#To plot through ggplot an MA plot with color coded significant genes and GOIs
library( dplyr )
library( ggplot2)
library( stringr)
GFPinj_vs_GFP$M <-   log2( GFPinj_vs_GFP$baseMeanB + 1 ) - log2( GFPinj_vs_GFP$baseMeanA + 1 )
GFPinj_vs_GFP$A <- ( log2( GFPinj_vs_GFP$baseMeanB + 1 ) + log2( GFPinj_vs_GFP$baseMeanA + 1 ) ) / 2


ggplot( GFPinj_vs_GFP ) +
  geom_point( aes( A, M ), col="grey")+
  labs( x = "Mean log2 normalized counts" , y ="log2FoldChange")+
  geom_point( data=filter( GFPinj_vs_GFP, padj < 0.05 ), aes( A, M ), col="skyblue" )+
  geom_point( data=filter( GFPinj_vs_GFP,  id %in% Canonical ), aes( A, M ), shape=21, col="red" )+
  geom_point( data=filter( GFPinj_vs_GFP, id %in% Tenocyte10X ), aes( A, M ), shape=21, col="black")+
  geom_point( data=filter( GFPinj_vs_GFP, id %in% Sheath10X ), aes( A, M ), shape=21, col="brown")+
  geom_point( data=filter( GFPinj_vs_GFP, id %in% Collagen ), aes( A, M ), shape=21, col="purple")



#Plot Differential Expression Volcano plots
plot(GFPinj_vs_GFP$log2FoldChange, -log10(GFPinj_vs_GFP$padj))

#To plot through ggplot a DE Volcano plot
GFPinj_vs_GFP$log2FC <-   log2( GFPinj_vs_GFP$baseMeanB + 1 ) - log2( GFPinj_vs_GFP$baseMeanA + 1 )
GFPinj_vs_GFP$log10padj <- ( -log10( GFPinj_vs_GFP$padj)  )

ggplot( GFPinj_vs_GFP ) +
  geom_point( aes( log2FC, log10padj ), col="grey" )+
  labs( x = "log2FoldChange" , y = "-log10padj")+
  geom_point( data=filter( GFPinj_vs_GFP, padj < 0.05 ), aes( log2FC, log10padj ), col="skyblue" )+
  geom_point( data=filter( GFPinj_vs_GFP,  id %in% Canonical ), aes( log2FC, log10padj ), shape=21, col="red" )+
  geom_point( data=filter( GFPinj_vs_GFP, id %in% Tenocyte10X ), aes( log2FC, log10padj ), shape=21, col="black")+
  geom_point( data=filter( GFPinj_vs_GFP, id %in% Sheath10X ), aes( log2FC, log10padj ), shape=21, col="brown")+
  geom_point( data=filter( GFPinj_vs_GFP, id %in% Collagen ), aes( log2FC, log10padj ), shape=21, col="purple")


