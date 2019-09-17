#set working directory
setwd("/Volumes/sequence/harvey/git/3populationbulkRNAseq/Dual_GFPScx")
getwd()

#Load metadata table and store it
DualGFPuninj_metadata<-read.table("DualGFPuninj.metadata.tsv", header=TRUE)
DualGFPuninj_metadata

#metadata data summary
summary(DualGFPuninj_metadata)

#Load DESeq library
library("DESeq")

#Define the count table??

#Load in the read count data
cds<-newCountDataSetFromHTSeqCount(DualGFPuninj_metadata, directory="/Volumes/sequence/harvey/git/3populationbulkRNAseq/Dual_GFPScx")
cds

#Normalizing Data
cds<-estimateSizeFactors(cds)
sizeFactors(cds)

#Estimating Dispersion; parametricDispersion failed so had to do local fitType
cds<-estimateDispersions(cds, fitType = c("parametric"), sharingMode= c("gene-est-only"))
plotDispEsts(cds)
head(fData(cds))


#Adjust Variance for PCA
vsd<-varianceStabilizingTransformation(cds)
vsd

#Make PCA plot with DESeq
plotPCA(vsd)


#Heatmap of vsd top # genes (specify in subset) using normalized counts, parametric fit
cdsFullBlind = estimateDispersions( cds, fitType = c("parametric"), method = "blind" )
vsdFull = varianceStabilizingTransformation( cdsFullBlind )
library("RColorBrewer")
library("gplots")
select500 = order(rowMeans(exprs(vsdFull)), decreasing=TRUE)[1:500]
select1K = order(rowMeans(exprs(vsdFull)), decreasing=TRUE)[1:1000]
hmcol = colorRampPalette(c("blue","yellow"))(12)
z<-heatmap.2(exprs(vsdFull)[select500,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
a<-heatmap.2(exprs(vsdFull)[select1K,], scale="row", col = hmcol, trace="none", margin=c(10, 6))

#create dataframe of counts data for calling out gois for subsequent heatmap generations
df <- data.frame( counts( cds ) )
df$id <- 1:nrow(df)
df
#goi specified lists
goi_names <- c( "Tppp3","Pdgfra","Pdgfa","Pdgfrl","Scx","Mkx","Egr1","Bgn","Ctgf","Dcn", "Fmod", "Fn1", "Tnc", "Tnmd", "Col1a1","Col1a2","Col5a1","Col6a1","Col12a1")
goi_heatmapDualUp <- c("Thbs1","Thbs2","Thbs4","Col5a1","Fmod","Fn1","Col1a1","Col1a2","Col12a1","Bgn","Ctgf")
goi_heatmapDualUpunmatched<- c("Col5a2","Col5a3","Col6a1","Col6a2","Col6a3","Fbn1","Col5a1","Col5a3","Col14a1","Col3a1","Dcn")
goi_notinheatmap<- c("Scx","Mkx","Egr1","Pdgfrl","Pdgfa")
goi_heatmap<-c("Thbs1","Thbs2","Thbs4","Col5a1","Fmod","Fn1","Col1a1","Col1a2","Col12a1","Bgn","Ctgf","Scx","Mkx","Egr1","Pdgfrl","Pdgfa")
#grab rows of interest (roi) with specified gois
roi <- df[ rownames( df ) %in% goi_names, "id" ]
roi
roi2 <- df[ rownames( df ) %in% goi_heatmapDualUp, "id" ]
roi2
roi3 <- df[ rownames( df ) %in% goi_heatmapDualUpunmatched, "id" ]
roi3
roi4 <- df[ rownames( df ) %in% goi_notinheatmap, "id" ]
roi4
roi5 <- df[ rownames( df ) %in% goi_heatmap, "id" ]
roi5
sapply( goi_names, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_heatmapDualUp, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_heatmapDualUpunmatched, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_notinheatmap, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_heatmap, function(x) grep( x, rownames( counts( cds ) ) ) )
heatmap.2(exprs(vsdFull)[roi,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi2,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi3,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi4,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi5,], scale="row", col = hmcol, trace="none", margin=c(10, 6))













