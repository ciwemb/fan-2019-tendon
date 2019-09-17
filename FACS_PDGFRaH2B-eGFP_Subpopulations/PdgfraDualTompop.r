#set working directory
setwd("/Volumes/sequence/harvey/git/3populationbulkRNAseq")

#Load metadata table and store it
TendonDualTom_metadata<-read.table("TendonDualTom.metadata.tsv", header=TRUE)
TendonDualTom_metadata

#metadata data summary
summary(TendonDualTom_metadata)

#Load DESeq library
library("DESeq")

#Define the count table??

#Load in the read count data
cds<-newCountDataSetFromHTSeqCount(TendonDualTom_metadata, directory="/Volumes/sequence/harvey/git/3populationbulkRNAseq")
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
##Make matrix of top # of genes and writeout as a .csv file
write.csv((exprs(vsdFull)[select500,])[rev(z$rowInd), z$colInd], file = "/Volumes/sequence/harvey/git/3populationbulkRNAseq/Heatmaps_parametric_DualTom/top500.csv")
write.csv((exprs(vsdFull)[select1K,])[rev(a$rowInd), a$colInd], file = "/Volumes/sequence/harvey/git/3populationbulkRNAseq/Heatmaps_parametric_DualTom/top1K.csv")

dists = dist( t( exprs(vsdFull) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(condition, libType, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))



#create dataframe of counts data for calling out gois for subsequent heatmap generations
df <- data.frame( counts( cds ) )
df$id <- 1:nrow(df)
df
#goi specified lists
goi_names <- c( "Tppp3","Pdgfra","Pdgfa","Pdgfrl","Scx","Mkx","Egr1","Bgn","Ctgf","Dcn", "Fmod", "Fn1", "Tnc", "Tnmd", "Col1a1","Col1a2","Col5a1","Col6a1","Col12a1")
goi_names_CanonClust4 <- c("Tppp3","Pdgfra","Pdgfa","Pdgfrl","Scx","Mkx","Egr1","Bgn","Ctgf","Dcn", "Fmod", "Fn1", "Tnc", "Tnmd", "Col1a1","Col1a2","Col5a1","Col6a1","Col12a1","Abibp","Adamts2","Aebp1","Angptl7","Apod","Ccdc3","Ccdc80","Cfh","Chad","Chrdl1","Cilp2","Clec11a","Clec3b","Clu","Col11a1","Col12a1","Col1a1","Col1a2","Comp","CrispId2","Cst3","Ctsk","CtsI","Dcn","Ecm2","Ergic3","Fbln1","Fibin","Fkbp7","Fkbp9","Fmod","Fndc1","Fxyd6","Gas1","Gm13305","Gm42418","Fpm6b","Gpx3","Hes1","Hpgd","Ifngr1","Ifg1","Igfbp6","Il11ra1","Itgbl1","Itm2c","Kera","Klf9","Laptm4a","Lhfp","Lum","Map1lc3a","Meg3","Mfge8","Mgp","Mt2","Mt3","Mxra8","Myoc","Ndufa4l2","Nucb2","Ogn","Olfml3","P4ha1","Pcolce","Pdgfrl","Pik3r1","Plpp1","Plxdc2","Prelp","Prg4","Ptgis","Rarres2","Rcn2","Rcn3","Rora","Selm","Serpinf1","Serping1","Sgk1","Sod3","Sparc","Thbs1","Thbs2","Thbs4","Timp2","Tnmd","Tob1","Tsc22d3","Wisp2")
goi_names_cluster4 <- c("Abibp","Adamts2","Aebp1","Angptl7","Apod","Ccdc3","Ccdc80","Cfh","Chad","Chrdl1","Cilp2","Clec11a","Clec3b","Clu","Col11a1","Col12a1","Col1a1","Col1a2","Comp","CrispId2","Cst3","Ctsk","CtsI","Dcn","Ecm2","Ergic3","Fbln1","Fibin","Fkbp7","Fkbp9","Fmod","Fndc1","Fxyd6","Gas1","Gm13305","Gm42418","Fpm6b","Gpx3","Hes1","Hpgd","Ifngr1","Ifg1","Igfbp6","Il11ra1","Itgbl1","Itm2c","Kera","Klf9","Laptm4a","Lhfp","Lum","Map1lc3a","Meg3","Mfge8","Mgp","Mt2","Mt3","Mxra8","Myoc","Ndufa4l2","Nucb2","Ogn","Olfml3","P4ha1","Pcolce","Pdgfrl","Pik3r1","Plpp1","Plxdc2","Prelp","Prg4","Ptgis","Rarres2","Rcn2","Rcn3","Rora","Selm","Serpinf1","Serping1","Sgk1","Sod3","Sparc","Thbs1","Thbs2","Thbs4","Timp2","Tnmd","Tob1","Tsc22d3","Wisp2")
goi_names_cluster2 <- c("150015010Rik","4930523C07Rik","Adk","Ahnak2","Anxa1","Anxa2","Anxa4","Anxa5","Anxa8","Aqp1","Arl6ip5","Aspn","Atp6ap2","Avpi1","Bgn","Capg","Ccnd1","Cd9","Cdh13","Cilp","Cilp2","Cltb","Col11a1","Col15a1","Col6a1","Col6a2","Col6a3","Crip1","Crip2","Cryab","Cyb5r3","Dap","Dpysl3","Dstn","Ecm1","Fermt2","Fhl1","Fn1","Fxyd5","Fxyd6","Glud1","Gm9780","Htra1","Igfbp6","Igfbp7","Lgals3","Lox","Mgll","Mgst3","Msn","Mt3","Myadm","Nbl1","Ndufa4l2","Npdc1","Nrep","Nupr1","Pcolce2","Pdlim4","Phlda3","Pkm","Plec","Plod2","Plp2","Plpp1","Pls3","Prelp","Prkcdbp","Prnp","Prr13","Prrx1","Ptgis","Ptms","Rgcc","Rnh1","Rora","S100a10","S100a4","Serpine2","Slurp1","Sod3","Spp1","Timp1","Tm4sf1","Tmbim1","Tnfrsf12a","Tob1","Tpm4","Tppp3","Trim47","Tspan3","Tspo","Tubb2a","Tubb4b","Txn","Uqcc2","Vim","Vimp")
goi_names_cluster24 <- c("Abibp","Adamts2","Aebp1","Angptl7","Apod","Ccdc3","Ccdc80","Cfh","Chad","Chrdl1","Cilp2","Clec11a","Clec3b","Clu","Col11a1","Col12a1","Col1a1","Col1a2","Comp","CrispId2","Cst3","Ctsk","CtsI","Dcn","Ecm2","Ergic3","Fbln1","Fibin","Fkbp7","Fkbp9","Fmod","Fndc1","Fxyd6","Gas1","Gm13305","Gm42418","Fpm6b","Gpx3","Hes1","Hpgd","Ifngr1","Ifg1","Igfbp6","Il11ra1","Itgbl1","Itm2c","Kera","Klf9","Laptm4a","Lhfp","Lum","Map1lc3a","Meg3","Mfge8","Mgp","Mt2","Mt3","Mxra8","Myoc","Ndufa4l2","Nucb2","Ogn","Olfml3","P4ha1","Pcolce","Pdgfrl","Pik3r1","Plpp1","Plxdc2","Prelp","Prg4","Ptgis","Rarres2","Rcn2","Rcn3","Rora","Selm","Serpinf1","Serping1","Sgk1","Sod3","Sparc","Thbs1","Thbs2","Thbs4","Timp2","Tnmd","Tob1","Tsc22d3","Wisp2","150015010Rik","4930523C07Rik","Adk","Ahnak2","Anxa1","Anxa2","Anxa4","Anxa5","Anxa8","Aqp1","Arl6ip5","Aspn","Atp6ap2","Avpi1","Bgn","Capg","Ccnd1","Cd9","Cdh13","Cilp","Cilp2","Cltb","Col11a1","Col15a1","Col6a1","Col6a2","Col6a3","Crip1","Crip2","Cryab","Cyb5r3","Dap","Dpysl3","Dstn","Ecm1","Fermt2","Fhl1","Fn1","Fxyd5","Fxyd6","Glud1","Gm9780","Htra1","Igfbp6","Igfbp7","Lgals3","Lox","Mgll","Mgst3","Msn","Mt3","Myadm","Nbl1","Ndufa4l2","Npdc1","Nrep","Nupr1","Pcolce2","Pdlim4","Phlda3","Pkm","Plec","Plod2","Plp2","Plpp1","Pls3","Prelp","Prkcdbp","Prnp","Prr13","Prrx1","Ptgis","Ptms","Rgcc","Rnh1","Rora","S100a10","S100a4","Serpine2","Slurp1","Sod3","Spp1","Timp1","Tm4sf1","Tmbim1","Tnfrsf12a","Tob1","Tpm4","Tppp3","Trim47","Tspan3","Tspo","Tubb2a","Tubb4b","Txn","Uqcc2","Vim","Vimp")                        
goi_names_cluster3 <- c("Ace","Ackr3","Adamts5","Aebp1","AI607873","Anxa3","Basp1","Bicc1","C1ra","C1s1","C3","C4b","Ccl11","Ccl7","Cd248","Cd34","Celf2","Col14a1","Col3a1","Col4a1","Col4a2","Col5a1","Col5a2","Col5a3","Col6a3","Cxcl1","Cxcl12","Cxcl13","Cxcl14","Cygb","Ddr2","Dkk2","Dpep1","Dpt","Ebf1","Efemp1","Egr1","Enpp2","Entpd2","Errfi1","Fbln2","Fbn1","Fgl2","Fosb","Fst","Fstl1","Gdf10","Gfpt2","Gpc3","Gstm1","Hk2","Htra3","Ifi204","Ifi205","Igfbp4","Irf1","Itih5","Itm2a","Jun","Jund","Lamc1","Loxl1","Lpl","Ly6a","Marcks","Mfap5","Mmp2","Mmp3","Ndn","Nfib","Nid1","Nrp1","Pdgfra","Phlda1","Pi16","Plac8","Plpp3","Procr","Prss23","Ptgs2","Ptx3","Ramp2","Rnase4","S100a16","Scara5","Serping1","Sfrp2","Smoc2","Smpd3","Srpx","Steap4","Sult5a1","Thbd","Tmsb10","Tnfaip2","Tnfaip6","Ugdh","Wisp2","Zfp36","Zfp36l1")
#goi merged list of canonical and FAP markers
goi_merged<- c(goi_names,goi_top500FAP)
goi_merged


#grab rows of interest (roi) with specified gois
roi <- df[ rownames( df ) %in% goi_names, "id" ]
roi
roi2 <- df[ rownames( df ) %in% goi_names_cluster4, "id" ]
roi2

roi3 <- df[ rownames( df ) %in% goi_names_cluster2, "id" ]
roi3

roi4 <- df[ rownames( df ) %in% goi_names_cluster24, "id" ]
roi4

roi5 <- df[ rownames( df ) %in% goi_names_CanonClust4, "id" ]
roi5

roi6 <- df[ rownames( df ) %in% goi_names_cluster3, "id" ]
roi6

roi13 <- df[ rownames( df ) %in% goi_merged, "id" ]
roi13

sapply( goi_names, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_names_cluster4, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_names_cluster2, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_names_cluster24, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_names_CanonClust4, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_names_cluster3, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_top500FAP, function(x) grep( x, rownames( counts( cds ) ) ) )
sapply( goi_merged, function(x) grep( x, rownames( counts( cds ) ) ) )
x<-heatmap.2(exprs(vsdFull)[roi,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi2,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi3,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi4,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi5,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi6,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi7,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
heatmap.2(exprs(vsdFull)[roi13,], scale="row", col = hmcol, trace="none", margin=c(10, 6))

#Extract gene ids by heatmap row cluster pattern
(exprs(vsdFull)[roi,])[rev(x$rowInd), x$colInd]


#call in top500 FAP dataset
top500FAPs_df<- read.csv(file="/Volumes/sequence/harvey/git/3populationbulkRNAseq/GSE89633_DESeq2_ExonCounts_WTFAPs.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
top500FAPs_df
#subset to just grab top500 FAP signature
top500FAPs_df[,1]
goi_top500FAP<- c(top500FAPs_df[,1])
goi_top500FAP

roi7 <- df[ rownames( df ) %in% goi_top500FAP, "id" ]
roi7




# gois from state1 pseudotime trajectory
state1<- read.csv(file ="/Volumes/sequence/harvey/git/3populationbulkRNAseq/Pseudotime_states/State1.csv", 
                  header = FALSE, sep = ",", stringsAsFactors = FALSE) 
state1
state1[,1]
goi_state1<- c(state1[,1])
goi_state1[1:10]
roi8 <- df[ rownames( df ) %in% goi_state1[1:10], "id" ]
roi8
sapply( goi_state1, function(x) grep( x, rownames( counts( cds ) ) ) )
heatmap.2(exprs(vsdFull)[roi8,], scale="row", col = hmcol, trace="none", margin=c(10, 6))

# gois from state2 pseudotime trajectory
state2<- read.csv(file ="/Volumes/sequence/harvey/git/3populationbulkRNAseq/Pseudotime_states/State2.csv", 
                  header = FALSE, sep = ",", stringsAsFactors = FALSE) 
state2
state2[,1]
goi_state2<- c(state2[,1])
goi_state2
roi9 <- df[ rownames( df ) %in% goi_state2, "id" ]
roi9
sapply( goi_state2, function(x) grep( x, rownames( counts( cds ) ) ) )
heatmap.2(exprs(vsdFull)[roi9,], scale="row", col = hmcol, trace="none", margin=c(10, 6))

# gois from state3 pseudotime trajectory
state3<- read.csv(file ="/Volumes/sequence/harvey/git/3populationbulkRNAseq/Pseudotime_states/State3.csv", 
                  header = FALSE, sep = ",", stringsAsFactors = FALSE) 
state3
state3[,1]
goi_state3<- c(state3[,1])
goi_state3
roi10 <- df[ rownames( df ) %in% goi_state3, "id" ]
roi10
sapply( goi_state3, function(x) grep( x, rownames( counts( cds ) ) ) )
heatmap.2(exprs(vsdFull)[roi10,], scale="row", col = hmcol, trace="none", margin=c(10, 6))
# gois from state4 pseudotime trajectory
state4<- read.csv(file ="/Volumes/sequence/harvey/git/3populationbulkRNAseq/Pseudotime_states/State4.csv", 
                  header = FALSE, sep = ",", stringsAsFactors = FALSE) 
state4
state4[,1]
goi_state4<- c(state4[,1])
goi_state4
roi11 <- df[ rownames( df ) %in% goi_state4, "id" ]
roi11
sapply( goi_state4, function(x) grep( x, rownames( counts( cds ) ) ) )
heatmap.2(exprs(vsdFull)[roi11,], scale="row", col = hmcol, trace="none", margin=c(10, 6))

# gois from state5 pseudotime trajectory
state5<- read.csv(file ="/Volumes/sequence/harvey/git/3populationbulkRNAseq/Pseudotime_states/State5.csv", 
                  header = FALSE, sep = ",", stringsAsFactors = FALSE) 
state5
state5[,1]
goi_state5<- c(state5[,1])
goi_state5
roi12 <- df[ rownames( df ) %in% goi_state5, "id" ]
roi12
sapply( goi_state5, function(x) grep( x, rownames( counts( cds ) ) ) )
heatmap.2(exprs(vsdFull)[roi12,], scale="row", col = hmcol, trace="none", margin=c(10, 6))

