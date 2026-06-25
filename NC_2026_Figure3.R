library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(bvt)
library(DESeq2)
library(edgeR)
library(lubridate)
library(GEOquery)

#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 3
#   Supplementary figure panel 2
#
# Special note - Uses r package bvt for plotting
# Bioconductor visualization tools (bvt) is available at https://github.com/ZachHunter/bvt
# See readme.txt for more information.
#============================#



##############################
# File dependencies
# Loading required files for analysis
##############################

#Set Input/Output directories
dataDir<-"~/Path/to/Directory/Data/"
outputDir<-"~/Path/to/Directory/"

#Loading eSets with relevant phenotype and feature data
# ExpressionSet of gene level TpM data from salmon.
load(file.path(dataDir,"studyTpM.RData"))
# SeqExpressionSet of Batch adjusted count data from STAR
load(file.path(dataDir,"studyCounts.RData"))
# variance stabilized transformation data from DESeq2 generated in Figure 1 code
load(file.path(dataDir,"studyVST.RData"))


##############################
# Custom Functions
# Includes style definitions and global settings
##############################

# Custom theme to make sure subtype colors remain the same across figures
WMSubtypeColor<-c(3,1,2,4:9)
WMHDSubtypeColor<-c(4,5,3,1,2,6:9)
SubtypeTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[WMSubtypeColor],fill=npDefaultTheme$plotColors$fill[WMSubtypeColor] ))
HDWMTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[WMHDSubtypeColor],fill=npDefaultTheme$plotColors$fill[WMHDSubtypeColor] ))

# Sample IDS for WM patients only
WMOnly<-sampleNames(studyCounts)[grep("WM", sampleNames(studyCounts))]

# Vector of ENSG with at least some expression in at least 20 patients
# Counts used to determine power to detect differences, TpM used for biological relevance
WMExpressed<-featureNames(studyCounts)[rowSums(counts(studyCounts)[,WMOnly]>10)>20]
WMExpressed<-WMExpressed[rowSums(exprs(studyTpM)[WMExpressed,WMOnly] > 1) > 20]

#WM Gene Sets for EScore and Subtype from NMF analysis in Figure 2
EScore_GeneSet<-c("LYZ", "FCN1", "VCAN", "S100A9", "S100A8", "CD14", "SERPINA1", "CST3", "CPVL",
                  "CD93", "CD163", "CSF3R", "CD36", "TGFBI", "C5AR1", "FPR1", "ANXA1", "LTF", "CD300E",
                  "FGR", "FCER1G", "CD4", "CSF1R", "CXCL8", "FCGR3A", "GAS7", "S100A12")

BCL_GeneSet<-c("AC092868.1", "KLF14", "GRIK3", "SGSM1", "ACTG2", "USP44", "MYOCD", "PDE7B",
               "ACE", "BMP4", "L1CAM", "GNAO1", "LINC00494", "MYO1E", "HTR3A", "TNR", "UST",
               "ERICH3", "TMEM56", "LIMA1", "PLD5", "ACKR3", "ALPL", "GABBR1", "AC136475.3",
               "LAMC3", "TRIM2", "DSCAML1")

PCL_GeneSet<-c("WNK2", "PRDM5", "LRRC31", "CADPS2", "RAB3B", "SNAI2", "CDHR1", "MIR99AHG",
               "SIX4", "FAM110C", "CCDC198", "WNT5A", "PIEZO2", "NOL4", "EHF", "MOCS1",
               "TBX15", "NPTX1", "PURPL", "IL17RB", "ANKRD29", "AC079466.1", "FAM19A5",
               "GPER1", "LINC00887", "AC023796.2", "CXXC4", "MYO6")

#calculating the mean of scaled probes per gene
probeToScaledGene<-function(ged, gene, gch="Gene Symbol") {
  probes<-rownames(fData(ged))[fData(ged)[,gch] ==gene]
  if(length(probes)>1) {
    exprs(ged)[probes,] %>% t() %>% scale() %>% rowMeans()
  } else {
    scale(exprs(ged)[probes,])
  }
}

#calculating gene signature scores
scoreSig<- function(ged, sigGenes, gch="Gene Symbol") {
  map(sigGenes[sigGenes %in% fData(ged)[,gch]], function(x) probeToScaledGene( ged,x,gch=gch) ) %>%
    data.frame() %>% rowSums()
}



##############################
# Figure 3A
# BM by EScore and Subtype
##############################

pdf(file=file.path(outputDir,"Figures/Figure3/F3A_BM_EScoreSubtype.pdf"), width = 7, height = 6)

genePlot(pData(studyTpM)[WMOnly,"bm"],
         group=pData(studyTpM)[WMOnly,"EScoreEL"],
         subgroup=pData(studyTpM)[WMOnly,"SimpleSubtype"],
         axisText=c("","%"),
         legendSize=1.2, main="",
         groupLabSize=1.3,
         axisLabelSize=1.5,
         yAxisLabSize=1.15,
         ylab="WM LPC Bone Marrow Involvement",
         theme=SubtypeTheme,
         RSOverride=TRUE)
dev.off()

EarlyID<-sampleNames(studyTpM)[pData(studyTpM)$EScoreEL=="Early EScore" & ! is.na(pData(studyTpM)$EScoreEL)]
LateID<-sampleNames(studyTpM)[pData(studyTpM)$EScoreEL=="Late EScore" & ! is.na(pData(studyTpM)$EScoreEL)]

#Testing impact on BM
p1<-wilcox.test(pData(studyTpM)[EarlyID[pData(studyTpM)[EarlyID,"SimpleSubtype"]=="BCL"],"bm"],
                pData(studyTpM)[EarlyID[pData(studyTpM)[EarlyID,"SimpleSubtype"]=="PCL"],"bm"])$p.value
p2<-wilcox.test(pData(studyTpM)[LateID[pData(studyTpM)[LateID,"SimpleSubtype"]=="BCL"],"bm"],
                pData(studyTpM)[LateID[pData(studyTpM)[LateID,"SimpleSubtype"]=="PCL"],"bm"])$p.value
p3<-wilcox.test(pData(studyTpM)[WMOnly,"bm"]~pData(studyTpM)[WMOnly,"EScoreEL"])$p.value

p.adjust(c(p1,p2,p3))


##############################
# Figure 3B
# Subtype prevalence at clinical staging
##############################

preTherapyMarrow<-interval(mdy(pData(studyTpM)$datebmbx),mdy(pData(studyTpM)$datetx1)) / weeks(1) < 12

a<-xtabs(~ preTherapyMarrow + pData(studyTpM)$SimpleSubtype)[,3:5]
fisher.test(a)
a<-as.data.frame(round(a/rowSums(a)*100,1))

pdf(file=file.path(outputDir,"Figures/Figure3/F3B_Subtype12W.pdf"), width = 7, height = 6)
genePlot(a$Freq,
         group=a$preTherapyMarrow,
         subgroup=a$pData.studyTpM..SimpleSubtype,
         plotType="bar",
         axisText=c("","%"),
         groupLabels=c("> 12 Weeks ", "< 12 Weeks"),
         main="",
         groupLabSize=1.4,
         yAxisLabSize=1.3,
         theme=SubtypeTheme,
         RSOverride=TRUE)
dev.off()


##############################
# Figure 3C
# Early/Late EScore prevalence at clinical staging
##############################

a<-xtabs(~ preTherapyMarrow + pData(studyTpM)$EScoreEL)[,3:4]
fisher.test(a)
a<-as.data.frame(round(a/rowSums(a)*100,1))

pdf(file=file.path(outputDir,"Figures/Figure3/F3C_EScore12W.png"), width = 7, height = 6)
genePlot(a$Freq,
         group=a$preTherapyMarrow,
         subgroup=a$pData.studyTpM..EScoreEL,
         plotType="bar",
         axisText=c("","%"),
         groupLabSize=1.4,
         yAxisLabSize=1.3,
         groupLabels=c("> 12 Weeks ", "< 12 Weeks"),
         main="",
         RSOverride=TRUE)
dev.off()


##############################
# Figure 3D
# Analysis of IgM MGUS and WM GEO data for validation
# Includes downloading and pre-processing steps
# as well as EScore/Subtype gene signature calculation
##############################

#Downloading Trojani et al. data from GEO
ATMGUS<-getGEO("GSE171739")

#Preprocessing data and selecting the CD19 selected data from compatibility with this study
ATMGUS<-ATMGUS[[1]]
exprs(ATMGUS)<-normalizeBetweenArrays(exprs(ATMGUS))
MGUSvWM<-ATMGUS[,pData(ATMGUS)[,"diagnosis:ch1"]!="CTRL" & pData(ATMGUS)[,"cell type:ch1"]=="CD19"]

#Calculating Signature Scores
MGUSvWM_EScore<-scoreSig(MGUSvWM,EScore_GeneSet)
MGUSvWM_BCL<-scoreSig(MGUSvWM,BCL_GeneSet)
MGUSvWM_PCL<-scoreSig(MGUSvWM,PCL_GeneSet)

# Kmeans clustering into BCL and PCL
set.seed(265)
ImputedSubtype<-kmeans(MGUSvWM_BCL-MGUSvWM_PCL,2,nstart = 20)$cluster
names(ImputedSubtype)<-names(MGUSvWM_PCL)

# Assigning apriori MGUS status and updating labels
ImputedSubtype[pData(MGUSvWM)[names(ImputedSubtype),"diagnosis:ch1"]=="IgMMGUS"]<-"IgM MGUS"
ImputedSubtype[ImputedSubtype==1]<-"PCL"
ImputedSubtype[ImputedSubtype==2]<-"BCL"
ImputedSubtype<-factor(ImputedSubtype, levels=c("IgM MGUS", "BCL","PCL"))
pData(MGUSvWM)$Subtype<-ImputedSubtype

WMSigData<- data.frame(
  Subtype=MGUSvWM_BCL-MGUSvWM_PCL,
  EScore=-MGUSvWM_EScore)

# Ploting data
pdf(file=file.path(outputDir,"Figures/Figure3/F3D_GEO_Validation.pdf"), width = 8, height = 6)

geneScatter(
  WMSigData,
  color=pData(MGUSvWM)[,"diagnosis:ch1"],
  shape=ImputedSubtype,
  legend=c("Sample Type", "Imputed Subtype"),
  legendSize=1,
  main="",
  pointSize=1.2,
  legendSpacing=.5,
  ylab="EScore Signature Score",
  xlab="Subtype Signature Score",
  theme=SubtypeTheme,
  RSOverride=TRUE)
dev.off()


##############################
# Figure 3E
# EScore signature differences between IgM MGUS and WM
# (Validation using Trojani GEO data)
##############################

pdf(file=file.path(outputDir,"Figures/Figure3/F3E_GEO_EScore.pdf"), width = 4, height = 6)
genePlot(-MGUSvWM_EScore,
         group=ImputedSubtype!="IgM MGUS",
         ylab="EScore",
         groupLabels=c("IgM MGUS", "WM"),
         main="",
         showCalc=TRUE,
         theme=SubtypeTheme,
         groupLabSize=1.4,
         yLabelSize=1.3,
         RSOverride=TRUE)
dev.off()


##############################
# Figure 3F
# Subtype Score density distribution in WM and IgM MGUS
# (Validation using Trojani GEO data)
##############################

pdf(file=file.path(outputDir,"Figures/Figure3/F3F_GEO_Density.pdf"), width = 8, height = 4)
genePlot(MGUSvWM_BCL-MGUSvWM_PCL,
         group=ImputedSubtype,
         plotType="density",
         legendSize=1,
         legendSpacing=.5,
         drawRug=T,
         main="",
         theme=SubtypeTheme,
         RSOverride=TRUE)
dev.off()

pairwise.wilcox.test(WMSigData$Subtype, ImputedSubtype)

##############################
# Supplemental Figure 2A
# Signature genes by subtype and CXCR4 status
##############################

pdf(file=file.path(outputDir,"Figures/SFigure2/SF2A_SubtypeCXCR4.pdf"), width = 9, height = 5)
genePlot(
  studyTpM[,WMOnly],
  gene = c("WNK2","DUSP22","GPER1","CABLES1"),
  group="CXCR4",
  highlight="SimpleSubtype",
  logScale=2,
  legendSize=1.1,
  guides=TRUE,
  lWidth=2,
  plotColors=list(lines=setAlpha("black",.8),fill="white"),
  errorBarLineType=1,
  pointSize=.65,
  theme=SubtypeTheme,
  ylab="Transcripts per Million + 1 (Log 2)",
  groupLabelSpacing=1.2,
  subgroupLabSize=.75,
  main="",
  RSOverride=TRUE)
dev.off()

##############################
# Supplementary Figure 2B
# Scatter plot of WNK2 and KLF14 by WM Subtype
# Note that we are using VST here to help with relative scaling
# Neither gene is expressed in B-cells or in the BM enviroment generally
##############################

pdf(file=file.path(outputDir,"Figures/SFigure2/SF2B_WNK2_KLF14.pdf"), width = 8, height = 6)
geneScatter(
  studyVST,
  c("KLF14","WNK2"),
  color="SimpleSubtype",
  size=round(pData(studyTpM)$EScore,3),
  legend=c("WM Subtype","EScore"),
  legendSize=1.2,
  main="",
  pointSize=1,
  sizeScale=2.3,
  legendSpacing=.3,
  xlab="KLF14 (VST)",
  ylab="WNK2 (VST)",
  theme=SubtypeTheme,
  yAxisLabSize=1.1,
  RSOverride=TRUE)
dev.off()

##############################
# Supplementary Figure 2C
# Scatter plot of GRIK3 and PRDM5 by WM Subtype
# Note that we are using VST here to help with relative scaling
# Neither gene is expressed in B-cells or in the BM environment generally
##############################

pdf(file=file.path(outputDir,"Figures/SFigure2/SF2C_PRDM5_GRIK3.pdf"), width = 8, height = 6)
geneScatter(
  studyVST,
  c("GRIK3","PRDM5"),
  color="SimpleSubtype",
  size=round(pData(studyTpM)$EScore,3),
  legend=c("WM Subtype","EScore"),
  legendSize=1.2,
  main="",
  pointSize=1,
  sizeScale=2.3,
  legendSpacing=.3,
  xlab="GRIK3 (VST)",
  ylab="PRDM5 (VST)",
  theme=SubtypeTheme,
  yAxisLabSize=1.1,
  RSOverride=TRUE)
dev.off()

##############################
# Supplementary Figure 2D
# Identification of Early WM Samples in Trojani et al data set
##############################

TrojaniSubtype<-rep("PCL",nrow(WMSigData))
TrojaniSubtype[MGUSvWM_BCL-MGUSvWM_PCL>0]<-"BCL"
#Note EScore is plotted is defined here as -MGUSvWM_EScore
#This is to make it increase with disease progression like the RNASeq DPT
TrojaniSubtype[MGUSvWM_EScore > 20]<-"Early WM"
TrojaniSubtype<-factor(TrojaniSubtype, levels=c("Early WM","BCL","PCL"))

pdf(file=file.path(outputDir,"Figures/SFigure2/SF2D_Trojani_Subtype.pdf"), width = 8, height = 6)
geneScatter(
  WMSigData,
  color=TrojaniSubtype,
  shape=TrojaniSubtype,
  legend="WM Subtype",
  legendSize=1,
  main="",
  pointSize=1.2,
  legendSpacing=.5,
  ylab="EScore Signature Score",
  xlab="Subtype Signature Score",
  theme=SubtypeTheme,
  yAxisLabSize=1.1,
  RSOverride=TRUE)
dev.off()

##############################
# Supplementary Figure 2E
# Documenting intermediate levels of BCL, PCL and Subtype gene signature
# scores in Early WM compared with the mature BCL and PCL subtypes shown
# this time in the Trojani et al data set.
##############################

pdf(file=file.path(outputDir,"Figures/SFigure2/SF2E_Trojani_SubtypeScoreLevels.pdf"), width = 8, height = 6)
genePlot(data.frame(MGUSvWM_BCL,MGUSvWM_PCL,MGUSvWM_BCL-MGUSvWM_PCL),
         group=TrojaniSubtype,
         theme=SubtypeTheme,
         ylab="Scaled Score",
         main="",legend=TRUE,
         groupLabSize=1.4,
         axisLabelSize=1.5,
         yAxisLabSize=1.15,
         legendSize=1.3,
         groupLabels=c("BCL Score","PCL Score", "Subtype Score"),
         RSOverride=TRUE)
dev.off()

# Pairwise testing by subtype
pairwise.wilcox.test(MGUSvWM_BCL,TrojaniSubtype)
pairwise.wilcox.test(MGUSvWM_PCL,TrojaniSubtype)
pairwise.wilcox.test(MGUSvWM_BCL-MGUSvWM_PCL,TrojaniSubtype)
