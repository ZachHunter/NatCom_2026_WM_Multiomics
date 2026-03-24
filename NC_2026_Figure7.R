library(dplyr)
library(tidyr)
library(purrr)
library(bvt)
library(edgeR)
library(destiny)
library(DESeq2)


#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 7
#
# Special note - Uses r package bvt for plotting
# Bioconductor visualization tools (bvt) is available at https://github.com/ZachHunter/bvt
# See readme.txt for more information.
#============================#



##############################
# File dependencies
# Loading required files for analysis
##############################

# Set Input/Output directories
dataDir<-"~/Desktop/DesktopData/Papers/Multiomics/Data/"
outputDir<-"~/Desktop/DesktopData/Papers/Multiomics"

# Loading eSets with relevant phenotype and feature data
# ExpressionSet of gene level TpM data from salmon.
load(file.path(dataDir,"studyTpM.RData"))
# SeqExpressionSet of Batch adjusted count data from STAR
load(file.path(dataDir,"studyCounts.RData"))
# variance stabilized transformation data from DESeq2 generated in Figure 1 code
load(file.path(dataDir,"studyVST.RData"))

# loading the dmap model derived previously
load(file=file.path(dataDir,"studyDmap.RData"))

# All of the final patient mutations in MAF format. Again, very large
load("~/Desktop/DesktopData/CurrentProjects/300/WES/MAF.RData")


##############################
# Custom Functions
# Includes style definitions and global settings
##############################

# Creating custom theme to make sure subtype colors remain the same across figures
WMSubtypeColor<-c(3,1,2,4:9)
EWMSubtypeColor<-c(3,8,1,2,9,5:7)
WMHDSubtypeColor<-c(4,5,3,1,2,6:9)
EWMHDSubtypeColor<-c(4,5,3,8,1,2,9,6:7)
SubtypeTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[WMSubtypeColor],fill=npDefaultTheme$plotColors$fill[WMSubtypeColor] ))
HDWMTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[WMHDSubtypeColor],fill=npDefaultTheme$plotColors$fill[WMHDSubtypeColor] ))
EHDWMTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[EWMHDSubtypeColor],fill=npDefaultTheme$plotColors$fill[EWMHDSubtypeColor] ))
ESubtypeTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[EWMSubtypeColor],fill=npDefaultTheme$plotColors$fill[EWMSubtypeColor] ))

# Sample IDS for WM patients only
WMOnly<-sampleNames(studyCounts)[grep("WM", sampleNames(studyCounts))]

# Vector of ENSG with at least some expression in at least 20 patients
# Counts used to determine power to detect differences, TpM used for biological relevance
WMExpressed<-featureNames(studyCounts)[rowSums(counts(studyCounts)[,WMOnly]>10)>20]
WMExpressed<-WMExpressed[rowSums(exprs(studyTpM)[WMExpressed,WMOnly] > 1) > 20]



##############################
# Figure 7A
# PreB, Tcell, Myeloid, and Stem cell Gene Expression by EScore Level
##############################

# Pre-Pro B-Cell Genes
pdf(file=file.path(outputDir,"Figures/Figure7/F7A_PreB_EScore.pdf"), width = 5, height = 8)
genePlot(studyTpM,c("RAG1","MYB","IGLL1"),
         group="EScore_5W",
         logScale=2,
         legend=TRUE,
         ylab="Transcripts per Million +1 (Log 2)",
         main="",
         legendSize=1.1,
         RSOverride=TRUE)
dev.off()

# T-cell Genes
pdf(file=file.path(outputDir,"Figures/Figure7/F7A_TCell_EScore.pdf"), width = 5, height = 8)
genePlot(studyTpM,c("CD3E","CD8A","CD4"),
         group="EScore_5W",
         logScale=2,
         legend=TRUE,
         ylab="Transcripts per Million +1 (Log 2)",
         main="",
         legendSize=1.1,
         RSOverride=TRUE)
dev.off()

# Myeloid Genes
pdf(file=file.path(outputDir,"Figures/Figure7/F7A_Myeloid_EScore.pdf"), width = 5, height = 8)
genePlot(studyTpM,c("CD33","FCGR3A","CD14"),
         group="EScore_5W",
         logScale=2,
         legend=TRUE,
         ylab="Transcripts per Million +1 (Log 2)",
         main="",
         legendSize=1.1,
         RSOverride=TRUE)
dev.off()

# Stem cell Genes
pdf(file=file.path(outputDir,"Figures/Figure7/F7A_StemCell_EScore.pdf"), width = 5, height = 8)
genePlot(studyTpM,c("CD34","KIT","PROM1"),
         group="EScore_5W",
         logScale=2,
         legend=TRUE,
         ylab="Transcripts per Million +1 (Log 2)",
         main="",
         legendSize=1.1,
         RSOverride=TRUE)
dev.off()


##############################
# Figure 7B
# EScore Gene Expression of XBP1 and PRDM1
##############################

pdf(file=file.path(outputDir,"Figures/Figure7/F7B_XBP1_PRDM1_EScore.pdf"), width = 7, height = 6)
genePlot(studyTpM,
         c("PRDM1","XBP1"),
         group="EScore_5W",
         logScale=2,
         theme=SubtypeTheme,
         plotColors=list(fill="lightgrey"),
         highlight="SimpleSubtype",
         legendSize=1.1,
         ylab="Transcripts per Million + 1 (Log 2 Scale)",
         main="",
         RSOverride=TRUE)
dev.off()


##############################
# Figure 7C
# BM vs EScore and MYD88 VAF
# Includes Generating MYD88 VAF Data
##############################

# Preprocessing MYD88 VAF data
MYD88<-MAF@data %>%
  dplyr::filter(Hugo_Symbol=="MYD88") %>%
  group_by(Tumor_Sample_Barcode) %>%
  dplyr::summarize(VAF=max(Raw_VAF,na.rm=T)) %>%
  dplyr::mutate(ID=gsub("_Tumor","",Tumor_Sample_Barcode)) %>%
  dplyr::select(VAF, ID) %>%
  dplyr::filter(! is.na(VAF), ID %in% sampleNames(studyTpM)) %>%
  as.data.frame()
rownames(MYD88) <- MYD88$ID

MYD88f10 <- MYD88[MYD88$VAF >= .4 & MYD88$VAF <= .6,]
HighP<-MYD88$ID %in% MYD88f10$ID
names(HighP)<-MYD88$ID

pdf(file=file.path(outputDir,"Figures/Figure7/F7C_BM_EScore.pdf"), width = 7, height = 7)
geneScatter(
  data.frame(
    EScore=pData(studyTpM)[MYD88$ID,"EScore"],
    BM=pData(studyTpM)[MYD88$ID,"bm"]),
  size=round(MYD88$VAF,3),
  color=HighP,
  legend=c("MYD88 VAF=.4-.6","MYD88 VAF"),
  legendSize=1.1, trendline="color",
  plotColors=list(fill=basicTheme$plotColors$points),
  axisText=list(x=c("",""),y=c("","%")),
  logScale=c(F,F),main="",
  minorTick=4,logAdjustment=0,rotateY=T,
  RSOverride=TRUE)
dev.off()


##############################
# Figure 7D
# MYD88 VAF by EScore - High Purity Only
##############################

pdf(file=file.path(outputDir,"Figures/Figure7/F7D_MYD88VAF_HighP_EScore.pdf"), width = 5, height = 7)
genePlot(MYD88f10$VAF, group=pData(studyTpM)[MYD88f10$ID,"EScore_5W"],
         ylab="MYD88 VAF", main="",RSOverride=TRUE)
dev.off()


##############################
# Figure 7E
# BM by EScore EL - High Purity Only
##############################

pdf(file=file.path(outputDir,"Figures/Figure7/F7E_BM_HighP_EScore.pdf"), width = 4, height = 7)
genePlot(pData(studyTpM)[MYD88f10$ID,"bm"], group=pData(studyTpM)[MYD88f10$ID,"EScore_EL"],
         ylab="WM LPC BM Involvement (%)", axisText=c("","%"),
         main="",RSOverride=TRUE)
dev.off()


##############################
# Figure 7F
# S100A9 by EScore and MYD88 VAF
##############################

pdf(file=file.path(outputDir,"Figures/Figure7/F7F_S100A9_EScore.pdf"), width = 5, height = 7)
genePlot(studyTpM[,MYD88$ID], "S100A9",
         group="EScore_5W",
         subgroup=HighP,
         logScale=2,
         ylab="Transcrips per Million (Log 2 Scale)",
         logAdjustment=0,
         main="",RSOverride=TRUE)
dev.off()




