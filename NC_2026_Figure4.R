library(dplyr)
library(tidyr)
library(purrr)
library(bvt)
library(edgeR)
library(DESeq2)



#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 4
#   Supplemental Figure panel 2
#   BCL vs. PCL DGE data
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
dataDir<-"~/Path/to/Multiomics/Data/"
outputDir<-"~/Path/to/Multiomics"

#Loading eSets with relevant phenotype and feature data
# ExpressionSet of gene level TpM data from salmon.
load(file.path(dataDir,"studyTpM.RData"))
# SeqExpressionSet of Batch adjusted count data from STAR
load(file.path(dataDir,"studyCounts.RData"))
# variance stabilized transformation data from DESeq2 generated in Figure 1 code
load(file.path(dataDir,"studyVST.RData"))

# Immune cell expression data from Human Protein Atlas
# https://www.proteinatlas.org/humanproteome/single+cell/immune+cell/data#immune_cells_hpa
Immune <- read.table(file.path(dataDir,"rna_immune_cell.tsv"), sep="\t", header=TRUE)


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



##############################
# Figure 4A
# BCL Subtype Asssociated Genes
##############################

BCLtoPlot<-c("KLF14","GRIK3", "SGSM1","GNAO1", "ACTG2")

pdf(file=file.path(outputDir,"Figures/Figure4/F4A_BCLGenes.pdf"), width = 6, height = 4)
genePlot(studyTpM,
         BCLtoPlot,
         group="SimpleSubtype", plotType="bar",
         theme=HDWMTheme,
         legend=TRUE, legendSize=1.1, ylab="Transcripts per Million (TpM)",
         RSOverride=TRUE, main="")
dev.off()


##############################
# Figure 4B
# PCL Subtype Asssociated Genes
##############################

PCLtoPlot<-c("WNK2","PRDM5", "GPER1","IL17RB", "FAM110C")

pdf(file=file.path(outputDir,"Figures/Figure4/F4B_PCLGenes.pdf"), width = 6, height = 4)
genePlot(studyTpM,
         PCLtoPlot,
         group="SimpleSubtype", plotType="bar",
         theme=HDWMTheme,
         legend=TRUE, legendSize=1.1, ylab="Transcripts per Million (TpM)",
         RSOverride=TRUE, main="")
dev.off()


##############################
# Figure 4C & 4D
# KLF14 and WNK2 TMA Path
# No R analysis needed
##############################


##############################
# Supplemental Data
# BCL vs PCL DGE Data
# Included in the supplemental data package
##############################

#Models and contrast
subtypeMod<-model.matrix(~ 0 + factor(SimpleSubtype) + EScore + factor(CXCR4) + gender+ agebmbx, data=pData(studyCounts)[WMOnly,])
colnames(subtypeMod)<-c("EWM","BCL","PCL","EScore","CXCR4","Sex","Age")
subtypeCon<-makeContrasts(BCL-PCL, levels=subtypeMod)

v<-voom(
  calcNormFactors(
    DGEList(
      counts(studyCounts)[WMExpressed,WMOnly],
      genes=fData(studyCounts)[WMExpressed,c("Chr","GeneSymbol")]
    )
  ),design=subtypeMod)
fit<-lmFit(v, design = subtypeMod)
fit<-contrasts.fit(fit, subtypeCon)
fit<-eBayes(fit)

# Recording top hist and removing unnecessary list structures
BCLvPCL<-topTable(fit,p.value =0.05,number=Inf, lfc=.8)
BCLvPCL$Chr<-unlist(BCLvPCL$Chr)
BCLvPCL$GeneSymbol<-unlist(BCLvPCL$GeneSymbol)
write.table(BCLvPCL,sep="\t", quote=FALSE, file=file.path(outputDir,"DGE/BCL_vs_PCL_DGE.tsv"))


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
# Supplemental Figure 2B
# BCL Signature genes expression in HD immune cells populations
##############################

Immune_BCL<-Immune[Immune$Gene.name %in% BCLtoPlot,]

pdf(file=file.path(outputDir,"Figures/SFigure2/SF2B_BCLImmune.pdf"), width = 8, height = 5)
genePlot(Immune_BCL$nTPM,
         group=Immune_BCL$Gene.name,
         subgroup=Immune_BCL$Immune.cell,
         plotType="bar",
         plotColors=list(fill=rainbow(19)),
         main="",
         legendSize=1,
         ylab="Normalized Transcripts per Million (nTpM)",
         RSOverride=TRUE)
dev.off()


##############################
# Supplemental Figure 2C
# BCL Signature genes expression in HD immune cells populations
##############################

Immune_PCL<-Immune[Immune$Gene.name %in% PCLtoPlot,]

pdf(file=file.path(outputDir,"Figures/SFigure2/SF2C_PCLImmune.pdf"), width = 8, height = 5)
genePlot(Immune_PCL$nTPM,
         group=Immune_PCL$Gene.name,
         subgroup=Immune_PCL$Immune.cell,
         plotType="bar",
         plotColors=list(fill=rainbow(19)),
         main="",
         legendSize=1,
         ylab="Normalized Transcripts per Million (nTpM)",
         RSOverride=TRUE)
dev.off()



