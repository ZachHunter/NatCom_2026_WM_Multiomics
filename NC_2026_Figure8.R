library(dplyr)
library(tidyr)
library(purrr)
library(bvt)
library(edgeR)
library(maftools)

#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 8
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
dataDir<-"~/Path/to/Directory/Data/"
outputDir<-"~/Path/to/Directory/"

# Loading eSets with relevant phenotype and feature data
# ExpressionSet of gene level TpM data from salmon.
load(file.path(dataDir,"studyTpM.RData"))
# SeqExpressionSet of Batch adjusted count data from STAR
load(file.path(dataDir,"studyCounts.RData"))
# variance stabilized transformation data from DESeq2 generated in Figure 1 code
load(file.path(dataDir,"studyVST.RData"))

# All of the final patient mutations in MAF format.
load("~/Path/to/Directory/WES/MAF.RData")


##############################
# Custom Functions
# Includes style definitions and global settings
##############################

# Sample IDS for WM patients only
WMOnly<-sampleNames(studyCounts)[grep("WM", sampleNames(studyCounts))]

# Vector of ENSG with at least some expression in at least 20 patients
# Counts used to determine power to detect differences, TpM used for biological relevance
WMExpressed<-featureNames(studyCounts)[rowSums(counts(studyCounts)[,WMOnly]>10)>20]
WMExpressed<-WMExpressed[rowSums(exprs(studyTpM)[WMExpressed,WMOnly] > 1) > 20]

# Preprocessing MYD88 VAF Data
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
HighP<-factor(MYD88$ID %in% MYD88f10$ID, labels=c("MYD8 VAF<0.4","MYD88 VAF≥0.4"))
names(HighP)<-MYD88$ID



##############################
# Figure 8A
# CD4, CD33 and CD34 Expression by Early/Late EScore and MYD88 VAF Purity
##############################

# This is combining the Early/Late EScore with Purity into a single factor for testing
VAF_EL_Fact<-factor(paste0(pData(studyTpM)[MYD88$ID,"EScore_EL"],HighP))

# Plotting CD4
pdf(file=file.path(outputDir,"Figures/Figure8/F8A_CD4_MYD88VAF_EScore.pdf"), width = 6, height = 7)
genePlot(studyTpM[,MYD88$ID], "CD4",
         group="EScore_EL",
         subgroup=HighP,
         logScale=2,
         groupLabSize=1.3,
         axisLabelSize=1.2,
         ylab="Transcripts per Million +1 (Log 2 Scale)",
         main="", RSOverride=TRUE)
dev.off()
pairwise.wilcox.test(exprs(studyTpM)[fData(studyTpM)$GeneSymbol=="CD4",MYD88$ID],VAF_EL_Fact)

# Plotting CD33
pdf(file=file.path(outputDir,"Figures/Figure8/F8A_CD33_MYD88VAF_EScore.pdf"),width = 6, height = 7)
genePlot(studyTpM[,MYD88$ID], "CD33",
         group="EScore_EL",
         subgroup=HighP,
         logScale=2,
         ylab="Transcripts per Million +1 (Log 2 Scale)",
         groupLabSize=1.3,
         axisLabelSize=1.2,
         main="", RSOverride=TRUE)
dev.off()
pairwise.wilcox.test(exprs(studyTpM)[fData(studyTpM)$GeneSymbol=="CD33",MYD88$ID],VAF_EL_Fact)

# Plotting CD34
pdf(file=file.path(outputDir,"Figures/Figure8/F8A_CD34_MYD88VAF_EScore.pdf"), width = 6, height = 7)
genePlot(studyTpM[,MYD88$ID], "CD34",
         group="EScore_EL",
         subgroup=HighP,
         logScale=2,
         ylab="Transcripts per Million +1 (Log 2 Scale)",
         groupLabSize=1.3,
         axisLabelSize=1.2,
         main="", RSOverride=TRUE)
dev.off()
pairwise.wilcox.test(exprs(studyTpM)[fData(studyTpM)$GeneSymbol=="CD34",MYD88$ID],VAF_EL_Fact)


##############################
# Figure 8B
# Path Figure Panel - No R Code
##############################
