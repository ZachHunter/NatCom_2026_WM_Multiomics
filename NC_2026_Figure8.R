library(dplyr)
library(tidyr)
library(purrr)
library(bvt)
library(edgeR)
library(maftools)
library(emmeans)
library(lmtest)
library(sandwich)
library(kableExtra)

#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 8
#   Supplementary table 11
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
VAF_5W_Fact<-factor(paste0(pData(studyTpM)[MYD88$ID,"EScore_5W"],HighP))

# Plotting CD4
pdf(file=file.path(outputDir,"Figures/Figure8/F8A_CD4_MYD88VAF_EScore.pdf"), width = 6, height = 7)
genePlot(studyTpM[,MYD88$ID], "CD4",
         group="EScore_5W",
         subgroup=HighP,
         logScale=2,
         groupLabSize=1.3,
         axisLabelSize=1.2,
         ylab="Transcripts per Million +1 (Log 2 Scale)",
         main="", RSOverride=TRUE)
dev.off()

# Plotting CD33
pdf(file=file.path(outputDir,"Figures/Figure8/F8A_CD33_MYD88VAF_EScore.pdf"),width = 6, height = 7)
genePlot(studyTpM[,MYD88$ID], "CD33",
         group="EScore_5W",
         subgroup=HighP,
         logScale=2,
         ylab="Transcripts per Million +1 (Log 2 Scale)",
         groupLabSize=1.3,
         axisLabelSize=1.2,
         main="", RSOverride=TRUE)
dev.off()

# Plotting CD34
pdf(file=file.path(outputDir,"Figures/Figure8/F8A_CD34_MYD88VAF_EScore.pdf"), width = 6, height = 7)
genePlot(studyTpM[,MYD88$ID], "CD34",
         group="EScore_5W",
         subgroup=HighP,
         logScale=2,
         ylab="Transcripts per Million +1 (Log 2 Scale)",
         groupLabSize=1.3,
         axisLabelSize=1.2,
         main="", RSOverride=TRUE)
dev.off()

##############################
# Supplementary table 11
# Statistical analysis of the impact of MYD88 VAF
# on EScore Associated gene expression
##############################

# All EScore genes explicitly discussed in manuscript
genesSets<-c("S100A9","CXCL1","CXCL8","CXCL12","CD38","SDC1","CD19","MS4A1","RAG1","MYB","IGLL1","CD3E","CD8A","CD4","CD33","FCGR3A","CD14","CD34","KIT","PROM1")

# First we will contrast samples with 0.4 ≤ MYD88 VAF ≤ 0.6 vs all others at each EScore level for each gene
geneStats<-lapply(genesSets, function(g) {
  a<-map_dbl(paste0("ESL",1:5), function(x) {
    IDS<-MYD88$ID[pData(studyTpM)[MYD88$ID,"EScore_5W"]==x]
    wilcox.test(exprs(studyTpM)[fData(studyTpM)$GeneSymbol==g,IDS] ~ factor(HighP[IDS]))$p.value
  }) %>%
    p.adjust(method = "fdr")
  names(a)<-paste0("ESL",1:5)
  a
})
names(geneStats)<-genesSets

as.data.frame(geneStats) %>%
  t() %>%
  kbl(caption = "FDR adjusted Wilcoxan tests of gene expression by EScore Level stratified by 0.4 ≤ MYD88 VAF ≤ 0.6") %>%
  kable_classic(full_width=FALSE)  %>%
  column_spec(column=1, italic = TRUE) %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST11_EScoreGenes_HighPurity_by_EScore.pdf"), width=6)

# As these genes are typically log linear with EScore, we can use linear models
# contrasting high purity status with continuous EScore or EScore levels useing
# estimated marginal means testing.

# Making continuous EScore models
modelList<-lapply(genesSets, function(g) {
  dat<-data.frame(
    L2TpM=log(exprs(studyTpM)[fData(studyTpM)$GeneSymbol==g,names(HighP)] +1 ,2),
    HighP,
    EScore=pData(studyTpM)[names(HighP),"EScore"]
  )
  lm(L2TpM ~ HighP + EScore, data=dat)
})

# Making discrete EScore models with ESL1-5
modelListLevel<-lapply(genesSets, function(g) {
  dat<-data.frame(
    L2TpM=log(exprs(studyTpM)[fData(studyTpM)$GeneSymbol==g,names(HighP)] +1 ,2),
    HighP,
    ESL=pData(studyTpM)[names(HighP),"EScore_5W"]
  )
  lm(L2TpM ~ HighP + ESL, data=dat)
})

#quick check to review models and make sure everything is OK
lapply(modelList, function(x) summary(x))
lapply(modelListLevel, function(x) summary(x))

# Lets look at the estimated marginal means results per gene
EMM_Genes<-map_dbl(modelList, function(m) {
  if(bptest(m)$p.value<0.05) {
    Mvcov<-vcovHC(m, type = "HC3")
    test(emmeans(m, pairwise ~ HighP, vcov. = Mvcov)$contrasts)[[6]]
  } else {
    test(emmeans(m, pairwise ~ HighP)$contrasts)[[6]]
  }
}) %>% p.adjust("fdr")
names(EMM_Genes)<-genesSets

# Same thing but using the discrete EScore model
EMM_Genes_ESL<-map_dbl(modelListLevel, function(m) {
  if(bptest(m)$p.value<0.05) {
    Mvcov<-vcovHC(m, type = "HC3")
    test(emmeans(m, pairwise ~ HighP, vcov. = Mvcov)$contrasts)[[6]]
  } else {
    test(emmeans(m, pairwise ~ HighP)$contrasts)[[6]]
  }
}) %>% p.adjust("fdr")

names(EMM_Genes_ESL)<-genesSets
data.frame(EScore=EMM_Genes,ESLevel=EMM_Genes_ESL) %>%
  kbl(caption = "FDR adjusted estimated marginal means testing of log 2 gene expression by 0.4 ≤ MYD88 VAF ≤ 0.6 and EScore or EScore Level (ESLevel)") %>%
  kable_classic(full_width=FALSE)  %>%
  column_spec(column=1, italic = TRUE) %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST11_EMM_HighPurity_by_EScore.pdf"), width=6)


# Now we will perform the reciprocal test using continuous MYD88 VAF and Early vs. Late EScore
modelEL_VAF<-lapply(genesSets, function(g) {
  dat<-data.frame(
    L2TpM=log(exprs(studyTpM)[fData(studyTpM)$GeneSymbol==g,MYD88$ID] +1 ,2),
    VAF=MYD88$VAF,
    ESL=pData(studyTpM)[MYD88$ID,"EScore_EL"]
  )
  lm(L2TpM ~ ESL + VAF, data=dat)
})

#quick check to make sure these models look OK
lapply(modelEL_VAF, function(x) summary(x))

EMM_Genes_EL_VAFL<-map_dbl(modelEL_VAF, function(m) {
  if(bptest(m)$p.value<0.05) {
    Mvcov<-vcovHC(m, type = "HC3")
    test(emmeans(m, pairwise ~ ESL, vcov. = Mvcov)$contrasts)[[6]]
  } else {
    test(emmeans(m, pairwise ~ ESL)$contrasts)[[6]]
  }
}) %>% p.adjust("fdr")

# The least significant gene for EScore was CD20 (MS4A1) and was still p <0.0007)
max(EMM_Genes_EL_VAFL)
genesSets[which(EMM_Genes_EL_VAFL==max(EMM_Genes_EL_VAFL))]


##############################
# Figure 8B
# Path Figure Panel - No R Code
##############################
       
