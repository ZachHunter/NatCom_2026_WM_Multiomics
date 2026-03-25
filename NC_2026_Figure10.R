library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(lubridate)
library(survival)
library(survminer)
library(kableExtra)
library(edgeR)

#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 10
#
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

# mSigDb data prepared for camera/fry analysis indexed for WMExpressed
load(file.path(dataDir,"study_mSig.RData"))


##############################
# Custom Functions
# Includes style definitions and global settings
##############################

WMOnly<-sampleNames(studyCounts)[grep("WM", sampleNames(studyCounts))]

WMExpressed<-featureNames(studyCounts)[rowSums(counts(studyCounts)[,WMOnly]>10)>20]
WMExpressed<-WMExpressed[rowSums(exprs(studyTpM)[WMExpressed,WMOnly] > 1) > 20]



##############################
# Figure 10A PFS in BTKi based Therapies by Early/Late EScore
# Includes generation of other relevant PFS data sets
##############################

# Collecting necessary clinical data
ClinDat<-pData(studyTpM)[WMOnly,]%>%
  mutate(ID=WMOnly)%>%
  dplyr::filter(!is.na(tx1)) %>%
  dplyr::select(ID, tx1, pfs1, pdtx1, responsetx1) %>%
  mutate(
    Subtype=factor(pData(studyTpM)[ID,"ExpandedSubtype"]),
    EScore=pData(studyTpM)[ID,"EScore_EL"],
    EScore5W=factor(pData(studyTpM)[ID,"EScore_5W"] %in% c("ESL4","ESL5"), labels=c("ESL1-3","ESL4-5"))
  )

# Subsetting data into BTKi, Ibrutinib monotherapy, and Proteasome Inhibitor base therapies
RxType<-ClinDat$tx1
ClinDatIB<-filter(ClinDat, RxType=="IB")

RxType[RxType %in% c("BDR" ,"CARD","IDR")]<-"Proteasome"
RxType[RxType %in% c("BENDA" ,"BENDAR","FR", "CPR", "CDR", "BR")]<-"Chemo"
RxType[RxType %in% c("OFA" ,"R","R-dex")]<-"CD20mab"
RxType[RxType %in% c("IB" ,"IB/ULO","ZANUBRUTINIB","IB/ULO + Zanubrutinib")]<-"BTKi"
RxType[RxType %in% c("RAD001" ,"IVEN")]<-"Other"
ClinDat<-mutate(ClinDat,RxType=RxType)
ClinDatBTKi<-filter(ClinDat, RxType=="BTKi")
ClinDatProt<-filter(ClinDat, RxType=="Proteasome")

# Plotting KM PFS Curves with risk table
pdf(file=file.path(outputDir,"Figures/Figure10/F10A_BTKi_EScore_EL_Surv.pdf"), width = 9, height = 6)
ggsurvplot(
  surv_fit(Surv(ClinDatBTKi$pfs1, event=ClinDatBTKi$pdtx1) ~ EScore, data=ClinDatBTKi),
  conf.int = FALSE,
  pval = TRUE,
  ylab="Progression Free Survival",
  risk.table=T,
  xlab="Time (Months)",
  risk.table.height=.35,
  title="BTKi Based Therapies"
)
dev.off()


##############################
# Figure 10B PFS in BTKi based Therapies by Late EScore (ESL4/5)
# Includes BTK monotherapy data that is discussed in text but not shown.
##############################

ClinDatBTKi$EScore<-ClinDatBTKi$EScore5W

# Plotting KM PFS Curves with risk table
pdf(file=file.path(outputDir,"Figures/Figure10/F10B_BTKi_LateEScore_Surv.pdf"), width = 9, height = 6)
ggsurvplot(
  surv_fit(Surv(ClinDatBTKi$pfs1, event=ClinDatBTKi$pdtx1) ~ EScore, data=ClinDatBTKi),
  conf.int = FALSE,
  pval = TRUE,
  ylab="Progression Free Survival",
  risk.table=T,
  xlab="Time (Months)",
  risk.table.height=.35,
  title="BTKi Based Therapies"
)
dev.off()

ClinDatIB$EScore<-ClinDatIB$EScore5W
ggsurvplot(
  surv_fit(Surv(ClinDatIB$pfs1, event=ClinDatIB$pdtx1) ~ EScore, data=ClinDatIB),
  conf.int = FALSE,
  pval = TRUE,
  ylab="Progression Free Survival",
  risk.table=T,
  xlab="Time (Months)",
  risk.table.height=.35,
  title="Ibrutinib Monotherapy"
)

##############################
# Figure 10C - BCR pathway up regulation with EScore Leveles 4/5
##############################

# defining model and contrast for GSEA analysis
EScoreMod<-model.matrix(~ 0 + factor(EScore_5W %in% c("ESL4","ESL5")) + factor(ExpandedSubtype) + factor(CXCR4) + gender + agebmbx, data=pData(studyCounts)[WMOnly,])
colnames(EScoreMod)<-c("Early","Late","EBCL","EPCL","PCL","EarlyWM","CXCR4","Sex","Age")
EScoreCon<-makeContrasts(Late-Early, levels=EScoreMod)

# voom/limma calculations
v<-voom(
  calcNormFactors(
    DGEList(counts=counts(studyCounts)[WMExpressed,WMOnly],genes=fData(studyCounts)[WMExpressed,c("Chr","GeneSymbol")])
  ),
  design=EScoreMod
)

fit<-lmFit(v,design = EScoreMod)
fit<-contrasts.fit(fit,EScoreCon)
fit<-eBayes(fit)

# Camera finds a ton of results. Limiting it to set of interest
camera(v,
       study_mSig[[3]],
       design = EScoreMod,
       contrast = EScoreCon)["SIG_BCR_SIGNALING_PATHWAY",]

# barcode plot of the result
pdf(file=file.path(outputDir,"Figures/Figure10/F10C_EScore_CameraIndex_BCRPath.pdf"), width = 6, height = 3)
barcodeplot(fit$t[,1],
            study_mSig[[3]]$SIG_BCR_SIGNALING_PATHWAY,
            design=EScoreMod,
            contrast=EScoreCon,
            main="Late vs Early EScore: SIG_BCR_SIGNALING_PATHWAY",
            xlab="moderated t-statistic")
dev.off()


##############################
# Figure 10D - PFS in Proteasome based therapies with EScore Levels 4/5
##############################

ClinDatProt$EScore<-ClinDatProt$EScore5W

pdf(file=file.path(outputDir,"Figures/Figure10/F10D_Proteasome_LateEScore_Surv.pdf"), width = 9, height = 6)
ggsurvplot(
  surv_fit(Surv(ClinDatProt$pfs1, event=ClinDatProt$pdtx1) ~ EScore, data=ClinDatProt),
  conf.int = FALSE,
  pval = TRUE,
  ylab="Progression Free Survival",
  risk.table=T,
  xlab="Time (Months)",
  risk.table.height=.35,
  title="Proteasome Inhibitor Based Therapies"
)
dev.off()

