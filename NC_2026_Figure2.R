library(dplyr)
library(tidyr)
library(purrr)
library(bvt)
library(survival)
library(survminer)
library(destiny)
library(broom)
library(NMF)
library(pheatmap)
library(lubridate)
library(kableExtra)

#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 2
#   Supplemental figure panel 1
#   Supplemental tables 4-5
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
dataDir<-"~/Desktop/DesktopData/Papers/Multiomics/Data/"
outputDir<-"~/Desktop/DesktopData/Papers/Multiomics"

#Loading eSets with relevant phenotype and feature data
# ExpressionSet of gene level TpM data from salmon.
load(file.path(dataDir,"studyTpM.RData"))
# SeqExpressionSet of Batch adjusted count data from STAR
load(file.path(dataDir,"studyCounts.RData"))
# variance stabilized transformation data from DESeq2 generated in Figure 1 code
load(file.path(dataDir,"studyVST.RData"))

#mSigDb data prepared for camera/fry analysis indexed for WMExpressed
load(file.path(dataDir,"study_mSig.RData"))

#All of the final patient mutations in MAF format and imported into maftools
load(file.path(dataDir,"MAF.RData"))

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


# Custom themes to make sure subtype colors remain the same across figures
WMSubtypeColor<-c(3,1,2,4:9)
WMHDSubtypeColor<-c(4,5,3,1,2,6:9)
SubtypeTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[WMSubtypeColor],fill=npDefaultTheme$plotColors$fill[WMSubtypeColor] ))
HDWMTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[WMHDSubtypeColor],fill=npDefaultTheme$plotColors$fill[WMHDSubtypeColor] ))

#Testing clinical associations between groups
ClinTest<-function(groupA,groupB) {
  continuousClin<-c("agebmbx","wbc","hb","igm","igg","iga","alymph","aneut", "bm" )
  categoricalClin<-c("b2m" ,"ipssbmbx", "familialType","Multiclonal","adenopathy","splenomegaly",  "cd5","cd10","cd23", "gender", "sxstatus" )

  cClin<-data.frame(
    medianA=map_dbl(continuousClin,function(x) median(pData(studyTpM)[groupA,x],na.rm=TRUE)),
    minA=map_dbl(continuousClin,function(x) min(pData(studyTpM)[groupA,x],na.rm=TRUE)),
    maxA=map_dbl(continuousClin,function(x) max(pData(studyTpM)[groupA,x],na.rm=TRUE)),
    medianB=map_dbl(continuousClin,function(x) median(pData(studyTpM)[groupB,x],na.rm=TRUE)),
    minB=map_dbl(continuousClin,function(x) min(pData(studyTpM)[groupB,x],na.rm=TRUE)),
    maxB=map_dbl(continuousClin,function(x) max(pData(studyTpM)[groupB,x],na.rm=TRUE)),
    p.value=map_dbl(continuousClin,function(x) wilcox.test(pData(studyTpM)[groupA,x],pData(studyTpM)[groupB,x])$p.value),
    adj.p.value=p.adjust(map_dbl(continuousClin,function(x) wilcox.test(pData(studyTpM)[groupA,x],pData(studyTpM)[groupB,x])$p.value),method="holm"),
    row.names = continuousClin)

  catClin<-data.frame(
    groupA=c(
      paste0(sum(pData(studyTpM)[groupA,"b2m"] > 3 ,na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupA,"b2m"]))," (",round(sum(pData(studyTpM)[groupA,"b2m"] >3,na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupA,"b2m"]))*100,1),"%)"),
      paste0(sum(pData(studyTpM)[groupA,"ipssbmbx"] == "High",na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupA,"ipssbmbx"]))," (",round(sum(pData(studyTpM)[groupA,"ipssbmbx"] == "High",na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupA,"ipssbmbx"]))*100,1),"%)"),
      paste0(sum(pData(studyTpM)[groupA,"familialType"] != "NH",na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupA,"familialType"]))," (",round(sum(pData(studyTpM)[groupA,"familialType"] != "NH",na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupA,"familialType"]))*100,1),"%)"),
      paste0(sum(pData(studyTpM)[groupA,"familialType"] == "WMH",na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupA,"familialType"]))," (",round(sum(pData(studyTpM)[groupA,"familialType"] == "WMH",na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupA,"familialType"]))*100,1),"%)"),
      map_chr(categoricalClin[4:11], function(x) {
        paste0(sum(pData(studyTpM)[groupA,x]==1,na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupA,x]))," (",round(sum(pData(studyTpM)[groupA,x] == 1,na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupA,x]))*100,1),"%)")
      })),
    groupB=c(
      paste0(sum(pData(studyTpM)[groupB,"b2m"] > 3 ,na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupB,"b2m"]))," (",round(sum(pData(studyTpM)[groupB,"b2m"] >3,na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupB,"b2m"]))*100,1),"%)"),
      paste0(sum(pData(studyTpM)[groupB,"ipssbmbx"] == "High",na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupB,"ipssbmbx"]))," (",round(sum(pData(studyTpM)[groupB,"ipssbmbx"] == "High",na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupB,"ipssbmbx"]))*100,1),"%)"),
      paste0(sum(pData(studyTpM)[groupB,"familialType"] != "NH",na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupB,"familialType"]))," (",round(sum(pData(studyTpM)[groupB,"familialType"] != "NH",na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupB,"familialType"]))*100,1),"%)"),
      paste0(sum(pData(studyTpM)[groupB,"familialType"] == "WMH",na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupB,"familialType"]))," (",round(sum(pData(studyTpM)[groupB,"familialType"] == "WMH",na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupB,"familialType"]))*100,1),"%)"),
      map_chr(categoricalClin[4:11], function(x) {
        paste0(sum(pData(studyTpM)[groupB,x]==1,na.rm=TRUE),"/",sum(!is.na(pData(studyTpM)[groupB,x]))," (",round(sum(pData(studyTpM)[groupB,x] == 1,na.rm=TRUE)/sum(!is.na(pData(studyTpM)[groupB,x]))*100,1),"%)")
      })),
    p.value=c(
      fisher.test(matrix(c(table(pData(studyTpM)[groupA,"b2m"] > 3), table(pData(studyTpM)[groupB,"b2m"] > 3)), nrow=2))$p.value,
      fisher.test(matrix(c(table(pData(studyTpM)[groupA,"ipssbmbx"] =="High"), table(pData(studyTpM)[groupB,"ipssbmbx"] =="High")), nrow=2))$p.value,
      fisher.test(matrix(c(table(pData(studyTpM)[groupA,"familialType"] !="NH"), table(pData(studyTpM)[groupB,"familialType"] !="NH")), nrow=2))$p.value,
      fisher.test(matrix(c(table(pData(studyTpM)[groupA,"familialType"] =="WMH"), table(pData(studyTpM)[groupB,"familialType"] =="WMH")), nrow=2))$p.value,
      map_dbl(categoricalClin[4:11], function(x) {
        fisher.test(matrix(c(sum(pData(studyTpM)[groupA,x] ==1,na.rm=TRUE),sum(pData(studyTpM)[groupA,x] !=1,na.rm=TRUE), sum(pData(studyTpM)[groupB,x] ==1,na.rm=TRUE),sum(pData(studyTpM)[groupB,x] !=1,na.rm=TRUE)), nrow=2))$p.value
      })),
    row.names=c("B2M>3", "High_IPSS","AnyFamilialHistory","WMFamilialHistory",categoricalClin[4:11])
  ) %>%
    mutate(adj.p.value=p.adjust(p.value, "fdr"))

  rownames(cClin)<-c("Age at Biopsy (Years)", "WBC (10E9/L)","Hgb (g/dL)","IgM (mg/dL)","IgG (mg/dL)", "IgA (mg/dL)", "ALymph (10E3/uL)", "ANeut (10E3/uL)", "LPC Intertrabecular Involvement (%)" )
  rownames(catClin)<-c("B2M (mg/L) > 3", "IPSS: High", "Family History of B-malignancies","Family History of WM","Multiclonal SPEP","Adenopathy (> 1.5cm)", "Splenomegaly (>15cm)","CD5+","CD10+","CD23+","Sex (male)","Symptomatic Disease")
  list(cClin,catClin)
}



##############################
# Figure 2A
# Non-Negative Matrix Factorization
##############################

NMFData<-assay(studyVST)[WMExpressed[order(rowVars(assay(studyVST)[WMExpressed,WMOnly]),decreasing = TRUE)[1:1000]],WMOnly]

set.seed(265)

NMFModel<-nmf(NMFData,3,nrun=20, method = "brunet")
NMFPlotting<-t(coef(NMFModel))
colnames(NMFPlotting)<-paste0("Metagene ",1:3)

pdf(file=file.path(outputDir,"Figures/Figure2/F2A_NMF_Heatmap.pdf"), width = 6, height = 5)
pheatmap(NMFPlotting,
         scale="none",
         clustering_method = "ward.D",
         annotation_row = data.frame(
           Subtype=pData(studyTpM)[WMOnly,"SimpleSubtype"],
           Symptomatic=factor(pData(studyTpM)[WMOnly,"sxstatus"],labels=c("Asymptomatic","Symptomatic")),
           row.names = WMOnly),
         show_rownames = FALSE,
         annotation_colors = list(
           Subtype=c("Early WM"=npDefaultTheme$plotColors$fill[3],
                     "BCL"=npDefaultTheme$plotColors$fill[1],
                     "PCL"=npDefaultTheme$plotColors$fill[2]),
           Symptomatic=c("Asymptomatic"=npDefaultTheme$plotColors$fill[4],"Symptomatic"=npDefaultTheme$plotColors$fill[5])
         ))
dev.off()


##############################
# Figure 2B
# Subtype correlation with Metagenes 2:3
##############################

pdf(file=file.path(outputDir,"Figures/Figure2/F2B_NMF_Subtype.pdf"), width = 6, height = 5)
pheatmap(NMFPlotting[,2:3],
         annotation_row = data.frame(Subtype=colData(studyVST)[WMOnly,"SimpleSubtype"]),
         clustering_method = "ward.D",
         show_rownames = FALSE,
         scale="row",
         annotation_colors = list(
           Subtype=c("Early WM"=npDefaultTheme$plotColors$fill[3],
                     "BCL"=npDefaultTheme$plotColors$fill[1],
                     "PCL"=npDefaultTheme$plotColors$fill[2]))
)
dev.off()

##############################
# Figure 2C
# 3D Rendering of Metagene Space
# Includes metagene gene signature extraction
# Note that genes are further filtered to ensure consistent FC direction and equal numbers between PCL and BCL
# Final gene lists can be seen in supplemental table 5
##############################

#Note that this generates an interactive RLG object
#Images captured manually
geneScatter(NMFPlotting,
            useRgl=TRUE,
            size=pData(studyTpM)[WMOnly,"bm"],
            color=pData(studyTpM)[WMOnly,"SimpleSubtype"],
            theme=SubtypeTheme,legend=c("WM Subtype","BM (%)"),
            pointSize=.75,
            main="",
            RSOverride=TRUE)

#Extract Features from NMF Model
#Note that metagene 1 throws an error with the default Kim method
NMFFeatures<-extractFeatures(NMFModel)
NMFFeatures[[1]]<-extractFeatures(NMFModel,"max")[[1]]
NMFFeatures<-map(NMFFeatures,function(g) rownames(NMFData[g,]))
names(NMFFeatures)<-c("EScore","BCL","PCL")

maxL<-map_dbl(NMFFeatures,length) %>% max()
FeatureSet=data.frame(rep("",maxL),rep("",maxL),rep("",maxL))
for(i in 1:3){
  FeatureSet[seq(length(NMFFeatures[[i]])),i]<-unlist(fData(studyTpM)[NMFFeatures[[i]],"GeneSymbol"])
}
colnames(FeatureSet)<-names(NMFFeatures)
kbl(FeatureSet) %>%
  kable_paper(full_width=FALSE) %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST5_CandidateSignatureGenes.pdf"))

#Just for convenience as ENSGs are easier to track
FeatureSetENST<-FeatureSet
for(i in 1:3){
  FeatureSetENST[seq(length(NMFFeatures[[i]])),i]<-unlist(NMFFeatures[[i]])
}


##############################
# Figure 2D
# Diffusion map of to 1K high variance genes
# Includes diffusion pseudotime analysis
# EScore introduced as study variable
##############################

set.seed(92820)
dmapData<-t(assay(studyVST)[WMExpressed[order(rowVars(assay(studyVST)[WMExpressed,WMOnly]),decreasing = TRUE)[1:1000]],WMOnly])
dmap<-DiffusionMap(dmapData)
dmapDpt<-DPT(dmap)

pdf(file=file.path(outputDir,"Figures/Figure2/F2D_DMAP_DPT.pdf"), width = 6, height = 5)
geneScatter(
  data.frame(DC1=dmapDpt$DC1,DC2=dmapDpt$DC2),
  color=colData(studyVST)[WMOnly,"SimpleSubtype"],
  size=round(dmapDpt$dpt,4),
  legend=c("Subtype","DPT"),
  theme=SubtypeTheme,
  pointSize=.75,
  legendSize=1.1, main="",labelSize=1,
  RSOverride=TRUE)
dev.off()

#Record EScore
FactorTemplate<-as.character(pData(studyTpM)$SimpleSubtype)
names(FactorTemplate)<-sampleNames(studyTpM)
FactorTemplate[WMOnly]<-as.character(factor(dmapDpt$DC1<0, labels=c("Late EScore","Early EScore")))
EScoreEL<-factor(FactorTemplate, levels=c("HDPB","HDMB","Early EScore","Late EScore"))
pData(studyTpM)$EScoreEL<-EScoreEL
EScore<-rep(NA,length(FactorTemplate))
names(EScore)<-names(FactorTemplate)
EScore[WMOnly]<-dmapDpt$dpt
pData(studyTpM)$EScore<-EScore
pData(studyCounts)<-pData(studyTpM)


##############################
# Supplementary Figures 1A-D
# PCA of 1K High Variance Genes by WM Subtype and EScore
# Correlation of EScore with DC1, PC1, and Metagene 1
# Early/Late EScore in DPT space analysis
# Demonstration how NMF metagenes are restricted to a plane in 3D space
##############################


#Making Supplementary Figure 1A
#Note that this is made using a screenshot of RGL output
pca1<-prcomp(scale(dmapData))
geneScatter(pca1$x[,1:3],
            color=pData(studyTpM)[WMOnly,"SimpleSubtype"],
            size=round(pData(studyTpM)[WMOnly,"EScore"],4),
            pointSize=.75,
            theme=SubtypeTheme,
            useRgl=TRUE,
            main="",
            legend=c("WM Subtype","EScore"),
            RSOverride=TRUE)


# Making Supplementary Figure 1B
pdf(file=file.path(outputDir,"Figures/SFigure1/SF1B_EScore_vs_NMF.pdf"), width = 5, height = 5)

a<-cor.test(pData(studyTpM)[WMOnly,"EScore"],coef(NMFModel)[1,])
geneScatter(data.frame(EScore=pData(studyTpM)[WMOnly,"EScore"],Metagene1=coef(NMFModel)[1,]),
            trendline=TRUE,
            main="",
            ylab="Metagene 1",
            sub=paste0("Pearson's r = ",round(a$estimate,4),"; p.value ",format.pval(a$p.value)),
            RSOverride=TRUE)
dev.off()

pdf(file=file.path(outputDir,"Figures/SFigure1/SF1B_EScore_vs_PCA.pdf"), width = 5, height = 5)
a<-cor.test(pData(studyTpM)[WMOnly,"EScore"],pca1$x[,1])
geneScatter(data.frame(EScore=pData(studyTpM)[WMOnly,"EScore"],PCA1=pca1$x[,1]),
            trendline=TRUE,
            main="",
            ylab="Principal Component 1",
            sub=paste0("Pearson's r = ",round(a$estimate,4),"; p.value ",format.pval(a$p.value)),
            RSOverride=TRUE)
dev.off()

pdf(file=file.path(outputDir,"Figures/SFigure1/SF1B_EScore_vs_DMAP.pdf"), width = 5, height = 5)
a<-cor.test(pData(studyTpM)[WMOnly,"EScore"],dmapDpt$DC1)
geneScatter(data.frame(EScore=pData(studyTpM)[WMOnly,"EScore"],DC1=dmapDpt$DC1),
            trendline=TRUE,
            main="",
            ylab="Diffusion Component 1",
            sub=paste0("Pearson's r = ",round(a$estimate,4),"; p.value ",format.pval(a$p.value)),
            RSOverride=TRUE)
dev.off()


# Making Supplementary Figure 1C - Early / Late EScore
pdf(file=file.path(outputDir,"Figures/SFigure1/SF1C_DPT_EL.pdf"), width = 6, height = 5)

geneScatter(data.frame(DPT1=dmapDpt$DPT1,DPT2=dmapDpt$DPT2),
            color=pData(studyTpM)[WMOnly,"EScore_EL"],
            legendSize=1.1,
            size=round(pData(studyTpM)[WMOnly,"EScore"],4),
            pointSize=.75,
            main="",
            legend=c("Early vs. Late","EScore/DPT"),
            RSOverride=TRUE)
dev.off()


#Supplemental Figure 1D
#Note that this generates an interactive RLG object
#Images captured manually
geneScatter(NMFPlotting,
            useRgl=TRUE,
            size=pData(studyTpM)[WMOnly,"bm"],
            color=pData(studyTpM)[WMOnly,"SimpleSubtype"],
            theme=SubtypeTheme,legend=c("WM Subtype","BM (%)"),
            pointSize=.75,
            main="",
            RSOverride=TRUE)


##############################
# Figure 2E
# Time to first therapy from biopsy by EScore
##############################

EScoreSuvData<-data.frame(
  pData(studyTpM)[WMOnly,c("TTFT","DTFT","treated")],
  EScore<-factor(pData(studyTpM)[WMOnly,"EScoreEL"])
)

pdf(file=file.path(outputDir,"Figures/Figure2/F2E_TTFT_EScore.pdf"), width = 9, height = 6)
ggsurvplot(
  surv_fit(
    Surv(EScoreSuvData$TTFT,EScoreSuvData$treated)~EScoreSuvData$EScore,
    data=EScoreSuvData),
  pval = TRUE,
  ylab="Time to First Therapy from Study Biopsy",
  risk.table=T,
  xlab="Time (Years)",
  conf.int = T)
dev.off()

##############################
# Figure 2F
# Time to first therapy from diagnosis by EScore
##############################

pdf(file=file.path(outputDir,"Figures/Figure2/F2F_DTFT_EScore.pdf"), width = 9, height = 6)
ggsurvplot(
  surv_fit(
    Surv(EScoreSuvData$DTFT,EScoreSuvData$treated)~EScoreSuvData$EScore,
    data=EScoreSuvData),
  pval = TRUE,
  ylab="Time to First Therapy From Diagnosis",
  risk.table=T,
  xlab="Time (Years)",
  conf.int = T)
dev.off()

##############################
# Figure 2G
# Cox PH time to first therapy from biopsy by EScore + sex + age + Early WM status
##############################

a<-coxph(
  Surv(pData(studyTpM)[WMOnly,"TTFT"], pData(studyTpM)[WMOnly,"treated"]) ~
    agebmbx + gender + factor(EScoreEL) + factor(SimpleSubtype=="Early WM",labels=c("BCL/PCL","Early WM")),
  data=pData(studyTpM)[WMOnly,]) %>%
  tidy()
a[,1]<-c("Age at Biopsy (Years)","Sex", "Late EScore","Early WM")
colnames(a)[1]<-""
a %>%
  kbl(caption="Time to First Therapy from Study Biopsy (Cox PH)") %>%
  kable_paper(full_width = FALSE) %>%
  column_spec(1,italic = TRUE) %>%
  column_spec(2,color=ifelse(a$estimate<0,"red","black")) %>%
  as_image(file = file.path(outputDir,"Figures/Figure2/F2G_COXPH_TTFT.pdf"),width = 6)


##############################
# Figure 2H
# Cox PH time to first therapy from diagnosis by EScore + sex + age + Early WM status
##############################

b<-coxph(
  Surv(pData(studyTpM)[WMOnly,"DTFT"], pData(studyTpM)[WMOnly,"treated"]) ~
    agewm + gender + factor(EScoreEL) + factor(SimpleSubtype=="Early WM",labels=c("BCL/PCL","Early WM")),
  data=pData(studyTpM)[WMOnly,]) %>%
  tidy()
b[,1]<-c("Age at Diagnosis (Years)","Sex", "Late EScore" ,"Early WM")
colnames(b)[1]<-""
b %>%
  kbl(caption="Time to First Therapy from Diagnosis (Cox PH)") %>%
  kable_paper(full_width = FALSE) %>%
  column_spec(1,italic = TRUE) %>%
  column_spec(2,color=ifelse(b$estimate<0,"red","black")) %>%
  as_image(file = file.path(outputDir,"Figures/Figure2/F2H_COXPH_DTFT.pdf"),width = 6)

##############################
# Supplemental Table 4
# Clinical associates between Early and Late EScore
##############################

EarlyES<-WMOnly[pData(studyTpM)[WMOnly,"EScore_EL"]=="Early"]
LateES<-WMOnly[pData(studyTpM)[WMOnly,"EScore_EL"]=="Late"]

ELClin<-ClinTest(EarlyES,LateES)

ELClin[[1]] %>%
  kbl(col.names = c(rep(c("median", "min", "max"),2),"p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  add_header_above(c("","Early EScore (N=148)"=3,"Late EScore (N=101)"=3,"","")) %>%
  row_spec(which(ELClin[[1]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST4_EScore_Early_vs_Late_Clin1.pdf"),width = 6)

ELClin[[2]] %>%
  kbl(col.names = c("Early EScore (N=148)","Late EScore (N=101)","p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  row_spec(which(ELClin[[2]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST4_EScore_Early_vs_Late_Clin2.pdf"),width = 6)


save(dmapDpt,file=file.path(dataDir,"studyDmap.RData"))
save(NMFModel,file=file.path(dataDir,"NMF.RData"))
save(studyCounts, file=file.path(dataDir,"studyCounts.RData"))
save(studyTpM, file=file.path(dataDir,"studyTpM.RData"))
