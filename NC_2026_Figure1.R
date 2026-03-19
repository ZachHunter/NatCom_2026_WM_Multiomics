library(dplyr)
library(tidyr)
library(purrr)
library(bvt)
library(maftools)
library(pheatmap)
library(survminer)
library(NMF)
library(edgeR)
library(DESeq2)
library(kableExtra)

#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 1
#   Supplemental Tables 1-2
#
# Special note - Uses r package bvt for ploting
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

#mSigDb data prepared for camera/fry analysis indexed for WMExpressed
load(file.path(dataDir,"study_mSig.RData"))

#All of the final patient mutations in MAF format and imported into maftools
load("~/Desktop/DesktopData/CurrentProjects/300/WES/MAF.RData")

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

#Testing for clinical associations between groups
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
#Figure 1A
# Maftools Oncoplot from WES
##############################

#Filtering to include MYD88 mutated samples only and restrict genes to of interest
MYD88tsb<-as.character(unique(unlist(MAF@data[MAF@data$SYMBOL=="MYD88" & MAF@data$HGVSp_Short=="p.L273P","Tumor_Sample_Barcode"])))
MYD88MAF<-subsetMaf(MAF, tsb=MYD88tsb)
pdf(file=file.path(outputDir,"Figures/Figure1/F1A_Oncoplot.pdf"), width = 10, height = 5)
oncoplot(MYD88MAF, genes=c("MYD88","CD79B","CXCR4","ARID1A","BIRC3","EP300","BTG2","H1-4","TRAF2","CTBP1","NOTCH1","TP53","NDUFA7", "NOTCH2", "TNFAIP3"))
dev.off()

##############################
# Figure 1B
# Stratifying samples on the CXCR4 mutant gene signature
##############################


#Deriving CXCR4 DEG Signature

CXCR4mod<-model.matrix(~ factor(CXCR4) + agebmbx + gender, data=pData(studyCounts)[WMOnly,])
v<-voom(
  calcNormFactors(
    DGEList(counts=counts(studyCounts)[WMExpressed,WMOnly])
    ),
  design = CXCR4mod)
fit<-lmFit(v, design=CXCR4mod)
fit<-treat(fit)
CXCR4GeneSig<-topTreat(fit, coef=2, p.value=0.01, lfc=1, number=1000)

#Creating VST data for downstream analysis.
detemp<-DESeqDataSetFromMatrix(countData = counts(studyCounts), colData=pData(studyCounts), design=~CXCR4)
studyVST<-varianceStabilizingTransformation(detemp, blind = TRUE)

#Subset and format data for clustering
CXCR4SigVST<-assay(studyVST)[rownames(CXCR4GeneSig),WMOnly] %>%
  t() %>%
  scale()

#Clustering and visualization
CXCR4Cluster<-kmeans(CXCR4SigVST, centers = 2,nstart = 20)$cluster


CXCR4prcomp<-prcomp(CXCR4SigVST)
pdf(file=file.path(outputDir,"Figures/Figure1/F1B_CXCR4Cluster.pdf"), width = 6, height = 4)

geneScatter(CXCR4prcomp$x[,1:2],
            shape = factor(CXCR4Cluster, labels=c("Cluster 1", "Cluster 2")),
            color = colData(studyVST)[WMOnly,"CXCR4"],
            legend=c("CXCR4 Status","Cluster Assignment"),
            legendSize=1,
            lWidth=1.5,
            main="",axisLabelSize=1.3, yAxisLabSize=1.3,RSOverride=TRUE)
genePlot(CXCR4prcomp$x[,1:2],plotType="density", drawPoints=F, add=TRUE, legendSize=0.00001,RSOverride=TRUE)
dev.off()

##############################
# Figure 1C
# Heatmap of top PC1 genes
# Includes creation of supplemental table 1
##############################

#Selecting data for heatmap of top 50 PC1 genes
CXCR4HeatGenes<-CXCR4prcomp$rotation[order(abs(CXCR4prcomp$rotation[,1]),decreasing = TRUE),1][1:50]

#Note we are capturing the return value here so we can use cuttree to define the groups later
pdf(file=file.path(outputDir,"Figures/Figure1/F1C_SubtypeHeatmap.pdf"), width = 5, height = 4)
a<-pheatmap(scale(t(assay(studyVST)[names(CXCR4HeatGenes),WMOnly])),
            clustering_method = "ward.D",
            cutree_rows = 3, cutree_cols = 2,
            show_colnames = FALSE, show_rownames = FALSE,
            annotation_row = data.frame(
              CXCR4=factor(colData(studyVST)[WMOnly,"CXCR4"]),
              `Del Chr6q`=factor(colData(studyVST)[WMOnly,"chr6qdel"], labels=c("Deleted","Intact")),
              row.names = WMOnly),RSOverride=TRUE, scale="none")
dev.off()

#Assigning the new groups
ssub<-as.character(pData(studyTpM)$CXCR4)
SimpleSub<-factor(cutree(a$tree_row,3),labels=c("BCL","PCL","Early WM"))
names(ssub)<-sampleNames(studyTpM)
ssub[names(SimpleSub)]<-as.character(SimpleSub)
ssub<-factor(ssub,levels=c("HDPB","HDMB","Early WM","BCL","PCL"))

pData(studyCounts)$SimpleSubtype<-ssub
pData(studyTpM)$SimpleSubtype<-ssub
colData(studyVST)$SimpleSubtype<-ssub

#Formatting gene list for supplementary table 1
exprs(studyTpM)[names(CXCR4HeatGenes),] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Subtype=pData(studyTpM)$SimpleSubtype) %>%
  group_by(Subtype) %>%
  dplyr::summarise_all(median) %>%
  pivot_longer(cols=-1,names_to = "ID",values_to = "Median TpM") %>%
  pivot_wider(id_cols = ID,names_from = Subtype,values_from = "Median TpM") %>%
  mutate("Symbol"=fData(studyTpM)[ID,"GeneSymbol"]) %>%
  select(ID,Symbol,"HD CD19+CD27-"=HDPB,"HD CD19+CD27+"=HDMB,"Early WM","BCL","PCL") %>%
  kbl(digits = 2,caption="Genes used in subtype identification by hierarchical clustering") %>%
  kable_classic() %>%
  column_spec(1,italic =TRUE) %>%
  as_image(file = file.path(outputDir,"Figures/STables/Subtyped_Cluster_Genes.pdf"),width = 6)


##############################
# Figure 1D
# GSEA analysis with Camera of BCL vs. PCL
# Highlighting B-cell / Plasma cell comparisons in results
##############################

#Modeling BCL vs PCL using CXCR4, sex, and age of biopsy as covariates
subtypeMod<-model.matrix(~ 0 + factor(SimpleSubtype) + factor(CXCR4) + gender + agebmbx, data=pData(studyCounts)[WMOnly,])
colnames(subtypeMod)<-c( "EarlyWM","BCL", "PCL","CXCR4","Gender","Age")
BCLvPCLcon<-makeContrasts(PCL-BCL, levels=subtypeMod)

v<-voom(
  calcNormFactors(
    DGEList(counts=counts(studyCounts)[WMExpressed,WMOnly], genes=fData(studyCounts)[WMExpressed,c("Chr","GeneSymbol")])
  ), design=subtypeMod
)


camera(v, study_mSig[[3]], design = subtypeMod, contrast = BCLvPCLcon)[c(4,6,8),] %>%
  head(3) %>%
  kbl(digits=c(0,0,20,20,20)) %>%
  kable_classic(full_width=FALSE)  %>%
  column_spec(1,italic = TRUE) %>%
  as_image(file = file.path(outputDir,"Figures/Figure1/F1D_PCLvBCL_GSEA1.pdf"),width = 6)

camera(v, study_mSig[[8]], design = subtypeMod, contrast = BCLvPCLcon) %>%
  head(7) %>%
  kbl(digits=c(0,0,20,20,20)) %>%
  kable_classic(full_width=FALSE)  %>%
  column_spec(1,italic = TRUE) %>%
  as_image(file = file.path(outputDir,"Figures/Figure1/F1D_PCLvBCL_GSEA2.pdf"),width = 6)




##############################
# Figure 1E
# Time to first therapy analysis by WM Subtype
# Includes supplemental table 2  - Clinical Characterisics of Early WM
##############################


#Time to first therapy analysis by subtype

subtypeSuvData<-data.frame(
  pData(studyTpM)[WMOnly,c("TTFT","treated")],
  Subtype=pData(studyTpM)[WMOnly,"SimpleSubtype"],
  Subtyped=factor(
    pData(studyTpM)[WMOnly,"SimpleSubtype"]=="Early WM",
    labels=c("BCL/PCL","Early WM")
    ))

pdf(file=file.path(outputDir,"Figures/Figure1/F1D_TTFT_Subtype.pdf"), width = 9, height = 6)
ggsurvplot(
  surv_fit(
    Surv(subtypeSuvData$TTFT,subtypeSuvData$treated)~Subtype,
    data=subtypeSuvData),
  pval = TRUE,
  ylab="Time to First Therapy",
  risk.table=T,
  xlab="Time (Years)",
  conf.int = T)
dev.off()



##############
#Subtype clinical comparisons

EWM<-WMOnly[pData(studyTpM)[WMOnly,"SimpleSubtype"]=="Early WM"]
SubtypedWM<-WMOnly[pData(studyTpM)[WMOnly,"SimpleSubtype"]!="Early WM"]

BCL<-WMOnly[pData(studyTpM)[WMOnly,"SimpleSubtype"]=="BCL"]
PCL<-WMOnly[pData(studyTpM)[WMOnly,"SimpleSubtype"]=="PCL"]

EWMvSubtyped<-ClinTest(EWM, SubtypedWM)

EWMvSubtyped[[1]] %>%
  kbl(col.names = c(rep(c("median", "min", "max"),2),"p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  add_header_above(c("","Early WM (N=86)"=3,"Subtyped WM (N=163)"=3,"","")) %>%
  row_spec(which(EWMvSubtyped[[1]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/STables/EWM_vs_Subtyped_Clin1.pdf"),width = 6)

EWMvSubtyped[[2]] %>%
  kbl(col.names = c("Early WM","Subtyped WM","p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  row_spec(which(EWMvSubtyped[[2]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/STables/EWM_vs_Subtyped_Clin2.pdf"),width = 6)


BCLvPCL<-ClinTest(BCL,PCL)
BCLvPCL[[1]] %>%
  kbl(col.names = c(rep(c("median", "min", "max"),2),"p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  add_header_above(c("","BCL Subtyped WM (N=96)"=3,"PCL Subtyped WM (N=67)"=3,"","")) %>%
  row_spec(which(BCLvPCL[[1]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST2_BCL_vs_PCL_Clin1.png"),width = 6)

BCLvPCL[[2]] %>%
  kbl(col.names = c("BCL Subtyped WM (N=96)","PCL Subtyped WM (N=67)","p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  row_spec(which(BCLvPCL[[2]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST2_BCL_vs_PCL_Clin2.png"),width = 6)


ClinTest(EWM,PCL)
ClinTest(EWM,BCL)

rowData(studyVST)<-fData(studyTpM)[rownames(assay(studyVST)),]
save(studyVST, file=file.path(dataDir,"studyVST.RData"))
save(studyCounts, file=file.path(dataDir,"studyCounts.RData"))
save(studyTpM, file=file.path(dataDir,"studyTpM.RData"))
