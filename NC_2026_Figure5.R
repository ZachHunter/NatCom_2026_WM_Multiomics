library(dplyr)
library(tidyr)
library(purrr)
library(bvt)
library(edgeR)
library(broom)
library(tibble)
library(lubridate)
library(survival)
library(survminer)
library(kableExtra)
library(DESeq2)
library(maftools)
library(umap)

#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 5
#   Supplemental Figure panel 3
#   EBCL and EPCL DGE data
#   Supplemental Tables 6-7
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

# Location of GATK CNV analysis files
CNADir<-"~/Path/to/Directory/WES/"

#Loading eSets with relevant phenotype and feature data
# ExpressionSet of gene level TpM data from salmon.
load(file.path(dataDir,"studyTpM.RData"))
# SeqExpressionSet of Batch adjusted count data from STAR
load(file.path(dataDir,"studyCounts.RData"))
# variance stabilized transformation data from DESeq2 generated in Figure 1 code
load(file.path(dataDir,"studyVST.RData"))

# mSigDb data prepared for camera/fry analysis indexed for WMExpressed
load(file.path(dataDir,"study_mSig.RData"))

# All of the final patient mutations in MAF format.
load("~/Path/to/Directory/WES/MAF.RData")

# High Quality % CpG Methylation data from MethylKit
load(file.path(dataDir,"studyCpG.RData"))


##############################
# Custom Functions
# Includes style definitions and global settings
##############################

#Creating custom theme to make sure subtype colors remain the same across figures
WMSubtypeColor<-c(3,1,2,4:9)
EWMSubtypeColor<-c(3,8,1,2,9,5:7)
WMHDSubtypeColor<-c(4,5,3,1,2,6:9)
EWMHDSubtypeColor<-c(4,5,3,8,1,2,9,6:7)
SubtypeTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[WMSubtypeColor],fill=npDefaultTheme$plotColors$fill[WMSubtypeColor] ))
HDWMTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[WMHDSubtypeColor],fill=npDefaultTheme$plotColors$fill[WMHDSubtypeColor] ))
EHDWMTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[EWMHDSubtypeColor],fill=npDefaultTheme$plotColors$fill[EWMHDSubtypeColor] ))
ESubtypeTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[EWMSubtypeColor],fill=npDefaultTheme$plotColors$fill[EWMSubtypeColor] ))

# Initial guesses for the Extreme BCL and Extreme PCL subsets
EBCL<-c("WM005", "WM024", "WM072", "WM185", "WM228", "WM248", "WM252", "WM270", "WM308", "WM350", "WM354", "WM357", "WM045")
EPCL<-c("WM025", "WM029", "WM035", "WM044", "WM076", "WM111", "WM162", "WM204", "WM244", "WM132", "WM314", "WM367")

# Sample IDS for WM patients only
WMOnly<-sampleNames(studyCounts)[grep("WM", sampleNames(studyCounts))]

# Vector of ENSG with at least some expression in at least 20 patients
# Counts used to determine power to detect differences, TpM used for biological relevance
WMExpressed<-featureNames(studyCounts)[rowSums(counts(studyCounts)[,WMOnly]>10)>20]
WMExpressed<-WMExpressed[rowSums(exprs(studyTpM)[WMExpressed,WMOnly] > 1) > 20]

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
# Figure 5A
# Forest plot of Mutations associated with BCL and PCL
# Also includes NOTCH pathway enrichment analysis mentioned in text
##############################

WES_Samples<-gsub("_Tumor","",as.character(MAF@variants.per.sample$Tumor_Sample_Barcode))

# Defining our populations of interest
BCLsamps<-WMOnly[pData(studyTpM)[WMOnly,"SimpleSubtype"]=="BCL" & WMOnly %in% WES_Samples]
PCLsamps<-WMOnly[pData(studyTpM)[WMOnly,"SimpleSubtype"]=="PCL" & WMOnly %in% WES_Samples]
BCLmaf<-subsetMaf(MAF,tsb=paste0(BCLsamps,"_Tumor"),)
PCLmaf<-subsetMaf(MAF,tsb=paste0(PCLsamps,"_Tumor"))

BCLpath<-maftools::pathways(BCLmaf)
PCLpath<-maftools::pathways(PCLmaf)

# First we will conduct the pathway enrichment analysis referenced in the manuscript
# Note we are really just interested in NOTCH based on EP300 and NOTCH mutation rates
data.frame(Pathway=BCLpath$Pathway,
           'BCL'=map_chr(BCLpath$Pathway, function(p) paste0(BCLpath$Mutated_samples[BCLpath$Pathway==p],"/",length(BCLsamps)," (",round(BCLpath$Mutated_samples[BCLpath$Pathway==p]/length(BCLsamps)*100,1),"%)")),
           'PCL'=map_chr(BCLpath$Pathway, function(p) paste0(PCLpath$Mutated_samples[PCLpath$Pathway==p],"/",length(PCLsamps)," (",round(PCLpath$Mutated_samples[PCLpath$Pathway==p]/length(PCLsamps)*100,1),"%)")),
           p.value=p.adjust(map_dbl(BCLpath$Pathway,function(p) fisher.test(matrix(c(BCLpath$Mutated_samples[BCLpath$Pathway==p],length(BCLsamps)-BCLpath$Mutated_samples[BCLpath$Pathway==p],PCLpath$Mutated_samples[PCLpath$Pathway==p],length(PCLsamps)-PCLpath$Mutated_samples[PCLpath$Pathway==p]),nrow=2))$p.value),method="none")) %>%
  knitr::kable(caption="Oncogenic Pathway Enrichment by Subtype")

pdf(file=file.path(outputDir,"Figures/Figure5/F5A_BCP_PCL_Forestplot.pdf"), width = 7, height = 4.5)
forestPlot(
  mafCompare(BCLmaf,PCLmaf,
             m1Name = "BCL",
             m2Name = "PCL"))
dev.off()


##############################
# Supplemental Figures 3A
# Positional GSEA based on chr6q deletions status in BCL
##############################

# Subsetting our population by WM BCL Subtype
studyBCL<-studyCounts[,pData(studyCounts)$SimpleSubtype=="BCL"]

# Setting up model and contrast
chr6qBCLMod<-model.matrix(~ 0 + chr6qdel + EScore +gender + agebmbx, data=pData(studyBCL))
colnames(chr6qBCLMod)<-c("del","intact","EScore","sex","age")
chr6qBCLCon<-makeContrasts(del-intact,levels=chr6qBCLMod)

# Running voom
vBCL<-voom(
  calcNormFactors(
    DGEList(counts=counts(studyBCL)[WMExpressed,rownames(chr6qBCLMod)])
    ),
  design = chr6qBCLMod)

# Camera positional enrichment for BCL
BCL6qEnrichment<-camera(vBCL, study_mSig[[2]], design = chr6qBCLMod, contrast =chr6qBCLCon ) %>% rownames_to_column(var = "Cytoband") %>%
  dplyr::filter(FDR<0.05, grepl("chr6q", Cytoband))


# Since we plan to filter on unique findings, need to do the same calculation for PCL
# Subsetting our population by WM PCL Subtype
studyPCL<-studyCounts[,pData(studyCounts)$SimpleSubtype=="PCL"]

# Setting up model and contrast
chr6qPCLMod<-model.matrix(~ 0 + chr6qdel + EScore + gender + agebmbx, data=pData(studyPCL))
colnames(chr6qPCLMod)<-c("del","intact","EScore","sex","age")
chr6qPCLCon<-makeContrasts(del-intact,levels=chr6qPCLMod)

# Running voom
vPCL<-voom(
  calcNormFactors(
    DGEList(
      counts=counts(studyPCL)[WMExpressed,rownames(chr6qPCLMod)])
  ),
  design = chr6qPCLMod)

# Camera positional enrichment for PCL
PCL6qEnrichment<-camera(vPCL, study_mSig[[2]], design = chr6qPCLMod, contrast =chr6qPCLCon ) %>%
  rownames_to_column(var = "Cytoband") %>%
  dplyr::filter(FDR<0.05, grepl("chr6q", Cytoband))

#Final table
BCL6qEnrichment %>%
  kbl() %>%
  kable_paper(full_width=FALSE) %>%
  row_spec(which(BCL6qEnrichment$Cytoband %in% PCL6qEnrichment$Cytoband), background = "lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/sFigure3/sF3A_chr6qSig_BCL.pdf"),width = 6)


##############################
# Supplemental Figures 3B
# Positional GSEA based on chr6q deletions status in PCL
##############################

#Final table
PCL6qEnrichment %>%
  kbl() %>%
  kable_paper(full_width=FALSE) %>%
  row_spec(which(PCL6qEnrichment$Cytoband %in% BCL6qEnrichment$Cytoband), background = "lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/sFigure3/sF3B_chr6qSig_PCL.pdf"),width=6)


##############################
# Figure 5B
# Comparison of relative chr6q coverage between BCL and PCL in 6q deleted patients
##############################

# Defining our study populations
BCL6q<-sampleNames(studyTpM[,pData(studyTpM)$chr6qdel=="del" & pData(studyTpM)$SimpleSubtype=="BCL"])
PCL6q<-sampleNames(studyTpM[,pData(studyTpM)$chr6qdel=="del" & pData(studyTpM)$SimpleSubtype=="PCL"])

#Removing NAs and other artifacts
BCL6q<-BCL6q[grep("WM",BCL6q)]
PCL6q<-PCL6q[grep("WM",PCL6q)]

#loading called segment files from GATK CNV analysis (BCL)
BCL6qData<-vector(mode="list",length=length(BCL6q))
names(BCL6qData)<-BCL6q
for(d in BCL6q) {
  cFile<-list.files(path=CNADir,pattern = paste0(d,".called.seg"), recursive = T,full.names = T)
  BCL6qData[[d]]<-read.table(cFile, header=T, comment.char = "@") %>%
    filter(CONTIG=="chr6", END>=59800000, NUM_POINTS_COPY_RATIO>4)
}

#loading called segment files from GATK CNV analysis (PCL)
PCL6qData<-vector(mode="list",length=length(PCL6q))
names(PCL6qData)<-PCL6q
for(d in PCL6q) {
  cFile<-list.files(path=CNADir,pattern = paste0(d,".called.seg"), recursive = T,full.names = T)
  PCL6qData[[d]]<-read.table(cFile, header=T, comment.char = "@") %>%
    filter(CONTIG=="chr6", END>=59800000, NUM_POINTS_COPY_RATIO>4)
}

#BP start postition on Chr6 for analysis
cPos<-59800000
#step size in bp
StepSize<-10000
#Number of steps to cover the remainder of chr6q
nTotal<-trunc((170805979-59800000)/StepSize)

#Making sure to draw each point in the middle of the window
posLoc<-map_dbl(seq(nTotal), function(x) cPos+StepSize*(x-1) + StepSize/2)

#Calculating mean coverage over each window for BCL samples
aggBCLDat<-vector(mode="list",length=nTotal)
for(i in seq(nTotal)){
  aggBCLDat[[i]]<-map_dbl(BCL6qData, function(p) {
    tDat<-filter(p,START<=cPos+StepSize*(i-1) & END>=cPos+StepSize*(i-1))
    if(nrow(tDat)==1){
      as.numeric(tDat$MEAN_LOG2_COPY_RATIO)
    } else {
      apply(tDat,1,function(x) {
        if(as.numeric(x[3])>cPos+StepSize*i) {
          as.numeric(x[5])*((cPos+StepSize*i-as.numeric(x[2]))/StepSize)
        } else {
          as.numeric(x[5])*(as.numeric(x[3])-cPos-StepSize*(i-1))/StepSize
        }
      }) %>% unlist() %>% sum()
    }
  })
}

#Calculating mean coverage over each window for PCL samples
aggPCLDat<-vector(mode="list",length=nTotal)
for(i in seq(nTotal)){
  aggPCLDat[[i]]<-map_dbl(PCL6qData, function(p) {
    tDat<-filter(p,START<=cPos+StepSize*(i-1) & END>=cPos+StepSize*(i-1))
    if(nrow(tDat)==1){
      as.numeric(tDat$MEAN_LOG2_COPY_RATIO)
    } else {
      apply(tDat,1,function(x) {
        if(as.numeric(x[3])>cPos+StepSize*i) {
          as.numeric(x[5])*((cPos+StepSize*i-as.numeric(x[2]))/StepSize)
        } else {
          as.numeric(x[5])*(as.numeric(x[3])-cPos-StepSize*(i-1))/StepSize
        }
      }) %>% unlist() %>% sum()
    }
  })
}

chr6qPval<-map_dbl(seq(nTotal), function(i) wilcox.test(aggBCLDat[[i]],aggPCLDat[[i]])$p.value)

#Data prep complete. Starting to plot.
pdf(file=file.path(outputDir,"Figures/Figure5/F5B_Chr6qCov.pdf"), width = 6, height = 4)

#Making the plotting environment
geneScatter(BCL6qData[[1]],
            genes =c("START","MEAN_LOG2_COPY_RATIO"),
            xLim=c(59800000,170805979),yLim=c(-1.5,.5),
            pointSize=0.000001, guides=F,
            ylab="Mean Log 2 Copy Ratio",
            xlab="Chromosome 6q Position (BP)",
            rotateLabels=T,main="",
            RSOverride=TRUE)

abline(h=0,lwd=1.5,col="navy")

#Drawing the smoothed splines of coverage with confidence intervals for BCL Data
cSpar<-.2
y1<-smooth.spline(
  x=posLoc,
  y=map_dbl(aggBCLDat, function(x) mean(x)+2*sd(x)/sqrt(length(BCL6q))),
  spar = cSpar)
y2<-smooth.spline(
  x=posLoc,
  y=map_dbl(aggBCLDat, function(x) mean(x)-2*sd(x)/sqrt(length(BCL6q))),
  spar = cSpar)

polygon(x=c(y1$x,rev(y2$x)),y=c(y1$y,rev(y2$y)), col=setAlpha("#E41A1C99"),lwd=0.1)
m1<-smooth.spline(
  posLoc,
  map_dbl(aggBCLDat,mean),
  spar = cSpar)

points(m1,type="l",lwd=2.5, col=setAlpha("#E41A1C99",.5))

#Drawing the smoothed splines of coverage with confidence intervals for PCL Data
y1<-smooth.spline(
  x=posLoc,
  y=map_dbl(aggPCLDat, function(x) mean(x)+2*sd(x)/sqrt(length(PCL6q))),
  spar = cSpar)
y2<-smooth.spline(
  x=posLoc,
  y=map_dbl(aggPCLDat, function(x) mean(x)-2*sd(x)/sqrt(length(PCL6q))),
  spar = cSpar)

polygon(x=c(y1$x,rev(y2$x)),y=c(y1$y,rev(y2$y)), col=setAlpha("#377EB899"),lwd=0.1)
m1<-smooth.spline(
  posLoc,
  map_dbl(aggPCLDat,mean),
  spar = cSpar)

points(m1,type="l",lwd=2.5, col=setAlpha("#377EB899",.5))

#Drawing fdr adjusted p-values as a rug plot
rug(x=posLoc[p.adjust(chr6qPval,method = "fdr")>0.05], col=setAlpha("grey"))
rug(x=posLoc[p.adjust(chr6qPval,method = "fdr")<=0.05], col=setAlpha("red"))
legend("topright",legend=c("BCL Subtyped WM", "PCL Subtyped WM"), fill = c(setAlpha("#E41A1C99",.75),setAlpha("#377EB899",.75)),bty = "n")
dev.off()


##############################
# Figure 5C
# Umap of high quality CpG data
##############################

set.seed(481)
umapCPG<-umap(t(scale(studyCpG)))

# preparing annotation
SelectionType<-rep("CD19+", ncol(studyCpG))
CellType<-rep("IgM MM", ncol(studyCpG))

SelectionType[grep("PC", colnames(studyCpG))]<-"CD138+"
CellType[grep("HDBMPC", colnames(studyCpG))] <- "HD PC"
CellType[grep("HDMB", colnames(studyCpG))] <- "HD MB"
CellType[grep("HDPB", colnames(studyCpG))] <- "HD PB"
CellType[grep("WM", colnames(studyCpG))] <- as.character(pData(studyTpM)[gsub("PC","",colnames(studyCpG)[grep("WM", colnames(studyCpG))]),"SimpleSubtype"])

SelectionType<-factor(SelectionType, levels=c("CD19+", "CD138+"))
CellType<-factor(CellType, levels=c("HD PB", "HD MB", "HD PC", "Early WM", "BCL", "PCL" ,"IgM MM"))

# Plotting data
pdf(file=file.path(outputDir,"Figures/Figure5/F5C_CpGumap.pdf"), width = 6, height = 4)
geneScatter(umapCPG$layout,
            color=CellType,
            shape=SelectionType,
            legend=c("Cell Type", "Selection"),
            plotColors=list(points=npDefaultTheme$plotColors$points[c(4,5,7,3,1,2,8)]),
            legendSize=1.1, pointSize=1,
            xlab="UMAP 1", ylab="UMAP 2",
            main="", RSOverride=TRUE)
dev.off()


##############################
# Figure 5D
# Umap of top 500 high variance genes
##############################

set.seed(338)
top500<-order(rowVars(assay(studyVST)[,WMOnly]),decreasing=TRUE)[1:1000]
a<-umap(scale(t(assay(studyVST)[top500,WMOnly])))

# Samples clusters identified at the extremes of the plot listed in variables
# EBCL and EPCL. These are starting guesses from exploratory analysis.
# We will try to refine these observations and test for significant associations

# Making a new subtype factor including the extreme observations
testTypes<-as.character(pData(studyTpM)[WMOnly,"SimpleSubtype"])
names(testTypes)<-WMOnly
testTypes[EBCL]<-"EBCL"
testTypes[EPCL]<-"EPCL"
extremeTypes<-factor(testTypes,levels=c("Early WM","EBCL","BCL","PCL","EPCL"))

#Next we perform DGE analysis to see if they form distinct clusters and if the group membership should be further refined
extremeModel<-model.matrix(~0 + extremeTypes + EScore + factor(CXCR4) + gender + agebmbx, data=pData(studyCounts)[WMOnly,])
colnames(extremeModel)<-c("EarlyWM","BCL","PCL","BCE","PCE","EScore","CXCR4","gender","age")
BCEcon<-makeContrasts(BCE-BCL, levels=extremeModel)
PCEcon<-makeContrasts(PCE-PCL, levels=extremeModel)

v<-voom(
  calcNormFactors(
    DGEList(
      counts=counts(studyCounts)[WMExpressed,WMOnly],
      genes=fData(studyCounts)[WMExpressed,c("Chr","GeneSymbol")]
    )
  )
  ,design=extremeModel)

#using treat instead of eBayes here to prioritize lfc differences for clustering
fit<-lmFit(v,design = extremeModel)
BCEfit<-contrasts.fit(fit,BCEcon)
PCEfit<-contrasts.fit(fit,PCEcon)
BCEfit<-treat(BCEfit)
PCEfit<-treat(PCEfit)
BCEtop<-topTreat(BCEfit, lfc=1, p.value=0.05, number=Inf)
PCEtop<-topTreat(PCEfit, lfc=1, p.value=0.05, number=Inf)

#principal component analysis within PCL and BCL, respectively
BCLextremeSort<-prcomp(
  scale(
    t(assay(studyVST)[rownames(BCEtop),colData(studyVST)$SimpleSubtype=="BCL"])))

PCLextremeSort<-prcomp(
  scale(
    t(assay(studyVST)[rownames(PCEtop),colData(studyVST)$SimpleSubtype=="PCL"])))

geneScatter(
  BCLextremeSort$x[,1:2],
  color=extremeTypes[rownames(BCLextremeSort$x)],
  size=pData(studyTpM)[rownames(BCLextremeSort$x),"EScore"])

geneScatter(
  PCLextremeSort$x[,1:2],
  color=extremeTypes[rownames(PCLextremeSort$x)],
  size=pData(studyTpM)[rownames(PCLextremeSort$x),"EScore"])

#The analysis above has us reevaluate our BCE and PCE assignments
extremeTypesfinal<-as.character(pData(studyCounts)$SimpleSubtype)
names(extremeTypesfinal)<-sampleNames(studyCounts)
extremeTypesfinal[names(extremeTypes)]<-as.character(extremeTypes)
extremeTypesfinal[c("WM252","WM045")]<-"BCL"
extremeTypesfinal[c("WM254","WM040")]<-"EBCL"
extremeTypesfinal[c("WM132","WM204", "WM367","WM244")]<-"PCL"
extremeTypesfinal[c("WM034","WM163","WM114")]<-"EPCL"

#Soring Expanded Subtype in the eSets
pData(studyCounts)$ExpandedSubtype<-factor(extremeTypesfinal, levels=c("HDPB","HDMB","Early WM","EBCL", "BCL", "PCL", "EPCL"))
pData(studyTpM)$ExpandedSubtype<-factor(extremeTypesfinal, levels=c("HDPB","HDMB","Early WM","EBCL", "BCL", "PCL", "EPCL"))

# plotting final assignments in initial umap
pdf(file=file.path(outputDir,"Figures/Figure5/F5E_ExtremeSubtype_umap.pdf"), width = 7, height = 6)
geneScatter(a$layout,
            color=pData(studyCounts)[rownames(a$layout),"ExpandedSubtype"],
            size=round(pData(studyTpM)[WMOnly,"EScore"],4),
            legend=c("Subtype", "EScore"), main="",
            xlab="UMAP 1", ylab="UMAP 2", legendSize=1.1,
            theme=ESubtypeTheme,
            pointSize=.3, RSOverride=TRUE)
dev.off()


##############################
# Figure 5D
# Supplementary tables 6-7
##############################

# Running clinical association testing between EBCL and BCL as well as EPCL and PCL
BCEvBCL_Clin<-ClinTest(names(extremeTypesfinal)[extremeTypesfinal=="EBCL"],names(extremeTypesfinal)[extremeTypesfinal=="BCL"])
PCEvPCL_Clin<-ClinTest(names(extremeTypesfinal)[extremeTypesfinal=="EPCL"],names(extremeTypesfinal)[extremeTypesfinal=="PCL"])

BCEvBCL_Clin[[1]] %>%
  kbl(col.names = c(rep(c("median", "min", "max"),2),"p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  add_header_above(c("","Extreme BCL (N=15)"=3,"BCL (N=81)"=3,"","")) %>%
  row_spec(which(BCEvBCL_Clin[[1]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST6_EBCL_vs_BCL_Clin1.png"),width = 6)

BCEvBCL_Clin[[2]] %>%
  kbl(col.names = c("Extreme BCL","BCL","p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  row_spec(which(BCEvBCL_Clin[[2]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST6_EBCL_vs_BCL_Clin2.png"),width = 6)

PCEvPCL_Clin[[1]] %>%
  kbl(col.names = c(rep(c("median", "min", "max"),2),"p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  add_header_above(c("","Extreme PCL (N=15)"=3,"PCL (N=52)"=3,"","")) %>%
  row_spec(which(PCEvPCL_Clin[[1]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/STables/ST7_EPCL_vs_PCL_Clin1.png"),width = 6)

PCEvPCL_Clin[[2]] %>%
  kbl(col.names = c("Extreme PCL","PCL","p.value","adj.p.value")) %>%
  kable_paper() %>%
  column_spec(1,italic = T) %>%
  row_spec(which(PCEvPCL_Clin[[2]]$adj.p.value<0.05),background="lightgrey") %>%
  as_image(file = file.path(outputDir,"Figures/SFigure3/ST7_EPCL_vs_PCL_Clin2.png"),width = 6)


##############################
# Supplemental Figure S3E
# linear model of BCL patient BM by EScore, CXCR4 status, gender, and EBCL status
# Includes MAF analysis mentioned in text but not plotted
##############################

#Modeling BM to explore BCE vs BCL differences
bmBceFit<-lm(bm ~ EScore + CXCR4 + gender + factor(extremeTypesfinal[extremeTypesfinal %in% c("BCL","EBCL")]), data=pData(studyCounts)[names(extremeTypesfinal)[extremeTypesfinal %in% c("BCL","EBCL")],])
summary(bmBceFit)

tidy(bmBceFit) %>%
  mutate(term=c("(Intercept)","EScore","CXCR4: Mutated","Sex: Male","Subtype: EBCL")) %>%
  kbl(caption="Linear Model of WM LPC BM Involment in BCL and EBCL Subtyped Samples",digits = 5) %>%
  kable_classic() %>%
  as_image(file=file.path(outputDir,"Figures/SFigure3/SF3E_BM_lm_EBCL.pdf"),width = 6)

# These are not plotted but are mentioned in the text
BCEmaf<-subsetMaf(MAF,tsb=paste0(BCLsamps[extremeTypesfinal[BCLsamps]=="EBCL"],"_Tumor"),)
BCLmaf<-subsetMaf(MAF,tsb=paste0(BCLsamps[extremeTypesfinal[BCLsamps]=="BCL"],"_Tumor"))
forestPlot(mafCompare(BCEmaf,BCLmaf,m1Name = "EBCL", m2Name = "BCL"))

PCLmaf<-subsetMaf(MAF,tsb=paste0(PCLsamps[extremeTypesfinal[PCLsamps]=="PCL"],"_Tumor"))
PCEmaf<-subsetMaf(MAF,tsb=paste0(PCLsamps[extremeTypesfinal[PCLsamps]=="EPCL"],"_Tumor"),)
#note: none found which throws an error
#forestPlot(mafCompare(PCEmaf,PCLmaf,m1Name = "EPCL", m2Name = "PCL"))


##############################
# Figure 5E
# Expression of PRDM1, XBP1, and SDC1 (CD138) by Expanded Subtype
##############################

pdf(file=file.path(outputDir,"Figures/Figure5/F5E_ExtremeSubtype_Genes.pdf"), width = 7, height = 6)
genePlot(studyTpM,c("XBP1","PRDM1","SDC1"),
         group="ExpandedSubtype",
         legend=T, logScale=2,
         legendSize=1.1,
         pointSize=.5,
         theme=EHDWMTheme,
         main="", ylab="Transcripts per Million + 1 (Log 2)",
         RSOverride=TRUE)
dev.off()


##############################
# Supplemental Data
# DGE analyses of EBCL and EPCL compared with BCL and PCL, respectively
##############################

#Finally we will finish off this section with a DGE analysis of EBCL and EPCL
eMod<-model.matrix(~ 0 + factor(ExpandedSubtype) + EScore + factor(CXCR4) + chr6qdel +gender + agebmbx, data=pData(studyCounts)[WMOnly,])
colnames(eMod)<-c("EWM","EBCL","BCL","PCL","EPCL","EScore","CXCR4","chr6q","Sex","Age")
ebclCon<-makeContrasts(EBCL-BCL, levels=eMod)
epclCon<-makeContrasts(EPCL-PCL, levels=eMod)

v<-voom(
  calcNormFactors(
    DGEList(
      counts=counts(studyCounts)[WMExpressed,rownames(eMod)],
      genes = fData(studyCounts)[WMExpressed,c("Chr","GeneSymbol")])),
  design = eMod)
efit<-lmFit(v, design = eMod)
ebclFit<-contrasts.fit(efit,ebclCon)
epclFit<-contrasts.fit(efit,epclCon)
ebclFit<-eBayes(ebclFit)
epclFit<-eBayes(epclFit)
EBCLvBCL<-topTable(ebclFit,p.value =0.05,number=Inf, lfc=.8)
EPCLvPCL<-topTable(epclFit,p.value =0.05,number=Inf, lfc=.8)
EPCLvPCL$Chr<-as.character(unlist(EPCLvPCL$Chr))
EPCLvPCL$GeneSymbol<-as.character(unlist(EPCLvPCL$GeneSymbol))
EBCLvBCL$Chr<-as.character(unlist(EBCLvBCL$Chr))
EBCLvBCL$GeneSymbol<-as.character(unlist(EBCLvBCL$GeneSymbol))
write.table(EBCLvBCL,sep="\t", quote=FALSE, file=file.path(outputDir,"DGE/EBCL_vs_BCL_DGE.tsv"))
write.table(EPCLvPCL,sep="\t", quote=FALSE, file=file.path(outputDir,"DGE/EPCL_vs_PCL_DGE.tsv"))


##############################
# Figure 5F
# PFS Following Proteasome Inhibitor Therapies by Expanded Subtype
##############################

ClinDat<-pData(studyTpM)[WMOnly,]%>%
  mutate(ID=WMOnly)%>%
  dplyr::filter(!is.na(tx1)) %>%
  dplyr::select(ID, tx1, pfs1, pdtx1, responsetx1) %>%
  mutate(Subtype=factor(pData(studyTpM)[ID,"ExpandedSubtype"]), EScore=pData(studyTpM)[ID,"EScoreEL"])
RxType<-ClinDat$tx1

RxType[RxType %in% c("BDR" ,"CARD","IDR")]<-"Proteasome"
RxType[RxType %in% c("BENDA" ,"BENDAR","FR", "CPR", "CDR", "BR")]<-"Chemo"
RxType[RxType %in% c("OFA" ,"R","R-dex")]<-"CD20mab"
RxType[RxType %in% c("IB" ,"IB/ULO","ZANUBRUTINIB","IB/ULO + Zanubrutinib")]<-"BTKi"
RxType[RxType %in% c("RAD001" ,"IVEN")]<-"Other"
ClinDat<-mutate(ClinDat,RxType=RxType)
ClinDatProt<-filter(ClinDat, RxType=="Proteasome")

pdf(file=file.path(outputDir,"Figures/Figure3/F3F_ProtPFSSubtype.pdf"), width = 9, height = 6)
ggsurvplot(
  surv_fit(Surv(ClinDatProt$pfs1, event=ClinDatProt$pdtx1) ~ Subtype, data=ClinDatProt),
  conf.int = FALSE,
  pval = TRUE,
  ylab="Progression Free Survival",
  risk.table=T,
  xlab="Time (Months)",
  risk.table.height=.35,
  title="Proteasome (BDR ,CARD,IDR)",
  palette =ESubtypeTheme$plotColors$fill)

dev.off()


save(studyCounts, file=file.path(dataDir,"studyCounts.RData"))
save(studyTpM, file=file.path(dataDir,"studyTpM.RData"))
