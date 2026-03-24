library(dplyr)
library(purrr)
library(DESeq2)
library(edgeR)
library(bvt)
library(gam)
library(destiny)
library(ggplot2)
library(kableExtra)


#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 6
#   Supplemental Figure panel 4
#   Early/Late EScore DGE data
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

# mSigDb data prepared for camera/fry analysis indexed for WMExpressed
load(file.path(dataDir,"study_mSig.RData"))

# Diffusion Map with Diffusion Pseudotime Analysis from Figure Panel 2
load(file=file.path(dataDir,"studyDmap.RData"))


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

WMOnly<-sampleNames(studyCounts)[grep("WM", sampleNames(studyCounts))]

WMExpressed<-featureNames(studyCounts)[rowSums(counts(studyCounts)[,WMOnly]>10)>20]
WMExpressed<-WMExpressed[rowSums(exprs(studyTpM)[WMExpressed,WMOnly] > 1) > 20]



##############################
# Figure 6A
# GSEA for Early/Late EScore
# Also generates Early/Late EScore DGE supplementary data
##############################

# Defining model and contrast.
# Covariates include expanded subtype, CXCR4, sex and age
EScoreMod<-model.matrix(~ 0 + EScore_EL + factor(ExpandedSubtype) + factor(CXCR4) + gender + agebmbx, data=pData(studyCounts)[WMOnly,])
colnames(EScoreMod)<-c("Early","Late","EBCL","BCL","PCL","EPCL","CXCR4","Sex","Age")
EScoreCon<-makeContrasts(Late-Early, levels=EScoreMod)

v<-voom(
  calcNormFactors(
    DGEList(counts=counts(studyCounts)[WMExpressed,WMOnly],genes=fData(studyCounts)[WMExpressed,c("Chr","GeneSymbol")])
  ),
  design=EScoreMod
)

fit<-lmFit(v,design = EScoreMod)
fit<-contrasts.fit(fit,EScoreCon)
fit<-eBayes(fit)

EScoreDGE<-topTable(fit,number=Inf,lfc=.8,p.value=.05)
EScoreDGE$Chr<-unlist(EScoreDGE$Chr)
EScoreDGE$GeneSymbol<-unlist(EScoreDGE$GeneSymbol)
write.table(EScoreDGE, file=file.path(outputDir,"DGE/Early_vs_Late_EScore.tsv"),quote=F,sep="\t")

# Camera HALLMARK analysis
camera(v, study_mSig[[1]], design = EScoreMod, contrast = EScoreCon) %>%
  dplyr::filter(FDR<0.05) %>%
  head() %>%
  kbl() %>% kable_classic() %>%
  as_image(file = file.path(outputDir,"Figures/Figure6/F6A_Camera_Halmark_EScore.pdf"), width=6)


##############################
# Figure 6B
# Barcode plot for HALLMARK Inflamtoary Response with EScore
##############################

pdf(file=file.path(outputDir,"Figures/Figure6/F6B_EScore_CameraIndex_Inflamatory.pdf"), width = 8, height = 4)
barcodeplot(fit$t[,1],
            study_mSig[[1]]$HALLMARK_INFLAMMATORY_RESPONSE,
            design=EScoreMod,
            contrast=EScoreCon,
            main="Late vs Early EScore: HALLMARK - Inflamatory Response",
            xlab = "moderated t-statistic")
dev.off()


##############################
# Supplemental Figure 4A
# S100A9 Correlation with EScore
##############################

# Calculating the stats ahead of plotting
a<-geneScatter(
  data.frame(
    EScore=colData(studyVST)[WMOnly,"EScore"],
    S100A9=assay(studyVST)[featureNames(studyTpM)[fData(studyTpM)$GeneSymbol=="S100A9"],WMOnly]),
  trendline=TRUE)

pdf(file=file.path(outputDir,"Figures/SFigure4/SF4A_S100A9_by_EScore.pdf"), width = 8, height = 7)
geneScatter(
  data.frame(
    EScore=colData(studyVST)[WMOnly,"EScore"],
    S100A9=assay(studyVST)[featureNames(studyTpM)[fData(studyTpM)$GeneSymbol=="S100A9"],WMOnly]),
  trendline=TRUE,
  color=colData(studyVST)[WMOnly,"SimpleSubtype"],
  size=colData(studyVST)[WMOnly,"bm"],
  legend=c("WM Subtype","WM BM Involvement (%)"),
  legendSize=1.1,
  RSOverride=TRUE,
  logScale=c(FALSE,2),
  ylab="S100A9 (VST Expression)",
  sub=paste0("Pearson's r=",round(a$stats$cor$estimate,4),"; p-value=",format.pval(a$stats$cor$p.value)),
  main="")
dev.off()


##############################
# Supplemental Figure 4B
# GAM Analysis of EScore
# Includes GAM DGE data generation
##############################

# Filter for expressed genes not already picked up by traditional DGE
toTest<-WMExpressed[! WMExpressed %in% rownames(EScoreDGE)]
GenesToTest<-assay(studyVST)[toTest,WMOnly]

# Prepare vectors for storing output
DPTgamTested<-vector(mode="list", length=nrow(GenesToTest))
lmVgam<-vector(mode="list", length=nrow(GenesToTest))

# Prepare testing data
testFrame<-data.frame(
  t(GenesToTest),
  pData(studyTpM)[WMOnly,c("CXCR4","SimpleSubtype","gender","EScore","agebmbx")])
testFrame$gender <- factor(testFrame$gender,labels=c("Female","Male"))
testFrame$SimpleSubtype <- factor(testFrame$SimpleSubtype)
testFrame$CXCR4<-factor(testFrame$CXCR4)

# Run gam and lm analysis
for (g in seq(nrow(GenesToTest))) {
  GeneX<-testFrame[,g]
  DPTgamTested[[g]]<-gam(GeneX ~ s(EScore) + SimpleSubtype + CXCR4 + gender + agebmbx, data=testFrame)
  compLM<-lm(GeneX ~ EScore + SimpleSubtype + CXCR4 + gender + agebmbx, data=testFrame)
  atest<-anova(compLM,DPTgamTested[[g]])
  lmVgam[[g]]<-atest
}

# Collecting p-values and adjusting for FDR
paramP<-map_dbl(DPTgamTested, function(x) summary(x)[4][[1]][1,5]) %>% p.adjust(method="BH")
NparamP<-map_dbl(DPTgamTested, function(x) summary(x)[3][[1]][2,3]) %>% p.adjust(method="BH")

# Restoring Gene ID labels
names(lmVgam)<-toTest
names(paramP)<-toTest
names(NparamP)<-toTest
names(DPTgamTested)<-toTest

# Restricting to significant reults and sorting them in order of p-value
gamResults<-fData(studyTpM)[rownames(GenesToTest)[NparamP<0.01],"GeneSymbol"]
gamResults<-gamResults[order(NparamP[names(gamResults)])]

# Location in gamResults list for MYD8 and IL17RA (6 and 34)
MYD88loc<-which(gamResults=="MYD88")
IL17RAloc<-which(gamResults=="IL17RA")

# Plotting IL17RA comparing GAM and LM models
pdf(file=file.path(outputDir,"Figures/SFigure4/SF4B_GAM_IL17RA.pdf"), width = 8, height = 7)
data.frame(EScore=pData(studyTpM)[colnames(GenesToTest),c("EScore")],
           IL17RA=assay(studyVST)[names(gamResults)[IL17RAloc],colnames(GenesToTest)],
           Subtype=pData(studyTpM[,colnames(GenesToTest)])$SimpleSubtype,
           EScore_Level=pData(studyTpM[,colnames(GenesToTest)])$EScore_EL) %>%
  ggplot(aes(x=EScore,y=IL17RA )) +
  geom_point(aes(color=EScore_Level, shape=Subtype)) +
  geom_smooth(method = "gam") +
  geom_smooth(method="lm", color="red") +
  ggtitle("VST IL17RA Expression with DPT Shows Non-Linear Effects",
          subtitle=paste0("Model RSS: ",round(lmVgam[[names(gamResults)[IL17RAloc]]][1,2],3), " (linear); ",round(lmVgam[[names(gamResults)[IL17RAloc]]][2,2],3)," (Spline); p-value: ", formatC(lmVgam[[names(gamResults)[IL17RAloc]]][2,6], format = "e", digits = 3)))
dev.off()

# Plotting MYD88 comparing GAM and LM models
pdf(file=file.path(outputDir,"Figures/SFigure4/SF4B_GAM_MYD88.pdf"), width = 8, height = 7)
data.frame(EScore=pData(studyTpM)[colnames(GenesToTest),c("EScore")],
           MYD88=assay(studyVST)[names(gamResults)[MYD88loc],colnames(GenesToTest)],
           Subtype=pData(studyTpM[,colnames(GenesToTest)])$SimpleSubtype,
           EScore_Level=pData(studyTpM[,colnames(GenesToTest)])$EScore_EL) %>%
  ggplot(aes(x=EScore,y=MYD88 )) +
  geom_point(aes(color=EScore_Level, shape=Subtype)) +
  geom_smooth(method = "gam") +
  geom_smooth(method="lm", color="red") +
  ggtitle("VST MYD88 Expression with DPT Shows Non-Linear Effects",
          subtitle=paste0("Model RSS: ",round(lmVgam[[names(gamResults)[MYD88loc]]][1,2],3), " (linear); ",round(lmVgam[[names(gamResults)[MYD88loc]]][2,2],3)," (Spline); p-value: ", formatC(lmVgam[[names(gamResults)[MYD88loc]]][2,6], format = "e", digits = 3)))
dev.off()

# Writing top 500 GAM results not found in limma DGE
# Note: List is restricted to 500 is the quality of the hits does taper off
top500<-names(gamResults)[1:500]
gamResultTable<-data.frame(
  Symbol=unlist(fData(studyTpM)[top500,"GeneSymbol"]),
  EScore.p.val=format.pval(paramP[top500]),
  LinearRSS=map_dbl(top500, function(x) round(lmVgam[[x]][1,2],3)),
  gamRSS=map_dbl(top500, function(x) round(lmVgam[[x]][2,2],3)),
  LinearVsGamPval=format.pval(map_dbl(top500, function(x) lmVgam[[x]][2,6]))
)

write.table(gamResultTable,sep="\t",quote=FALSE,file=file.path(outputDir,"DGE/GAM_EScore_results.tsv"))


##############################
# Supplemental Figure 4C
# Derivation of the 5 EScore levels
##############################

# Since EScore is fairly continuous, its hard to know how best group values
# The DPT space captures many of the features seen in DC1/2 while collapsing BCL/PCL differences
# Combining DPT1/2 with EScore allows k-means to select cut offs that reflect both spaces
# We will then use this clustering solutoin to select hard EScore cut offs.
EScore_KM<-kmeans(
  data.frame(
    dmapDpt$DPT1,
    dmapDpt$DPT2,
    dmapDpt$dpt),
  5, nstart = 20)$cluster

# Now we want to sort ths KM clusters into accending order
temp<-data.frame(
  EScore=pData(studyTpM)[WMOnly,"EScore"],
  Level=EScore_KM) %>%
  group_by(Level) %>%
  summarize(EScore=mean(EScore))

# Assigning all samples to EScore Level 1 (ESL1)
EScore_5W<-rep("ESL1",length(WMOnly))
names(EScore_5W)<-WMOnly

# Updating assignments for ELS2-5
for(i in 2:5){
  EScore_5W[EScore_KM==temp[order(temp$EScore),]$Level[i]]<-paste0("ESL",i)
}
EScore_5W<-factor(EScore_5W,levels=c("ESL1","ESL2","ESL3","ESL4","ESL5"))

# Plotting initial clustering solution
# Note that this mostly, but not perfectly corresponds to EScore
pdf(file=file.path(outputDir,"Figures/SFigure4/SF4C_EScore_Cutoffs.pdf"), width = 8, height = 7)
geneScatter(sort(pData(studyTpM)[WMOnly,"EScore"],decreasing = F),
            color=EScore_5W[order(pData(studyTpM)[WMOnly,"EScore"],decreasing = F)],
            type="p",
            ylab="WM Evolutionary Score (EScore)",
            main="",xlab="WM Samples",
            legend="K-means Level",legendSize=1.1,
            RSOverride=TRUE)

# Assigning hard cut offs fro EScore based on results
EScore_Cutoffs<-c(-.02,0.1519637,0.2292154,0.3361360,0.4287816,.56)
abline(h=EScore_Cutoffs,lty=2,col="darkgrey")

text(labels=paste0("ESL",1:5),
     x=25,
     y=map_dbl(1:5, function(i) mean(EScore_Cutoffs[c(i,i+1)])),
     font=2,col=npDefaultTheme$plotColors$fill[1:5])
dev.off()


# Updating EScore_5W Using Hard Cut Offs
EScore_5W[pData(studyTpM)[WMOnly,"EScore"]<=EScore_Cutoffs[2]]<-"ESL1"
EScore_5W[pData(studyTpM)[WMOnly,"EScore"]>EScore_Cutoffs[2] & pData(studyTpM)[WMOnly,"EScore"]<=EScore_Cutoffs[3]]<-"ESL2"
EScore_5W[pData(studyTpM)[WMOnly,"EScore"]>EScore_Cutoffs[3] & pData(studyTpM)[WMOnly,"EScore"]<EScore_Cutoffs[4]]<-"ESL3"
EScore_5W[pData(studyTpM)[WMOnly,"EScore"]>=EScore_Cutoffs[4] & pData(studyTpM)[WMOnly,"EScore"]<=EScore_Cutoffs[5]] <-"ESL4"


##############################
# Supplemental Figure 4D
# Diffusion Map with EScore Level, Subtype and BM
##############################

pdf(file=file.path(outputDir,"Figures/SFigure4/SF4D_dmap_EScore_Levels_Subtype.pdf"), width = 9, height = 7)
geneScatter(data.frame(DC1=dmapDpt$DC1,DC2=dmapDpt$DC2),
            color=EScore_5W,size=pData(studyTpM)[WMOnly,"bm"],
            shape=pData(studyTpM)[WMOnly,"SimpleSubtype"],
            legend=c("EScore Level","WM Subtype","WM BM Involvement (%)"),
            main="",
            legendSize=1.1,
            RSOverride=TRUE)
dev.off()


##############################
# Figure 6C
# Path Analysis. No R Code
##############################


##############################
# Figure 6D
# Chemokines by EScore Level
##############################

pdf(file=file.path(outputDir,"Figures/Figure6/F6D_Chemokines_by_EScore.pdf"), width = 8, height = 7)
genePlot(studyTpM[,WMOnly],
         c("CXCL1","CXCL8","CXCL12"),
         group=EScore_5W,
         logScale=2,
         legend=TRUE,
         legendSize=1.1,
         ylab="Transcripts per Million +1 (Log 2)",
         RSOverride=TRUE,
         main="")
dev.off()


##############################
# Figure 6E
# B-cell and Plasma cell Markers by EScore Level
##############################

pdf(file=file.path(outputDir,"Figures/Figure6/F6E_CD19CD20_by_EScore.pdf"), width = 7, height = 6)
genePlot(studyTpM[,WMOnly],
         c("CD19","MS4A1"),
         group=EScore_5W,
         logScale=2,
         legend=TRUE,
         legendSize=1.1,
         ylab="Transcripts per Million +1 (Log 2)",
         RSOverride=TRUE,
         main="")
dev.off()

pdf(file=file.path(outputDir,"Figures/Figure6/F6E_CD38CD138_by_EScore.pdf"), width = 6, height = 6)
genePlot(studyTpM[,WMOnly],
         c("CD38","SDC1"),
         group=EScore_5W,
         logScale=2,
         legend=TRUE,
         legendSize=0.001,
         ylab="Transcripts per Million +1 (Log 2)",
         RSOverride=TRUE,
         main="")
dev.off()


#Saving EScore_5W to study data
temp<-rep(NA,length(sampleNames(studyTpM)))
names(temp)<-sampleNames(studyTpM)
temp[WMOnly]<-as.character(EScore_5W)
temp<-factor(temp,levels=c("ESL1","ESL2","ESL3","ESL4","ESL5"))
pData(studyTpM)$EScore_5W<-temp
pData(studyCounts)$EScore_5W<-temp
colData(studyVST)$EScore_5W<-temp

save(studyCounts, file=file.path(dataDir,"studyCounts.RData"))
save(studyTpM, file=file.path(dataDir,"studyTpM.RData"))
save(studyVST, file=file.path(dataDir,"studyVST.RData"))
