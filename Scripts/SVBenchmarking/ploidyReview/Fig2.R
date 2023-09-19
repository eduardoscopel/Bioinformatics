library(data.table)
library(ggplot2)
library(cowplot)

df <- fread("~/Documents/GitHub/eduardo/ploidy_review/data/supTab1.tsv")

############# Scer #################
Scer <- df[df$Species == "Scer",]
ScerHet <- Scer[Scer$Heterozygous == "TRUE",]
ScerHom <- Scer[Scer$Heterozygous == "FALSE",]

############# VarScan2 #############
ScerVsSens <- (sum(Scer$VarScanTP))/(sum(Scer$VarScanTP)+sum(Scer$VarScanFN))
ScerVsSpec <- (sum(Scer$VarsScanTN))/(sum(Scer$VarScanFP)+sum(Scer$VarsScanTN))
ScerVsFDR <- (sum(Scer$VarScanFP))/(sum(Scer$VarScanTP)+sum(Scer$VarScanFP))

ScerHetVsSens <- (sum(ScerHet$VarScanTP))/(sum(ScerHet$VarScanTP)+sum(ScerHet$VarScanFN))
ScerHetVsSpec <- (sum(ScerHet$VarsScanTN))/(sum(ScerHet$VarScanFP)+sum(ScerHet$VarsScanTN))
ScerHetVsFDR <- (sum(ScerHet$VarScanFP))/(sum(ScerHet$VarScanTP)+sum(ScerHet$VarScanFP))

ScerHomVsSens <- (sum(ScerHom$VarScanTP))/(sum(ScerHom$VarScanTP)+sum(ScerHom$VarScanFN))
ScerHomVsSpec <- (sum(ScerHom$VarsScanTN))/(sum(ScerHom$VarScanFP)+sum(ScerHom$VarsScanTN))
ScerHomVsFDR <- (sum(ScerHom$VarScanFP))/(sum(ScerHom$VarScanTP)+sum(ScerHom$VarScanFP))

ScVS <- data.frame(Tool = "VarScan", Het = c("Heterozygous","Homozygous"), Sensitivity = c(ScerHetVsSens, ScerHomVsSens), Specificity = c(ScerHetVsSpec, ScerHomVsSpec), FDR = c(ScerHetVsFDR, ScerHomVsFDR))

############# Control-Freec2 #############
ScerCFSens <- (sum(Scer$CFTP))/(sum(Scer$CFTP)+sum(Scer$CFFN))
ScerCFSpec <- (sum(Scer$CFTN))/(sum(Scer$CFFP)+sum(Scer$CFTN))
ScerCFFDR <- (sum(Scer$CFFP))/(sum(Scer$CFTP)+sum(Scer$CFFP))

ScerHetCFSens <- (sum(ScerHet$CFTP))/(sum(ScerHet$CFTP)+sum(ScerHet$CFFN))
ScerHetCFSpec <- (sum(ScerHet$CFTN))/(sum(ScerHet$CFFP)+sum(ScerHet$CFTN))
ScerHetCFFDR <- (sum(ScerHet$CFFP))/(sum(ScerHet$CFTP)+sum(ScerHet$CFFP))

ScerHomCFSens <- (sum(ScerHom$CFTP))/(sum(ScerHom$CFTP)+sum(ScerHom$CFFN))
ScerHomCFSpec <- (sum(ScerHom$CFTN))/(sum(ScerHom$CFFP)+sum(ScerHom$CFTN))
ScerHomCFFDR <- (sum(ScerHom$CFFP))/(sum(ScerHom$CFTP)+sum(ScerHom$CFFP))

ScCF <- data.frame(Tool = "Control-Freec", Het = c("Heterozygous","Homozygous"), Sensitivity = c(ScerHetCFSens, ScerHomCFSens), Specificity = c(ScerHetCFSpec, ScerHomCFSpec), FDR = c(ScerHetCFFDR, ScerHomCFFDR))

############# nQuire #############
ScernQSens <- (sum(Scer$nQuireTP))/(sum(Scer$nQuireTP)+sum(Scer$nQuireFN))
ScernQSpec <- (sum(Scer$nQuireTN))/(sum(Scer$nQuireFP)+sum(Scer$nQuireTN))
ScernQFDR <- (sum(Scer$nQuireFP))/(sum(Scer$nQuireTP)+sum(Scer$nQuireFP))

ScerHetnQSens <- (sum(ScerHet$nQuireTP))/(sum(ScerHet$nQuireTP)+sum(ScerHet$nQuireFN))
ScerHetnQSpec <- (sum(ScerHet$nQuireTN))/(sum(ScerHet$nQuireFP)+sum(ScerHet$nQuireTN))
ScerHetnQFDR <- (sum(ScerHet$nQuireFP))/(sum(ScerHet$nQuireTP)+sum(ScerHet$nQuireFP))

ScerHomnQSens <- (sum(ScerHom$nQuireTP))/(sum(ScerHom$nQuireTP)+sum(ScerHom$nQuireFN))
ScerHomnQSpec <- (sum(ScerHom$nQuireTN))/(sum(ScerHom$nQuireFP)+sum(ScerHom$nQuireTN))
ScerHomnQFDR <- (sum(ScerHom$nQuireFP))/(sum(ScerHom$nQuireTP)+sum(ScerHom$nQuireFP))

ScnQuire <- data.frame(Tool = "nQuire", Het = c("Heterozygous","Homozygous"), Sensitivity = c(ScerHetnQSens, ScerHomnQSens), Specificity = c(ScerHetnQSpec, ScerHomnQSpec), FDR = c(ScerHetnQFDR, ScerHomnQFDR))

############# ploidyNGS #############
ScerPlNGSSens <- (sum(Scer$ploidyNGSTP))/(sum(Scer$ploidyNGSTP)+sum(Scer$ploidyNGSFN))
ScerPlNGSSpec <- (sum(Scer$ploidyNGSTN))/(sum(Scer$ploidyNGSFP)+sum(Scer$ploidyNGSTN))
ScerPlNGSFDR <- (sum(Scer$ploidyNGSFP))/(sum(Scer$ploidyNGSTP)+sum(Scer$ploidyNGSFP))

ScerHetPlNGSSens <- (sum(ScerHet$ploidyNGSTP))/(sum(ScerHet$ploidyNGSTP)+sum(ScerHet$ploidyNGSFN))
ScerHetPlNGSSpec <- (sum(ScerHet$ploidyNGSTN))/(sum(ScerHet$ploidyNGSFP)+sum(ScerHet$ploidyNGSTN))
ScerHetPlNGSFDR <- (sum(ScerHet$ploidyNGSFP))/(sum(ScerHet$ploidyNGSTP)+sum(ScerHet$ploidyNGSFP))

ScerHomPlNGSSens <- (sum(ScerHom$ploidyNGSTP))/(sum(ScerHom$ploidyNGSTP)+sum(ScerHom$ploidyNGSFN))
ScerHomPlNGSSpec <- (sum(ScerHom$ploidyNGSTN))/(sum(ScerHom$ploidyNGSFP)+sum(ScerHom$ploidyNGSTN))
ScerHomPlNGSFDR <- (sum(ScerHom$ploidyNGSFP))/(sum(ScerHom$ploidyNGSTP)+sum(ScerHom$ploidyNGSFP))

ScPloidyNGS <- data.frame(Tool = "ploidyNGS", Het = c("Heterozygous","Homozygous"), Sensitivity = c(ScerHetPlNGSSens, ScerHomPlNGSSens), Specificity = c(ScerHetPlNGSSpec, ScerHomPlNGSSpec), FDR = c(ScerHetPlNGSFDR, ScerHomPlNGSFDR))

############# vcf2alleleplot #############
ScerVcfSens <- (sum(Scer$vcf2alleleplotTP))/(sum(Scer$vcf2alleleplotTP)+sum(Scer$vcf2alleleplotFN))
ScerVcfSpec <- (sum(Scer$vcf2alleleplotTN))/(sum(Scer$vcf2alleleplotFP)+sum(Scer$vcf2alleleplotTN))
ScerVcfFDR <- (sum(Scer$vcf2alleleplotFP))/(sum(Scer$vcf2alleleplotTP)+sum(Scer$vcf2alleleplotFP))

ScerHetVcfSens <- (sum(ScerHet$vcf2alleleplotTP))/(sum(ScerHet$vcf2alleleplotTP)+sum(ScerHet$vcf2alleleplotFN))
ScerHetVcfSpec <- (sum(ScerHet$vcf2alleleplotTN))/(sum(ScerHet$vcf2alleleplotFP)+sum(ScerHet$vcf2alleleplotTN))
ScerHetVcfFDR <- (sum(ScerHet$vcf2alleleplotFP))/(sum(ScerHet$vcf2alleleplotTP)+sum(ScerHet$vcf2alleleplotFP))

ScerHomVcfSens <- (sum(ScerHom$vcf2alleleplotTP))/(sum(ScerHom$vcf2alleleplotTP)+sum(ScerHom$vcf2alleleplotFN))
ScerHomVcfSpec <- (sum(ScerHom$vcf2alleleplotTN))/(sum(ScerHom$vcf2alleleplotFP)+sum(ScerHom$vcf2alleleplotTN))
ScerHomVcfFDR <- (sum(ScerHom$vcf2alleleplotFP))/(sum(ScerHom$vcf2alleleplotTP)+sum(ScerHom$vcf2alleleplotFP))

ScVcf <- data.frame(Tool = "vcf2alleleplot", Het = c("Heterozygous","Homozygous"), Sensitivity = c(ScerHetVcfSens, ScerHomVcfSens), Specificity = c(ScerHetVcfSpec, ScerHomVcfSpec), FDR = c(ScerHetVcfFDR, ScerHomVcfFDR))

############# Ymap #############
ScerYmapSens <- (sum(Scer$YmapTP))/(sum(Scer$YmapTP)+sum(Scer$YmapFN))
ScerYmapSpec <- (sum(Scer$YmapTN))/(sum(Scer$YmapFP)+sum(Scer$YmapTN))
ScerYmapFDR <- (sum(Scer$YmapFP))/(sum(Scer$YmapTP)+sum(Scer$YmapFP))

ScerHetYmapSens <- (sum(ScerHet$YmapTP))/(sum(ScerHet$YmapTP)+sum(ScerHet$YmapFN))
ScerHetYmapSpec <- (sum(ScerHet$YmapTN))/(sum(ScerHet$YmapFP)+sum(ScerHet$YmapTN))
ScerHetYmapFDR <- (sum(ScerHet$YmapFP))/(sum(ScerHet$YmapTP)+sum(ScerHet$YmapFP))

ScerHomYmapSens <- (sum(ScerHom$YmapTP))/(sum(ScerHom$YmapTP)+sum(ScerHom$YmapFN))
ScerHomYmapSpec <- (sum(ScerHom$YmapTN))/(sum(ScerHom$YmapFP)+sum(ScerHom$YmapTN))
ScerHomYmapFDR <- (sum(ScerHom$YmapFP))/(sum(ScerHom$YmapTP)+sum(ScerHom$YmapFP))

ScYmap <- data.frame(Tool = "Ymap", Het = c("Heterozygous","Homozygous"), Sensitivity = c(ScerHetYmapSens, ScerHomYmapSens), Specificity = c(ScerHetYmapSpec, ScerHomYmapSpec), FDR = c(ScerHetYmapFDR, ScerHomYmapFDR))



ScerDF <- rbind(ScVS, ScCF, ScnQuire, ScPloidyNGS, ScVcf, ScYmap)
ScerDF
mScerDF <- cbind(ScerDF[,1:2], reshape2::melt(ScerDF[3:5]))
mScerDF

sensP <- ggplot(ScerDF, aes(Tool, Sensitivity, fill = Het)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
  scale_fill_manual(values = c("#a6cee3","#b2df8a"))+
  scale_y_continuous(n.breaks = 10, expand = c(0,0.01))+
  theme_bw() 

specP <- ggplot(ScerDF, aes(Tool, Specificity, fill = Het)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
  scale_fill_manual(values = c("#a6cee3","#b2df8a"))+
  scale_y_continuous(n.breaks = 10, expand = c(0,0.01))+
  theme_bw() 

FDRP <- ggplot(ScerDF, aes(Tool, FDR, fill = Het)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
  scale_fill_manual(values = c("#a6cee3","#b2df8a"))+
  scale_y_continuous(n.breaks = 10, expand = c(0,0.01))+
  theme_bw() 

prowHet <- plot_grid(
  sensP + theme(legend.position = "none"),
  specP + theme(legend.position = "none"), 
  FDRP + theme(legend.position = "none"), 
  nrow = 1, ncol = 3)
legendHet<- get_legend(sensP +
                       theme(legend.position = "bottom"))
plot_grid(prowHet, legendHet, ncol = 1, rel_heights = c(1.5,.1))

ggsave("~/Documents/GitHub/eduardo/ploidy_review/data/Fig2B.png",width = 15, height = 4, units = "in", limitsize = FALSE, dpi = 150)


############# Calb #################
Calb <- df[df$Species == "Calb",]

CalbVsSens <- (sum(Calb$VarScanTP))/(sum(Calb$VarScanTP)+sum(Calb$VarScanFN))
CalbVsSpec <- (sum(Calb$VarsScanTN))/(sum(Calb$VarScanFP)+sum(Calb$VarsScanTN))
CalbVsFDR <- (sum(Calb$VarScanFP))/(sum(Calb$VarScanTP)+sum(Calb$VarScanFP))

CaVs <- data.frame(Tool = "VarScan2", Species = "Calb", Sensitivity = CalbVsSens, Specificity = CalbVsSpec, FDR = CalbVsFDR)
ScVsO <- data.frame(Tool = "VarScan2", Species = "Scer", Sensitivity = ScerVsSens, Specificity = ScerVsSpec, FDR = ScerVsFDR)

############# Control-Freec2 #############
CalbCFSens <- (sum(Calb$CFTP))/(sum(Calb$CFTP)+sum(Calb$CFFN))
CalbCFSpec <- (sum(Calb$CFTN))/(sum(Calb$CFFP)+sum(Calb$CFTN))
CalbCFFDR <- (sum(Calb$CFFP))/(sum(Calb$CFTP)+sum(Calb$CFFP))

CaCF <- data.frame(Tool = "Control-Freec2", Species = "Calb", Sensitivity = CalbCFSens, Specificity = CalbCFSpec, FDR = CalbCFFDR)
ScCFO <- data.frame(Tool = "Control-Freec2", Species = "Scer", Sensitivity = ScerCFSens, Specificity = ScerCFSpec, FDR = ScerCFFDR)

############# nQuire #############
CalbnQSens <- (sum(Calb$nQuireTP))/(sum(Calb$nQuireTP)+sum(Calb$nQuireFN))
CalbnQSpec <- (sum(Calb$nQuireTN))/(sum(Calb$nQuireFP)+sum(Calb$nQuireTN))
CalbnQFDR <- (sum(Calb$nQuireFP))/(sum(Calb$nQuireTP)+sum(Calb$nQuireFP))

CanQ <- data.frame(Tool = "nQuire", Species = "Calb", Sensitivity = CalbnQSens, Specificity = CalbnQSpec, FDR = CalbnQFDR)
ScnQO <- data.frame(Tool = "nQuire", Species = "Scer", Sensitivity = ScernQSens, Specificity = ScernQSpec, FDR = ScernQFDR)

############# ploidyNGS #############
CalbPlNGSSens <- (sum(Calb$ploidyNGSTP))/(sum(Calb$ploidyNGSTP)+sum(Calb$ploidyNGSFN))
CalbPlNGSSpec <- (sum(Calb$ploidyNGSTN))/(sum(Calb$ploidyNGSFP)+sum(Calb$ploidyNGSTN))
CalbPlNGSFDR <- (sum(Calb$ploidyNGSFP))/(sum(Calb$ploidyNGSTP)+sum(Calb$ploidyNGSFP))

CaPlNGS <- data.frame(Tool = "ploidyNGS", Species = "Calb", Sensitivity = CalbPlNGSSens, Specificity = CalbPlNGSSpec, FDR = CalbPlNGSFDR)
ScPlNGSO <- data.frame(Tool = "ploidyNGS", Species = "Scer", Sensitivity = ScerPlNGSSens, Specificity = ScerPlNGSSpec, FDR = ScerPlNGSFDR)


############# vcf2alleleplot #############
CalbVcfSens <- (sum(Calb$vcf2alleleplotTP))/(sum(Calb$vcf2alleleplotTP)+sum(Calb$vcf2alleleplotFN))
CalbVcfSpec <- (sum(Calb$vcf2alleleplotTN))/(sum(Calb$vcf2alleleplotFP)+sum(Calb$vcf2alleleplotTN))
CalbVcfFDR <- (sum(Calb$vcf2alleleplotFP))/(sum(Calb$vcf2alleleplotTP)+sum(Calb$vcf2alleleplotFP))

CaVcf <- data.frame(Tool = "vcf2alleleplot", Species = "Calb", Sensitivity = CalbVcfSens, Specificity = CalbVcfSpec, FDR = CalbVcfFDR)
ScVcfO <- data.frame(Tool = "vcf2alleleplot", Species = "Scer", Sensitivity = ScerVcfSens, Specificity = ScerVcfSpec, FDR = ScerVcfFDR)

############# Ymap #############
CalbYmapSens <- (sum(Calb$YmapTP))/(sum(Calb$YmapTP)+sum(Calb$YmapFN))
CalbYmapSpec <- (sum(Calb$YmapTN))/(sum(Calb$YmapFP)+sum(Calb$YmapTN))
CalbYmapFDR <- (sum(Calb$YmapFP))/(sum(Calb$YmapTP)+sum(Calb$YmapFP))

CaYmap <- data.frame(Tool = "Ymap", Species = "Calb", Sensitivity = CalbYmapSens, Specificity = CalbYmapSpec, FDR = CalbYmapFDR)
ScYmapO <- data.frame(Tool = "Ymap", Species = "Scer", Sensitivity = ScerYmapSens, Specificity = ScerYmapSpec, FDR = ScerYmapFDR)

CaDF <- rbind(CaVs, CaCF, CanQ, CaPlNGS, CaVcf, CaYmap)
ScDFO <- rbind(ScVsO, ScCFO, ScnQO, ScPlNGSO, ScVcfO, ScYmapO)

OverallDF <- rbind(CaDF, ScDFO)

OsensP <- ggplot(OverallDF, aes(Tool, Sensitivity, fill = Species)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
  scale_fill_manual(values = c("#1b9e77","#d95f02"))+
  scale_y_continuous(n.breaks = 10, expand = c(0,0.01))+
  theme_bw()
OsensP
OspecP <- ggplot(OverallDF, aes(Tool, Specificity, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
  scale_fill_manual(values = c("#1b9e77","#d95f02"))+
  scale_y_continuous(n.breaks = 10, expand = c(0,0.01))+
  theme_bw()
OspecP
OFDRP <- ggplot(OverallDF, aes(Tool, FDR, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
  scale_fill_manual(values = c("#1b9e77","#d95f02"))+
  scale_y_continuous(n.breaks = 10, expand = c(0,0.01))+
  theme_bw()
OFDRP

prow <- plot_grid(
  OsensP + theme(legend.position = "none"),
  OspecP + theme(legend.position = "none"), 
  OFDRP + theme(legend.position = "none"), 
  nrow = 1, ncol = 3)
legend <- get_legend(OsensP +
                       theme(legend.position = "bottom"))
plot_grid(prow, legend, ncol = 1, rel_heights = c(1.5,.1))

ggsave("~/Documents/GitHub/eduardo/ploidy_review/data/Fig2A.png",width = 15, height = 4, units = "in", limitsize = FALSE, dpi = 150)
