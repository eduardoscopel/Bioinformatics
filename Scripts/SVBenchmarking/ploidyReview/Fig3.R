library(data.table)
library(ggplot2)
library(cowplot)

df <- fread("~/Documents/GitHub/eduardo/ploidy_review/Manuscript/SupTable3.txt")

############# Scer #################
Scer <- df[df$Species == "Scer",]
ScerHet <- Scer[Scer$Heterozygous == "TRUE",]
ScerHom <- Scer[Scer$Heterozygous == "FALSE",]

############# Calb #################
Calb <- df[df$Species == "Calb",]

############# CF #############
CFscore <- sum(Scer$`Control-Freec` == Scer$TruePloidy)/nrow(Scer)
CFHetScore <- sum(ScerHet$`Control-Freec` == ScerHet$TruePloidy)/nrow(ScerHet)
CFHomScore <- sum(ScerHom$`Control-Freec` == ScerHom$TruePloidy)/nrow(ScerHom)
CFdataFrame <- data.frame(Tool = "Control-Freec", 
                          Het = c("Heterozygous", "Homozygous", "All"),
                          Score = c(CFHetScore, CFHomScore, CFscore))
############# nQuire #############
NQscore <- sum(Scer$nQuire == Scer$TruePloidy)/nrow(Scer)
NQHetScore <- sum(ScerHet$nQuire == ScerHet$TruePloidy)/nrow(ScerHet)
NQHomScore <- sum(ScerHom$nQuire == ScerHom$TruePloidy)/nrow(ScerHom)
NQdataFrame <- data.frame(Tool = "nQuire", 
                          Het = c("Heterozygous", "Homozygous", "All"),
                          Score = c(NQHetScore, NQHomScore, NQscore))

############# ploidyNGS #############
PLscore <- sum(Scer$ploidyNGS == Scer$TruePloidy)/nrow(Scer)
PLHetScore <- sum(ScerHet$ploidyNGS == ScerHet$TruePloidy)/nrow(ScerHet)
PLHomScore <- sum(ScerHom$ploidyNGS == ScerHom$TruePloidy)/nrow(ScerHom)
PLdataFrame <- data.frame(Tool = "PloidyNGS", 
                          Het = c("Heterozygous", "Homozygous", "All"),
                          Score = c(PLHetScore, PLHomScore, PLscore))

############# vcf2alleleplot #############
VCFscore <- sum(Scer$vcf2alleleplot == Scer$TruePloidy)/nrow(Scer)
VCFHetScore <- sum(ScerHet$vcf2alleleplot == ScerHet$TruePloidy)/nrow(ScerHet)
VCFHomScore <- sum(ScerHom$vcf2alleleplot == ScerHom$TruePloidy)/nrow(ScerHom)
VCFdataFrame <- data.frame(Tool = "vcf2alleleplot", 
                          Het = c("Heterozygous", "Homozygous", "All"),
                          Score = c(VCFHetScore, VCFHomScore, VCFscore))


#### Plot #####
ScerDF <- rbind(CFdataFrame, NQdataFrame, PLdataFrame, VCFdataFrame)
ggplot(ScerDF, aes(Tool, Score, fill = Het)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
  scale_fill_manual(values = c( "black", "#a6cee3","#b2df8a"))+
  scale_y_continuous(n.breaks = 10, expand = c(0,0.01))+
  ylab("Proportion of correct ploidy calls") +
  theme_bw() 
ggsave("~/Documents/GitHub/eduardo/ploidy_review/Manuscript/Figures/Fig3A.png",width = 5, height = 5, units = "in", limitsize = FALSE, dpi = 300)
