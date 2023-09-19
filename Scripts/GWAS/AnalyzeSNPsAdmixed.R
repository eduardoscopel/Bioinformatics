library(ggplot2)
library(data.table)
library(dplyr)

#GWASlist <- list.files("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/GWASoutput/",pattern = "ChrG.GWAS.Results.csv")
GWASlist <- list.files("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/",pattern = "ChrG.GWAS.Results.csv")
#simpleMeff <- 212602 # for regular
simpleMeff <- 112501 # for admixed
adjustedAlpha <- 0.05/simpleMeff



plotGWAS <- c()
counter <- 1
sigSNPs <- data.frame()
for(file in GWASlist){
  #temp <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/GWASoutput/",file, sep="")) # for regular 
  temp <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/",file, sep="")) # for admixed
  temp$SimpleM_Adj_P <- temp$P.value*simpleMeff
  temp$logP <- -log10(temp$P.value)
  temp$size <- ifelse(temp$SimpleM_Adj_P < 0.05, 2, 1)
  temp$shape <- ifelse(temp$SimpleM_Adj_P < 0.05, 19, 1)
  temp$model <- gsub(".*GAPIT.(.+).ChrG.*", "\\1",file)
  tempSigSNPs <- temp[temp$SimpleM_Adj_P < 0.05,]
  temp_cum <- temp %>%
    group_by(Chromosome) %>%
    summarise(max_pos = max(Position)) %>%
    mutate(bp_add = lag(cumsum(max_pos), default = 0)) %>%
    select(Chromosome, bp_add)
  temp_data <- temp %>%
    inner_join(temp_cum, by = "Chromosome") %>%
    mutate(bp_cum = Position + bp_add)
  axis_set <- temp_data %>%
    group_by(Chromosome) %>%
    summarize(center = mean(bp_cum))
  if(counter == 1){
    plotGWAS[[counter]] <- ggplot(temp_data, aes(x = bp_cum, y = logP, color = as.factor(Chromosome))) +
      geom_point(shape = temp_data$shape, alpha =0.8, size = temp_data$size) + 
      scale_x_continuous(expand = c(0.02,0.02), label = axis_set$Chromosome, breaks = axis_set$center)+
      scale_y_continuous(expand = c(0.02,0.02))+
      scale_color_manual(values = rep(c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854'),4)) + 
      geom_hline(yintercept = -log10(adjustedAlpha), linetype = "dashed")+
      labs(x= NULL, 
           y = expression(-log[10](p)))+
      theme_bw()+
      theme(legend.position = 'none',
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            plot.margin = unit(c(0.2,0.1,0,0.2),"in"))
  }
  else{
    plotGWAS[[counter]] <- ggplot(temp_data, aes(x = bp_cum, y = logP, color = as.factor(Chromosome))) +
      geom_point(shape = temp_data$shape, alpha =0.8, size = temp_data$size) + 
      scale_x_continuous(expand = c(0.02,0.02), label = axis_set$Chromosome, breaks = axis_set$center)+
      scale_y_continuous(expand = c(0.02,0.02))+
      scale_color_manual(values = rep(c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854'),4)) + 
      geom_hline(yintercept = -log10(adjustedAlpha), linetype = "dashed")+
      labs(x= NULL, 
           y = expression(-log[10](p)))+
      theme_bw()+
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.position = 'none',
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            plot.margin = unit(c(0.2,0.1,0,0.2),"in"))
  }
  counter <- counter +1
  sigSNPs <- rbind(sigSNPs,tempSigSNPs)
}
#plotGWAS2 <- plotGWAS[c(5,2,6,3,1)]
#GWASgrid <- plot_grid(plotlist = plotGWAS2,ncol = 1, labels = c("(A)","(B)","(C)","(D)","(E)"), label_size = 12, vjust = 1, hjust = -0.1)
#ggsave("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/GWASoutput/GWASplotv1.png",
#       plot = GWASgrid,device = "png",
#       width = 12,
#       height = 12,
#       limitsize = FALSE)

sigSNPs <- sigSNPs[sigSNPs$model != "GLM" & sigSNPs$model != "SUPER",]
sigSNPs$SNP <- gsub("10:", "X:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("11:", "XI:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("12:", "XII:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("13:", "XIII:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("14:", "XIV:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("15:", "XV:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("16:", "XVI:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("1:", "I:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("2:", "II:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("3:", "III:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("4:", "IV:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("5:", "V:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("6:", "VI:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("7:", "VII:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("8:", "VIII:",sigSNPs$SNP)
sigSNPs$SNP <- gsub("9:", "IX:",sigSNPs$SNP)
uniSNPs <- unique(sigSNPs$SNP)
#write.table(x = uniSNPs, 
#            file = "~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/uniSigSnpsAdmix.txt",
#            quote = FALSE,
#            col.names = FALSE,
#            row.names = FALSE)

#EffSNPs <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/GWASoutput/sigSNPsEff.txt",
#                      col.names = c("SNP","Ref","Alt","INFO"))
EffSNPs <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/AdmixSigSNPsEff.txt",
                      col.names = c("SNP","Ref","Alt","INFO")) # for admixed


sigSNPs <- sigSNPs[,c(1,4,5,11,12,15)]
INFO <- data.frame(do.call('rbind', strsplit(as.character(EffSNPs$INFO),'|',fixed=TRUE)))
INFO <- apply(INFO, 2, function(x) gsub("^$|^ $", NA, x))

geneList <- data.frame()
for(i in 1:nrow(INFO)){
  if(INFO[i,2] == "synonymous_variant" | INFO[i,2] == "missense_variant"){
    geneList[i,c(1,2,3,4,5)] <- INFO[i,c(2,3,4,10,15)]
  }
  else{
    upIndices <- grep(pattern = "upstream_gene_variant",x = INFO[i,])
    distIndices <- upIndices + 13
    if(min(as.numeric(INFO[i,distIndices])) <= 750){
      minIndex <- which.min(as.numeric(INFO[i,distIndices]))
      distIndex <- 15 + (minIndex-1)*15
      mutIndex <- distIndex - 5
      geneIndex <- distIndex - 11
      effIndex <- distIndex - 12
      typeIndex <- distIndex - 13
      geneList[i,c(1,2,3,4,5)] <- INFO[i,c(typeIndex, effIndex, geneIndex, mutIndex, distIndex)]
    }
    else{
      minIndex <- which.min(as.numeric(INFO[i,seq(from=15, to=150,by = 15)]))
      distIndex <- 15 + (minIndex-1)*15
      mutIndex <- distIndex - 5
      geneIndex <- distIndex - 11
      effIndex <- distIndex - 12
      typeIndex <- distIndex - 13
      geneList[i,c(1,2,3,4,5)] <- INFO[i,c(typeIndex, effIndex, geneIndex, mutIndex, distIndex)]
    }
  }
}
#geneList[6,3] <- "YGL182C"

EffSNPs <- EffSNPs[,c(1,2,3)]
EffSNPs <- cbind(EffSNPs, geneList)
colnames(EffSNPs)[c(4:8)] <- c("type","effect","gene","region","distToGene")
finalSNPtab <- merge(EffSNPs, sigSNPs, by = "SNP")
finalSNPtab <- finalSNPtab[order(finalSNPtab$SNP,finalSNPtab$P.value),]

finalSNPtab$funct <- rep(NA,35)
finalSNPtab$CIN <- rep(NA,35)

finalSNPtab$funct[1] <- "transmembrane transporter"
finalSNPtab$CIN[1] <- "None"
finalSNPtab$funct[2] <- "sphingolipid biosynthesis"
finalSNPtab$CIN[2] <- "None"
finalSNPtab$funct[c(3:7)] <- "cell cycle"
finalSNPtab$CIN[c(3:7)] <- "Increased (HI)"
finalSNPtab$funct[8] <- "recombination regulation"
finalSNPtab$CIN[8] <- "None"
finalSNPtab$funct[9] <- "NPC component"
finalSNPtab$CIN[9] <- "None"
finalSNPtab$funct[10] <- "nuclear pore assembly"
finalSNPtab$CIN[10] <- "None"
finalSNPtab$funct[11] <- "DNA replication"
finalSNPtab$CIN[11] <- "Chr rearrangements (ts, rf)"
finalSNPtab$funct[12] <- "Protein kinase"
finalSNPtab$CIN[12] <- "None"
finalSNPtab$funct[c(13,14)] <- "mRNA 5' capping"
finalSNPtab$CIN[c(13,14)] <- "Increased (HI)"
finalSNPtab$funct[c(15,16)] <- "Unknown function"
finalSNPtab$CIN[c(15,16)] <- "None"
finalSNPtab$funct[17] <- "mRNA binding"
finalSNPtab$CIN[17] <- "None"
finalSNPtab$funct[18] <- "Dubious"
finalSNPtab$CIN[18] <- "None"
finalSNPtab$funct[19] <- "Dubious"
finalSNPtab$CIN[19] <- "None"
finalSNPtab$funct[20] <- "catabolism of amino acids"
finalSNPtab$CIN[20] <- "None"
finalSNPtab$funct[21] <- "Dubious"
finalSNPtab$CIN[21] <- "None"
finalSNPtab$funct[22] <- "DNA repair"
finalSNPtab$CIN[22] <- "Chr rearrangements/loss (null, overXP)"
finalSNPtab$funct[23] <- "cytoskeleton organization"
finalSNPtab$CIN[23] <- "None"
finalSNPtab$funct[24] <- "Putative"
finalSNPtab$CIN[24] <- "None"
finalSNPtab$funct[25] <- "cell cycle"
finalSNPtab$CIN[25] <- "Colony sectoring, chr loss (overXP)"
finalSNPtab$funct[26] <- "Putative"
finalSNPtab$CIN[26] <- "None"
finalSNPtab$funct[27] <- "mitochondrial genome maintenance"
finalSNPtab$CIN[27] <- "None"
finalSNPtab$funct[28] <- "peroxisome biogenesis"
finalSNPtab$CIN[28] <- "None"
finalSNPtab$funct[29] <- "endosomal recycling"
finalSNPtab$CIN[29] <- "None"
finalSNPtab$funct[c(30:34)] <- "macroautophagy and reticulophagy"
finalSNPtab$CIN[c(30:34)] <- "Increased (HI)"
finalSNPtab$funct[35] <- "Putative"
finalSNPtab$CIN[35] <- "None"


finalTabForLatex <- data.frame(SNP = finalSNPtab$SNP,
                               Type = finalSNPtab$type,
                               Gene = finalSNPtab$gene,
                               Function = finalSNPtab$funct,
                               CIN = finalSNPtab$CIN,
                               p.value = finalSNPtab$P.value,
                               Model = finalSNPtab$model,
                               MAF = finalSNPtab$maf)
finalTabForLatex$color <- ifelse(finalTabForLatex$Type == "missense_variant", "firebrick1",
                           ifelse(finalTabForLatex$Type == "synonymous_variant", "darkgreen","gray57"))
finalTabForLatex$TypeAbb <- ifelse(finalTabForLatex$Type == "upstream_gene_variant", "Intergenic",
                                ifelse(finalTabForLatex$Type == "missense_variant", "Mis","Syn"))

finalTabForLatex2 <- finalTabForLatex[,c(1, 10, 3,4,5,6,7,8)]
AdmixSNPs <- finalTabForLatex2$SNP

write.table(paste(as.character(finalTabForLatex2[,1]),
                  as.character(finalTabForLatex2[,2]),
                  as.character(finalTabForLatex2[,3]),
                  as.character(finalTabForLatex2[,4]),
                  as.character(finalTabForLatex2[,5]),
                  as.character(finalTabForLatex2[,6]),
                  as.character(finalTabForLatex2[,7]),
                  as.character(round(finalTabForLatex2[,8], digits= 3)),
                  sep = " & "),file = "~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/TableC4T3.txt",quote = FALSE,sep = "",row.names = FALSE,col.names = FALSE)

#dipGainOnly <- fread("~/Documents/GitHub/eduardo/Dissertation/ScopelThesis/STs/C4ST2.txt", data.table = FALSE,blank.lines.skip = TRUE)
AdmixedOnly <- fread("~/Documents/GitHub/eduardo/Dissertation/ScopelThesis/STs/C4ST3.txt", data.table = FALSE,blank.lines.skip = TRUE)
#myY <- dipGainOnly[, c(20, 44, 45, 16, 27:42)]
myY <- AdmixedOnly[, c(1, 44, 45, 16, 27:42)]
colnames(myY)[1:4] <- c("taxa", "ChrG", "clade", "specAneuploidy")
nrow(myY)

aneuploidyDF <- reshape2::melt(myY,id.vars = c("taxa","ChrG","clade", "specAneuploidy"), variable.name = "chr",value.name = "an")
aneuploidyDF[is.na(aneuploidyDF$an),]$an <- 0
aneuploidyDF$an <- as.factor(aneuploidyDF$an)

#ggplot(aneuploidyDF, aes(x=chr, y=taxa, fill = an)) + 
#  geom_tile() +
#  scale_fill_steps2(low = "white",  mid="#ffffbf",high = "#d01c8b", midpoint = 1)

#sigSNPsVCF <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/SNPeff/sigSNPs.recode.vcf",data.table = FALSE)
sigSNPsVCF <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/AnnoSigSNPsAdmix.vcf",data.table = FALSE)

hitsTab <- finalTabForLatex[!duplicated(finalTabForLatex$SNP),]
temp <- sigSNPsVCF[sigSNPsVCF$ID == finalTabForLatex[1,1],]
strainHits <- data.frame(row.names = colnames(temp)[10:length(temp)])
strainHits$taxa <- rownames(strainHits)

count = 2
for(SNP in hitsTab$SNP){
  temp <- sigSNPsVCF[sigSNPsVCF$ID == SNP,]
  strainHits$genotype <- as.character(temp[,10:length(temp)])
  strainHits$genotype <- ifelse(strainHits$genotype == "0/0", 0,
                      ifelse(strainHits$genotype == "1/1", 2,
                             ifelse(strainHits$genotype == "0/1", 1,NA)))
  #colnames(a)[count] <- as.character(SNPeff[SNPeff$variable == i, ]$geneNames)
  colnames(strainHits)[count] <- as.character(hitsTab[hitsTab$SNP == SNP,"Gene"])
  count = count + 1
}



mergedMSH <- reshape2::melt(strainHits)
mergedMSH <- merge(mergedMSH, aneuploidyDF, by = "taxa")
#mergedMSH <- merge(meltedSH, myY, by = "taxa")
#colnames(mergedYa)[2] <- "geneNames"
colnames(mergedMSH)[2] <- "Gene"
colnames(mergedMSH)[3] <- "Genotype"
mergedMSH <- merge(mergedMSH, hitsTab, by = "Gene")
mergedMSH <- merge(mergedMSH, NJorder, by = "taxa")

#orderedMSH <- mergedMSH[order(mergedMSH$ChrG, decreasing = TRUE),] # order by Chromosome Gain
#orderedMSH$taxa <- factor(orderedMSH$taxa, levels = unique(orderedMSH$taxa))

#cladeOrderedMSH <- mergedMSH[order(mergedMSH$clade,mergedMSH$ChrG, decreasing = TRUE),]
#cladeOrderedMSH$taxa <- factor(cladeOrderedMSH$taxa, levels = unique(cladeOrderedMSH$taxa))

NJOrderedMSH <- mergedMSH[order(mergedMSH$order,mergedMSH$ChrG, decreasing = TRUE),]
NJOrderedMSH$taxa <- factor(NJOrderedMSH$taxa, levels = unique(NJOrderedMSH$taxa))


cladeDF <- NJOrderedMSH[NJOrderedMSH$Gene == "YME2",c(1,5,17)]
cladeDF <- cladeDF[!duplicated(cladeDF$taxa),]
minVec <- c()
maxVec <- c()
count <- 1
for(i in c("Ale Beer", "French Dairy", "Mixed Origin", "Sake", "Wine")){
  minVec[count] <- 903 - min(cladeDF[cladeDF$clade == i,"order"]) + 0.45
  maxVec[count] <- 903 - max(cladeDF[cladeDF$clade == i,"order"]) - 0.45
  count <- count + 1
} 
cladeDF

we1 <- "#800080"
ol2 <- "#808000"
bb3 <- "#4B0082"
mo4 <- "#228B22"
fd5 <- "skyblue1"
ab6 <- "#bf812d"
mb7 <- "black"
mo8 <- "#FF1493"
ma9 <- "#8B4513" 
fgh10 <- "#FF8C00"
ab11 <- "#ffff99"
wac12 <- "#D2B48C"
apw13 <- "#9400D3"
ch14 <- "#98FB98"
ch15 <- "#3CB371"
ch16 <- "#006400"
t17 <- "#008080"
fea18 <- "#00FFFF"
m19 <-  "#00BFFF"
ch20 <- "#2F4F4F"
e21 <- "#0000FF"
fer22 <- "#00008B"
nao23 <- "#BDB76B"
ai24 <- "#FF8C00"
s25 <- "#FF0000"
af26 <- "#8B0000"
israel <- "darkblue"
ecuador1 <- "grey26"
clinical <- "coral2"
siberia <- "darkslategray4"
ethiopia <- "yellowgreen"
taiwan2 <- "springgreen3"
malaysian <- "gray70"
saf <- "lightsalmon"

cladecolors <- c(ab11, apw13, ab6, ol2, af26, ai24, bb3, ch14, clinical, ecuador1, e21, ethiopia, fea18, fer22, fd5, fgh10, 
                 israel, malaysian, mo4, ma9, mo8, "black", nao23, s25, siberia, saf, t17, taiwan2, wac12, we1)
cladeDF$color <- ifelse(cladeDF$clade == "Wine", we1,
                        ifelse(cladeDF$clade == "West African Cocoa", wac12,
                               ifelse(cladeDF$clade == "Taiwan2", taiwan2,
                                      ifelse(cladeDF$clade == "Taiwan", t17,
                                             ifelse(cladeDF$clade == "South Africa", saf,
                                                    ifelse(cladeDF$clade == "Siberia", siberia,
                                                           ifelse(cladeDF$clade == "Sake", s25,
                                                                  ifelse(cladeDF$clade == "North American Oak", nao23,
                                                                         ifelse(cladeDF$clade == "Mosaic", "black",
                                                                                ifelse(cladeDF$clade == "Mixed Origin", mo8,
                                                                                       ifelse(cladeDF$clade == "Mexican Agave", ma9,
                                                                                              ifelse(cladeDF$clade == "Med Oak", mo4,
                                                                                                     ifelse(cladeDF$clade == "Malaysian", malaysian,
                                                                                                            ifelse(cladeDF$clade == "Israel", israel,
                                                                                                                   ifelse(cladeDF$clade == "French Guiana Human", fgh10,
                                                                                                                          ifelse(cladeDF$clade == "French Dairy", fd5,
                                                                                                                                 ifelse(cladeDF$clade == "Far East Russian", fer22,
                                                                                                                                        ifelse(cladeDF$clade == "Far East Asian", fea18,
                                                                                                                                               ifelse(cladeDF$clade == "Ethiopia", ethiopia,
                                                                                                                                                      ifelse(cladeDF$clade == "Ecuadorean", e21,
                                                                                                                                                             ifelse(cladeDF$clade == "Ecuador1", ecuador1,
                                                                                                                                                                    ifelse(cladeDF$clade == "Clinical", clinical,
                                                                                                                                                                           ifelse(cladeDF$clade == "China", ch14,
                                                                                                                                                                                  ifelse(cladeDF$clade == "Brazilian Bioethanol", bb3,
                                                                                                                                                                                         ifelse(cladeDF$clade == "Asian Islands", ai24,
                                                                                                                                                                                                ifelse(cladeDF$clade == "Asian Ferm.", af26,
                                                                                                                                                                                                       ifelse(cladeDF$clade == "Alpechin", ol2,
                                                                                                                                                                                                              ifelse(cladeDF$clade == "Ale Beer", ab6,
                                                                                                                                                                                                                     ifelse(cladeDF$clade == "African Wine", apw13,
                                                                                                                                                                                                                            ifelse(cladeDF$clade == "African Beer", ab11,NA))))))))))))))))))))))))))))))


NJOrderedMSH$Genotype <- as.factor(NJOrderedMSH$Genotype)
NJOrderedMSH$order <- as.factor(NJOrderedMSH$order)
hitsTab <- hitsTab[order(hitsTab$Type, hitsTab$p.value),]
GeneOrder <- unique(hitsTab$Gene)
SNPdistrPlot <- ggplot(NJOrderedMSH, aes(x = factor(Gene, level = GeneOrder), taxa, fill = Genotype)) + 
  geom_tile(aes(width = 0.95), show.legend = TRUE) + 
  scale_fill_manual(values = c("gray","#ffffbf","#d01c8b"),na.value = "white", name = "Allele", labels = c("Ref","AltHet","AltHom","NA"))+
  #scale_fill_steps2(low = "#2c7bb6",  mid="#ffffbf",high = "#d01c8b", midpoint = 1,na.value = "white", name = "Allele", labels = c("Ref",""))+
  geom_hline(yintercept = c(maxVec[1],minVec[1]), size =0.5,color = ab6) +
  geom_hline(yintercept = c(maxVec[2],minVec[2]), size =0.5,color = fd5) +
  geom_hline(yintercept = c(maxVec[3],minVec[3]), size =0.5,color = mo8) +
  geom_hline(yintercept = c(maxVec[4],minVec[4]), size =0.5,color = s25) +
  theme_minimal()+
  theme(axis.text.y = element_text(color = cladeDF$color,size = 1), 
        axis.text.x = element_text(color = hitsTab$color,size = 6,angle = 90, hjust = 1, vjust = 0.5), 
        axis.title = element_blank(), axis.ticks.x = element_line(size = 1),
        legend.text = element_text(size=6)) + 
  geom_tile(aes(x = 0.5, y=taxa, color = -ChrG), width = 0, size = 4, height = 0.9, show.legend = FALSE) + 
  geom_tile(aes(x = length(unique(NJOrderedMSH$Gene))+0.5, y=taxa, color = -ChrG), width = 0, size = 4, height = 0.9, show.legend = FALSE) +
  scale_color_gradient(low="black",high = "white")+
  scale_x_discrete(expand=c(0,0))

ggsave("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/AdmixSNPdistrPerClade.png", 
       plot = SNPdistrPlot,
       width = 6, 
       height = 4, 
       units = "in",
       limitsize = FALSE, 
       dpi = 600)

p1 <- ggplot(orderedMSH, aes(Gene, taxa, fill = Genotype)) + 
  geom_tile(aes(width = 0.95), show.legend = FALSE) + 
  #scale_fill_distiller(palette = "Accent") + 
  scale_fill_steps2(low = "#2c7bb6",  mid="#ffffbf",high = "#d01c8b", midpoint = 1)+
  theme_minimal()+
  theme(axis.text.y = element_blank(), axis.text.x = element_text(color = hitsTab$color,size = 12, angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank(), axis.ticks.x = element_line(size = 1)) + 
  geom_tile(aes(x = 0.5, y=taxa, color = -ChrG), width = 0, size = 2, height = 0.9, show.legend = FALSE) + 
  geom_tile(aes(x = length(unique(orderedMSH$Gene))+0.5, y=taxa, color = -ChrG), width = 0, size = 2, height = 0.9, show.legend = FALSE) +
  scale_color_gradient(low="black",high = "white")+
  geom_hline(yintercept = c(1,164), size = 1, color = "black") + 
  scale_x_discrete(expand=c(0,0))

StrainSNPclade <- merge(strainHits, myY, by= "taxa")
admixVec <- c(as.numeric(table(StrainSNPclade$ChrG)[2]), sum(table(StrainSNPclade$ChrG)))
pList <- c()
pListB <- c()
# FUN26
table(StrainSNPclade[StrainSNPclade$FUN26 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$FUN26 == 2 | StrainSNPclade$FUN26 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$FUN26 == 2 | StrainSNPclade$FUN26 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$FUN26 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$FUN26 == 0,"ChrG"])[2])
# Contingency Table
FUN26df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FUN26df, alternative = "greater")
pX <- as.numeric(fisher.test(FUN26df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)


# SUR2
table(StrainSNPclade[StrainSNPclade$SUR2 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$SUR2 == 2 | StrainSNPclade$SUR2 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$SUR2 == 2 | StrainSNPclade$SUR2 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$SUR2 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$SUR2 == 0,"ChrG"])[2])
# Contingency Table
SUR2df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(SUR2df, alternative = "greater")
pX <- as.numeric(fisher.test(SUR2df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# MBP1
table(StrainSNPclade[StrainSNPclade$MBP1 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$MBP1 == 2 | StrainSNPclade$MBP1 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$MBP1 == 2 | StrainSNPclade$MBP1 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$MBP1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$MBP1 == 0,"ChrG"])[2])
# Contingency Table
MBP1df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(MBP1df, alternative = "greater")
pX <- as.numeric(fisher.test(MBP1df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)


# HED1
table(StrainSNPclade[StrainSNPclade$HED1 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$HED1 == 2 | StrainSNPclade$HED1 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$HED1 == 2 | StrainSNPclade$HED1 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$HED1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$HED1 == 0,"ChrG"])[2])
# Contingency Table
HED1df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(HED1df, alternative = "greater")
pX <- as.numeric(fisher.test(HED1df, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+EuMut), n = admixVec,alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# NUP42
table(StrainSNPclade[StrainSNPclade$NUP42 != 0,"ChrG"])
EuMut <- 0
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$NUP42 == 2 | StrainSNPclade$NUP42 == 1,"ChrG"])[1])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$NUP42 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$NUP42 == 0,"ChrG"])[2])
# Contingency Table
NUP42df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(NUP42df[1,1]/NUP42df[2,1])*100
sum(NUP42df[,1])
(NUP42df[1,2]/sum(NUP42df[,2]))*100
sum(NUP42df[,2])
fisher.test(NUP42df, alternative = "greater")
pX <- as.numeric(fisher.test(NUP42df, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+EuMut), n = admixVec,alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# RTN2
table(StrainSNPclade[StrainSNPclade$RTN2 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$RTN2 == 2 | StrainSNPclade$RTN2 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$RTN2 == 2 | StrainSNPclade$RTN2 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$RTN2 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$RTN2 == 0,"ChrG"])[2])
# Contingency Table
RTN2df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(RTN2df[1,1]/sum(RTN2df[,1]))
sum(RTN2df[,1])
(RTN2df[1,2]/sum(RTN2df[,2]))*100
sum(RTN2df[,2])
fisher.test(RTN2df, alternative = "greater")
pX <- as.numeric(fisher.test(RTN2df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]
paste("\textit{RTN2} & ", 
      round(100*(RTN2df[1,1]/sum(RTN2df[,1])),digits = 0), '% (', sum(RTN2df[,1]), ") & ", 
      round(100*(RTN2df[1,2]/sum(RTN2df[,2])), digits = 0), '% (', sum(RTN2df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# PRI1
table(StrainSNPclade[StrainSNPclade$PRI1 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$PRI1 == 2 | StrainSNPclade$PRI1 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$PRI1 == 2 | StrainSNPclade$PRI1 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$PRI1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$PRI1 == 0,"ChrG"])[2])
# Contingency Table
PRI1df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(PRI1df[1,1]/sum(PRI1df[,1]))*100
sum(PRI1df[,1])
(PRI1df[1,2]/sum(PRI1df[,2]))*100
sum(PRI1df[,2])
fisher.test(PRI1df, alternative = "greater")
pX <- as.numeric(fisher.test(PRI1df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{PRI1} & ", 
      round(100*(PRI1df[1,1]/sum(PRI1df[,1])),digits = 0), '% (', sum(PRI1df[,1]), ") & ", 
      round(100*(PRI1df[1,2]/sum(PRI1df[,2])), digits = 0), '% (', sum(PRI1df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# TOS3
table(StrainSNPclade[StrainSNPclade$TOS3 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$TOS3 == 2 | StrainSNPclade$TOS3 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$TOS3 == 2 | StrainSNPclade$TOS3 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$TOS3 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$TOS3 == 0,"ChrG"])[2])
# Contingency Table
TOS3df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(TOS3df[1,1]/sum(TOS3df[,1]))*100
sum(TOS3df[,1])
(TOS3df[1,2]/sum(TOS3df[,2]))*100
sum(TOS3df[,2])
fisher.test(TOS3df, alternative = "greater")
pX <- as.numeric(fisher.test(TOS3df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{TOS3} & ", 
      round(100*(TOS3df[1,1]/sum(TOS3df[,1])),digits = 0), '% (', sum(TOS3df[,1]), ") & ", 
      round(100*(TOS3df[1,2]/sum(TOS3df[,2])), digits = 0), '% (', sum(TOS3df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# CEG1
table(StrainSNPclade[StrainSNPclade$CEG1 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$CEG1 == 2 | StrainSNPclade$CEG1 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$CEG1 == 2 | StrainSNPclade$CEG1 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$CEG1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$CEG1 == 0,"ChrG"])[2])
# Contingency Table
CEG1df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(CEG1df[1,1]/sum(CEG1df[,1]))*100
sum(CEG1df[,1])
(CEG1df[1,2]/sum(CEG1df[,2]))*100
sum(CEG1df[,2])
fisher.test(CEG1df, alternative = "greater")
pX <- as.numeric(fisher.test(CEG1df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{CEG1} & ", 
      round(100*(CEG1df[1,1]/sum(CEG1df[,1])),digits = 0), '% (', sum(CEG1df[,1]), ") & ", 
      round(100*(CEG1df[1,2]/sum(CEG1df[,2])), digits = 0), '% (', sum(CEG1df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# AIM18
table(StrainSNPclade[StrainSNPclade$AIM18 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$AIM18 == 2 | StrainSNPclade$AIM18 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$AIM18 == 2 | StrainSNPclade$AIM18 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$AIM18 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$AIM18 == 0,"ChrG"])[2])
# Contingency Table
AIM18df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(AIM18df[1,1]/sum(AIM18df[,1]))*100
sum(AIM18df[,1])
(AIM18df[1,2]/sum(AIM18df[,2]))*100
sum(AIM18df[,2])
fisher.test(AIM18df, alternative = "greater")
pX <- as.numeric(fisher.test(AIM18df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{AIM18} & ", 
      round(100*(AIM18df[1,1]/sum(AIM18df[,1])),digits = 0), '% (', sum(AIM18df[,1]), ") & ", 
      round(100*(AIM18df[1,2]/sum(AIM18df[,2])), digits = 0), '% (', sum(AIM18df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# JSN1
table(StrainSNPclade[StrainSNPclade$JSN1 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$JSN1 == 2 | StrainSNPclade$JSN1 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$JSN1 == 2 | StrainSNPclade$JSN1 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$JSN1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$JSN1 == 0,"ChrG"])[2])
# Contingency Table
JSN1df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(JSN1df[1,1]/sum(JSN1df[,1]))*100
sum(JSN1df[,1])
(JSN1df[1,2]/sum(JSN1df[,2]))*100
sum(JSN1df[,2])
fisher.test(JSN1df, alternative = "greater")
pX <- as.numeric(fisher.test(JSN1df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{JSN1} & ", 
      round(100*(JSN1df[1,1]/sum(JSN1df[,1])),digits = 0), '% (', sum(JSN1df[,1]), ") & ", 
      round(100*(JSN1df[1,2]/sum(JSN1df[,2])), digits = 0), '% (', sum(JSN1df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# YJL185C
table(StrainSNPclade[StrainSNPclade$YJL185C != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YJL185C == 2 | StrainSNPclade$YJL185C == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YJL185C == 2 | StrainSNPclade$YJL185C == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YJL185C == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YJL185C == 0,"ChrG"])[2])
# Contingency Table
YJL185Cdf <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(YJL185Cdf[1,1]/sum(YJL185Cdf[,1]))*100
sum(YJL185Cdf[,1])
(YJL185Cdf[1,2]/sum(YJL185Cdf[,2]))*100
sum(YJL185Cdf[,2])
fisher.test(YJL185Cdf, alternative = "greater")
pX <- as.numeric(fisher.test(YJL185Cdf, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{YJL185C} & ", 
      round(100*(YJL185Cdf[1,1]/sum(YJL185Cdf[,1])),digits = 0), '% (', sum(YJL185Cdf[,1]), ") & ", 
      round(100*(YJL185Cdf[1,2]/sum(YJL185Cdf[,2])), digits = 0), '% (', sum(YJL185Cdf[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)


# YKL165C-A
table(StrainSNPclade[StrainSNPclade$`YKL165C-A` != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$`YKL165C-A` == 2 | StrainSNPclade$`YKL165C-A` == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$`YKL165C-A` == 2 | StrainSNPclade$`YKL165C-A` == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$`YKL165C-A` == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$`YKL165C-A` == 0,"ChrG"])[2])
# Contingency Table
YKL165CAdf <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(YKL165CAdf[1,1]/sum(YKL165CAdf[,1]))*100
sum(YKL165CAdf[,1])
(YKL165CAdf[1,2]/sum(YKL165CAdf[,2]))*100
sum(YKL165CAdf[,2])
fisher.test(YKL165CAdf, alternative = "greater")
pX <- as.numeric(fisher.test(YKL165CAdf, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{YKL165C-A} & ", 
      round(100*(YKL165CAdf[1,1]/sum(YKL165CAdf[,1])),digits = 0), '% (', sum(YKL165CAdf[,1]), ") & ", 
      round(100*(YKL165CAdf[1,2]/sum(YKL165CAdf[,2])), digits = 0), '% (', sum(YKL165CAdf[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)


# SRY1
table(StrainSNPclade[StrainSNPclade$SRY1 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$SRY1 == 2 | StrainSNPclade$SRY1 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$SRY1 == 2 | StrainSNPclade$SRY1 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$SRY1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$SRY1 == 0,"ChrG"])[2])
# Contingency Table
SRY1df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(SRY1df[1,1]/sum(SRY1df[,1]))*100
sum(SRY1df[,1])
(SRY1df[1,2]/sum(SRY1df[,2]))*100
sum(SRY1df[,2])
fisher.test(SRY1df, alternative = "greater")
pX <- as.numeric(fisher.test(SRY1df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{SRY1} & ", 
      round(100*(SRY1df[1,1]/sum(SRY1df[,1])),digits = 0), '% (', sum(SRY1df[,1]), ") & ", 
      round(100*(SRY1df[1,2]/sum(SRY1df[,2])), digits = 0), '% (', sum(SRY1df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# YKL202W
table(StrainSNPclade[StrainSNPclade$YKL202W != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YKL202W == 2 | StrainSNPclade$YKL202W == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YKL202W == 2 | StrainSNPclade$YKL202W == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YKL202W == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YKL202W == 0,"ChrG"])[2])
# Contingency Table
YKL202Wdf <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(YKL202Wdf[1,1]/sum(YKL202Wdf[,1]))*100
sum(YKL202Wdf[,1])
(YKL202Wdf[1,2]/sum(YKL202Wdf[,2]))*100
sum(YKL202Wdf[,2])
fisher.test(YKL202Wdf, alternative = "greater")
pX <- as.numeric(fisher.test(YKL202Wdf, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{YKL202W} & ", 
      round(100*(YKL202Wdf[1,1]/sum(YKL202Wdf[,1])),digits = 0), '% (', sum(YKL202Wdf[,1]), ") & ", 
      round(100*(YKL202Wdf[1,2]/sum(YKL202Wdf[,2])), digits = 0), '% (', sum(YKL202Wdf[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# RAD5
table(StrainSNPclade[StrainSNPclade$RAD5 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$RAD5 == 2 | StrainSNPclade$RAD5 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$RAD5 == 2 | StrainSNPclade$RAD5 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$RAD5 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$RAD5 == 0,"ChrG"])[2])
# Contingency Table
RAD5df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(RAD5df[1,1]/sum(RAD5df[,1]))*100
sum(RAD5df[,1])
(RAD5df[1,2]/sum(RAD5df[,2]))*100
sum(RAD5df[,2])
fisher.test(RAD5df, alternative = "greater")
pX <- as.numeric(fisher.test(RAD5df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{RAD5} & ", 
      round(100*(RAD5df[1,1]/sum(RAD5df[,1])),digits = 0), '% (', sum(RAD5df[,1]), ") & ", 
      round(100*(RAD5df[1,2]/sum(RAD5df[,2])), digits = 0), '% (', sum(RAD5df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# ACF2
table(StrainSNPclade[StrainSNPclade$ACF2 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$ACF2 == 2 | StrainSNPclade$ACF2 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$ACF2 == 2 | StrainSNPclade$ACF2 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$ACF2 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$ACF2 == 0,"ChrG"])[2])
# Contingency Table
ACF2df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(ACF2df[1,1]/sum(ACF2df[,1]))*100
sum(ACF2df[,1])
(ACF2df[1,2]/sum(ACF2df[,2]))*100
sum(ACF2df[,2])
fisher.test(ACF2df, alternative = "greater")
pX <- as.numeric(fisher.test(ACF2df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{ACF2} & ", 
      round(100*(ACF2df[1,1]/sum(ACF2df[,1])),digits = 0), '% (', sum(ACF2df[,1]), ") & ", 
      round(100*(ACF2df[1,2]/sum(ACF2df[,2])), digits = 0), '% (', sum(ACF2df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# YLR257W
table(StrainSNPclade[StrainSNPclade$YLR257W != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YLR257W == 2 | StrainSNPclade$YLR257W == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YLR257W == 2 | StrainSNPclade$YLR257W == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YLR257W == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YLR257W == 0,"ChrG"])[2])
# Contingency Table
YLR257Wdf <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(YLR257Wdf[1,1]/sum(YLR257Wdf[,1]))*100
sum(YLR257Wdf[,1])
(YLR257Wdf[1,2]/sum(YLR257Wdf[,2]))*100
sum(YLR257Wdf[,2])
fisher.test(YLR257Wdf, alternative = "greater")
pX <- as.numeric(fisher.test(YLR257Wdf, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{YLR257W} & ", 
      round(100*(YLR257Wdf[1,1]/sum(YLR257Wdf[,1])),digits = 0), '% (', sum(YLR257Wdf[,1]), ") & ", 
      round(100*(YLR257Wdf[1,2]/sum(YLR257Wdf[,2])), digits = 0), '% (', sum(YLR257Wdf[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# CLN1
table(StrainSNPclade[StrainSNPclade$CLN1 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$CLN1 == 2 | StrainSNPclade$CLN1 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$CLN1 == 2 | StrainSNPclade$CLN1 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$CLN1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$CLN1 == 0,"ChrG"])[2])
# Contingency Table
CLN1df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(CLN1df[1,1]/sum(CLN1df[,1]))*100
sum(CLN1df[,1])
(CLN1df[1,2]/sum(CLN1df[,2]))*100
sum(CLN1df[,2])
fisher.test(CLN1df, alternative = "greater")
pX <- as.numeric(fisher.test(CLN1df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{CLN1} & ", 
      round(100*(CLN1df[1,1]/sum(CLN1df[,1])),digits = 0), '% (', sum(CLN1df[,1]), ") & ", 
      round(100*(CLN1df[1,2]/sum(CLN1df[,2])), digits = 0), '% (', sum(CLN1df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)


# YMR230W-A
table(StrainSNPclade[StrainSNPclade$`YMR230W-A` != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$`YMR230W-A` == 2 | StrainSNPclade$`YMR230W-A` == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$`YMR230W-A` == 2 | StrainSNPclade$`YMR230W-A` == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$`YMR230W-A` == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$`YMR230W-A` == 0,"ChrG"])[2])
# Contingency Table
YMR230WAdf <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(YMR230WAdf[1,1]/sum(YMR230WAdf[,1]))*100
sum(YMR230WAdf[,1])
(YMR230WAdf[1,2]/sum(YMR230WAdf[,2]))*100
sum(YMR230WAdf[,2])
fisher.test(YMR230WAdf, alternative = "greater")
pX <- as.numeric(fisher.test(YMR230WAdf, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{YMR230W-A} & ", 
      round(100*(YMR230WAdf[1,1]/sum(YMR230WAdf[,1])),digits = 0), '% (', sum(YMR230WAdf[,1]), ") & ", 
      round(100*(YMR230WAdf[1,2]/sum(YMR230WAdf[,2])), digits = 0), '% (', sum(YMR230WAdf[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# YME2
table(StrainSNPclade[StrainSNPclade$YME2 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YME2 == 2 | StrainSNPclade$YME2 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YME2 == 2 | StrainSNPclade$YME2 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YME2 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YME2 == 0,"ChrG"])[2])
# Contingency Table
YME2df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(YME2df[1,1]/sum(YME2df[,1]))*100
sum(YME2df[,1])
(YME2df[1,2]/sum(YME2df[,2]))*100
sum(YME2df[,2])
fisher.test(YME2df, alternative = "greater")
pX <- as.numeric(fisher.test(YME2df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{YME2} & ", 
      round(100*(YME2df[1,1]/sum(YME2df[,1])),digits = 0), '% (', sum(YME2df[,1]), ") & ", 
      round(100*(YME2df[1,2]/sum(YME2df[,2])), digits = 0), '% (', sum(YME2df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# PEX15
table(StrainSNPclade[StrainSNPclade$PEX15 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$PEX15 == 2 | StrainSNPclade$PEX15 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$PEX15 == 2 | StrainSNPclade$PEX15 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$PEX15 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$PEX15 == 0,"ChrG"])[2])
# Contingency Table
PEX15df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(PEX15df[1,1]/sum(PEX15df[,1]))*100
sum(PEX15df[,1])
(PEX15df[1,2]/sum(PEX15df[,2]))*100
sum(PEX15df[,2])
fisher.test(PEX15df, alternative = "greater")
pX <- as.numeric(fisher.test(PEX15df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{PEX15} & ", 
      round(100*(PEX15df[1,1]/sum(PEX15df[,1])),digits = 0), '% (', sum(PEX15df[,1]), ") & ", 
      round(100*(PEX15df[1,2]/sum(PEX15df[,2])), digits = 0), '% (', sum(PEX15df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# RTT10
table(StrainSNPclade[StrainSNPclade$RTT10 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$RTT10 == 2 | StrainSNPclade$RTT10 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$RTT10 == 2 | StrainSNPclade$RTT10 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$RTT10 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$RTT10 == 0,"ChrG"])[2])
# Contingency Table
RTT10df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(RTT10df[1,1]/sum(RTT10df[,1]))*100
sum(RTT10df[,1])
(RTT10df[1,2]/sum(RTT10df[,2]))*100
sum(RTT10df[,2])
fisher.test(RTT10df, alternative = "greater")
pX <- as.numeric(fisher.test(RTT10df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{RTT10} & ", 
      round(100*(RTT10df[1,1]/sum(RTT10df[,1])),digits = 0), '% (', sum(RTT10df[,1]), ") & ", 
      round(100*(RTT10df[1,2]/sum(RTT10df[,2])), digits = 0), '% (', sum(RTT10df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# SEC23
table(StrainSNPclade[StrainSNPclade$SEC23 != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$SEC23 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$SEC23 == 0,"ChrG"])[2])
# Contingency Table
SEC23df <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(SEC23df[1,1]/sum(SEC23df[,1]))*100
sum(SEC23df[,1])
(SEC23df[1,2]/sum(SEC23df[,2]))*100
sum(SEC23df[,2])
fisher.test(SEC23df, alternative = "greater")
pX <- as.numeric(fisher.test(SEC23df, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{SEC23} & ", 
      round(100*(SEC23df[1,1]/sum(SEC23df[,1])),digits = 0), '% (', sum(SEC23df[,1]), ") & ", 
      round(100*(SEC23df[1,2]/sum(SEC23df[,2])), digits = 0), '% (', sum(SEC23df[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

# YPR202W
table(StrainSNPclade[StrainSNPclade$YPR202W != 0,"ChrG"])
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YPR202W == 2 | StrainSNPclade$YPR202W == 1,"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$YPR202W == 2 | StrainSNPclade$YPR202W == 1,"ChrG"])[2])
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YPR202W == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$YPR202W == 0,"ChrG"])[2])
# Contingency Table
YPR202Wdf <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
(YPR202Wdf[1,1]/sum(YPR202Wdf[,1]))*100
sum(YPR202Wdf[,1])
(YPR202Wdf[1,2]/sum(YPR202Wdf[,2]))*100
sum(YPR202Wdf[,2])
fisher.test(YPR202Wdf, alternative = "greater")
pX <- as.numeric(fisher.test(YPR202Wdf, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")
paste("\textit{YPR202W} & ", 
      round(100*(YPR202Wdf[1,1]/sum(YPR202Wdf[,1])),digits = 0), '% (', sum(YPR202Wdf[,1]), ") & ", 
      round(100*(YPR202Wdf[1,2]/sum(YPR202Wdf[,2])), digits = 0), '% (', sum(YPR202Wdf[,2]), ") & & ", 
      round(as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3]),digits = 6),
      sep="")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = admixVec[1]/admixVec[2],alternative = "greater")[3])
pListB <- append(pListB, pB)

### FDR adjustment
FDRList <- mt.rawp2adjp(pList, proc = "BH", alpha = 0.05)
adjP <- FDRList$adjp
index <- FDRList$index
round(adjP[order(index),],digits = 6)


### FDR adjustment
FDRList <- mt.rawp2adjp(pListB, proc = "BH", alpha = 0.05)
adjP <- FDRList$adjp
index <- FDRList$index
round(adjP[order(index),],digits = 6)




