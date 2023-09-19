library(ggplot2)
library(data.table)
library(dplyr)

 GWASlist <- list.files("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/",pattern = "ChrG.GWAS.Results.csv")

simpleMeff <- 212602 # for regular
#simpleMeff <- 112501 # for admixed
adjustedAlpha <- 0.05/simpleMeff



plotGWAS <- c()
counter <- 1
sigSNPs <- data.frame()
for(file in GWASlist){
  temp <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/",file, sep="")) # for regular 
  #temp <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/",file, sep="")) # for admixed
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
            axis.text = element_text(size=12),
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
      theme(axis.text.y = element_text(size=12),
            axis.title.x = element_blank(), 
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
plotGWAS2 <- plotGWAS[c(5,2,6,3,1)]
GWASgrid <- plot_grid(plotlist = plotGWAS2,ncol = 1, labels = c("(A)","(B)","(C)","(D)","(E)"), label_size = 12, vjust = 1, hjust = -0.1)
ggsave("~/Documents/Github/eduardo/Dissertation/ScopelThesis/figures/C4F6v1.png",
       plot = GWASgrid,device = "png",
       width = 12,
       height = 12,
       limitsize = FALSE)


ggsave("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/GWASplot.png",
       plot = GWASgrid,device = "png",
       width = 12,
       height = 12,
       limitsize = FALSE)

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
#            file = "~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/uniSigSnps.txt",
#            quote = FALSE,
#            col.names = FALSE,
#            row.names = FALSE)

EffSNPs <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/sigSNPsEff.txt",
                      col.names = c("SNP","Ref","Alt","INFO"), na.strings = c("","NA"))
#EffSNPs <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/AdmixSigSNPsEff.txt",
#                      col.names = c("SNP","Ref","Alt","INFO")) # for admixed


sigSNPs <- sigSNPs[,c(1,4,5,11,12,15)]
INFO <- data.frame(do.call('rbind', strsplit(as.character(EffSNPs$INFO),'|',fixed = TRUE)))
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
      minIndex <- which.min(as.numeric(INFO[i,seq(from=15, to=195,by = 15)]))
      distIndex <- 15 + (minIndex-1)*15
      mutIndex <- distIndex - 5
      geneIndex <- distIndex - 11
      effIndex <- distIndex - 12
      typeIndex <- distIndex - 13
      geneList[i,c(1,2,3,4,5)] <- INFO[i,c(typeIndex, effIndex, geneIndex, mutIndex, distIndex)]
    }
  }
}
geneList[6,3] <- "YGL182C_u"



EffSNPs <- EffSNPs[,c(1,2,3)]
EffSNPs <- cbind(EffSNPs, geneList)

colnames(EffSNPs)[c(4:8)] <- c("type","effect","gene","region","distToGene")
finalSNPtab <- merge(EffSNPs, sigSNPs, by = "SNP")
finalSNPtab <- finalSNPtab[order(finalSNPtab$SNP,finalSNPtab$P.value),]

finalSNPtab$funct <- rep(NA,26)
finalSNPtab$CIN <- rep(NA,26)
finalSNPtab$funct[1] <- "mRNA binding"
finalSNPtab$CIN[1] <- "Increased (null)"
finalSNPtab$funct[2] <- "assembly of the NPC"
finalSNPtab$CIN[2] <- "Chr loss (overXP)"
finalSNPtab$funct[3] <- "sphingolipid catabolism"
finalSNPtab$CIN[3] <- "None"
finalSNPtab$funct[4] <- "chromatin silencing"
finalSNPtab$CIN[4] <- "None"
finalSNPtab$funct[c(5,6)] <- "protein glycosylation"
finalSNPtab$CIN[c(5,6)] <- "None"
finalSNPtab$funct[c(7,8)] <- "Dubious"
finalSNPtab$CIN[c(7,8)] <- "None"
finalSNPtab$funct[c(9:11)] <- "Dubious"
finalSNPtab$CIN[c(9:11)] <- "None"
finalSNPtab$funct[c(12:15)] <- "HIstone Regulation"
finalSNPtab$CIN[c(12:15)] <- "Chr loss (null)"
finalSNPtab$funct[16] <- "Dubious"
finalSNPtab$CIN[16] <- "None"
finalSNPtab$funct[c(17,18)] <- "nuclear transport"
finalSNPtab$CIN[c(17,18)] <- "Chr loss (overXp)"
finalSNPtab$funct[c(19,20)] <- "RNA catabolism"
finalSNPtab$CIN[c(19,20)] <- "None"
finalSNPtab$funct[21] <- "chitin biosynthesis and septation"
finalSNPtab$CIN[21] <- "Chr loss (overXp)"
finalSNPtab$funct[22] <- "cytoskeleton organization"
finalSNPtab$CIN[22] <- "Increased (HI)"
finalSNPtab$funct[c(23,24)] <- "Unknown"
finalSNPtab$CIN[c(23,24)] <- "None"
finalSNPtab$funct[25] <- "Thiosulfate sulfurtransferase"
finalSNPtab$CIN[25] <- "None"
finalSNPtab$funct[26] <- "macroautophagy and reticulophagy"
finalSNPtab$CIN[26] <- "Increased (null)"


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
AllSNPs <- finalTabForLatex2$SNP
intersect(AllSNPs, AdmixSNPs)
write.table(paste(as.character(finalTabForLatex2[,1]),
      as.character(finalTabForLatex2[,2]),
      as.character(finalTabForLatex2[,3]),
      as.character(finalTabForLatex2[,4]),
      as.character(finalTabForLatex2[,5]),
      as.character(finalTabForLatex2[,6]),
      as.character(finalTabForLatex2[,7]),
      as.character(round(finalTabForLatex2[,8], digits= 3)),
      sep = " & "),file = "~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/TableC4T1.txt",quote = FALSE,sep = "",row.names = FALSE,col.names = FALSE)


dipGainOnly <- fread("~/Documents/GitHub/eduardo/Dissertation/ScopelThesis/STs/C4ST2.txt", data.table = FALSE,blank.lines.skip = TRUE)
myY <- dipGainOnly[, c(20, 44, 45, 16, 27:42)]
colnames(myY)[1:4] <- c("taxa", "ChrG", "clade", "specAneuploidy")
nrow(myY)

aneuploidyDF <- reshape2::melt(myY,id.vars = c("taxa","ChrG","clade", "specAneuploidy"), variable.name = "chr",value.name = "an")
aneuploidyDF[is.na(aneuploidyDF$an),]$an <- 0
aneuploidyDF$an <- as.factor(aneuploidyDF$an)

#ggplot(aneuploidyDF, aes(x=chr, y=taxa, fill = an)) + 
#  geom_tile() +
#  scale_fill_steps2(low = "white",  mid="#ffffbf",high = "#d01c8b", midpoint = 1)

sigSNPsVCF <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/AnnoSigSNPs.vcf",data.table = FALSE)

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


cladeDF <- NJOrderedMSH[NJOrderedMSH$Gene == "NUP100",c(1,5,18)]
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
hitsTab$Type <- gsub("downstream_gene_variant","zdownstream_gene_variant",hitsTab$Type)
hitsTab <- hitsTab[order(hitsTab$Type, hitsTab$p.value),]
GeneOrder <- hitsTab$Gene
SNPdistrPlot <- ggplot(NJOrderedMSH, aes(x = factor(Gene, level = GeneOrder), taxa, fill = Genotype)) + 
  geom_tile(aes(width = 1), show.legend = TRUE) + 
  scale_fill_manual(values = c("gray","#ffffbf","#d01c8b"),na.value = "white", name = "Allele", labels = c("Ref","AltHet","AltHom","NA"))+
  #scale_fill_steps2(low = "#2c7bb6",  mid="#ffffbf",high = "#d01c8b", midpoint = 1,na.value = "white", name = "Allele", labels = c("Ref",""))+
  geom_hline(yintercept = c(maxVec[1],minVec[1]), size =0.5,color = ab6) +
  geom_hline(yintercept = c(maxVec[2],minVec[2]), size =0.5,color = fd5) +
  geom_hline(yintercept = c(maxVec[3],minVec[3]), size =0.5,color = mo8) +
  geom_hline(yintercept = c(maxVec[4],minVec[4]), size =0.5,color = s25) +
  geom_tile(aes(x = 0, y=taxa, color = -ChrG), width = 0.5, size = 4, height = 0.9, show.legend = FALSE) + 
  geom_tile(aes(x = length(unique(NJOrderedMSH$Gene))+1, y=taxa, color = -ChrG), width = 0.5, size = 4, height = 0.9, show.legend = FALSE) +
  scale_color_gradient(low="black",high = "white")+
  scale_x_discrete(expand=c(0,0))+
  theme_minimal()+
  theme(axis.text.y = element_text(color = cladeDF$color,size = 1.5), 
        axis.text.x = element_text(color = hitsTab$color,size = 14,angle = 90, hjust = 1, vjust = 0.5), 
        axis.title = element_blank(), axis.ticks.x = element_line(size = 1),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        panel.background = element_rect(fill = "white", color = "white")) 
  
ggsave("/Users/es47540/Documents/GitHub/eduardo/Dissertation/ScopelThesis/figures/C4F7v6.png", 
       plot = SNPdistrPlot,
       width = 9, 
       height = 12, 
       units = "in",
       limitsize = FALSE, 
       dpi = 600)

ggsave("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newSNPdistrPerCladev2.png", 
       plot = SNPdistrPlot,
       width = 9, 
       height = 12, 
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
table(StrainSNPclade[StrainSNPclade$HIR3 == 1,"clade"])
table(StrainSNPclade[StrainSNPclade$HIR3 == 2,"clade"])

pList <- c()
pListB <- c()
### Fisher's Exact Test for each variant by clade
# HIR3
table(StrainSNPclade[StrainSNPclade$HIR3 != 0,"clade"])
# Sake and HIR3+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$HIR3 == 2 | StrainSNPclade$HIR3 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$HIR3 == 2 | StrainSNPclade$HIR3 == 1),"ChrG"])[2])
# Sake and HIR3-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$HIR3 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$HIR3 == 0,"ChrG"])[2])
# Contingency Table
SakeHIR3 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(SakeHIR3, alternative = "greater")
pX <- as.numeric(fisher.test(SakeHIR3, alternative = "greater")[1])
pList <- append(pList, pX)
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# NUP100
table(StrainSNPclade[StrainSNPclade$NUP100 != 0,"clade"])
# Wine and NUP100+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$NUP100 == 2 | StrainSNPclade$NUP100 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$NUP100 == 2 | StrainSNPclade$NUP100 == 1),"ChrG"])[2])
# Wine and NUP100-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$NUP100 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$NUP100 == 0,"ChrG"])[2])
# Contingency Table
WineNUP100 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(WineNUP100, alternative = "greater")
pX <- as.numeric(fisher.test(WineNUP100, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# MNT3
table(StrainSNPclade[StrainSNPclade$MNT3 != 0,"clade"])
# Ale Beer and MNT3+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & (StrainSNPclade$MNT3 == 2 | StrainSNPclade$MNT3 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & (StrainSNPclade$MNT3 == 2 | StrainSNPclade$MNT3 == 1),"ChrG"])[2])
# Ale Beer and MNT3-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & StrainSNPclade$MNT3 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & StrainSNPclade$MNT3 == 0,"ChrG"])[2])
# Contingency Table
AleBeerMNT3 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AleBeerMNT3, alternative = "greater")
pX <- as.numeric(fisher.test(AleBeerMNT3, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# Mixed Origin and MNT3+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & (StrainSNPclade$MNT3 == 2 | StrainSNPclade$MNT3 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & (StrainSNPclade$MNT3 == 2 | StrainSNPclade$MNT3 == 1),"ChrG"])[2])
# Ale Beer and MNT3-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & StrainSNPclade$MNT3 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & StrainSNPclade$MNT3 == 0,"ChrG"])[2])
# Contingency Table
MixedOriginMNT3 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(MixedOriginMNT3, alternative = "greater")
pX <- as.numeric(fisher.test(MixedOriginMNT3, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# Admixed and MNT3+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mosaic" & (StrainSNPclade$MNT3 == 2 | StrainSNPclade$MNT3 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mosaic" & (StrainSNPclade$MNT3 == 2 | StrainSNPclade$MNT3 == 1),"ChrG"])[2])
# Ale Beer and MNT3-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mosaic" & StrainSNPclade$MNT3 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mosaic" & StrainSNPclade$MNT3 == 0,"ChrG"])[2])
# Contingency Table
MosaicMNT3 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(MosaicMNT3, alternative = "greater")
pX <- as.numeric(fisher.test(MosaicMNT3, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)

# UTP20
table(StrainSNPclade[StrainSNPclade$UTP20 != 0,"clade"])
# Mixed Origin and UTP20+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & (StrainSNPclade$UTP20 == 2 | StrainSNPclade$UTP20 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & (StrainSNPclade$UTP20 == 2 | StrainSNPclade$UTP20 == 1),"ChrG"])[2])
# Mixed Origin and UTP20-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & StrainSNPclade$UTP20 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & StrainSNPclade$UTP20 == 0,"ChrG"])[2])
# Contingency Table
MixedOriginUTP20 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(MixedOriginUTP20, alternative = "greater")
pX <- as.numeric(fisher.test(MixedOriginUTP20, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# French Dairy and UTP20+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$UTP20 == 2 | StrainSNPclade$UTP20 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$UTP20 == 2 | StrainSNPclade$UTP20 == 1),"ChrG"])[2])
# French Dairy and UTP20-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$UTP20 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$UTP20 == 0,"ChrG"])[2])
# Contingency Table
FrenchDairyUTP20 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FrenchDairyUTP20, alternative = "greater")
pX <- as.numeric(fisher.test(FrenchDairyUTP20, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# SEC23
table(StrainSNPclade[StrainSNPclade$SEC23 != 0,"clade"])
# Ale Beer and SEC23+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[2])
# Ale Beer and SEC23-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & StrainSNPclade$SEC23 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & StrainSNPclade$SEC23 == 0,"ChrG"])[2])
# Contingency Table
AleBeerSEC23 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AleBeerSEC23, alternative = "greater")
pX <- as.numeric(fisher.test(AleBeerSEC23, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# Admixed and SEC23+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mosaic" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mosaic" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[2])
# Admixed and SEC23-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mosaic" & StrainSNPclade$SEC23 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mosaic" & StrainSNPclade$SEC23 == 0,"ChrG"])[2])
# Contingency Table
AdmixedSEC23 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AdmixedSEC23, alternative = "greater")
pX <- as.numeric(fisher.test(AdmixedSEC23, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# Wine and SEC23+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[2])
# Wine and SEC23-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$SEC23 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$SEC23 == 0,"ChrG"])[2])
# Contingency Table
WineSEC23 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(WineSEC23, alternative = "greater")
pX <- as.numeric(fisher.test(WineSEC23, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)

# French Dairy and SEC23+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[2])
# Wine and SEC23-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$SEC23 == 0,"ChrG"])[1])
AnRef <-0
# Contingency Table
FDSEC23 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FDSEC23, alternative = "greater")
pX <- as.numeric(fisher.test(WineSEC23, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)

# African Beer and SEC23+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "African Beer" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "African Beer" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[2])
# Wine and SEC23-
EuRef <- 0
AnRef <- 0
# Contingency Table
AfBeerSEC23 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FDSEC23, alternative = "greater")
pX <- as.numeric(fisher.test(WineSEC23, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# Ethiopia and SEC23+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ethiopia" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ethiopia" & (StrainSNPclade$SEC23 == 2 | StrainSNPclade$SEC23 == 1),"ChrG"])[2])
# Ale Beer and SEC23-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ethiopia" & StrainSNPclade$SEC23 == 0,"ChrG"])[1])
AnRef <- 0
# Contingency Table
EthiopiaSEC23 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AleBeerSEC23, alternative = "greater")
pX <- as.numeric(fisher.test(AleBeerSEC23, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef))
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)



# ASF2
table(StrainSNPclade[StrainSNPclade$ASF2 != 0,"clade"])
# Mixed Origin and ASF2+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & (StrainSNPclade$ASF2 == 2 | StrainSNPclade$ASF2 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & (StrainSNPclade$ASF2 == 2 | StrainSNPclade$ASF2 == 1),"ChrG"])[2])
# Mixed Origin and ASF2-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & StrainSNPclade$ASF2 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Mixed Origin" & StrainSNPclade$ASF2 == 0,"ChrG"])[2])
# Contingency Table
MixedOriginASF2 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(MixedOriginASF2, alternative = "greater")
pX <- as.numeric(fisher.test(MixedOriginASF2, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# Sake and ASF2+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$ASF2 == 2 | StrainSNPclade$ASF2 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$ASF2 == 2 | StrainSNPclade$ASF2 == 1),"ChrG"])[2])
# Sake and ASF2-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$ASF2 == 0,"ChrG"])[1])
AnRef <- 0
# Contingency Table
SakeASF2 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(SakeASF2, alternative = "greater")
pX <- as.numeric(fisher.test(SakeASF2, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)

# BNI4
table(StrainSNPclade[StrainSNPclade$BNI4 != 0,"clade"])
# Asian Ferm. and BNI4+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Asian Ferm." & (StrainSNPclade$BNI4 == 2 | StrainSNPclade$BNI4 == 1),"ChrG"])[1])
AnMut <- 0
# Asian Ferm. and BNI4-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Asian Ferm." & StrainSNPclade$BNI4 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Asian Ferm." & StrainSNPclade$BNI4 == 0,"ChrG"])[2])
# Contingency Table
AsianFermBNI4 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AsianFermBNI4, alternative = "less")
pX <- as.numeric(fisher.test(AsianFermBNI4, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# French Dairy and BNI4+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$BNI4 == 2 | StrainSNPclade$BNI4 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$BNI4 == 2 | StrainSNPclade$BNI4 == 1),"ChrG"])[2])
# French Dairy and BNI4-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$BNI4 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$BNI4 == 0,"ChrG"])[2])
# Contingency Table
FrenchDairyBNI4 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FrenchDairyBNI4, alternative = "less")
#pX <- as.numeric(fisher.test(FrenchDairyBNI4, alternative = "greater")[1])
#pList <- append(pList, pX)
# Sake and BNI4+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$BNI4 == 2 | StrainSNPclade$BNI4 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$BNI4 == 2 | StrainSNPclade$BNI4 == 1),"ChrG"])[2])
# Sake and BNI4-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$BNI4 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$BNI4 == 0,"ChrG"])[2])
# Contingency Table
SakeBNI4 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(SakeBNI4, alternative = "greater")
pX <- as.numeric(fisher.test(SakeBNI4, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)



# YGL182C
table(StrainSNPclade[StrainSNPclade$YGL182C != 0,"clade"])
# Sake and YGL182C+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$YGL182C == 2 | StrainSNPclade$YGL182C == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$YGL182C == 2 | StrainSNPclade$YGL182C == 1),"ChrG"])[2])
# Sake and YGL182C-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$YGL182C == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$YGL182C == 0,"ChrG"])[2])
# Contingency Table
SakeYGL182C <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(SakeYGL182C, alternative = "greater")
pX <- as.numeric(fisher.test(SakeYGL182C, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)

# YGL182C_u
table(StrainSNPclade[StrainSNPclade$YGL182C_u != 0,"clade"])
# Sake and YGL182C+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$YGL182C_u == 2 | StrainSNPclade$YGL182C_u == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$YGL182C_u == 2 | StrainSNPclade$YGL182C_u == 1),"ChrG"])[2])
# Sake and YGL182C-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$YGL182C_u == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$YGL182C_u == 0,"ChrG"])[2])
# Contingency Table
SakeYGL182C_u <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(SakeYGL182C, alternative = "greater")
pX <- as.numeric(fisher.test(SakeYGL182C, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# Wine and RDL1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$RDL1 == 2 | StrainSNPclade$RDL1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$RDL1 == 2 | StrainSNPclade$RDL1 == 1),"ChrG"])[2])
# Wine and RDL1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$RDL1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$RDL1 == 0,"ChrG"])[2])
# Contingency Table
WineRDL1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FDRDL1, alternative = "greater")
pX <- as.numeric(fisher.test(WineRDL1, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# RDL1
table(StrainSNPclade[StrainSNPclade$RDL1 != 0,"clade"])
# French Dairy and RDL1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$RDL1 == 2 | StrainSNPclade$RDL1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$RDL1 == 2 | StrainSNPclade$RDL1 == 1),"ChrG"])[2])
# French Dairy and RDL1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$RDL1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$RDL1 == 0,"ChrG"])[2])
# Contingency Table
FDRDL1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FDRDL1, alternative = "greater")
pX <- as.numeric(fisher.test(FDRDL1, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# ARC35
table(StrainSNPclade[StrainSNPclade$ARC35 != 0,"clade"])
# Wine and ARC35+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$ARC35 == 2 | StrainSNPclade$ARC35 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$ARC35 == 2 | StrainSNPclade$ARC35 == 1),"ChrG"])[2])
# Wine and ARC35-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$ARC35 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$ARC35 == 0,"ChrG"])[2])
# Contingency Table
WineARC35 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(WineARC35, alternative = "greater")
pX <- as.numeric(fisher.test(WineARC35, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)


# ZPS1
table(StrainSNPclade[StrainSNPclade$ZPS1 != 0,"clade"])
# French Dairy and ZPS1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$ZPS1 == 2 | StrainSNPclade$ZPS1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$ZPS1 == 2 | StrainSNPclade$ZPS1 == 1),"ChrG"])[2])
# French Dairy and ZPS1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$ZPS1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$ZPS1 == 0,"ChrG"])[2])
# Contingency Table
FDZPS1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FDZPS1, alternative = "greater")
pX <- as.numeric(fisher.test(FDZPS1, alternative = "greater")[1])
pList <- append(pList, pX)
prop.test(x = c(AnMut, AnMut+AnRef), n = c(AnMut+EuMut, AnMut+AnRef + EuMut+EuRef),alternative = "greater")
binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")
pB <- as.numeric(binom.test(x = AnMut, n = AnMut+EuMut, p = (AnMut+AnRef)/(AnMut+AnRef + EuMut+EuRef),alternative = "greater")[3])
pListB <- append(pListB, pB)



# Wine and ZPS1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$ZPS1 == 2 | StrainSNPclade$ZPS1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$ZPS1 == 2 | StrainSNPclade$ZPS1 == 1),"ChrG"])[2])
# Wine and ZPS1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$ZPS1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$ZPS1 == 0,"ChrG"])[2])
# Contingency Table
WineZPS1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FDZPS1, alternative = "greater")
pX <- as.numeric(fisher.test(WineZPS1, alternative = "greater")[1])
pList <- append(pList, pX)


### FDR adjustment
FDRList <- mt.rawp2adjp(pList, proc = "BH", alpha = 0.05)
adjP <- FDRList$adjp
index <- FDRList$index
round(adjP[order(index),],digits = 4)

### FDR adjustment
FDRList <- mt.rawp2adjp(pListB, proc = "BH", alpha = 0.05)
adjP <- FDRList$adjp
index <- FDRList$index
round(adjP[order(index),],digits = 6)


Model <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/GAPIT.Blink.ChrG.GWAS.Results.csv", sep=""))
Model <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/GAPIT.FarmCPU.ChrG.GWAS.Results.csv", sep=""))
Model <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/GAPIT.MLM.ChrG.GWAS.Results.csv", sep=""))
Model <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/GAPIT.MLMM.ChrG.GWAS.Results.csv", sep=""))
Model <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/newGWASoutput/GAPIT.CMLM.ChrG.GWAS.Results.csv", sep=""))

Model <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/GAPIT.FarmCPU.ChrG.GWAS.Results.csv", sep=""))
Model <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/GAPIT.MLM.ChrG.GWAS.Results.csv", sep=""))
Model <- fread(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newAdmixRunQC40Miss10MAF2.0noRecomb/GWASoutput/GAPIT.CMLM.ChrG.GWAS.Results.csv", sep=""))


head(Model[order(`FDR_Adjusted_P-values`)],20)

### SSD1
ChrIV <- Model[Model$Chromosome == 4,]
SSD1 <- ChrIV[ChrIV$Position >= 1040643 & ChrIV$Position <= 1054326,]
SSD1[order(`FDR_Adjusted_P-values`)]

### SLX8
ChrV <- Model[Model$Chromosome == 5,]
SLX8 <- ChrV[ChrV$Position >= 395348 & ChrV$Position <= 396172,]
SLX8[order(`FDR_Adjusted_P-values`)]

### SLX5
ChrIV <- Model[Model$Chromosome == 4,]
SLX5 <- ChrIV[ChrIV$Position >= 429067 & ChrIV$Position <= 430926,]
SLX5[order(`FDR_Adjusted_P-values`)]

### ULP2
ChrIX <- Model[Model$Chromosome == 9,]
ULP2 <- ChrIX[ChrIX$Position >= 292633 & ChrIX$Position <= 295737,]
ULP2[order(`FDR_Adjusted_P-values`)]

### IPL1
ChrXVI <- Model[Model$Chromosome == 16,]
IPL1 <- ChrXVI[ChrXVI$Position >= 156490 & ChrXVI$Position <= 157593,]
IPL1[order(`FDR_Adjusted_P-values`)]


############# OLD
# KAP123
table(StrainSNPclade[StrainSNPclade$KAP123 != 0,"clade"])
# Ale Beer and KAP123+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & (StrainSNPclade$KAP123 == 2 | StrainSNPclade$KAP123 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & (StrainSNPclade$KAP123 == 2 | StrainSNPclade$KAP123 == 1),"ChrG"])[2])
# Ale Beer and KAP123-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & StrainSNPclade$KAP123 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & StrainSNPclade$KAP123 == 0,"ChrG"])[2])
# Contingency Table
AleBeerKAP123 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AleBeerKAP123, alternative = "greater")

# RTC1
table(StrainSNPclade[StrainSNPclade$RTC1 != 0,"clade"])
# Ale Beer and RTC1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & (StrainSNPclade$RTC1 == 2 | StrainSNPclade$RTC1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & (StrainSNPclade$RTC1 == 2 | StrainSNPclade$RTC1 == 1),"ChrG"])[2])
# Ale Beer and KAP123-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Ale Beer" & StrainSNPclade$RTC1 == 0,"ChrG"])[1])
AnRef <- 0
# Contingency Table
AleBeerRTC1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AleBeerRTC1, alternative = "greater")

# Clinical+Admixed and RTC1+
EuMut <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & (StrainSNPclade$RTC1 == 2 | StrainSNPclade$RTC1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & (StrainSNPclade$RTC1 == 2 | StrainSNPclade$RTC1 == 1),"ChrG"])[2])
# Ale Beer and KAP123-
EuRef <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & StrainSNPclade$RTC1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & StrainSNPclade$RTC1 == 0,"ChrG"])[2])
# Contingency Table
AdmixedRTC1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AdmixedRTC1, alternative = "greater")

# TPM2
table(StrainSNPclade[StrainSNPclade$TPM2 != 0,"clade"])
# French Dairy and TPM2+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$TPM2 == 2 | StrainSNPclade$TPM2 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$TPM2 == 2 | StrainSNPclade$TPM2 == 1),"ChrG"])[2])
# French Dairy and TPM2-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$TPM2 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$TPM2 == 0,"ChrG"])[2])
# Contingency Table
FrenchDairyTPM2 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FrenchDairyTPM2, alternative = "greater")

# Alpechin and TPM2+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Alpechin" & (StrainSNPclade$TPM2 == 2 | StrainSNPclade$TPM2 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Alpechin" & (StrainSNPclade$TPM2 == 2 | StrainSNPclade$TPM2 == 1),"ChrG"])[2])
# Alpechin and TPM2-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Alpechin" & StrainSNPclade$TPM2 == 0,"ChrG"])[1])
AnRef <- 0
# Contingency Table
AlpechinTPM2 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AlpechinTPM2, alternative = "greater")


# ROK1 (no more than 5 strains per lineage with the alternative)
table(StrainSNPclade[StrainSNPclade$ROK1 != 0,"clade"])
# African Wine and ROK1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "African Wine" & (StrainSNPclade$ROK1 == 2 | StrainSNPclade$ROK1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "African Wine" & (StrainSNPclade$ROK1 == 2 | StrainSNPclade$ROK1 == 1),"ChrG"])[2])
# African Wine and ROK1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "African Wine" & StrainSNPclade$ROK1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "African Wine" & StrainSNPclade$ROK1 == 0,"ChrG"])[2])
# Contingency Table
AfricanWineROK1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AfricanWineROK1, alternative = "greater")

# ROK1 (no more than 5 strains per lineage with the alternative)
table(StrainSNPclade[StrainSNPclade$ROK1 != 0,"clade"])
# Clinical + Admixed and ROK1+
EuMut <- 0
AnMut <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & (StrainSNPclade$ROK1 == 2 | StrainSNPclade$ROK1 == 1),"ChrG"])[1])
# Clinical + Admixed and ROK1-
EuRef <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & StrainSNPclade$ROK1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & StrainSNPclade$ROK1 == 0,"ChrG"])[2])
# Contingency Table
AdmixedClinicalROK1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AdmixedClinicalROK1, alternative = "greater")

# LPD1
table(StrainSNPclade[StrainSNPclade$LPD1 != 0,"clade"])
# Wine and LPD1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$LPD1 == 2 | StrainSNPclade$LPD1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$LPD1 == 2 | StrainSNPclade$LPD1 == 1),"ChrG"])[2])
# Wine and LPD1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$LPD1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$LPD1 == 0,"ChrG"])[2])
# Contingency Table
WineLPD1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(WineLPD1, alternative = "greater")


# SGS1
table(StrainSNPclade[StrainSNPclade$SGS1 != 0,"clade"])
# Brazilian Bioethanol and SGS1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Brazilian Bioethanol" & (StrainSNPclade$SGS1 == 2 | StrainSNPclade$SGS1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Brazilian Bioethanol" & (StrainSNPclade$SGS1 == 2 | StrainSNPclade$SGS1 == 1),"ChrG"])[2])
# Brazilian Bioethanol and SGS1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Brazilian Bioethanol" & StrainSNPclade$SGS1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Brazilian Bioethanol" & StrainSNPclade$SGS1 == 0,"ChrG"])[2])
# Contingency Table
BBSGS1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(BBSGS1, alternative = "greater")

# Admixed + Clinical and SGS1+
EuMut <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & (StrainSNPclade$SGS1 == 2 | StrainSNPclade$SGS1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & (StrainSNPclade$SGS1 == 2 | StrainSNPclade$SGS1 == 1),"ChrG"])[2])
# Admixed + Clinical and SGS1-
EuRef <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & StrainSNPclade$SGS1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[(StrainSNPclade$clade == "Clinical" | StrainSNPclade$clade == "Mosaic") & StrainSNPclade$SGS1 == 0,"ChrG"])[2])
# Contingency Table
AdmixedSGS1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(AdmixedSGS1, alternative = "greater")


# TPN1
table(StrainSNPclade[StrainSNPclade$TPN1 != 0,"clade"])
# Sake and TPN1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$TPN1 == 2 | StrainSNPclade$TPN1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & (StrainSNPclade$TPN1 == 2 | StrainSNPclade$TPN1 == 1),"ChrG"])[2])
# Sake and TPN1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$TPN1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Sake" & StrainSNPclade$TPN1 == 0,"ChrG"])[2])
# Contingency Table
SakeTPN1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(SakeTPN1, alternative = "greater")
pX <- as.numeric(fisher.test(SakeTPN1, alternative = "greater")[1])
pList <- append(pList, pX)

# RFM1
table(StrainSNPclade[StrainSNPclade$RFM1 != 0,"clade"])
# Wine and RFM1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$RFM1 == 2 | StrainSNPclade$RFM1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$RFM1 == 2 | StrainSNPclade$RFM1 == 1),"ChrG"])[2])
# Wine and RFM1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$RFM1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$RFM1 == 0,"ChrG"])[2])
# Contingency Table
WineRFM1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(WineRFM1, alternative = "greater")
pX <- as.numeric(fisher.test(WineRFM1, alternative = "greater")[1])
pList <- append(pList, pX)
# French Dairy and RFM1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$RFM1 == 2 | StrainSNPclade$RFM1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$RFM1 == 2 | StrainSNPclade$RFM1 == 1),"ChrG"])[2])
# French Dairy and RFM1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$RFM1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$RFM1 == 0,"ChrG"])[2])
# Contingency Table
FrenchDairyRFM1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FrenchDairyRFM1, alternative = "greater")
pX <- as.numeric(fisher.test(FrenchDairyRFM1, alternative = "greater")[1])
pList <- append(pList, pX)

# DBP6
table(StrainSNPclade[StrainSNPclade$DBP6 != 0,"clade"])
# Wine and DBP6+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$DBP6 == 2 | StrainSNPclade$DBP6 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$DBP6 == 2 | StrainSNPclade$DBP6 == 1),"ChrG"])[2])
# Wine and DBP6-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$DBP6 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$DBP6 == 0,"ChrG"])[2])
# Contingency Table
WineDBP6 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(WineDBP6, alternative = "greater")
pX <- as.numeric(fisher.test(WineDBP6, alternative = "greater")[1])
pList <- append(pList, pX)

# HPF1
table(StrainSNPclade[StrainSNPclade$HPF1 != 0,"clade"])
# French Dairy and HPF1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$HPF1 == 2 | StrainSNPclade$HPF1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & (StrainSNPclade$HPF1 == 2 | StrainSNPclade$HPF1 == 1),"ChrG"])[2])
# French Dairy and HPF1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$HPF1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "French Dairy" & StrainSNPclade$HPF1 == 0,"ChrG"])[2])
# Contingency Table
FrenchDairyHPF1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(FrenchDairyHPF1, alternative = "greater")
pX <- as.numeric(fisher.test(FrenchDairyHPF1, alternative = "greater")[1])
pList <- append(pList, pX)

# DPL1
table(StrainSNPclade[StrainSNPclade$DPL1 != 0,"clade"])
# French Dairy and DPL1+
EuMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$DPL1 == 2 | StrainSNPclade$DPL1 == 1),"ChrG"])[1])
AnMut <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & (StrainSNPclade$DPL1 == 2 | StrainSNPclade$DPL1 == 1),"ChrG"])[2])
# Wine and DPL1-
EuRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$DPL1 == 0,"ChrG"])[1])
AnRef <- as.numeric(table(StrainSNPclade[StrainSNPclade$clade == "Wine" & StrainSNPclade$DPL1 == 0,"ChrG"])[2])
# Contingency Table
WineDPL1 <- matrix(c(AnMut, EuMut, AnRef, EuRef), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Mut","Ref")))
fisher.test(WineDPL1, alternative = "greater")
pX <- as.numeric(fisher.test(WineDPL1, alternative = "greater")[1])
pList <- append(pList, pX)

# DPL1
