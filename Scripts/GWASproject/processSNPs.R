library(data.table)

df <- fread("~/Documents/Papers_Scopel/GWAS/sc_table_most_recent.txt", data.table = FALSE,blank.lines.skip = TRUE)
nrow(df)
newdf <- data.frame(matrix(NA, nrow = nrow(df),ncol = 16))
chr <- seq(1,16)
chr <- paste("c", chr, sep = "")
colnames(newdf) <- chr
df <- cbind(df, newdf)
head(df)

first <- strsplit(df$aneuploidy_RD[388],split = ";")
second <- strsplit(first[[1]], split = "@")
third <- strsplit(second[[3]], split = ",")

for(i in as.numeric(rownames(df))){
  first <- strsplit(df$aneuploidy_RD[i], split = ";")
  print(first[[1]])
  if(!(first %in% c(NA, "No", "LOWCOV"))){
    for(j in first){
      second <- strsplit(j, split = "@")
      for(k in second){
        third <- strsplit(k, split = ",")
        type <- third[[1]]
        print(third)
        for(w in third[[2]]){
          column <- as.numeric(w) + 26
          df[i,column]
          df[i,column] <- type
          df[i,column]
        }
      }
    }
  }
}

df

df <- df[df$reference == "Fay_et_al_2018" |
           df$reference == "Gallone_et_al_2016" |
           df$reference == "Han_et_al_2021" |
           df$reference == "Our_Data" |
           df$reference == "Peter_et_al_2018" |
           df$reference == "Pontes_et_al_2020" |
           df$reference == "Zhu_et_al_2016",]
nrow(df)

df$ploidyCons <- ifelse(is.na(df$ploidyFACS),df$ploidyBAF,
                        ifelse(!is.na(df$ploidyFACS) & df$ploidyBAF == "Unknown", df$ploidyFACS,
                               ifelse(df$ploidyFACS == df$ploidyBAF, df$ploidyBAF,
                                      ifelse(df$ploidyFACS != df$ploidyBAF, df$ploidyBAF,"Unknown"))))
df <- df[complete.cases(df$aneuploidy_binary),]
nrow(df)
# Remove low coverage sequences
df <- df[df$aneuploidy_binary != "LOWCOV" | 
           df$avg_cov >= 50,]
nrow(df)
# Remove strains derived from a single spore
df <- df[df$monosporic_derivative != "Yes",]
nrow(df)
# Remove highly contaminated (>10%) strains
df <- df[df$Contaminated < 0.1 |
           df$Contaminated == "No",]
nrow(df)
df$aneuploidy_num <- ifelse(df$aneuploidy_binary == "Yes", 1, 0)
diploids <- df[df$ploidyCons == 2,]
nrow(diploids)
table(diploids$aneuploidy_type, useNA = "ifany")
table(diploids$ecology_category, diploids$aneuploidy_type)
dipGainOnly <- diploids[diploids$aneuploidy_type == "Gain" | diploids$aneuploidy_type == "No",]
nrow(dipGainOnly)
dipGainOnly$lineage <- ifelse(dipGainOnly$population_old %in% c("Europe/wine", "Wine","wine/European","1._Wine/European_(subclade_4)","1._Wine/European_" ,"1._Wine/European_(subclade_1)","1._Wine/European_(subclade_3)","1._Wine/European_(subclade_2)","Wine - Main"),"Wine",
                              ifelse(dipGainOnly$population_old %in% c("Beer2","11._Ale_beer_","Beer 2"), "Ale Beer",
                                     ifelse(dipGainOnly$population_old %in% c("Asian","24._Asian_islands_","Philippines"),"Asian",
                                            ifelse(dipGainOnly$population_old %in% c("Mosaic","Mosaic Honey Wine","Mosaic 3 (Baijiu)","Mosaic 4 (South Africa)","Mosaic 5 (South Africa)","M2._Mosaic_region_2","M3._Mosaic_region_3","M3._Mosaic_region_3_","M1._Mosaic_region_1","M1._Mosaic_region_1_","M2._Mosaic_region_2_"), "Mosaic",
                                                   ifelse(dipGainOnly$population_old %in% c("Mixed 1", "8._Mixed_origin_"), "Mixed Origin",
                                                          ifelse(dipGainOnly$population_old %in% c("Huangjiu","Baijiu","Qingkejiu","Mantou 7","Daqu/Baijiu","26._Asian_fermentation_"), "Asian Ferm.",
                                                                 ifelse(dipGainOnly$population_old %in% c("African Honey Wine","African_Palm_Wine","13._African_palm_wine_"), "African Wine",
                                                                        ifelse(dipGainOnly$population_old %in% c("Mauritius/South Africa"),"South Africa",
                                                                               ifelse(dipGainOnly$population_old %in% c("West African Beer","South African Beer","6._African_beer_"), "African Beer",
                                                                                      ifelse(dipGainOnly$population_old %in% c("7._Mosaic_beer_"), "Mosaic Beer",
                                                                                             ifelse(dipGainOnly$population_old %in% c("2._Alpechin_"), "Alpechin",
                                                                                                    ifelse(dipGainOnly$population_old %in% c("4._Mediterranean_oak_"), "Med Oak",
                                                                                                           ifelse(dipGainOnly$population_old %in% c("12._West_African_cocoa_"), "West African Cocoa",
                                                                                                                  ifelse(dipGainOnly$population_old %in% c("18._Far_East_Asia_"), "Far East Asian",
                                                                                                                         ifelse(dipGainOnly$population_old %in% c("25._Sake_","Sake"), "Sake",
                                                                                                                                ifelse(dipGainOnly$population_old %in% c("5._French_dairy_","Dairy"), "French Dairy",
                                                                                                                                       ifelse(dipGainOnly$population_old %in% c("10._French_Guiana_human_"), "French Guiana Human",
                                                                                                                                              ifelse(dipGainOnly$population_old %in% c("21._Ecuadorean_"),"Ecuadorean",
                                                                                                                                                     ifelse(dipGainOnly$population_old %in% c("17._Taiwanese_"), "Taiwan",
                                                                                                                                                            ifelse(dipGainOnly$population_old %in% c("14._CHNIII_","20._CHN_V_","16._CHNI","15._CHNII_"), "China",
                                                                                                                                                                   ifelse(dipGainOnly$population_old %in% c("9._Mexican_agave_", "3._Brazilian_bioethanol_","22._Far_East_Russian_","19._Malaysian_"), "Other",
                                                                                                                                                                          ifelse(dipGainOnly$population_old %in% c("23._North_American_oak_"),"North American Oak","NA"))))))))))))))))))))))

nrow(dipGainOnly)
dipGainOnly <- dipGainOnly[!duplicated(dipGainOnly$basename),]
nrow(dipGainOnly)
head(dipGainOnly)

myY <- dipGainOnly[, c(20, 44, 45, 16, 27:42)]
#myY <- data.frame(taxa = dipGainOnly$basename, ChrG = dipGainOnly$aneuploidy_num, clade = dipGainOnly$lineage, specAneuploidy = dipGainOnly$aneuploidy_RD)
myY <- myY[-c(255,257),]
colnames(myY)[1:4] <- c("taxa", "ChrG", "clade", "specAneuploidy")
myY$taxa <- gsub('_','',myY$taxa)
nrow(myY)
aneuploidyDF <- reshape2::melt(myY,id.vars = c("taxa","ChrG","clade", "specAneuploidy"), variable.name = "chr",value.name = "an")
aneuploidyDF[is.na(aneuploidyDF$an),]$an <- 0
ggplot(aneuploidyDF, aes(chr, taxa, fill = an)) + 
  geom_tile() + 
  scale_fill_manual(values = c("red","green","blue","white"))


myG <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/gain2n_gSNPs_QC4.hmp.txt", header = TRUE,data.table = FALSE)
header <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.2/newIDs.txt")
colnames(myG)[12:914] <- as.character(header$V1)


MLMM <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/Filter_MLMM.ChrG_GWAS_result.txt")
FarmCPU <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/Filter_FarmCPU.ChrG_GWAS_result.txt")
Blink <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/Filter_Blink.ChrG_GWAS_result.txt")

mergedSNPs <- rbind(MLMM,Blink)
mergedSNPs <- mergedSNPs[,c(1,10,11)]
mergedSNPs <- c(as.character(MLMM$SNP),as.character(FarmCPU$SNP), as.character(Blink$SNP))
mergedSNPs <- mergedSNPs[!duplicated(mergedSNPs)]
mergedSNPs[29] <- "14:692971"


geneNames <- c("HIR3", "RSM18", "RFM1", "ASI1", "CTL1", "HPF1", "NUP100",
               "RPS5", "YGL182C","ASF2","CDC25","SPO22","RFM1","MNT4","PCA1","PIB2","SRP40","MPM1","MFM1","JHD2","YPL083C","UBR2",
               "MNT3","CMP2","TRM732","RNH1","PRE4","YGR293C","DBP6")
SNPeff <- c("mis", "syn", "syn", "mis", "int", "int", "mis",
            "syn", "syn", "syn", "mis", "mis", "int", "mis",
            "int", "syn", "int", "syn", "mis", "int", "int",
            "syn", "syn", "mis", "syn", "int", "int", "int", "int")
SNPeff <- cbind(mergedSNPs,geneNames,SNPeff)
SNPeff <- as.data.frame(SNPeff)

mergedSNPs$color <- ifelse(mergedSNPs$effect == "mis", "firebrick1",
                          ifelse(SNPeff$SNPeff == "syn", "darkgreen","gray57"))
colnames(mergedSNPs)[1] <- "variable"

temp <- myG[myG$`rs`== mergedSNPs[7],]
#alleles <- temp$alleles
#ref <- strsplit(alleles, "/")[[1]][1]
#mut <- strsplit(alleles, "/")[[1]][2]
#ttemp <- t(temp[,-c(1:11)])
#colnames(ttemp)[1] <- mergedSNPs[1]

a <- data.frame(row.names = colnames(temp[12:length(temp)]))
rownames(a) <- gsub('_','',rownames(a))
a$taxa <- rownames(a)

count = 2
SNPeff <- mergedSNPs
for(i in SNPeff$variable){
  temp <- myG[myG$`rs`== i,]
  alleles <- temp$alleles
  ref <- strsplit(alleles, "/")[[1]][1]
  mut <- strsplit(alleles, "/")[[1]][2]
  a[,count] <- as.character(temp[12:length(temp)])
  a[,count] <- ifelse(a[,count] == ref, 0,
                  ifelse(a[,count] == mut, 2,
                         1))
  #colnames(a)[count] <- as.character(SNPeff[SNPeff$variable == i, ]$geneNames)
  colnames(a)[count] <- as.character(SNPeff[SNPeff$variable == i, ]$variable)
  count = count + 1
  }

#mergedYa <- merge(a, myY, by = "taxa")
#orderedYa <- mergedYa[order(mergedYa$ChrG,decreasing = TRUE),]
#for(col in 1:ncol(a)){
#  colnames(a)[col] <- sub("V", "", colnames(a)[col])
#}
melteda <- reshape2::melt(a)
mergedYa <- merge(melteda, aneuploidyDF, by = "taxa")
mergedYa <- merge(melteda, myY, by = "taxa")
#colnames(mergedYa)[2] <- "geneNames"
mergedYa <- merge(mergedYa, SNPeff, by = "variable")

orderedYa <- mergedYa[order(mergedYa$ChrG, decreasing = TRUE),] # order by Chromosome Gain
orderedYa$taxa <- factor(orderedYa$taxa, levels = unique(orderedYa$taxa))

for(c in unique(orderedYa$variable)){
  chr <- strsplit(as.character(c),":")[[1]][1]
  colnumber <- 6+as.numeric(chr)
  print(paste("this is SNP ", c), sep = "")
  print("Gain at any chr vs mutation type")
  tab1 <- table(orderedYa[orderedYa$variable == c,]$ChrG, orderedYa[orderedYa$variable == c,]$value)
  print(tab1)
  print(paste("Gain at chr ", chr, " vs mutation type"), sep = "")
  tab2 <- table(orderedYa[orderedYa$variable == c,colnumber], orderedYa[orderedYa$variable == c,]$value)
  print(tab2)
  RefGen <- tab1[2]
  MutGen <- tab1[4]+tab1[6]
  RefChr <- sum(tab2[,1])
  MutChr <- sum(tab2[,2:3])
  print(fisher.test(matrix(c(RefGen, MutGen, RefChr, MutChr),nrow = 2,dimnames = list(c("Ref","Mut"),c("Gen","C10")))))
}

cladeOrderedYa <- mergedYa[order(mergedYa$clade,mergedYa$ChrG, decreasing = TRUE),]
cladeOrderedYa$taxa <- factor(cladeOrderedYa$taxa, levels = unique(cladeOrderedYa$taxa))
cladeDF <- cladeOrderedYa[cladeOrderedYa$gene == "ASF2",c(2,5)]
cladeDF$color <- ifelse(cladeDF$clade == "Wine", "Purple",
                        ifelse(cladeDF$clade == "Sake", "Red", 
                               ifelse(cladeDF$clade == "Ale Beer", "orange",
                                      ifelse(cladeDF$clade == "French Dairy", "blue",
                                             ifelse(cladeDF$clade == "Mixed Origin", "brown",
                                                    ifelse(cladeDF$clade == "Mosaic", "pink",
                                                           ifelse(cladeDF$clade == "NA", "black",
                                                                  ifelse(cladeDF$clade %in% c("Asian Ferm.","Asian"),"forestgreen",
                                      "White"))))))))


p1 <- ggplot(orderedYa, aes(gene, taxa, fill = value)) + 
  geom_tile(aes(width = 0.95), show.legend = FALSE) + 
  #scale_fill_distiller(palette = "Accent") + 
  scale_fill_steps2(low = "#2c7bb6",  mid="#ffffbf",high = "#d01c8b", midpoint = 1)+
  theme_minimal()+
  theme(axis.text.y = element_blank(), axis.text.x = element_text(color = SNPeff$color,size = 12, angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank(), axis.ticks.x = element_line(size = 1)) + 
  geom_tile(aes(x = 0.5, y=taxa, color = -ChrG), width = 0, size = 2, height = 0.9, show.legend = FALSE) + 
  geom_tile(aes(x = length(unique(orderedYa$gene))+0.5, y=taxa, color = -ChrG), width = 0, size = 2, height = 0.9, show.legend = FALSE) +
  scale_color_gradient(low="black",high = "white")+
  geom_hline(yintercept = c(1,164), size = 1, color = "black") + 
  scale_x_discrete(expand=c(0,0))
ggsave("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF5.0noRecomb/viewSNPs.png", width = 12, height = 8, units = "in", limitsize = FALSE, dpi = 600)

ggplot(cladeOrderedYa, aes(gene, taxa, fill = value)) + 
  geom_tile(aes(width = 0.95), show.legend = FALSE) + 
  #scale_fill_distiller(palette = "Accent") + 
  scale_fill_steps2(low = "#2c7bb6",  mid="#ffffbf",high = "#d01c8b", midpoint = 1)+
  theme_minimal()+
  theme(axis.text.y = element_text(color = cladeDF$color,size = 3), axis.text.x = element_text(color = SNPeff$color,size = 12, angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank(), axis.ticks.x = element_line(size = 1)) + 
  geom_tile(aes(x = 0.5, y=taxa, color = -ChrG), width = 0, size = 2, height = 0.9, show.legend = FALSE) + 
  geom_tile(aes(x = length(unique(orderedYa$gene))+0.5, y=taxa, color = -ChrG), width = 0, size = 2, height = 0.9, show.legend = FALSE) +
  scale_color_gradient(low="black",high = "white")+
  scale_x_discrete(expand=c(0,0))
ggsave("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF5.0noRecomb/CladeviewSNPs.png", width = 12, height = 8, units = "in", limitsize = FALSE, dpi = 600)

an <- 164
eu <- 735

for(i in SNPeff$gene){
  tmp <- table(orderedYa[orderedYa$gene == i,"ChrG"],orderedYa[orderedYa$gene == i,"value"])  
  tmpMat <- matrix(c(tmp[1],tmp[2],(tmp[3]+tmp[5]),(tmp[4]+tmp[6])), nrow = 2, dimnames = list(c("Euploids","Aneuploids"),c("Ref","Alt")))
  FT <- fisher.test(tmpMat)
  print(paste(i,FT$estimate, FT$p.value,sep=" "))
}

KKQ8<-orderedYa[orderedYa$gene== "KKQ8",]
table(KKQ8[KKQ8$value == 2,]$clade)
table(KKQ8[KKQ8$value == 2,]$ChrG)

ASF2 <-orderedYa[orderedYa$gene== "ASF2",]
table(ASF2[ASF2$value == 2,]$clade)
table(ASF2[ASF2$value == 2,]$ChrG)

SFH1 <-orderedYa[orderedYa$gene== "SFH1",]
table(SFH1[SFH1$value == 2,]$clade)
table(SFH1[SFH1$value == 2,]$ChrG)


numSNPsdf <- data.frame(row.names = unique(orderedYa$taxa))
count2 = 1
for(i in unique(orderedYa$taxa)){
  numSNPsdf$taxa[count2] <- i
  numSNPsdf$ChrG[count2] <- myY[myY$taxa == i,"ChrG"]
  numSNPsdf$SNPs[count2] <- nrow(orderedYa[orderedYa$taxa == i & orderedYa$value !=0,])
  count2 = count2 + 1
}
model0 <- glm(ChrG ~ SNPs, family = "binomial", data = numSNPsdf)
numSNPsdf$ChrG <- as.factor(numSNPsdf$ChrG)
p <- ggplot(numSNPsdf, aes(x = ChrG, y = SNPs)) + geom_boxplot(outlier.alpha = 0)
p + geom_jitter(shape = 16,position = position_jitter(width = 0.1, height = 0.1),size = 0.5, alpha = 0.4, color = "forestgreen") + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0,10,1),name = "Number of alleles associated with Chromosome Gain") + 
  scale_x_discrete(labels = c("Euploids","Aneuploids"))+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 12))+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2)
ggsave("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/multSNPsBox.png", width = 8, height = 8, units = "in", limitsize = FALSE, dpi = 600)

library(ape)
library(treeio)
library(ggtree)
tree <- read.tree("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/RunQC40Miss10MAF1.0Recomb0.8.NJ.tree")
newdf <- data.frame(tree$tip.label, row.names = TRUE)
for(i in 1:length(tree$tip.label)){
  for(j in 1:length(dipGainOnly$basename)){
    if(tree$tip.label[i] == dipGainOnly$basename[j]){
      rownames(newdf)[i] <- dipGainOnly$basename[j]
      newdf$amplification[i] <- dipGainOnly$aneuploidy_binary[j]
      newdf$ecology[i] <- dipGainOnly$ecology_category[j]
      newdf$clade[i] <- dipGainOnly$lineage[j]
    }
  }
}

newdf[!newdf$clade %in% c("Wine","Sake","NA","Mosaic","Mixed Origin","French Dairy", "Asian", "Asian Ferm.","Ale Beer"),"clade"] <- "Other"


newdf[is.na(newdf$ecology),] <- "unknown"
amplification <- as.matrix(newdf)[,1]
eco <- as.matrix(newdf)[,2]
clade <- as.matrix(newdf)[,3]
ampdf <- data.frame(node = nodeid(tree, names(amplification)),
                    amplification = amplification)
traitsdf <- cbind(ampdf, eco, clade)
# Make node numeric
traitsdf$node <- as.numeric(traitsdf$node)
# Add traits to tree
y <- full_join(tree, traitsdf, by = 'node')

ggtree(y, size = 0.2, layout = "fan") + 
  geom_tree(aes(color = clade))+
  scale_color_manual(values = c(ecocolors,"black","yellow","pink"))
  scale_color_manual(values=c("yellow",(unique(cladeDF$color))))

  
  subset <- dipGainOnly[dipGainOnly$lineage %in% c("Wine"),]
  table(subset$lineage,subset$aneuploidy_binary)
  nrow(cladeOrderedYa[cladeOrderedYa$geneNames == "DBP6" & cladeOrderedYa$value == 2 & cladeOrderedYa$ChrG == 1,])
  fisher.test(matrix(c(5,38,6,264),nrow = 2, dimnames = list(c("WineSNP", "WineOv"),c("Aneuploid","Euploid"))))
  