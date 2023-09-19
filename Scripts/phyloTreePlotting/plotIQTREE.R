library(data.table)
library(ggtree)

library(ape)
library(treeio)
library(dplyr)
library(ggrepel)
library(cowplot)
library(ggnewscale)

petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt",header=TRUE)
nrow(petertab)
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="Euploid","No","Yes")
petertab$amp_response <- ifelse(petertab$Aneuploidy_type == "Euploid",NA,"Yes")
petertab$amp_binary <- ifelse(petertab$amp_response == "Yes", 1,0)
petertab$polyploidy <- ifelse(petertab$Ploidy == 3 | 
                                petertab$Ploidy == 4 |
                                petertab$Ploidy == 5, "poly",
                              ifelse(petertab$Ploidy == 1, "haploid",
                                     "diploid"))
petertab$MD <- ifelse(petertab$Reference == 2, "Yes","No")

for(i in 1:nrow(petertab)){
  if(!is.na(petertab$amp_response[i])){
    petertab$chr_amp[i] <- substr(petertab$Aneuploidies[i],6,nchar(petertab$Aneuploidies[i])-1)
  }
  else{
    petertab$chr_amp[i] <- NA
  }
}


# Remove HO deleted strains and strains with no clade, ploidy or heterozygosity information
nrow(petertab)
petertab1 <- petertab[petertab$`HO deletion` == "no",]
nrow(petertab1)
petertab2 <- petertab1[petertab1$MD == "No",]
nrow(petertab2)
petertab3 <- petertab2[complete.cases(petertab2$heterozygosity),]
nrow(petertab3)
petertab4 <- petertab3[complete.cases(petertab3$Clades),]
nrow(petertab4)
petertab5 <- petertab4[complete.cases(petertab4$Aneuploidies),]                                    # Remove strains with no aneuploidy information
nrow(petertab5)
petertab6 <- petertab5[petertab5$Aneuploidy_type == "Euploid" | 
                         petertab5$Aneuploidy_type == "Gain",]
nrow(petertab6)
# Remove strains from the Mosaic clade
petertab7 <- petertab6[petertab6$Clades != "M. Mosaic",]
nrow(petertab7)
#petertab8 <- petertab7[petertab7$Clades != "7. Mosaic beer",]
#nrow(petertab8)
petertab8 <- petertab7[petertab7$Ploidy != 1,]
nrow(petertab8)
petertab9 <- petertab8[petertab8$basename != "CBS1593",]
nrow(petertab9)
petertab10 <- petertab9[petertab9$basename != "CBS382",]
nrow(petertab10)
peter_clades <- petertab10
nrow(peter_clades)

set.seed(12345)
contStrains <- peter_clades[,.SD[sample(.N, min(2,.N))],by = Clades]
contStrains <- data.frame(ID = contStrains$basename, Clade = contStrains$Clades)

write.table(contStrains,"/Users/es47540/Documents/GitHub/eduardo/Contamination/tree/51strains.txt", row.names = FALSE, sep = "\t")


tree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/tree/IQ-Tree-621.contree")
rtree <- root(tree, c("EM14S01_3B","EN14S01","GE14S01_7B"))
newdf <- data.frame(rtree$tip.label, row.names = TRUE)

# Add important info to data frame
for(i in 1:length(rtree$tip.label)){
  for(j in 1:length(peter_clades$basename)){
    if(rtree$tip.label[i] == peter_clades$basename[j]){
      rownames(newdf)[i] <- peter_clades$basename[j]
      newdf$amplification[i] <- peter_clades$amp_response[j]
      newdf$ecology[i] <- peter_clades$`Ecological origins`[j]
      newdf$clade[i] <- peter_clades$Clades[j]
      newdf$chr[i] <- peter_clades$chr_amp[j]
      newdf$ssd1[i] <- peter_clades$SSD1_allele[j]
      newdf$ploidy[i] <- peter_clades$Ploidy[j]
      newdf$ssd1geno[i] <- peter_clades$ssd1geno[j]
      newdf$chrI[i] <- peter_clades$chrI[j]
      newdf$chrII[i] <- peter_clades$chrII[j]
      newdf$chrIII[i] <- peter_clades$chrIII[j]
      newdf$chrIV[i] <- peter_clades$chrIV[j]
      newdf$chrV[i] <- peter_clades$chrV[j]
      newdf$chrVI[i] <- peter_clades$chrVI[j]
      newdf$chrVII[i] <- peter_clades$chrVII[j]
      newdf$chrVIII[i] <- peter_clades$chrVIII[j]
      newdf$chrIX[i] <- peter_clades$chrIX[j]
      newdf$chrX[i] <- peter_clades$chrX[j]
      newdf$chrXI[i] <- peter_clades$chrXI[j]
      newdf$chrXII[i] <- peter_clades$chrXII[j]
      newdf$chrXIII[i] <- peter_clades$chrXIII[j]
      newdf$chrXIV[i] <- peter_clades$chrXIV[j]
      newdf$chrXV[i] <- peter_clades$chrXV[j]
      newdf$chrXVI[i] <- peter_clades$chrXVI[j]
    }
  }
}

# Create a matrix for each trait
amplification <- as.matrix(newdf)[,1]
eco <- as.matrix(newdf)[,2]
clade <- as.matrix(newdf)[,3]
chr_amp <- as.matrix(newdf)[,4]
ssd1 <- as.matrix(newdf)[,5]
ploidy <- as.matrix(newdf)[,6]
ssd1geno <- as.matrix(newdf)[,7]
chr <- as.matrix(newdf)[,8:23]

# Create a data frame to merge node ID with traits
ampdf <- data.frame(node = nodeid(rtree, names(amplification)),
                    amplification = amplification)
traitsdf <- cbind(ampdf, eco, clade, chr_amp,ssd1, ploidy,ssd1geno)
# Make node numeric
traitsdf$node <- as.numeric(traitsdf$node)
# Add traits to tree
y <- full_join(rtree, traitsdf, by = 'node')
# Group root nodes
EM14S01_3Bnode <- which(rownames(newdf) == "EM14S01_3B")
EN14S01node <- which(rownames(newdf) == "EN14S01")
GE14S01_7Bnode<- which(rownames(newdf) == "GE14S01_7B")
m<- MRCA(rtree, EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode)
y <- groupClade(y, m)
y <- full_join(y, traitsdf, by ='node')


darkred <- "#8B0000"
red <- "#FF0000"
darksalmon <- "#E9967A"
darkorange <- "#FF8C00"
darkgolden <- "#FF8C00"
darkkhaki <- "#BDB76B"
olive <- "#808000"
yellow <- "#ffed6f"
lightyellow <- "#ffff99"
darkgreen <- "#006400"
forestgreen <- "#228B22"
palegreen <- "#98FB98"
mediumgreen <- "#3CB371"
darkgray <- "#2F4F4F"
teal <- "#008080"
aqua <- "#00FFFF"
darkblue <- "#00008B"
blue <- "#0000FF"
lightblue <- "#00BFFF"
darkviolet <- "#9400D3"
indigo <- "#4B0082"
purple <- "#800080"
lightpink <- "#EE82EE"
mediumpink <- "#C71585"
pink <- "#FF1493"
brown <- "#8B4513"
tan <- "#D2B48C"

we1 <- purple
ol2 <- olive
bb3 <- "#14B5AF"
mo4 <- forestgreen
fd5 <- lightblue
ab6 <- "#bf812d"
mb7 <- "black"
mo8 <- pink
ma9 <- brown 
fgh10 <- darkorange
ab11 <- darkkhaki
wac12 <- tan
apw13 <- darkviolet
ch14 <- palegreen
ch15 <- mediumgreen
ch16 <- darkgreen
t17 <- teal
fea18 <- aqua
m19 <- lightblue
ch20 <- darkgray
e21 <- blue
fer22 <- darkblue
nao23 <- darkkhaki
ai24 <- darksalmon
s25 <- red
af26 <- darkred

cladecolors <- c(we1, fgh10,ab11,wac12,apw13, ch14, ch15, ch16, t17, fea18, m19, ol2, ch20, e21, fer22, nao23, ai24,s25,af26,bb3, mo4, fd5, ab6, mb7, mo8, ma9)

bak <- darksalmon
bee <- yellow
bio <- indigo
cid <- lightpink
dai <- lightblue
dist <- mediumpink
ferm <- darkgolden
flow <- palegreen
fruit <- mediumgreen
hum <- darkorange
humc <- "#D16900"
ind <- brown
ins <- aqua
lab <- darkblue
nat <- teal
palm <- darkviolet
prob <- tan
sak <- red
soil <- darkkhaki
tr <- forestgreen
unk <- "#000000"
wat <- blue
wine <- purple

ecocolors <- c(bak, bee, bio, cid, dai, dist, ferm, flow, fruit, hum, humc, ind, ins, lab, nat, palm, prob, sak, soil, tr, unk, wat, wine)
#ecocolors <- c(bak, bee, bio, cid, dai, dist, ferm, flow, fruit, hum, humc, ind, ins, nat, palm, prob, sak, soil, tr, unk, wat, wine)


ptree <- ggtree(y,size=0.2,layout="rectangular") + 
  geom_tree(aes(color=ssd1geno, linetype = group), show.legend = FALSE, size=0.5)+
  #scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3"))+
  scale_color_manual(values=c("#de77ae","#4d9221","#1f78b4"))+
  #scale_color_manual(values=c("#de77ae","#4CD000","#50B9FF"))+
  geom_treescale(x=0.015,y=30, width=0.01, offset = 5, linesize = 0.2, fontsize=3) + 
  #geom_point2(aes(subset = label == "100", label = label), size=1.5, alpha = 0.3, color = "black",shape =16, show.legend = FALSE) +
  labs(color="ssd1") +
  new_scale_color()+
  geom_tiplab(size=0.4,aes(color=clade),align = TRUE,linesize = 0.05, linetype = "dotted",show.legend = FALSE, x=0.0525) +
  scale_colour_manual(values = cladecolors)
  #guides(colour = guide_legend(override.aes = list(size=3, label = "")))

  
saketree <- viewClade(ptree, 1066)
#ptree$data[ptree$data$node %in% c(EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode), "x"] <- mean(ptree$data$x)
#ptree$data[ptree$data$node %in% c(622, 1032, 1241), "x"] <- 0.02
#ptree$data <- ptree$data %>% mutate(group = replace(group, label == "GE14S01_7B", 1))

ampdf <- data.frame(chr_amp = traitsdf$amplification)
rownames(ampdf) <- row.names(traitsdf)
ptree <- gheatmap(ptree, 
                  ampdf, 
                  offset = 0.001, 
                  width = 0.025, 
                  colnames_position = "top",
                  font.size = 2,
                  colnames_offset_y = 4) +
  scale_fill_manual(values = c("black"), na.translate = FALSE) +
  labs(fill = "Chr Amplification") + 
  scale_x_continuous(expand = c(0.001,0.001))+
  scale_y_continuous(expand = c(0.01,0.01))

for(i in (levels(ptree$data$clade))){
   if(i == "19. Malaysian" | 
      i == "18. Far East Asia" | 
      i == "16. CHNI" | 
      i == "17. Taiwanese"){
    ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                label = paste("                    ", i, sep = ""),
                                barsize = 0,
                                offset.text = 0.0028,
                                angle = 0,
                                fontsize = 1.6,
                                color = cladecolors[which(i == levels(ptree$data$clade))],
                                align = T)
   }
  else{
    ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                              ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                              label = i,
                              barsize = 0,
                              offset.text = 0.0028,
                              angle = 0,
                              fontsize = 1.6,
                              color = cladecolors[which(i == levels(ptree$data$clade))],
                              align = T)
  }
}
ptree <- ptree + theme(legend.position = 'none')
#ptree <- ptree + theme(legend.position = c(0.12,0.9), 
#              legend.key.height = unit(0.1, "lines"),
#              legend.key.width = unit(1, "lines"),
#              legend.text = element_text(size=7), 
#              legend.title = element_text(size=7),
#              legend.background = element_rect(color = NA)) +
#  guides(fill = guide_legend(ncol=1))

ptree <- ptree + geom_hilight(node = 963, fill = "black", alpha=0.1, extend = 0.00585)+ # French dairy + African Beer
  #geom_hilight(node = 994, fill = "steelblue", alpha=0.3, extend = 0.00945)+ # African beer
  #geom_hilight(node = 1191, fill = "steelblue", alpha=0.3, extend = 0.0124)+ # Ale beer
  geom_hilight(node = 629, fill = "black", alpha=0.1, extend = 0.0067)+ # Mixed origin
  geom_hilight(node = 1066, fill = "black", alpha=0.1, extend = 0.0144) + # sake
  geom_hilight(node = 633, fill = "black", alpha=0.2, extend = 0.0115) # Mosaic beer

ptree <- ptree + scale_x_continuous(limits = c(0, 0.06))+scale_y_continuous(expand = c(0.01,0.01))
ptree
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/Fig6A.noBS.eps",device=cairo_ps, width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)

ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/Fig6A.v2.png", width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)
ptree + scale_x_continuous(limits = c(0, 0.06))+scale_y_continuous(expand = c(0.01,0.01))
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/Fig6A.pdf", width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)


p2tree <- ggtree(y,size=0.2,layout="rectangular") + 
  geom_tree(aes(color=ssd1geno), show.legend = FALSE, size=0.4)+
  scale_color_manual(values=c("#de77ae","#4d9221","#1f78b4"))+
  geom_point2(aes(subset = label == "100", label = label), size=1.5, alpha = 0.5,color = "black") +
  new_scale_color()+
  geom_tiplab(size=2,aes(color=clade),align = TRUE,linesize = 0, show.legend = FALSE, x=0.0407) +
  scale_colour_manual(values = cladecolors)+
  geom_point2(aes(shape=amplification),size=1.5, color = 'black',na.rm = TRUE, x = 0.0417)+
  scale_shape_manual(values=15, na.translate = FALSE)

sakedf <- traitsdf[traitsdf$clade == "25. Sake",]
chrdf <- data.frame(chrI = sakedf$chrI , chrII = sakedf$chrII, chrIII = sakedf$chrIII, chrIV = sakedf$chrIV, chrV = sakedf$chrV, chrVI = sakedf$chrVI, chrVII = sakedf$chrVII, 
                    chrVIII = sakedf$chrVIII, chrIX = sakedf$chrIX, chrX = sakedf$chrX, chrXI = sakedf$chrXI, chrXII = sakedf$chrXII, chrXIII = sakedf$chrXIII, chrXIV = sakedf$chrXIV, 
                    chrXV = sakedf$chrXV, chrXVI = sakedf$chrXVI)
rownames(chrdf) <- row.names(sakedf)

p2tree <- gheatmap(p2tree, 
                  chrdf, 
                  offset = 0, 
                  width = 0.1, 
                  colnames_position = "top",
                  font.size = 1,
                  colnames_offset_y = 0,
                  colnames_angle = 45) +
  scale_fill_manual(values = c('#fee0d2','#fc9272','#de2d26'), na.translate = FALSE) +
  labs(fill = "Copy Number") + 
  scale_x_continuous(expand = c(0.001,0.001))+
  scale_y_continuous(expand = c(0.01,0.01))

p2tree

p2tree <- p2tree + theme(legend.position = c(0.038,0.85), 
                       legend.key.height = unit(0.1, "lines"),
                       legend.key.width = unit(1, "lines"),
                       legend.text = element_text(size=7), 
                       legend.title = element_text(size=7),
                       legend.background = element_rect(color = NA)) +
  guides(fill = guide_legend(ncol=1))

saketree <- viewClade(p2tree, 1066) + scale_x_continuous(limits = c(0, 0.2))+scale_y_continuous(expand = c(1,1))+
  geom_treescale(x=0.038,y=122, width=0.001, offset = 0.5, linesize = 0.2, fontsize=3)
saketree
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS1.png", width=10, height = 5, units="in",limitsize = FALSE, dpi = 600)







basetree <- ggtree(y,size=0.2,layout="rectangular") + 
  geom_tree(aes(color=ssd1geno), show.legend = FALSE, size=0.4)+
  scale_color_manual(values=c("#de77ae","#4d9221","#1f78b4"))+
  geom_treescale(width=0.005, offset = 0., linesize = 0.2, fontsize=3) + 
  labs(color="ssd1") +
  new_scale_color()+
  geom_tiplab(size=3,aes(color=clade),align = FALSE,linesize = 0, show.legend = FALSE, x = 0.15) +
  scale_colour_manual(values = cladecolors)+
  guides(colour = guide_legend(override.aes = list(size=3, label = ""))) + 
  geom_strip("RIB6007",
             "CBS2270",
             label = "25. Sake",
             barsize = 0,
             offset.text = 0,
             angle = 0,
             fontsize = 3,
             color = cladecolors[which("25. Sake" == levels(ptree$data$clade))],
             align = T)
basetree <- gheatmap(basetree, 
                  ampdf, 
                  offset = 0, 
                  width = 0.05, 
                  colnames_position = "top",
                  font.size = 1,
                  colnames_offset_y = 4,
                  colnames_angle = 45) +
  scale_fill_manual(values = c("black")) +
  theme(legend.position = "none")
basetree



saketree <- viewClade(basetree, 1066)+geom_treescale(x=0.14,y=135, width=0.01, offset = 0.2, linesize = 0.2, fontsize=3)
saketree

ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/629_iqtree_sake_amphm.png", width=3, height = 5, units="in",limitsize = FALSE, dpi = 600)


#ssd1df <- data.frame(ssd1 = traitsdf$ssd1geno)
#rownames(ssd1df) <- row.names(traitsdf)
#ptree <- ptree+ new_scale_fill()
#ptree <- gheatmap(ptree,
#                  ssd1df,
#                  offset = 0.0125,
#                  width = 0.05, 
#                  colnames_position = "top",
#                  font.size = 1,
#                  colnames_offset_y = 3,
#                  colnames_angle = 45)+
#  scale_fill_manual(values=c("#de77ae","#4d9221","#1f78b4"))+
#  labs(fill = "SSD1") + 
#  scale_x_continuous(expand = c(0.001,0.001))+
#  scale_y_continuous(expand = c(0.01,0.01))


#ecodf <- data.frame(eco = traitsdf$eco)
#rownames(ecodf) <- row.names(traitsdf)
#ptree <- ptree + new_scale_fill()
#ptree <- gheatmap(ptree,
#                  ecodf, 
#                  offset = 0.0175, 
#                  width = 0.05, 
#                  colnames_position = "top",
#                  font.size = 1,
#                  colnames_offset_y = 2.5,
#                  colnames_angle = 45) +
#  scale_fill_manual(values = c(ecocolors)) +
#  labs(fill = "ecology") + 
#  scale_x_continuous(expand = c(0.001,0.001))+
#  scale_y_continuous(expand = c(0.01,0.01))
