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

for(i in 1:nrow(petertab)){
  petertab$chrI[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i]),
                             ifelse(grepl("-1*1;", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])-1,
                                    ifelse(grepl("-1*1\n", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])-1,
                             ifelse(grepl("+1*1;", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+1,
                                    ifelse(grepl("+1*1\n", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+1,
                             ifelse(grepl("+2*1;", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+2,
                                    ifelse(grepl("+2*1\n", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+2,       
                             ifelse(grepl("+3*1;", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+3,
                                    ifelse(grepl("+3*1\n", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+3,
                             ifelse(grepl("+4*1;", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+4,
                                    ifelse(grepl("+4*1\n", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i])))))))))))
  petertab$chrII[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i]),
                              ifelse(grepl("-1*2", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])-1,
                             ifelse(grepl("+1*2", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+1,
                             ifelse(grepl("+2*2", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+2,
                             ifelse(grepl("+3*2", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+3,
                             ifelse(grepl("+4*2", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrIII[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i]),
                             ifelse(grepl("-1*3", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])-1,
                             ifelse(grepl("+1*3", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+1,
                             ifelse(grepl("+2*3", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+2,
                             ifelse(grepl("+3*3", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+3,
                             ifelse(grepl("+4*3", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrIV[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i]),
                              ifelse(grepl("-1*4", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])-1,
                              ifelse(grepl("+1*4", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+1,
                              ifelse(grepl("+2*4", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+2,
                              ifelse(grepl("+3*4", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+3,
                              ifelse(grepl("+4*4", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrV[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i]),
                             ifelse(grepl("-1*5", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])-1,
                             ifelse(grepl("+1*5", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+1,
                             ifelse(grepl("+2*5", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+2,
                             ifelse(grepl("+3*5", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+3,
                             ifelse(grepl("+4*5", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrVI[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i]),
                              ifelse(grepl("-1*6", petertab$Aneuploidies[i], fixed = TRUE), as.numeric(petertab$Ploidy[i])-1,
                              ifelse(grepl("+1*6", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                              ifelse(grepl("+2*6", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                              ifelse(grepl("+3*6", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                              ifelse(grepl("+4*6", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrVII[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                             ifelse(grepl("-1*7", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                             ifelse(grepl("+1*7", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                             ifelse(grepl("+2*7", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                             ifelse(grepl("+3*7", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                             ifelse(grepl("+4*7", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrVIII[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                              ifelse(grepl("-1*8", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                              ifelse(grepl("+1*8", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                              ifelse(grepl("+2*8", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                              ifelse(grepl("+3*8", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                              ifelse(grepl("+4*8", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrIX[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                             ifelse(grepl("-1*9", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                             ifelse(grepl("+1*9", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                             ifelse(grepl("+2*9", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                             ifelse(grepl("+3*9", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                             ifelse(grepl("+4*9", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrX[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                              ifelse(grepl("-1*10", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                              ifelse(grepl("+1*10", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                              ifelse(grepl("+2*10", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                              ifelse(grepl("+3*10", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                              ifelse(grepl("+4*10", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrXI[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                             ifelse(grepl("-1*11", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                             ifelse(grepl("+1*11", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                             ifelse(grepl("+2*11", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                             ifelse(grepl("+3*11", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                             ifelse(grepl("+4*11", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrXII[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                              ifelse(grepl("-1*12", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                              ifelse(grepl("+1*12", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                              ifelse(grepl("+2*12", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                              ifelse(grepl("+3*12", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                              ifelse(grepl("+4*12", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrXIII[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                             ifelse(grepl("-1*13", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                             ifelse(grepl("+1*13", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                             ifelse(grepl("+2*13", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                             ifelse(grepl("+3*13", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                             ifelse(grepl("+4*13", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrXIV[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                              ifelse(grepl("-1*14", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                              ifelse(grepl("+1*14", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                              ifelse(grepl("+2*14", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                              ifelse(grepl("+3*14", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                              ifelse(grepl("+4*14", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrXV[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                             ifelse(grepl("-1*15", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                             ifelse(grepl("+1*15", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                             ifelse(grepl("+2*15", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                             ifelse(grepl("+3*15", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                             ifelse(grepl("+4*15", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
  petertab$chrXVI[i] <- ifelse(grepl("euploid", petertab$Aneuploidies[i]),as.numeric(petertab$Ploidy[i]),
                              ifelse(grepl("-1*16", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])-1,
                              ifelse(grepl("+1*16", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+1,
                              ifelse(grepl("+2*16", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+2,
                              ifelse(grepl("+3*16", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+3,
                              ifelse(grepl("+4*16", petertab$Aneuploidies[i], fixed = TRUE),as.numeric(petertab$Ploidy[i])+4,petertab$Ploidy[i]))))))
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
peter_clades <- petertab7
nrow(peter_clades)


tree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/tree/IQ-Tree.629.WG.gene.tree")
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
      newdf$chrVIII[i] <- peter_clades$chrVII[j]
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
traitsdf <- cbind(ampdf, eco, clade, chr_amp,ssd1, ploidy,ssd1geno, chr)
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
bb3 <- indigo
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


ptree <- ggtree(y,size=0.2,layout="circular") + 
  geom_tree(aes(color=ssd1geno), show.legend = FALSE, size=0.2)+
  scale_color_manual(values=c("#de77ae","#4d9221","#1f78b4"))+
  #geom_tree(aes(color=clade), show.legend = FALSE)+
  #scale_color_manual(values=cladecolors)+
  geom_treescale(x=0.05,y=315, width=0.05, offset = -2, linesize = 0.2, fontsize=2) + 
  #geom_text2(aes(subset = node> 620 , label = label), size=1) +
  labs(color="ssd1") +
  new_scale_color()+
  geom_tiplab(size=0.4,aes(color=clade),align = TRUE,linesize = 0.1, show.legend = FALSE, x=0.1815) +
  #geom_tiplab(size=0.8,aes(color=clade), offset = 0.002, align = TRUE,linesize = 0.05, show.legend = FALSE) +
  scale_colour_manual(values = cladecolors)+
  guides(colour = guide_legend(override.aes = list(size=3, label = "")))

#ptree$data[ptree$data$node %in% c(EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode), "x"] <- mean(ptree$data$x)
#ptree$data[ptree$data$node %in% c(630, 1001, 1257), "x"] <- 0.08
#ptree$data <- ptree$data %>% mutate(group = replace(group, label == "GE14S01_7B", 1))

#ampdf <- data.frame(chr_amp = traitsdf$amplification)
#rownames(ampdf) <- row.names(traitsdf)
#ptree <- gheatmap(ptree, 
#                  ampdf, 
#                  offset = 0.02, 
#                  width = 0.05, 
#                  colnames_position = "top",
#                  font.size = 1,
#                  colnames_offset_y = 4,
#                  colnames_angle = 45) +
#  scale_fill_manual(values = c("black")) +
#  labs(fill = "amplification") + 
#  scale_x_continuous(expand = c(0.001,0.001))+
#  scale_y_continuous(expand = c(0.01,0.01))

chrdf <- data.frame(chrI = traitsdf$chrI , chrII = traitsdf$chrII, chrIII = traitsdf$chrIII, chrIV = traitsdf$chrIV, chrV = traitsdf$chrV, chrVI = traitsdf$chrVI, chrVII = traitsdf$chrVII, 
                    chrVIII = traitsdf$chrVIII, chrIX = traitsdf$chrIX, chrX = traitsdf$chrX, chrXI = traitsdf$chrXI, chrXII = traitsdf$chrXII, chrXIII = traitsdf$chrXIII, chrXIV = traitsdf$chrXIV, 
                    chrXV = traitsdf$chrXV, chrXVI = traitsdf$chrXVI)
rownames(chrdf) <- row.names(traitsdf)
ptree <- gheatmap(ptree, 
                  chrdf, 
                  offset = 0.025, 
                  width = 0.4, 
                  colnames_position = "top",
                  font.size = 0.5,
                  colnames_offset_y = 4,
                  colnames_angle = 0) +
  scale_fill_manual(values = c('#feedde','#fdd0a2','#fdae6b','#fd8d3c','#e6550d','#a63603')) +
  labs(fill = "amplification") + 
  scale_x_continuous(expand = c(0.001,0.001))+
  scale_y_continuous(expand = c(0.01,0.01))


for(i in (levels(ptree$data$clade))){
  if(i == "7. Mosaic beer"){
    ptree <- ptree + geom_strip("DBVPG6590",
                                "CBS7539",
                                label = i,
                                barsize = 0,
                                offset.text = 0.15,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(ptree$data$clade))],
                                align = T)
  }
  else if(i == "15. CHNII"){
    ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.16,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(ptree$data$clade))],
                                align = T)
  }
  else if(i == "3. Brazilian bioethanol"){
    ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.2,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(ptree$data$clade))],
                                align = T)
  }
  else if(i == "8. Mixed origin"){
    ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.15,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(ptree$data$clade))],
                                align = T)
  }
  else if(i == "9. Mexican agave"){
    ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.14,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(ptree$data$clade))],
                                align = T)
  }
  else if(i == "12. West African cocoa"){
    ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.13,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(ptree$data$clade))],
                                align = T)
  }
  else{
    ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.14,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(ptree$data$clade))],
                                align = T)
  }
}

ptree <- ptree + theme(legend.position = c(0.5,0.5), 
                       legend.key.height = unit(0.1, "lines"),
                       legend.key.width = unit(0.5, "lines"),
                       legend.text = element_text(size=3), 
                       legend.title = element_text(size=3),
                       legend.background = element_rect(color = NA)) +
  guides(fill = guide_legend(ncol=1))

ptree <- ptree + geom_hilight(node = 1037, fill = "steelblue", alpha=0.3, extend = 0.055)+
  geom_hilight(node = 1206, fill = "steelblue", alpha=0.3, extend = 0.046)+
  geom_hilight(node = 1193, fill = "steelblue", alpha=0.3, extend = 0.046)+
  geom_hilight(node = 636, fill = "steelblue", alpha=0.3, extend = 0.029)+
  geom_hilight(node = 1163, fill = "steelblue", alpha=0.3, extend = 0.036)
ptree <- ptree +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
ptree
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/629_iqtree_wg_circ_chrhm.png", width=5, height = 5, units="in",limitsize = FALSE, dpi = 600)








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



ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/629_iqtree_wg_chrheat.png", width=3, height = 8, units="in",limitsize = FALSE, dpi = 600)

basetree <- ggtree(y,size=0.2,layout="rectangular") + 
  geom_tree(aes(color=ssd1geno), show.legend = FALSE, size=0.4)+
  scale_color_manual(values=c("#de77ae","#4d9221","#1f78b4"))+
  geom_treescale(width=0.01, offset = 0., linesize = 0.2, fontsize=3) + 
  labs(color="ssd1") +
  new_scale_color()+
  geom_tiplab(size=3,aes(color=clade),align = FALSE,linesize = 0, show.legend = FALSE, x = 0.15) +
  scale_colour_manual(values = cladecolors)+
  guides(colour = guide_legend(override.aes = list(size=3, label = ""))) + 
  geom_strip("RIB1016",
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



saketree <- viewClade(basetree, 1037)+geom_treescale(x=0.14,y=135, width=0.01, offset = 0.2, linesize = 0.2, fontsize=3)
saketree

ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/629_iqtree_sake_amphm.png", width=3, height = 5, units="in",limitsize = FALSE, dpi = 600)
