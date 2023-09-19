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
petertab8 <- petertab7[petertab7$Ploidy != 1,]
nrow(petertab8)
petertab9 <- petertab8[petertab8$basename != "CBS1593",]
nrow(petertab9)
petertab10 <- petertab9[petertab9$basename != "CBS382",]
nrow(petertab10)

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


############################################################################################################################################
#################################### Figure S3B - rectangular ##############################################################################
############################################################################################################################################
peter_clades <- petertab10
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
traitsdf <- cbind(ampdf, eco, clade, chr_amp,ssd1, ploidy,ssd1geno, chr)

# Make node numeric
traitsdf$node <- as.numeric(traitsdf$node)
chrdf <- data.frame(`1` = traitsdf$chrI , `2` = traitsdf$chrII, `3` = traitsdf$chrIII, `4` = traitsdf$chrIV, `5` = traitsdf$chrV, `6` = traitsdf$chrVI, `7` = traitsdf$chrVII, 
                    `8` = traitsdf$chrVIII, `9` = traitsdf$chrIX, `10` = traitsdf$chrX, `11` = traitsdf$chrXI, `12` = traitsdf$chrXII, `13` = traitsdf$chrXIII, `14` = traitsdf$chrXIV, 
                    `15` = traitsdf$chrXV, `16` = traitsdf$chrXVI, check.names = FALSE)
rownames(chrdf) <- row.names(traitsdf)
# Add traits to tree
y <- full_join(rtree, traitsdf, by = 'node')
# Group root nodes
EM14S01_3Bnode <- which(rownames(newdf) == "EM14S01_3B")
EN14S01node <- which(rownames(newdf) == "EN14S01")
GE14S01_7Bnode<- which(rownames(newdf) == "GE14S01_7B")
m<- MRCA(rtree, EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode)
y <- groupClade(y, m)
y <- full_join(y, traitsdf, by ='node')

S3Arect <- ggtree(y,layout="rectangular",size=0) + 
  geom_tree(aes(linetype = group),show.legend = FALSE, size=0.5)+
  geom_treescale(x=0.018,y=165, width=0.01, offset = 5, linesize = 0.2, fontsize=3) + 
  geom_point2(aes(subset = label == "100", label = label), size=2, alpha = 0.5, color = "darkgreen",shape =16, show.legend = FALSE) +
  new_scale_color()+
  geom_tiplab(size=0.4,aes(color=clade),align = TRUE,linesize = 0.05, linetype = "dotted",show.legend = FALSE, x=0.0525) +
  scale_colour_manual(values = cladecolors)

S3Arect$data[S3Arect$data$node %in% c(EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode), "x"] <- mean(S3Arect$data$x)
S3Arect$data[S3Arect$data$node %in% c(622, 1031, 1241), "x"] <- 0.018
S3Arect$data <- S3Arect$data %>% mutate(group = replace(group, label == "GE14S01_7B", 1))

S3Arect <- gheatmap(S3Arect, 
                  chrdf, 
                  offset = 0.0015,
                  width = 0.07, 
                  colnames_position = "top",
                  font.size = 0.6,
                  colnames_offset_y = 1,
                  colnames_angle = 0) +
  scale_fill_manual(values = c('#fff5f0','#fc9272','#ef3b2c','#a50f15','black')) +
  labs(fill = "Chr copy number") + 
  scale_x_continuous(expand = c(0.001,0.001))+
  scale_y_continuous(expand = c(0.01,0.01)) +
  theme(legend.position = c(0.2,0.8), 
        legend.key.height = unit(0.2, "lines"),
        legend.key.width = unit(1.5, "lines"),
        legend.text = element_text(size=10), 
        legend.title = element_text(size=10),
        legend.background = element_rect(color = NA)) +
  guides(fill = guide_legend(ncol=1))
  #+ theme(legend.position = 'none')

for(i in (levels(S3Arect$data$clade))){
  if(i == "19. Malaysian" | 
     i == "16. CHNI"){
    S3Arect <- S3Arect + geom_strip(S3Arect$data[S3Arect$data$y == min(S3Arect$data[S3Arect$data$clade == i,]$y, na.rm = T),]$label,
                                    S3Arect$data[S3Arect$data$y == max(S3Arect$data[S3Arect$data$clade == i,]$y, na.rm = T),]$label,
                                label = paste("                    ", i, sep = ""),
                                barsize = 0,
                                offset.text = 0.0034,
                                angle = 0,
                                fontsize = 1.6,
                                color = cladecolors[which(i == levels(S3Arect$data$clade))],
                                align = T)
  }
  else{
    S3Arect <- S3Arect + geom_strip(S3Arect$data[S3Arect$data$y == min(S3Arect$data[S3Arect$data$clade == i,]$y, na.rm = T),]$label,
                                    S3Arect$data[S3Arect$data$y == max(S3Arect$data[S3Arect$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.0034,
                                angle = 0,
                                fontsize = 1.6,
                                color = cladecolors[which(i == levels(S3Arect$data$clade))],
                                align = T)
  }
}
S3Arect <- S3Arect +
  geom_hilight(node = 963, fill = "black", alpha=0.1, extend = 0.00675)+ # French dairy + African Beer
  geom_hilight(node = 629, fill = "black", alpha=0.1, extend = 0.0076)+ # Mixed origin
  geom_hilight(node = 1066, fill = "black", alpha=0.1, extend =0.0153) + # sake
  geom_hilight(node = 633, fill = "black", alpha=0.2, extend = 0.0124) + # Mosaic beer
  scale_x_continuous(limits = c(0.018, 0.059))



ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS3A-rect.png", width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS3A-rect.eps",device = cairo_ps, width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)

############################################################################################################################################
#################################### Figure S3B - rectangular ##############################################################################
############################################################################################################################################

tokeep <- fread("Documents/GitHub/eduardo/aneuploidy/peter/distmatrix/621strains/outfile.txt", header = FALSE)
tokeep2 <- fread("Documents/GitHub/eduardo/aneuploidy/peter/distmatrix/outfile-0.000007.txt", header = FALSE)

petertab11 <- petertab10[petertab10$basename %in% tokeep$V1]

nrow(petertab11)
peter_clades <- petertab11
nrow(peter_clades)

tree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/tree/IQ-Treee.453.noamb.WG.gene.tree")
rtree <- root(tree, c("EM14S01_3B","GE14S01_7B"))
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
traitsdf <- cbind(ampdf, eco, clade, chr_amp,ssd1, ploidy,ssd1geno, chr)
# Make node numeric
traitsdf$node <- as.numeric(traitsdf$node)
chrdf <- data.frame(`1` = traitsdf$chrI , `2` = traitsdf$chrII, `3` = traitsdf$chrIII, `4` = traitsdf$chrIV, `5` = traitsdf$chrV, `6` = traitsdf$chrVI, `7` = traitsdf$chrVII, 
                    `8` = traitsdf$chrVIII, `9` = traitsdf$chrIX, `10` = traitsdf$chrX, `11` = traitsdf$chrXI, `12` = traitsdf$chrXII, `13` = traitsdf$chrXIII, `14` = traitsdf$chrXIV, 
                    `15` = traitsdf$chrXV, `16` = traitsdf$chrXVI, check.names = FALSE)
rownames(chrdf) <- row.names(traitsdf)
# Add traits to tree
y <- full_join(rtree, traitsdf, by = 'node')
# Group root nodes
EM14S01_3Bnode <- which(rownames(newdf) == "EM14S01_3B")
GE14S01_7Bnode<- which(rownames(newdf) == "GE14S01_7B")
m<- MRCA(rtree, EM14S01_3Bnode,  GE14S01_7Bnode)
y <- groupClade(y, m)
y <- full_join(y, traitsdf, by ='node')




S3Brect <- ggtree(y,size=0,layout="rectangular") + 
  geom_tree(aes(linetype = group), show.legend = FALSE, size=0.5)+
  scale_color_manual(values=c("#de77ae","#4CD000","#50B9FF"))+
  geom_treescale(x=0.018,y=120, width=0.01, offset = 5, linesize = 0.2, fontsize=3) + 
  geom_point2(aes(subset = label == "100", label = label), size=1.5, alpha = 0.3, color = "black",shape =16, show.legend = FALSE) +
  labs(color="ssd1") +
  new_scale_color()+
  geom_tiplab(size=0.4,aes(color=clade),align = TRUE,linesize = 0.05, linetype = "dotted",show.legend = FALSE, x=0.0482) +
  scale_colour_manual(values = cladecolors)


S3Brect$data[S3Brect$data$node %in% c(EM14S01_3Bnode, GE14S01_7Bnode), "x"] <- mean(S3Brect$data$x)
S3Brect$data[S3Brect$data$node %in% c(454, 905), "x"] <- 0.018
S3Brect$data <- S3Brect$data %>% mutate(group = replace(group, label == "GE14S01_7B", 1))


S3Brect <- gheatmap(S3Brect, 
                    chrdf, 
                    offset = 0.0012,
                    width = 0.07, 
                    colnames_position = "top",
                    font.size = 0.6,
                    colnames_offset_y = 1,
                    colnames_angle = 0) +
  scale_fill_manual(values = c('#fff5f0','#fc9272','#ef3b2c','#a50f15','black')) +
  labs(fill = "amplification") + 
  scale_x_continuous(expand = c(0.001,0.001))+
  scale_y_continuous(expand = c(0.01,0.01))


for(i in (levels(S3Brect$data$clade))){
  if(i == "19. Malaysian" | 
     i == "16. CHNI"){
    S3Brect <- S3Brect + geom_strip(S3Brect$data[S3Brect$data$y == min(S3Brect$data[S3Brect$data$clade == i,]$y, na.rm = T),]$label,
                                    S3Brect$data[S3Brect$data$y == max(S3Brect$data[S3Brect$data$clade == i,]$y, na.rm = T),]$label,
                                    label = paste("                    ", i, sep = ""),
                                    barsize = 0,
                                    offset.text = 0.0034,
                                    angle = 0,
                                    fontsize = 1.6,
                                    color = cladecolors[which(i == levels(S3Brect$data$clade))],
                                    align = T)
  }
  else{
    S3Brect <- S3Brect + geom_strip(S3Brect$data[S3Brect$data$y == min(S3Brect$data[S3Brect$data$clade == i,]$y, na.rm = T),]$label,
                                    S3Brect$data[S3Brect$data$y == max(S3Brect$data[S3Brect$data$clade == i,]$y, na.rm = T),]$label,
                                    label = i,
                                    barsize = 0,
                                    offset.text = 0.0031,
                                    angle = 0,
                                    fontsize = 1.6,
                                    color = cladecolors[which(i == levels(S3Brect$data$clade))],
                                    align = T)
  }
}

S3Brect <- S3Brect + theme(legend.position = 'none')
#S3Brect <- S3Brect + theme(legend.position = c(.15,0.4), 
#                           legend.key.height = unit(0.1, "lines"),
#                           legend.key.width = unit(1, "lines"),
#                           legend.text = element_text(size=7), 
#                           legend.title = element_text(size=7),
#                           legend.background = element_rect(color = NA)) +
#  guides(fill = guide_legend(ncol=1))

S3Brect <- S3Brect + 
  #geom_hilight(node = 830, fill = "black", alpha=0.2, extend = 0.00607)+ # French dairy + African Beer
  geom_hilight(node = 460, fill = "black", alpha=0.1, extend = 0.0066)+ # Mixed origin
  geom_hilight(node = 748, fill = "black", alpha=0.1, extend = 0.01385)  # sake
  #geom_hilight(node = 463, fill = "black", alpha=0.2, extend = 0.0108) # Mosaic beer

S3Brect <- S3Brect + scale_x_continuous(limits = c(0.018, 0.055))+scale_y_continuous(expand = c(0.01,0.01))
ggarrange(S3Arect, S3Brect)
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS3-rect.eps",device = cairo_ps, width=15, height = 15, units="in",limitsize = FALSE, dpi = 600)

ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS3B-rect.png", width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS3B-rect.eps",device = cairo_ps, width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)



############################################################################################################################################
#################################### Figure S3A - circular #################################################################################
############################################################################################################################################


S3Acirc <- ggtree(y,size=0,layout="circular", branch.length = 'none') + 
  geom_tree(aes(linetype = group), show.legend = FALSE, size=0.2)+
  #geom_treescale(x=0.001,y=315, width=0.01, offset = -2, linesize = 0.2, fontsize=2) + 
  geom_point2(aes(subset = label == "100", label = label), size=1.5, alpha = 0.5, color = "darkgreen",shape =16, show.legend = FALSE) +
  geom_tiplab(size=1,aes(color=clade),show.legend = FALSE) +
  #geom_tiplab(size=0.8,aes(color=clade), offset = 0.002, align = TRUE,linesize = 0.05, show.legend = FALSE) +
  scale_colour_manual(values = cladecolors)

#S3Acirc$data[S3Acirc$data$node %in% c(EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode), "x"] <- mean(S3Acirc$data$x)
#S3Acirc$data[S3Acirc$data$node %in% c(622, 1031, 1241), "x"] <- 0.02
#S3Acirc$data <- S3Acirc$data %>% mutate(group = replace(group, label == "GE14S01_7B", 1))

S3Acirc <- gheatmap(S3Acirc, 
                  chrdf, 
                  offset = 2, 
                  width = 0.4, 
                  colnames_position = "top",
                  font.size = 0.5,
                  colnames_offset_y = 4,
                  colnames_angle = 0) +
  scale_fill_manual(values = c('#fff5f0','#fc9272','#ef3b2c','#a50f15','black')) +
  labs(fill = "amplification")
  #scale_x_continuous(expand = c(0.001,0.001))+
  #scale_y_continuous(expand = c(0.01,0.01))


for(i in (levels(S3Acirc$data$clade))){
  if(i == "7. Mosaic beer"){
    S3Acirc <- S3Acirc + geom_strip("CBS7957",
                                "CBS7539",
                                label = i,
                                barsize = 0,
                                offset.text = 0.04,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(S3Acirc$data$clade))],
                                align = T)
  }
  else if(i == "3. Brazilian bioethanol"){
    S3Acirc <- S3Acirc + geom_strip(S3Acirc$data[S3Acirc$data$y == min(S3Acirc$data[S3Acirc$data$clade == i,]$y, na.rm = T),]$label,
                                S3Acirc$data[S3Acirc$data$y == max(S3Acirc$data[S3Acirc$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.04,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(S3Acirc$data$clade))],
                                align = T)
  }
  else if(i == "9. Mexican agave"){
    S3Acirc <- S3Acirc + geom_strip(S3Acirc$data[S3Acirc$data$y == min(S3Acirc$data[S3Acirc$data$clade == i,]$y, na.rm = T),]$label,
                                S3Acirc$data[S3Acirc$data$y == max(S3Acirc$data[S3Acirc$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.04,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(S3Acirc$data$clade))],
                                align = T)
  }
  else{
    S3Acirc <- S3Acirc + geom_strip(S3Acirc$data[S3Acirc$data$y == min(S3Acirc$data[S3Acirc$data$clade == i,]$y, na.rm = T),]$label,
                                S3Acirc$data[S3Acirc$data$y == max(S3Acirc$data[S3Acirc$data$clade == i,]$y, na.rm = T),]$label,
                                label = i,
                                barsize = 0,
                                offset.text = 0.03,
                                angle = 0,
                                fontsize = 1,
                                color = cladecolors[which(i == levels(S3Acirc$data$clade))],
                                align = T)
  }
}

S3Acirc <- S3Acirc + theme(legend.position = c(0.5,0.4), 
                       legend.key.height = unit(0.1, "lines"),
                       legend.key.width = unit(0.5, "lines"),
                       legend.text = element_text(size=3), 
                       legend.title = element_text(size=3),
                       legend.background = element_rect(color = NA)) +
  guides(fill = guide_legend(ncol=1))

S3Arect <- S3Arect +
  geom_hilight(node = 963, fill = "black", alpha=0.1, extend = 0.00675)+ # French dairy + African Beer
  geom_hilight(node = 629, fill = "black", alpha=0.1, extend = 0.0076)+ # Mixed origin
  geom_hilight(node = 1066, fill = "black", alpha=0.1, extend =0.0153) + # sake
  geom_hilight(node = 633, fill = "black", alpha=0.2, extend = 0.0124)  # Mosaic beer

S3Acirc <- S3Acirc +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS3A-circ.png", width=5, height = 5, units="in",limitsize = FALSE, dpi = 600)



























