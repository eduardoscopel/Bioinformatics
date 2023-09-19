library(ggtree)
library(ape)
library(treeio)
library(data.table)
library(dplyr)
library(ggrepel)
library(cowplot)
library(ggnewscale)


### Read tables and match
petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /Scopel_etalSuppTable.txt",header=TRUE)
#sctab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /sc_table_most_recent.txt", header = TRUE)
nrow(petertab)
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="Euploid","No","Yes")
petertab$amp_response <- ifelse(petertab$Aneuploidy_type == "Euploid",NA,"Yes")
petertab$amp_binary <- ifelse(petertab$amp_response == "Yes", 1,0)
petertab$admixture <- ifelse(petertab$Clades == "M. Mosaic" |
                               petertab$Clades == "8. Mixed origin",
                             "Admixed","Non-admixed")
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
petertab8 <- petertab7[petertab7$Clades != "7. Mosaic beer",]
nrow(petertab8)
peter_clades <- petertab8
nrow(peter_clades)
tokeep <- fread("Documents/GitHub/eduardo/aneuploidy/peter/distmatrix/outfile-0.000007.txt", header = FALSE)
petertab9 <- petertab8[petertab8$basename %in% tokeep$V1]
nrow(petertab9)


tree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/tree/tree5.treefile")
ssd1tree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/ssd1/iq-tree-619.ssd1.prot.treefile")
raxml_file <- system.file("extdata/RAxML","/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/ssd1/RAxML.619.ssd1.tree", package='treeio')
#raxml <- read.raxml("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/ssd1/RAxML.619.ssd1.tree")
#rtree <- raxml
rtree <- root(tree, c("EM14S01_3B","GE14S01_7B"))
rtree <- root(ssd1tree, c("EM14S01_3B","EN14S01","GE14S01_7B"))
# Create a data frame with strain IDs as row names, and add traits as columns (clade, ecology, aneuploidy, ploidy)
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
    }
  }
}

# Create a matrix for each trait
amplification <- as.matrix(newdf)[,1]
eco <- as.matrix(newdf)[,2]
clade <- as.matrix(newdf)[,3]
chr <- as.matrix(newdf)[,4]
ssd1 <- as.matrix(newdf)[,5]
ploidy <- as.matrix(newdf)[,6]
# Create a data frame to merge node ID with traits
ampdf <- data.frame(node = nodeid(rtree, names(amplification)),
                    amplification = amplification)
traitsdf <- cbind(ampdf, eco, clade, chr,ssd1, ploidy)
# Make node numeric
traitsdf$node <- as.numeric(traitsdf$node)
# Add traits to tree
rtree <- full_join(rtree, traitsdf, by = 'node')
# Group root nodes
EM14S01_3Bnode <- which(rownames(newdf) == "EM14S01_3B")
#EN14S01node <- which(rownames(newdf) == "EN14S01")
GE14S01_7Bnode<- which(rownames(newdf) == "GE14S01_7B")
#m<- MRCA(rtree, EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode)
m<- MRCA(rtree, EM14S01_3Bnode, GE14S01_7Bnode)
y <- groupClade(rtree, m)
y <- full_join(y, traitsdf, by ='node')
y <- rtree
ecodf <- data.frame(traitsdf$eco)
rownames(ecodf) <- row.names(traitsdf)
ampdf <- data.frame(traitsdf$amplification)
rownames(ampdf) <- row.names(traitsdf)
test <- cbind(ecodf, ampdf)

darkred <- "#8B0000"
red <- "#FF0000"
darksalmon <- "#E9967A"
darkorange <- "#FF8C00"
darkgolden <- "#FF8C00"
darkkhaki <- "#BDB76B"
olive <- "#808000"
yellow <- "#FFFF00"
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
ab6 <- yellow
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

cladecolors <- c(we1, fgh10,ab11,wac12,apw13, ch14, ch15, ch16, t17, fea18, m19, ol2, ch20, e21, fer22, nao23, ai24,s25,af26,bb3, mo4, fd5, ab6, mo8, ma9)

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
#ecocolors <- c(bak, bee, bio, cid, dai, dist, ferm, flow, fruit, hum, humc, ind, ins, lab, nat, palm, sak, soil, tr, unk, wat, wine)

# Plot tree

ptree <- ggtree(y,size=0.2,layout="rectangular") + 
  geom_tree(aes(color=ssd1))+
  scale_color_manual(values=c("magenta","blue","gray"))+
  geom_treescale(x=0.001,y=150,width=0.0005, offset = 5, linesize = 1, fontsize=5) + 
  geom_nodelab(aes(label = label),size=1.5, hjust=1, vjust=0.5) +
  labs(color="ssd1") +
  new_scale_color()+
  geom_tiplab(size=0.8,aes(color=clade), offset = 0.002, align = TRUE,linesize = 0.05, show.legend = FALSE) +
  scale_colour_manual(values = c(cladecolors,"magenta","blue","gray"))+
  guides(colour = guide_legend(override.aes = list(size=3, label = "")))+
#  geom_tippoint(aes(shape=amplification),size=2, color = '#D2691E')+
  geom_point2(aes(shape=amplification),size=1, color = 'black',na.rm = TRUE, alpha = 0.6, x = 0.0495)+
  scale_shape_manual(values=16)+
   geom_tiplab(aes(label=chr),geom='text',size=1.5,color="black") 
ptree
#ptree + theme_classic()
# Change root position so it's pushed to the right of the tree
#ptree$data[ptree$data$node %in% c(EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode), "x"] <- mean(ptree$data$x)
#dftree[dftree$x < 0.02,"x"] <- mean(ptree$data$x)
#ptree$data[335,"x"] <- mean(ptree$data$x)
#ptree$data[336,"x"] <- mean(ptree$data$x)
#ptree$data[337,"x"] <- mean(ptree$data$x)
#ptree$data[620,"x"] <- mean(ptree$data$x)
#ptree$data[975,"x"] <- mean(ptree$data$x)
#ptree$data[1237,"x"] <- mean(ptree$data$x)

subtree <- ptree

tdf <- ptree$data
tdf <- tdf[!tdf$isTip,]
tdf$label <- as.numeric(tdf$label)
tdf <- tdf[tdf$label>80,]

ptree <- ptree + 
  geom_text2(data = tdf, aes(label = label), size=1, hjust=1)

for(i in (levels(ptree$data$clade))){
  #curr_color <- cladecolors[which(i == levels(ptree$data$clade))]
  #print(i)
  #print(curr_color)
  ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                              ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
                              label = i,
                              barsize = 0,
                              offset.text = 0.0023,
                              fontsize = 1.5,
                              color = cladecolors[which(i == levels(ptree$data$clade))],
                              align = T)
}
#ptree +theme_classic()

names(ecodf)[1] <- "eco"
names(test)[2] <- "chr_amp"

ptree <- gheatmap(ptree, 
         ecodf, 
         offset = 0.0055, 
         width = .005, 
         colnames = TRUE, 
         colnames_position = "top",
         colnames_offset_y = 1,
         font.size = 1.5,
         hjust=0.5) +
  scale_fill_manual(values = c(ecocolors)) +
  labs(fill = "ecology") + 
  scale_x_continuous(expand = c(0.001,0.001))+
  scale_y_continuous(expand = c(0.01,0.01))
ptree <- ptree + theme(legend.position = c(0.2, 0.5), legend.text = element_text(size = 20),legend.title = element_text(size = 20))
#ptree
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/rect_455_wg_tree.png", width=50, height = 60, units="cm",limitsize = FALSE)
ptree + theme_classic()
