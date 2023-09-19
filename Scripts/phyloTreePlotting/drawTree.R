library(ggtree)
library(ape)
library(treeio)
library(data.table)
library(dplyr)
library(ggrepel)
library(cowplot)
library(ggnewscale)

# Read strains table
petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt",header=TRUE)
#sctab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /sc_table_most_recent.txt", header = TRUE)

# Get important categories from my table into peter df
#for(i in 1:nrow(petertab)){
#  for(j in 1:nrow(sctab)){
#    if(!is.element(petertab$`Isolate name`[i],sctab$ID)){
#      petertab$heterozygosity[i] <- NA
#      petertab$new_eco[i] <- NA
#      petertab$basename[i] <- NA
#    }
#    else if(petertab$`Isolate name`[i] == sctab$ID[j]){
#      petertab$heterozygosity[i] <- sctab$heterozygosity[j]
#      petertab$new_eco[i] <- sctab$ecology_category[j]
#      petertab$basename[i] <- sctab$basename[j]
#    }
#  }
#}
petertab <- petertab[complete.cases(petertab$Clades),]
petertab <- petertab[complete.cases(petertab$basename),]

# Add a column with admixture information and remove admixed strains
petertab$admixture <- ifelse(petertab$Clades == "M. Mosaic" |
                               petertab$Clades == "8. Mixed origin" |
                               petertab$Clades == "7. Mosaic beer","Admixed","Non-admixed")
petertab <- petertab[petertab$admixture != "Admixed",]
# Remove non-diploids
petertab <- petertab[petertab$Ploidy == 2,]
# Remove highly heterozygous strains
petertab <- petertab[petertab$heterozygosity <= 0.003,]


nrow(petertab)
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="euploid",NA,"Yes")
petertab$amp_response <- ifelse(petertab$Aneuploidy_type =="Loss" | 
                                  petertab$Aneuploidy_type =="Euploid",NA,"Yes")

for(i in 1:nrow(petertab)){
  if(!is.na(petertab$amp_response[i])){
    petertab$chr_amp[i] <- substr(petertab$Aneuploidies[i],6,nchar(petertab$Aneuploidies[i])-1)
  }
  else{
    petertab$chr_amp[i] <- NA
  }
}

# Remove admixed strains
# peter$admixture <- ifelse(peter$pop_summary == "8._Mixed_origin_" |
#                            peter$pop_summary == "Mosaics" |
#                            peter$pop_summary == "7._Mosaic_beer_","Admixed","Non-admixed")
# nonadpeter <- peter[peter$admixture != "Admixed",]

# Read tree
#njtree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/tree/bionjtrees.suptree")
#alltree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/tree/bionjtrees")

tree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/tree/boot1.suptree")
rtree <- root(tree, c("EM14S01_3B","EN14S01","GE14S01_7B"))

#pall <- ggtree(alltree[1:50], branch.length = 'none') + facet_wrap(~.id, ncol =10)
#ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/alltrees.png", width=75, height = 75, units="cm",limitsize = FALSE)
#pdensi <- ggdensitree(alltree, alpha=.1, colour = "steelblue", branch.length = 'none') + geom_tiplab(size=0.5) 
#ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/densitree.png", width=30, height = 30, units="cm",limitsize = FALSE)

# Create a data frame with strain IDs as row names, and add traits as columns (clade, ecology, aneuploidy, ploidy)
newdf <- data.frame(rtree$tip.label, row.names = TRUE)

# Add important info to data frame
for(i in 1:length(rtree$tip.label)){
  for(j in 1:length(petertab$basename)){
    if(rtree$tip.label[i] == petertab$basename[j]){
      rownames(newdf)[i] <- petertab$basename[j]
      newdf$amplification[i] <- petertab$amp_response[j]
      newdf$clade[i] <- petertab$Clades[j]
      newdf$ecology[i] <- petertab$`Ecological origins`[j]
      newdf$new_eco[i] <- petertab$new_eco[j]
      newdf$chr[i] <- petertab$chr_amp[j]
      newdf$ssd1[i] <- petertab$SSD1_allele[j]
    }
  }
}

# Create a matrix for each trait
amplification <- as.matrix(newdf)[,1]
clade <- as.matrix(newdf)[,2]
eco <- as.matrix(newdf)[,3]
new_eco <- as.matrix(newdf)[,4]
chr <- as.matrix(newdf)[,5]
ssd1 <- as.matrix(newdf)[,6]
# Create a data frame to merge node ID with traits
ampdf <- data.frame(node = nodeid(rtree, names(amplification)),
                      amplification = amplification)
traitsdf <- cbind(ampdf, clade, eco, new_eco,chr,ssd1)
# Make node numeric
traitsdf$node <- as.numeric(traitsdf$node)
# Add traits to tree
rtree <- full_join(rtree, traitsdf, by = 'node')
# Group root nodes
EM14S01_3Bnode <- which(rownames(newdf) == "EM14S01_3B")
EN14S01node <- which(rownames(newdf) == "EN14S01")
GE14S01_7Bnode<- which(rownames(newdf) == "GE14S01_7B")
m<- MRCA(rtree, EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode)
y <- groupClade(rtree,m)
y <- full_join(y, traitsdf, by ='node')

ecodf <- data.frame(traitsdf$eco)
rownames(ecodf) <- row.names(traitsdf)
ampdf <- data.frame(traitsdf$amplification)
rownames(ampdf) <- row.names(traitsdf)
test <- cbind(ecodf, ampdf)

# Create df with node IDs for each clade
cladeIDs <- c()
nodeIDs <- c()
for(i in 1:length(levels(traitsdf$clade))){
  cladeIDs[i] <- levels(traitsdf$clade)[i]
  nodeIDs[i] <- as.numeric(MRCA(y,traitsdf[traitsdf$clade == levels(traitsdf$clade)[i],"node"]))
}
cl_node <- cbind(cladeIDs,nodeIDs)
cl_node <- as.data.frame(cl_node)
we <- "#6600CC"
fgh <- "#FFB266"
ab<-"#CC6600"
wac<-"#FF9933"
wpw<-"#FF8000"
ciii<-"#336600"
cii<-"#006600"
ci<-"#006633"
t<-"#009999"
fea<-"#00CCCC"
m<-"#00FFFF"
ol<-"#000066"
cv<-"#00994C"
e<-"#4C9900"
fer<-"#66CC00"
nao<-"#80FF00"
ai<-"#99FFFF"
s<-"#FF0000"
af<-"#FF9999"
bb<-"#66B2FF"
mo<-"#333300"
fd<-"#FF66FF"
ab<-"#CCCC00"
ma<-"#999900"

bak <- "#FFFF00"
bee <- "#CCCC00"
bio <- "#66B2FF"
cid <- "#A46CDB"
dai <- "#FF66FF"
dist <- "#0757A8"
ferm <- "#8692DF"
flow <- "#30C0A8"
fruit <- "#6D9F8B"
hum <- "#FFB266"
humc <- "#D16900"
ind <- "#0022FD"
ins <- "#EEF5D3"
lab <- "#03005C"
nat <- "#005C03"
palm <- "#FF8000"
prob <- "#AF6BAF"
sak <- "#FF0000"
soil <- "#ABAB3B"
tr <- "#333300"
unk <- "#000000"
wat <- "#A0C0BD"
wine <- "#6600CC"

cladecolors <- c(we,fgh,ab,wac,wpw,ciii,cii,ci,t,fea,m,ol,cv,e,fer,nao,ai,s,af,bb,mo,fd,ab,ma)
ecocolors <- c(bak, bee, bio, cid, dai, dist, ferm, flow, fruit, hum, humc, ind, ins, lab, nat, palm, prob, sak, soil, tr, unk, wat, wine)
#colors <- c("#890084","#FF8101","#FFCD01","#705B05","#D000AD","#22FF00","#169B02","#69BD5C","#98E28D","#8DE2CB","#456C02",
#           "#9FFF9F","#05A68B","#00B8E6","#356433","#82E4DD","#AE6F6D","#C71D17","#0F26F3","#597951","#95540A","#BAB405",
#           "#6C6B56","#DADA06","#FFCD01","#0F26F3","#B606DA","#BAB405","#0A020B","#8A7BD8","#38DB06","#1E7F01","#FC9004",
#           "#E8BF89","#B2CAF9","#0991AC","#004DFF","#015E2F","#B603EC","#F58974","gray","#FF0000","#615F28","#093A08",
#           "black","#B3F6FA","#7F049E")

# Plot tree

ptree <- ggtree(y,size=0.3,layout="rectangular") + 
  geom_tree(aes(color=ssd1))+
  scale_color_manual(values=c("magenta","blue","gray"))+
  geom_treescale(x=-0.001,y=150,width=0.01, offset = 50, linesize =0.4, fontsize=2) + 
  geom_nodelab(size=1) +
  labs(color="ssd1")



#ptree <- ggtree(y,size=0.3,layout="circular") + 
#  geom_tree(aes(color=ssd1))+
#  scale_color_manual(values=c("magenta","blue","gray"))+
#  geom_treescale(x=-0.001,y=150,width=0.01, offset = 50, linesize =0.4, fontsize=2) + 
#  geom_nodelab(size=1) +
#  labs(color="ssd1")

# Change root position so it's pushed to the right of the tree
ptree$data[ptree$data$node %in% c(EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode), "x"] <- mean(ptree$data$x)

subtree <- ptree
ptree <- ptree + 
  new_scale_color()+
  geom_tiplab(size=0.8,aes(color=clade), offset = 0.002, align = TRUE,linesize = 0) +
  scale_colour_manual(values = c(cladecolors,"magenta","blue","gray"))+
  guides(colour = guide_legend(override.aes = list(size=3, label = "")))
ptree <- ptree + 
  geom_tippoint(aes(shape=amplification),size=1, color = '#D2691E')
ptree <- ptree + 
  geom_tiplab(aes(label=chr),geom='text',size=1,color="#D2691E")

tdf <- ptree$data
tdf <- tdf[!tdf$isTip,]
tdf$label <- as.numeric(tdf$label)
tdf <- tdf[tdf$label>80,]

#ptree <- gheatmap(ptree, ecodf, offset = .005, width = .02, colnames = FALSE) + 
#  scale_fill_manual(values = ecocolors) + 
#  labs(fill = "ecology")
ptree <- ptree + 
  geom_text2(data = tdf, aes(label = label), size=1, hjust=1)

ptree <- gheatmap(ptree, test, offset = 0.006, width = .02, colnames = FALSE) +
  scale_fill_manual(values = c(ecocolors,"#D2691E")) + 
  labs(fill = "ecology")
ptree <- ptree + 
  theme(legend.position = "bottom")

#for(i in (unique(ptree$data$clade))){
#  ptree <- ptree + geom_strip(ptree$data[ptree$data$y == min(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
#                              ptree$data[ptree$data$y == max(ptree$data[ptree$data$clade == i,]$y, na.rm = T),]$label,
#                              barsize = 0,
#                              label = i,
#                              offset.text = 0.02,
#                              align = T)
#}
#ptree + geom_strip('CCY_21_4_119','Lib73', barsize=0, color = we, label = "Wine/EU", offset.text = 0.01)

#ptree <- ptree + theme(legend.position = 'none')
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/rect_tree.png", width=35, height = 25, units="cm",limitsize = FALSE)

### Plot Clades individually
rectree <- ggtree(y,size=0.3) + geom_treescale(x = 0.0332, y = 135, width=0.001, offset = 0.5, linesize =0.3, fontsize=3) + geom_nodelab(size=2,hjust=1)
rectree$data[rectree$data$node %in% c(EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode), "x"] <- mean(rectree$data$x)
rectree$data[rectree$data$node == 580,"x"] <- 0.018
rectree$data[rectree$data$node == 921,"x"] <- 0.018
rectree$data[rectree$data$node == 1157,"x"] <- 0.018
rectree$data[rectree$data$node == 580,"branch.length"] <- rectree$data[rectree$data$node == 580,"branch.length"] + 0.018
rectree$data[rectree$data$node == 921,"branch.length"] <- rectree$data[rectree$data$node == 921,"branch.length"] + 0.018
rectree$data[rectree$data$node == 1157,"branch.length"] <- rectree$data[rectree$data$node == 1157,"branch.length"] + 0.018
rectree$data[rectree$data$node == 580,"branch"] <- rectree$data[rectree$data$node == 580,"branch"] + 0.018
rectree$data[rectree$data$node == 921,"branch"] <- rectree$data[rectree$data$node == 921,"branch"] + 0.018
rectree$data[rectree$data$node == 1157,"branch"] <- rectree$data[rectree$data$node == 1157,"branch"] + 0.018
rectree <- rectree + theme(legend.position = 'none')

sake <- viewClade(rectree, MRCA(subtree, 'RIB6007','CBS2270')) 
sake <- sake + geom_tiplab(x = 0.0385,size=3,aes(color=clade), align=TRUE, linesize = 0) + 
  scale_colour_manual(values = cladecolors)
sake <- sake +  geom_tippoint(x=0.03605, aes(shape=amplification),size=2, color = '#D2691E') 
allcolors <- c(cladecolors,ecocolors)
sake <- sake + geom_point2(aes(subset = isTip, x=0.0405, color=eco),size=2,shape=15) + scale_color_manual(values = allcolors)
sake <- sake + geom_tiplab(x=0.0362,aes(label=chr),geom='text',size=3,color='#D2691E') 
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/sake2.png", width=15, height = 20, units="cm",limitsize = FALSE)

winepre <- rectree + xlim_tree(0.0461)

wine <- viewClade(winepre, MRCA(subtree, 'ULG84F88I90','Lib73'))
wine +theme_classic()

wine <- wine  + geom_tiplab(x = 0.0451,size=1,aes(color=clade), align=TRUE, linesize = 0) + 
  scale_colour_manual(values = cladecolors)
wine <- wine +  geom_tippoint(aes(shape=amplification),size=1, color = '#D2691E')
wine <- wine + geom_point2(aes(subset = isTip, x=0.0459, color=eco),size=1,shape=15) + scale_color_manual(values = allcolors)
wine <- wine + geom_tiplab(aes(label=chr),geom='text',size=1,color='#D2691E') 
wine <- wine +  geom_tippoint(x=0.045,aes(shape=amplification),size=1, color = '#D2691E')
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/wine.png", width=20, height = 30, units="cm",limitsize = FALSE)

wine + theme_classic()
+ geom_tiplab(size=1,aes(color=clade), offset = 0.0009) +
  geom_tippoint(aes(shape=amplification),size=0.5, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0014, color=eco),size=0.5,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=1,color="blue") +
  scale_color_manual(values = c("#890084","#FF8101","#FFCD01","#705B05","#D000AD","#22FF00","#169B02","#69BD5C","#98E28D","#8DE2CB","#456C02",
                                "#9FFF9F","#05A68B","#00B8E6","#356433","#82E4DD","#AE6F6D","#C71D17","#0F26F3","#597951","#95540A","#BAB405",
                                "#6C6B56","#DADA06","#FFCD01","#0F26F3","#B606DA","#BAB405","#0A020B","#8A7BD8","#38DB06","#1E7F01","#FC9004",
                                "#E8BF89","#B2CAF9","#0991AC","#004DFF","#015E2F","#B603EC","#F58974","gray","#FF0000","#615F28","#093A08",
                                "black","#B3F6FA","#7F049E"))
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/wine.png", width=25, height = 20, units="cm",limitsize = FALSE)

bioeth <- viewClade(subtree, MRCA(subtree, 'SA_9_3_VR1','RP_10_14'))
bioeth <- bioeth + geom_tiplab(size=5,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=2, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0051, color=eco),size=2,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=5,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/bioeth.png", width=15, height = 20, units="cm",limitsize = FALSE)

dairy <- viewClade(subtree, MRCA(subtree, 'CLIB649','CBS420'))
dairy <- dairy + geom_tiplab(size=5,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=2, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0051, color=eco),size=2,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=5,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/dairy.png", width=15, height = 20, units="cm",limitsize = FALSE)

afbeer <- viewClade(subtree, MRCA(subtree, 'N34_2_4_a','CBS4454'))
afbeer <- afbeer + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/afbeer.png", width=15, height = 20, units="cm",limitsize = FALSE)

moak <- viewClade(subtree, MRCA(subtree, 'STG_2','N19_1C'))
moak <- moak + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/moak.png", width=15, height = 20, units="cm",limitsize = FALSE)

ale <- viewClade(subtree, MRCA(subtree, 'CBS2165a','CBS6308'))
ale <- ale + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/ale.png", width=15, height = 20, units="cm",limitsize = FALSE)

fgh <- viewClade(subtree, MRCA(subtree, 'HE020','HE008'))
fgh <- fgh + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/fgh.png", width=15, height = 20, units="cm",limitsize = FALSE)

agave <- viewClade(subtree, MRCA(subtree, 'LCBG_Mosca3','LCBG_3D2'))
agave <- agave+ geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/agave", width=15, height = 20, units="cm",limitsize = FALSE)

cocoa <- viewClade(subtree, MRCA(subtree, 'MTF2554','MTF2550'))
cocoa <- cocoa + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/cocoa.png", width=15, height = 20, units="cm",limitsize = FALSE)


aferm <- viewClade(subtree, MRCA(subtree, 'CLQCA_10_620','K12'))
aferm <- aferm + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/aferm.png", width=15, height = 20, units="cm",limitsize = FALSE)

aisl <- viewClade(subtree, MRCA(subtree, 'CBS1593','CBS1576'))
aisl <- aisl + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/aisl.png", width=15, height = 20, units="cm",limitsize = FALSE)

naoak <- viewClade(subtree, MRCA(subtree, 'YPS128','YPS142'))
naoak <- naoak + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/naoak.png", width=15, height = 20, units="cm",limitsize = FALSE)

russ <- viewClade(subtree, MRCA(subtree, 'N22_00_5D','N163_01_5A'))
russ <- russ + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/russ.png", width=15, height = 20, units="cm",limitsize = FALSE)

afwine <- viewClade(subtree, MRCA(subtree, 'YJM1248','DJ74'))
afwine <- afwine + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/afwine.png", width=15, height = 20, units="cm",limitsize = FALSE)

equat <- viewClade(subtree, MRCA(subtree, 'YPS615','CLQCA_20_259'))
equat <- equat + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/equat.png", width=15, height = 20, units="cm",limitsize = FALSE)

chnV <- viewClade(subtree, MRCA(subtree, 'HN15','HN16'))
chnV <- chnV + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/chnV.png", width=15, height = 20, units="cm",limitsize = FALSE)


malaysian <- viewClade(subtree, MRCA(subtree, 'YJM1447','UWOPS03_459_1'))
malaysian <- malaysian + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/malaysian.png", width=15, height = 20, units="cm",limitsize = FALSE)


chnIII <- viewClade(subtree, MRCA(subtree, 'HN19','HN10'))
chnIII <- chnIII + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/chnIII.png", width=15, height = 20, units="cm",limitsize = FALSE)

feasia <- viewClade(subtree, MRCA(subtree, 'Ksc73_2D','Ksc40_8A'))
feasia <- feasia + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/feasia.png", width=15, height = 20, units="cm",limitsize = FALSE)

chnII <- viewClade(subtree, MRCA(subtree, 'SX1','SX3'))
chnII <- chnII + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/chnII.png", width=15, height = 20, units="cm",limitsize = FALSE)

chnI <- viewClade(subtree, MRCA(subtree, 'HN6'))
chnI <- chnI + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/chnI.png", width=15, height = 20, units="cm",limitsize = FALSE)

taiwan <- viewClade(subtree, MRCA(subtree, 'GE14S01_7B','EN14S01'))
taiwan <- taiwan + geom_tiplab(size=0.6,aes(color=clade), offset = 0.002) +
  geom_tippoint(aes(shape=amplification),size=0.4, color = 'blue') +
  geom_point2(aes(subset = isTip, x=x+0.0029, color=eco),size=0.4,shape=15) +
  geom_tiplab(aes(label=chr),geom='text',size=0.6,color="blue")
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/taiwan.png", width=15, height = 20, units="cm",limitsize = FALSE)


#########################################################################################################################################################

### Stuff that doesn't quite work yet
tdf <- as.data.frame(ptree$data)
for(i in levels(tdf$clade)){
  ctree <- viewClade(subtree, MRCA(subtree, min(tdf[tdf$clade == i & tdf$isTip == TRUE,12]),
                     max(tdf[tdf$clade == i & tdf$isTip == TRUE,12])))
  ctree <- ctree + geom_tiplab(size=5,aes(color=clade), offset = 0.0001) +
    geom_tippoint(aes(shape=amplification),size=2, color = '#FF0000')+
    geom_point2(aes(subset = isTip, x=x+0.0002, color=eco),size=2,shape=17)
  ggsave(paste0("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/",substr(i,1,2),".png",sep=""), width=15, height = 20, units="cm",limitsize = FALSE)
}

as.data.frame(ptree$data[ptree$data$clade == "23. North American oak" & ptree$data$isTip == TRUE,"y"])
ptree <- ptree + geom_strip('ULG84F88I90','Lib51', color = 'purple', label = "1. Wine/European")

for(i in 1:nrow(cl_node)){
  ptree <- ptree + geom_cladelabel(node = cl_node[i,2], label = cl_node[i,1])
}




ptree + geom_point(aes(x=0.0002, shape= ifelse((isTip == TRUE), "Yes",NA), color = `Ecological origins`), size = 1)



ptree + geom_text2(aes(subset=!isTip, label=label),hjust=-0.0001, size=1)


ptree

ptree + geom_text(aes(subset = !isTip, x=branch, label=label))

ptree + theme(legend.position = c(.4,.6), legend.key.size = unit(1,"line")) 
ptree
  
  geom_treescale(x=0.003, y=30, width=0.005,fontsize = 10,linesize = 0.5,offset=20) + 
  geom_tiplab(size=2,hjust = -0.2) +
  geom_rootedge(0)+
  geom_tippoint(shape = ifelse((aneuploidy=="Yes"),8,NA), size = 2, color='black')+
  theme(legend.position = c(.4,.6), legend.text = 'none', legend.key.size = unit(6,"line")) 
  guides(fill = guide_legend(title = "Legend Title",override.aes = list(size=4)))



#

treedata <- as.data.frame(ptree$data)
clade_names <- c()
clade_numbers <- c()
# Get clade numbers
for(i in levels(treedata$clade)){
  clade_names <- append(clade_names,i)
  clade_names <- clade_names[!is.na(clade_names)]
  temp <- unique(treedata[treedata$clade == i, 2])
  temp <- temp[!is.na(temp)]
  clade_numbers <- append(clade_numbers,MRCA(y,temp))
}
clade_numbers[8] <- 923
names(clade_numbers) <- clade_names
newtree <- groupClade(y,clade_numbers)
cols <- c(we = "#6600CC",
          fgh= "#FFB266",
          ab="#CC6600",
          wac="#FF9933",
          wpw="#FF8000",
          ciii="#336600",
          cii="#006600",
          ci="#006633",
          t="#009999",
          fea="#00CCCC",
          m="#00FFFF",
          a="#000066",
          cv="#00994C",
          e="#4C9900",
          fer="#66CC00",
          nao="#80FF00",
          ai="#99FFFF",
          s="#FF0000",
          af="#FF9999",
          bb="#66B2FF",
          mo="#333300",
          fd="#FF66FF",
          ab="#CCCC00",
          ma="#999900")
p <- ggtree(newtree, aes(color=group)) + 
  geom_tiplab(size=1.5,hjust = -0.2) +
  geom_treescale(x=0,y=1,width=0.002) +
  scale_color_manual(values=c(cols,"black"),na.value="black",name="Lineage",
                     breaks=c("we","fgh","ab","wac","wpw","ciii","cii","ci","t","fea","m","a","cv","e","fer","nao","ai","s","af","bb","mo",
                              "fd","ab","ma")) +
  guides(color = guide_legend(override.aes = list(size=5,shape=15))) + 
  theme_tree2(legend.position = c(.1,.88))
  
p
