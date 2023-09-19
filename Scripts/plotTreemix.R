library(RColorBrewer)
library(R.utils)
library(data.table)
library(ggplot2)
library(cowplot)
library(ape)
library(treeio)
library(ggtree)
source("/Users/es47540/Downloads/treemix-1.13/src/plotting_funcs.R")

setwd("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/corrected_pop_names/")
edges5="new_621-merged-snps.noN.5"
edges6="new_621-merged-snps.noN.6"
edges7="new_621-merged-snps.noN.7"
source("/Users/es47540/Documents/GitHub/eduardo/Scripts/treemix.R")

cladeColors5 <- c(fgh10, nao23,ol2, m19, ai24, fea18, ch20, mo4, we1, bb3, ch14, fer22, af26,wac12, ch16, e21, apw13, t17, ch15, ma9)
aneupCladeColors5 <- c(mb7, ab11, ab6, fd5, s25, mo8)

e5 <- plot.treemix(edges5,plot.nodes = FALSE, branch.colour = "black", clade_color = cladeColors5, aneup_clade_color = aneupCladeColors5)

cladeColors6 <- c(wac12, ma9, fgh10, fea18, t17, mo4, e21, fer22, ol2, ch15, ch16, we1, apw13, ch14, m19, ai24, bb3, nao23, ch20, af26)
aneupCladeColors6 <- c(fd5, s25,mb7,ab6, mo8, ab11)
e6 <- plot.treemix(edges6,plot.nodes = FALSE, branch.colour = "black", clade_color = cladeColors6, aneup_clade_color = aneupCladeColors6)

cladeColors7 <- c(ch20, mo4, af26, nao23, fer22, wac12, m19, t17, we1, ol2, ch15, ai24, e21, fgh10, ma9, bb3, ch16, apw13, fea18, ch14)
aneupCladeColors7 <- c(mb7, fd5, ab6, s25, ab11, mo8)
e7 <- plot.treemix(edges7,plot.nodes = FALSE, branch.colour = "black", clade_color = cladeColors7, aneup_clade_color = aneupCladeColors7)

lliketab <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/k500/s1/k500_llike.txt")
colnames(lliketab) <- c("M", "llike")
lliketab <- lliketab[lliketab$M < 11,]
llikeplot <- ggplot(lliketab) + 
  geom_point(aes(x=M, y= llike)) +
  geom_line(aes(x=M, y= llike)) +
  scale_x_continuous(breaks = c(0:11),expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(-10000,2000, by = 2000),limits = c(-10000,2000),expand = c(0.01,0.01)) +
  xlab("Migration events")+
  ylab("Log likelihood") + 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.border = element_blank(), axis.line.x = element_line(size=0.5),axis.line.y = element_line(size=0.5))
options(scipen=999)
llikeplot

suppFigS4_3 <- plot_grid(llikeplot, e5, e6, e7, nrow=2,ncol=2)
suppFigS4_3
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v3.png", width=12, height = 12, units="in",limitsize = FALSE, dpi = 600)
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v3.pdf", width=10, height = 10, units="in",limitsize = FALSE, dpi = 600)


prefix="new_621-merged-snps.noN"
par(mfrow=c(1,1))
for(edge in 0:3){
  plot_tree(cex=0.8, paste0(prefix,".",edge))
  title(paste(edge,"edges"))
}
plot_tree(cex=0.8, paste0(prefix,".",6))
write.table(unique(peter_clades$Clades),"/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/pop.list", row.names = FALSE)
for(edge in 5:7){
  plot_resid(stem=paste0(prefix,".",edge),pop_order = "pop.list")
}



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

cladeColors <- c(ab6,apw13,ab11,ol2,af26,ai24,bb3,ch16,ch15,ch14,ch20,e21,we1,fea18,fer22,fd5,fgh10,m19,mo4,ma9,mo8,mb7,nao23,s25,t17,wac12)


### Five migration events 

fiveEtree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/5e.tree")

newdf5 <- data.frame(fiveEtree$tip.label, row.names = TRUE)

newdf5$highAneuploidy5 <- ifelse(rownames(newdf5) == "Sake" |
                                 rownames(newdf5) == "Ale_beer" |
                                 rownames(newdf5) == "Mixed_origin" |
                                 rownames(newdf5) == "Mosaic_beer" |
                                 rownames(newdf5) == "French_dairy" |
                                 rownames(newdf5) == "African_beer", TRUE, FALSE)

highAneuploidy5 <- as.matrix(newdf5)[,1]
fiveTreedf <- data.frame(node = nodeid(fiveEtree, names(highAneuploidy5)),
                        clades = fiveEtree$tip.label,
                        highAneuploidy5 = newdf5$highAneuploidy5)

fiveTreedf$highAneuploidy5 <- ifelse(rownames(fiveTreedf) == "Sake" |
                                     rownames(fiveTreedf) == "Ale_beer" |
                                     rownames(fiveTreedf) == "Mixed_origin" |
                                     rownames(fiveTreedf) == "Mosaic_beer" |
                                     rownames(fiveTreedf) == "French_dairy" |
                                     rownames(fiveTreedf) == "African_beer", TRUE, FALSE)

fiveEtree <- full_join(fiveEtree, fiveTreedf, by = 'node')
fiveEdges <- ggtree(fiveEtree, layout = 'rectangular') + 
  geom_tiplab(size=2.8,align = TRUE, aes(color = clades), show.legend=FALSE, linesize = 0.25) + 
  scale_color_manual(values = c(ab6,apw13,ab11,ol2,af26,ai24,bb3,ch16,ch15,ch14,ch20,e21,we1,fea18,fer22,fd5,fgh10,m19,mo4,ma9,mo8,mb7,nao23,s25,t17,wac12))+
  geom_taxalink(taxa1 = 51, taxa2 = 'Ale_beer',color='red',size = 0.32/0.44,alpha = 0.32/0.44, arrow = arrow(length = unit(0.02,'npc'))) + # Sake/AF -> Mo
  geom_taxalink(taxa1 = 51, taxa2 = 'Mixed_origin', color = 'red', size = 0.32/0.44,alpha = 0.32/0.44, arrow = arrow(length = unit(0.02,'npc'))) + # Sake/AF -> AleB
  geom_taxalink(taxa1 = 51, taxa2 = 'Mosaic_beer', color = 'red', size = 0.13/0.44,alpha = 0.13/0.44, arrow = arrow(length = unit(0.02,'npc'))) + 
  geom_taxalink(taxa1 = 41,taxa2 = 'West_African_cocoa', color = 'red', size = 0.44/0.44,alpha = 0.44/0.44, arrow = arrow(length = unit(0.02,'npc'))) + # EUW/A -> WAc
  geom_taxalink(taxa1 = 'Alpechin',taxa2 = 'Mexican_agave', color = 'red', size = 0.4/0.44, alpha = 0.4/0.44, arrow = arrow(length = unit(0.02,'npc'))) + 
  geom_hilight(node = 23, fill = "black", alpha = 0.1, extend = 0.05)+ # sake
  geom_hilight(node = 11, fill = "black", alpha = 0.1, extend = 0.051)+ # mosaic beer
  geom_hilight(node = 10, fill = "black", alpha = 0.1, extend = 0.0395)+ # Ale beer
  geom_hilight(node = 7, fill = "black", alpha = 0.1, extend = 0.0485)+ # African beer
  geom_hilight(node = 6, fill = "black", alpha = 0.1, extend = 0.0488)+ # French dairy
  geom_hilight(node = 9, fill = "black", alpha = 0.1, extend = 0.0438)+ # Mixed origin
  scale_x_continuous(limits = c(0, 0.16))+
  theme_tree2()+
  labs(caption = "Genetic drift")

fiveEdges <- flip(fiveEdges, MRCA(fiveEdges, "Alpechin", "French_Guiana"), MRCA(fiveEdges, "Asian_fermentation", "CHN_III"))
fiveEdges

### Six migration events 


sixEtree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/6e.tree")

newdf <- data.frame(sixEtree$tip.label, row.names = TRUE)
newdf$highAneuploidy <- ifelse(rownames(newdf) == "Sake" |
                                     rownames(newdf) == "Ale_beer" |
                                     rownames(newdf) == "Mixed_origin" |
                                     rownames(newdf) == "Mosaic_beer" |
                                     rownames(newdf) == "French_dairy" |
                                     rownames(newdf) == "African_beer", TRUE, FALSE)

highAneuploidy <- as.matrix(newdf)[,1]
sixTreedf <- data.frame(node = nodeid(sixEtree, names(highAneuploidy)),
                        clades = sixEtree$tip.label,
                        highAneuploidy = newdf$highAneuploidy)

sixTreedf$highAneuploidy <- ifelse(rownames(sixTreedf) == "Sake" |
                                     rownames(sixTreedf) == "Ale_beer" |
                                     rownames(sixTreedf) == "Mixed_origin" |
                                     rownames(sixTreedf) == "Mosaic_beer" |
                                     rownames(sixTreedf) == "French_dairy" |
                                     rownames(sixTreedf) == "African_beer", TRUE, FALSE)

sixEtree <- full_join(sixEtree, sixTreedf, by = 'node')
sixEdges <- ggtree(sixEtree, layout = 'slanted') + 
  geom_tiplab(size=2.8,align = FALSE, aes(color = clades), show.legend=FALSE, linesize = 0.25) + 
  scale_color_manual(values = c(ab6,apw13,ab11,ol2,af26,ai24,bb3,ch16,ch15,ch14,ch20,e21,we1,fea18,fer22,fd5,fgh10,m19,mo4,ma9,mo8,mb7,nao23,s25,t17,wac12))+
  geom_taxalink(taxa1 = 50, taxa2 = 'Mixed_origin',color='red',size = 0.35/0.43,alpha = 0.35/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # Sake/af -> mo
  geom_taxalink(taxa1 = 36,taxa2 = 'West_African_cocoa', color = 'red', size = 0.43/0.43,alpha = 0.43/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # WE/al/bb -> wac
  geom_taxalink(taxa1 = 'Alpechin',taxa2 = 'Mexican_agave', color = 'red', size = 0.39/0.43,alpha = 0.39/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # al -> ma
  geom_taxalink(taxa1 = 50,taxa2 = 'Mosaic_beer', color = 'red', size = 0.16/0.43,alpha = 0.16/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # sake/af -> mb
  geom_taxalink(taxa1 = 50,taxa2 = 'Brazilian_bioethanol', color = 'red', size = 0.22/0.43, alpha = 0.22/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # sake/af -> bb
  geom_taxalink(taxa1 = 50,taxa2 = 'Ale_beer', color = 'red', size = 0.35/0.43, alpha = 0.35/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # sake/af -> bb
  geom_hilight(node = 21, fill = "black", alpha = 0.1, extend = 0.046)+ #sake
  geom_hilight(node = 9, fill = "black", alpha = 0.1, extend = 0.0335)+ # ale beer
  geom_hilight(node = 8, fill = "black", alpha = 0.1, extend = 0.0382)+ # mixed origin
  geom_hilight(node = 10, fill = "black", alpha = 0.1, extend = 0.046)+ # mosaic beer
  geom_hilight(node = 3, fill = "black", alpha = 0.1, extend = 0.046)+ # African beer
  geom_hilight(node = 2, fill = "black", alpha = 0.1, extend = 0.0465)+ # French dairy
  geom_curve(aes(x=0.08,y=11,xend=0.09, yend = 13),arrow = arrow(length = unit(0.5, "cm")))+
  scale_x_continuous(limits = c(0, 0.16))+
  theme_tree2()+
  labs(caption = "Genetic drift")

sixEdges <- flip(sixEdges, MRCA(sixEdges, "Ale_beer", "French_Guiana"), MRCA(sixEdges, "Asian_fermentation", "CHN_III"))
sixEdges <- flip(sixEdges, MRCA(sixEdges, "Ale_beer", "Mosaic_beer"), MRCA(sixEdges, "European_wine", "Brazilian_bioethanol"))
sixEdges


### seven migration events 


sevenEtree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/7e.tree")

newdf <- data.frame(sevenEtree$tip.label, row.names = TRUE)
newdf$highAneuploidy7 <- ifelse(rownames(newdf) == "Sake" |
                                 rownames(newdf) == "Ale_beer" |
                                 rownames(newdf) == "Mixed_origin" |
                                 rownames(newdf) == "Mosaic_beer" |
                                 rownames(newdf) == "French_dairy" |
                                 rownames(newdf) == "African_beer", TRUE, FALSE)

highAneuploidy7 <- as.matrix(newdf)[,1]
sevenTreedf <- data.frame(node = nodeid(sevenEtree, names(highAneuploidy7)),
                        clades = sevenEtree$tip.label,
                        highAneuploidy7 = newdf$highAneuploidy7)

sevenTreedf$highAneuploidy <- ifelse(rownames(sevenTreedf) == "Sake" |
                                     rownames(sevenTreedf) == "Ale_beer" |
                                     rownames(sevenTreedf) == "Mixed_origin" |
                                     rownames(sevenTreedf) == "Mosaic_beer" |
                                     rownames(sevenTreedf) == "French_dairy" |
                                     rownames(sevenTreedf) == "African_beer", TRUE, FALSE)

sevenEtree <- full_join(sevenEtree, sevenTreedf, by = 'node')
sevenEdges <- ggtree(sevenEtree, layout = 'rectangular') + 
  geom_tiplab(size=2.8,align = TRUE, aes(color = clades), show.legend=FALSE, linesize = 0.25) + 
  scale_color_manual(values = c(ab6,apw13,ab11,ol2,af26,ai24,bb3,ch16,ch15,ch14,ch20,e21,we1,fea18,fer22,fd5,fgh10,m19,mo4,ma9,mo8,mb7,nao23,s25,t17,wac12))+
  geom_taxalink(taxa1 = 'Alpechin',taxa2 = 'Mexican_agave', color = 'red', size = 0.39/0.43,alpha = 0.39/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # al -> ma
  geom_taxalink(taxa1 = 51, taxa2 = 'Mixed_origin',color='red',size = 0.34/0.43,alpha = 0.34/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # Sake/af -> mo
  geom_taxalink(taxa1 = 51,taxa2 = 'Ale_beer', color = 'red', size = 0.35/0.43, alpha = 0.35/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # sake/af -> bb
  geom_taxalink(taxa1 = 51,taxa2 = 'Mosaic_beer', color = 'red', size = 0.16/0.43,alpha = 0.16/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # sake/af -> mb
  geom_taxalink(taxa1 = 51,taxa2 = 'Brazilian_bioethanol', color = 'red', size = 0.22/0.43, alpha = 0.22/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # sake/af -> bb
  geom_taxalink(taxa1 = 'CHN_II',taxa2 = 'Far_East_Asia', color = 'red', size = 0.16/0.43,alpha = 0.16/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # CHNII -> FEA
  geom_taxalink(taxa1 = 39,taxa2 = 'West_African_cocoa', color = 'red', size = 0.43/0.43,alpha = 0.43/0.43, arrow = arrow(length = unit(0.02,'npc'))) + # WE/al/bb -> wac
  geom_hilight(node = 22, fill = "black", alpha = 0.1, extend = 0.046)+ #sake
  geom_hilight(node = 5, fill = "black", alpha = 0.1, extend = 0.0335)+ # ale beer
  geom_hilight(node = 4, fill = "black", alpha = 0.1, extend = 0.0382)+ # mixed origin
  geom_hilight(node = 6, fill = "black", alpha = 0.1, extend = 0.046)+ # mosaic beer
  geom_hilight(node = 2, fill = "black", alpha = 0.1, extend = 0.046)+ # African beer
  geom_hilight(node = 1, fill = "black", alpha = 0.1, extend = 0.0465)+ # French dairy
  scale_x_continuous(limits = c(0, 0.16))+
  theme_tree2()+
  labs(caption = "Genetic drift")

sevenEdges <- flip(sevenEdges, MRCA(sevenEdges, "European_wine", "French_Guiana"), MRCA(sevenEdges, "Asian_fermentation", "CHN_III"))
sevenEdges


suppFigS4_1 <- plot_grid(llikeplot, sixEdges, rel_widths = c(1,2))
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v1.png", width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v1.pdf", width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)


top_row <- plot_grid(llikeplot, fiveEdges, nrow=1, ncol=2, rel_widths = c(1,2))
bottom_row <- plot_grid(sixEdges, sevenEdges, nrow = 1, ncol =2)
suppFigS4_2 <- plot_grid(top_row, sixEdges,sevenEdges, nrow=3,ncol=1)
suppFigS4_2
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v2.png", width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v2.pdf", width=9, height = 12, units="in",limitsize = FALSE, dpi = 600)

suppFigS4_3 <- plot_grid(llikeplot, fiveEdges, sixEdges, sevenEdges, nrow=2,ncol=2)
suppFigS4_3
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v3.png", width=12, height = 12, units="in",limitsize = FALSE, dpi = 600)
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v3.pdf", width=10, height = 10, units="in",limitsize = FALSE, dpi = 600)
