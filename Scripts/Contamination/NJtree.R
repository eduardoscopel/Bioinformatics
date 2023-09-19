library(ape)
library(ggtree)
library(ggplot2)

alcat <- read.dna("fasta/chr4.cont_hh0.fa", format = "fasta")
f <- function(x) nj(dist.dna(x,model="F84"))
alcatree <- f(alcat)
set.seed(657)
alcatbs <- boot.phylo(alcatree, alcat, f, B=10,quiet=TRUE)
write.tree(alcatree,"test.tree")
alcatree$node.label <- alcatbs

ggtree(alcatree) + geom_tiplab() + geom_nodelab(size=2) + theme_classic()

ggsave("test.pdf", height = 8, width = 14, units = "in", limitsize = FALSE)

library(data.table)
library(dplyr)


alcatree <- root(alcatree, c("EN14S01","GE14S01_7B"))
alcatree$tip.label[c(52,53)] <- c("CLIB219.26-cont","CBS1479-hh0")
newdf <- data.frame(alcatree$tip.label, row.names = TRUE)

strains <- fread("Scer/51Strains.txt")
strains <- as.data.frame(strains)
strains <- strains[order(strains$ID),]
strains[52,] = c("CLIB219.26","cont")
strains[53,] = c("CBS1479","hh0")
newdf <- cbind(newdf, strains)
rownames(newdf)[c(52,53)] <- c("CLIB219.26-cont","CBS1479-hh0")
newdf$node <- nodeid(alcatree, alcatree$tip.label)
newtree <- full_join(alcatree, newdf, by = "node")

we1 <- "#800080"
ol2 <- "#808000"
bb3 <- "#14B5AF"
mo4 <- "#228B22"
fd5 <- "#00BFFF"
ab6 <- "#bf812d"
mb7 <- "black"
mo8 <- "#FF1493"
ma9 <- "#8B4513" 
fgh10 <- "#FF8C00"
ab11 <- "#BDB76B"
wac12 <- "#D2B48C"
apw13 <-  "#9400D3"
ch14 <- "#98FB98"
ch15 <- "#3CB371"
ch16 <- "#006400"
t17 <- "#008080"
fea18 <- "#00FFFF"
m19 <- "#4B0082"
ch20 <- "#2F4F4F"
e21 <- "#0000FF"
fer22 <- "#00008B"
nao23 <- "#FF8C00"
ai24 <- "#E9967A"
s25 <- "#FF0000"
af26 <-  "#8B0000"

cladecolors <- c(we1, fgh10,ab11,wac12,apw13, ch14, ch15, ch16, t17, fea18, m19, ol2, ch20, e21, fer22, nao23, ai24,s25,af26,bb3, mo4, fd5, ab6, mb7, mo8, ma9, "darkgray", "darkgray")


ggtree(newtree) + 
  geom_tiplab(aes(color = Clade), size = 3, show.legend = FALSE) + 
  scale_color_manual(values = cladecolors) +
  geom_nodelab(size=2, hjust = 0) + 
  geom_treescale(x= 0, y=10) 

p1 <- ggtree(newtree) + 
  geom_tiplab(aes(color = Clade), size = 3, show.legend = FALSE) + 
  scale_color_manual(values = cladecolors) +
  geom_nodelab(size=2, hjust = 0) + 
  geom_treescale(x= 0.008, y=1) + 
  geom_hilight(node = 52, fill = 'steelblue', alpha=.3, extend = 0.001) +
  geom_hilight(node = 53, fill = 'forestgreen', alpha=.3, extend = 0.001) +
  scale_x_break(breaks = c(0.001,0.008), scales = 'fixed', ticklabels = NULL) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())
p1

ggsave("test.pdf", height = 8, width = 14, units = "in", limitsize = FALSE)

p2 <- rotate(p1, 98)
p2

