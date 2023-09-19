library(ape)
library(ggtree)
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
ecocolors <- c(bee, bio, bak, humc, hum, ind, dai, ferm, flow, fruit, ind, ins, lab, palm, sak, ferm, tr, unk, wat, wine)

distMat <- read.table("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.8/gain2n_gSNPs_QC4_distance.txt")
rownames(distMat) <- distMat[,1]
distMat <- distMat[,2:904]
colnames(distMat) <- rownames(distMat)
dist.gene(distMat)
distMat2 <- as.dist(distMat)
tree <- nj(distMat2)

ggtree(tree, layout = "circular")

options(
  repos = c(
    zkamvar = "https://zkamvar.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)
install.packages("adegenet")
library(adegenet)
alignment <- read.PLINK("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.8/gain2n_gSNPs_QC4.raw")
tre <- nj(dist(as.matrix(alignment)))

library(adegenet)
library(ape)
library(parallel)

df <- read.PLINK("gain2n_gSNPs_QC4.raw", parallel = TRUE, n.cores = 6)

pca1 <- glPca(df, nf = 2, parallel = TRUE, n.cores = 6)

tre <- nj(dist(as.matrix(df)))

tree <- read.tree("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/RunQC40Miss10MAF1.0Recomb0.8.NJ.tree")
newdf <- data.frame(tree$tip.label, row.names = TRUE)
for(i in 1:length(tree$tip.label)){
  for(j in 1:length(dipGainOnly$basename)){
    if(tree$tip.label[i] == dipGainOnly$basename[j]){
      rownames(newdf)[i] <- dipGainOnly$basename[j]
      newdf$amplification[i] <- dipGainOnly$aneuploidy_binary[j]
      newdf$ecology[i] <- dipGainOnly$ecology_category[j]
      newdf$clade[i] <- dipGainOnly$
    }
  }
}
newdf[is.na(newdf$ecology),] <- "unknown"
amplification <- as.matrix(newdf)[,1]
eco <- as.matrix(newdf)[,2]
ampdf <- data.frame(node = nodeid(tree, names(amplification)),
                    amplification = amplification)
traitsdf <- cbind(ampdf, eco)
# Make node numeric
traitsdf$node <- as.numeric(traitsdf$node)
# Add traits to tree
y <- full_join(tree, traitsdf, by = 'node')

ggtree(y, size = 0.2, layout = "fan") + 
  geom_tree(aes(color = eco))+
  scale_color_manual(values=ecocolors)
