library(ape)
library(phytools)
library(data.table)

tree <- read.tree("IQ-tree-WG-MD-1k.tree")
par(xpd=T)							# allow text in the margins

strain_table <- fread("sc_table.txt",header=TRUE)
strain_table <- strain_table[strain_table$population_new == "BEER" | 
                               strain_table$population_new == "GRAPE.WINE" |
                               strain_table$population_new == "EU.OAK" |
                               strain_table$population_new == "NC.OAK" |
                               strain_table$population_new == "PA.OAK" |
                               strain_table$population_new == "TAIWAN",]
strain_table <- strain_table[strain_table$ID != "DAV1e",]

plot.phylo(tree, cex = 0.4, label.offset=0.001, x.lim=0.1, no.margin = TRUE)
rtree <- root(tree, outgroup = c("GE14S01-7B", "EN14S01", "EM14S01-3B"), resolve.root = TRUE)
plot.phylo(rtree, cex = 0.4, label.offset=0.001, x.lim=0.1)

getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

# popcolors is in the same order as pop: GRAPE.WINE, EU.OAK, NC.OAK, PA.OAK, BEER, TAIWAN
pop <- unique(strain_table$population_new)
popcolors <- c("#a6dba0","#762a83","#5aae61","#005a32","orange","#87CEFA")

# ecocolors is in the same order as eco: 
eco <- unique(strain_table$ecology_category)
eco[6] <- "unknown"
ecocolors <- c("#009900", "#C316DE", "#01FF2C","#FFF801","#28CB9A", "#000000","#6E5104",
               "#D7AB33","#FF8B04","#FF5004","#0471FF","#EFA3C2")

pew <- 2
ncoak <- 288
paoak <- 303
taiwan <- 317
euoak <- 260
wine <- 165
boulardii <- 235

ec<-rep("black",length(rtree$edge[,2]))	# default edge color = black
ew<-rep(1,length(rtree$edge[,2]))		# default edge width = 1
tc<-rep("black",length(rtree$tip))		# default tip color = black (use for lineage)
tc2<-rep("black",length(rtree$tip))		# default tip color = black (use for environment)
tc3<-rep("black",length(rtree$tip))		# default tip color = black (use for aneuploidy_binary)


### adding color to tc vector based on population on strain table
i=1 
for(x in strain_table$ID){
  tlab <- grep(x,rtree$tip.label)
  pc <- grep(strain_table$population_new[i],pop)
  tc[tlab] <- popcolors[pc]
  i=i+1
}

d <- getDescendants(rtree, node=ncoak)
for (x in d) {ec[rtree$edge[,2]==x]<-popcolors[3]; ew[rtree$edge[,2]==x]<-pew; tc2[x]<-popcolors[3] }
ec[rtree$edge[,2]==ncoak]<-popcolors[3]		# color the edge leading to this clade
ew[rtree$edge[,2]==ncoak]<-pew

d <- getDescendants(rtree, node=paoak)
for (x in d) {ec[rtree$edge[,2]==x]<-popcolors[4]; ew[rtree$edge[,2]==x]<-pew; tc2[x]<-popcolors[4] }
ec[rtree$edge[,2]==paoak]<-popcolors[3]		# color the edge leading to this clade
ew[rtree$edge[,2]==paoak]<-pew

d <- getDescendants(rtree, node=taiwan)
for (x in d) {ec[rtree$edge[,2]==x]<-popcolors[6]; ew[rtree$edge[,2]==x]<-pew; tc2[x]<-popcolors[6] }
ec[rtree$edge[,2]==taiwan]<-popcolors[6]		# color the edge leading to this clade
ew[rtree$edge[,2]==taiwan]<-pew

d <- getDescendants(rtree, node=euoak)
for (x in d) {ec[rtree$edge[,2]==x]<-popcolors[1]; ew[rtree$edge[,2]==x]<-pew; tc2[x]<-popcolors[1] }
ec[rtree$edge[,2]==euoak]<-popcolors[1]		# color the edge leading to this clade
ew[rtree$edge[,2]==euoak]<-pew

d <- getDescendants(rtree, node=wine)
for (x in d) {ec[rtree$edge[,2]==x]<-popcolors[2]; ew[rtree$edge[,2]==x]<-pew; tc2[x]<-popcolors[2] }
ec[rtree$edge[,2]==wine]<-popcolors[2]		# color the edge leading to this clade
ew[rtree$edge[,2]==wine]<-pew

png('aneuploidy-non-admixed-tree.png',2200,2200,res=200,pointsize=11)
plot.phylo(rtree,cex=0.5,edge.color=ec,tip.color=tc,edge.width=ew,font=2, label.offset=0.0002, x.lim=0.026, no.margin = TRUE)
boostratps <- as.numeric(rtree$node.label)
boostratps[boostratps<70] <- NA
nodelabels(boostratps, frame="none",cex=0.6,bg="white",adj = c(1.5,-0.4))
add.scale.bar(0,-3,cex=0.8)

j=1 
for(x in strain_table$ID){
  tlab <- grep(x,rtree$tip.label)
  ecc <- grep(strain_table$ecology_category[j],eco)
  tc2[tlab] <- ecocolors[ecc]
  j=j+1
}

tiplabels("",grep(ecocolors[1],tc2),pch=20,frame="none",col=ecocolors[1],adj=0.5001)
tiplabels("",grep(ecocolors[2],tc2),pch=20,frame="none",col=ecocolors[2],adj=0.5001)
tiplabels("",grep(ecocolors[3],tc2),pch=20,frame="none",col=ecocolors[3],adj=0.5001)
tiplabels("",grep(ecocolors[4],tc2),pch=20,frame="none",col=ecocolors[4],adj=0.5001)
tiplabels("",grep(ecocolors[5],tc2),pch=20,frame="none",col=ecocolors[5],adj=0.5001)
tiplabels("",grep(ecocolors[6],tc2),pch=20,frame="none",col=ecocolors[6],adj=0.5001)
tiplabels("",grep(ecocolors[7],tc2),pch=20,frame="none",col=ecocolors[7],adj=0.5001)
tiplabels("",grep(ecocolors[8],tc2),pch=20,frame="none",col=ecocolors[8],adj=0.5001)
tiplabels("",grep(ecocolors[9],tc2),pch=20,frame="none",col=ecocolors[9],adj=0.5001)
tiplabels("",grep(ecocolors[10],tc2),pch=20,frame="none",col=ecocolors[10],adj=0.5001)
tiplabels("",grep(ecocolors[11],tc2),pch=20,frame="none",col=ecocolors[11],adj=0.5001)
tiplabels("",grep(ecocolors[12],tc2),pch=20,frame="none",col=ecocolors[12],adj=0.5001)

anp <- unique(strain_table$aneuploidy_binary)
anpcolors <- c("Red", "Black","Blue")
k=1 
for(x in strain_table$ID){
  tlab <- grep(x,rtree$tip.label)
  antc <- grep(strain_table$aneuploidy_binary[k],anp)
  tc3[tlab] <- anpcolors[antc]
  k=k+1
}

#tiplabels("",grep(anpcolors[1],tc3),pch=4,frame="none",col=anpcolors[1],adj=0.5095)
tiplabels("",grep(anpcolors[2],tc3),pch=8,cex =0.8,frame="none",col=anpcolors[2],adj=0.5012)


lineagelty<-rep(1,6)
lineagelwd<-rep(2,6)
op <- par(cex=2)
legend(0.002,130,pop,col=popcolors,lty=lineagelty,lwd=lineagelwd,bty="n",cex=0.8,title="Populations",title.adj=0,text.col=popcolors,title.col="black")

#legend("left",eco,col=ecocolors,lty=ecolty,lwd=ecolwd,bty="n",cex=0.8,title="Environment",title.adj=0,text.col=ecocolors,title.col="black")
legend(0.002,90,eco,col=ecocolors,pch=20,bty="n",cex=0.8,title="Environment",title.adj=0,text.col=ecocolors,title.col="black")

legend(0.002,20,"Aneuploids",col="Black",pch=8,bty="n",cex=0.8,title.adj=0,text.col="Black",title.col="Black")

dev.off()