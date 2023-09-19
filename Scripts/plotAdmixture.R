library(data.table)
library(ggplot2)
petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt",header=TRUE)
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





### High aneuploidy clades

highAneuploidy <- peter_clades[peter_clades$Clades == "7. Mosaic beer" | 
                                 peter_clades$Clades == "8. Mixed origin" | 
                                 peter_clades$Clades == "5. French dairy" | 
                                 peter_clades$Clades == "6. African beer" | 
                                 peter_clades$Clades == "11. Ale beer" | 
                                 peter_clades$Clades == "25. Sake",]
highAneuploidy <- highAneuploidy[order(highAneuploidy$Clades),]
#HAdf <- paste("/work/dblab/escopel/Scer/1011peter/vcf/",highAneuploidy$basename, ".vcf.gz", sep = "")
#write.table(HAdf,"/Users/es47540//Documents/GitHub/eduardo/aneuploidy/peter/admixture/highAneuploidy/HA_149.txt", row.names = FALSE)

#treemix.clust <- data.frame(id = peter_clades$basename, id2 = peter_clades$basename, pop = peter_clades$Clades)
#write.table(treemix.clust, "/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/treemix.clust", row.names = FALSE)

logtab12345 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/highAneuploidy/s12345.txt")
logtab12345$seed <- "Run1"
logtab22222 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/highAneuploidy/s22222.txt")
logtab22222$seed <- "Run2"
logtab54321 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/highAneuploidy/s54321.txt")
logtab54321$seed <- "Run3"
logtab <- rbind(logtab12345, logtab22222, logtab54321)
colnames(logtab) <- c("K","CV","Loglikelihood","seed")
logtab <- logtab[logtab$K <=10,]
logplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= Loglikelihood,color = seed)) +
  geom_line(aes(x=K, y= Loglikelihood, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw() +
  theme(legend.position = 'none')
logplot
CVplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= CV,color = seed)) +
  geom_line(aes(x=K, y= CV, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw()+ 
  theme(legend.position = c(0.8,0.8))
CVplot

plot_grid(logplot, CVplot)

popTable <- data.frame(id = highAneuploidy$basename, pop = highAneuploidy$Clades)

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "steelblue"
)

ticklist <- c()
c <- 1
for(i in unique(popTable$pop)){
  ticklist[c] <- round(mean(as.integer(rownames(popTable[popTable$pop == i,]))))
  c= c+1
}


par(mfrow = c(3,1), mar = c(4, 4.1, 7,0))
for(i in 7:9){
  Qmat <- fread(paste("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/highAneuploidy/S22222/HA-149-merged-snps.", i, ".Q", sep=""))
  Qmat <- cbind(Qmat, popTable)
  barplot(t(as.matrix(Qmat[,1:i])), col=c25[1:i], cex.names = 0.2, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2, main = paste("K=",i),cex.main =2)
  axis(1, at=1:149, labels = Qmat$id, las = 2,cex.axis = 0.4)
  axis(3, at=ticklist, labels = unique(Qmat$pop), las = 2,cex.axis = 1)
}


#### Wine High aneuploidy rate clades
wineHighAneuploidy <- peter_clades[peter_clades$Clades == "1. Wine/European" | 
                                     peter_clades$Clades == "7. Mosaic beer" | 
                                     peter_clades$Clades == "8. Mixed origin" | 
                                     peter_clades$Clades == "5. French dairy" | 
                                     peter_clades$Clades == "6. African beer" | 
                                     peter_clades$Clades == "11. Ale beer" | 
                                     peter_clades$Clades == "25. Sake",]
wineHighAneuploidy <- wineHighAneuploidy[order(wineHighAneuploidy$Clades),]
#WHAdf <- paste("/work/dblab/escopel/Scer/1011peter/vcf/",wineHighAneuploidy$basename, ".vcf.gz", sep = "")
#write.table(WHAdf,"/Users/es47540//Documents/GitHub/eduardo/aneuploidy/peter/admixture/wineHighAneuploidy/WHA_429.txt", row.names = FALSE)


logtab12345 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/wineHighAneuploidy/s12345.txt")
logtab12345$seed <- "Run1"
logtab22222 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/wineHighAneuploidy/s22222.txt")
logtab22222$seed <- "Run2"
logtab54321 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/wineHighAneuploidy/s54321.txt")
logtab54321$seed <- "Run3"
logtab <- rbind(logtab12345, logtab22222, logtab54321)
colnames(logtab) <- c("K","CV","Loglikelihood","seed")
logtab <- logtab[logtab$K <=10,]
logplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= Loglikelihood,color = seed)) +
  geom_line(aes(x=K, y= Loglikelihood, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw() +
  theme(legend.position = 'none')
logplot
CVplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= CV,color = seed)) +
  geom_line(aes(x=K, y= CV, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw()+ 
  theme(legend.position = c(0.8,0.8))
CVplot

plot_grid(logplot, CVplot)

popTable <- data.frame(id = wineHighAneuploidy$basename, pop = wineHighAneuploidy$Clades)

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "steelblue"
)

ticklist <- c()
c <- 1
for(i in unique(popTable$pop)){
  ticklist[c] <- round(mean(as.integer(rownames(popTable[popTable$pop == i,]))))
  c= c+1
}


par(mfrow = c(3,1), mar = c(4, 4.1, 8,0))
for(i in 7:9){
  Qmat <- fread(paste("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/wineHighAneuploidy/S12345/WHA_429-merged-snps.", i, ".Q", sep=""))
  Qmat <- cbind(Qmat, popTable)
  barplot(t(as.matrix(Qmat[,1:i])), col=c25[1:i], cex.names = 0.2, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2, main = paste("K=",i),cex.main =2)
  axis(1, at=1:429, labels = Qmat$id, las = 2,cex.axis = 0.4)
  axis(3, at=ticklist, labels = unique(Qmat$pop), las = 2,cex.axis = 1)
}




### No thinning 
idList <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/noThinning/listStrains.txt")
colnames(idList) <- "id"
popTable <- data.frame(id = peter_clades$basename, pop = peter_clades$Clades)
popTable <- merge(idList, popTable, by = "id")

logtab12345 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/s12345.txt")
logtab12345$seed <- "Run1"
logtab22222 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/s22222.txt")
logtab22222$seed <- "Run2"
logtab54321 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/s54321.txt")
logtab54321$seed <- "Run3"
logtab <- rbind(logtab12345, logtab22222, logtab54321)

logplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= Loglikelihood,color = seed)) +
  geom_line(aes(x=K, y= Loglikelihood, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw()
logplot
CVplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= CV,color = seed)) +
  geom_line(aes(x=K, y= CV, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw()
CVplot

test2KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/noThinning/621-merged-snps.2.Q")
merged2KQ <- cbind(test2KQ, popTable)
merged2KQ <- merged2KQ[order(merged2KQ$pop),]
barplot(t(as.matrix(merged2KQ[,1:2])), col=c("steelblue3","orange"), cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged6KQ$id, las =2,cex.axis = 0.2)
axis(3, at=1:621, labels = merged6KQ$pop, las =2,cex.axis = 0.2)


test6KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/noThinning/621-merged-snps.6.Q")
merged6KQ <- cbind(test6KQ, popTable)
merged6KQ <- merged6KQ[order(merged6KQ$pop),]
barplot(t(as.matrix(merged6KQ[,1:6])), col=c("steelblue3","orange","#e41a1c","#5ddb3e","#acdb3e","royalblue2"), cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged6KQ$id, las =2,cex.axis = 0.2)
axis(3, at=1:621, labels = merged6KQ$pop, las =2,cex.axis = 0.2)

test10KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/noThinning/621-merged-snps.10.Q")
merged10KQ <- cbind(test10KQ, popTable)
merged10KQ <- merged10KQ[order(merged10KQ$pop),]
barplot(t(as.matrix(merged10KQ[,1:10])), col=c("steelblue3","orange","#e41a1c","#5ddb3e","#acdb3e","royalblue2","darkcyan","#762a83","#7db777","gray"), cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged10KQ$id, las =2,cex.axis = 0.2)
axis(3, at=1:621, labels = merged10KQ$pop, las =2,cex.axis = 0.2)


test20KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/noThinning/621-merged-snps.20.Q")
merged20KQ <- cbind(test20KQ, popTable)
merged20KQ <- merged20KQ[order(merged20KQ$pop),]
barplot(t(as.matrix(merged20KQ[,1:20])), col=rainbow(20), cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged20KQ$id, las =2,cex.axis = 0.2)
axis(3, at=1:621, labels = merged20KQ$pop, las =2,cex.axis = 0.2)



### 1kb window run 
logtab12345 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/s12345.txt")
logtab12345$seed <- "Run1"
logtab22222 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/s22222.txt")
logtab22222$seed <- "Run2"
logtab54321 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/s54321.txt")
logtab54321$seed <- "Run3"
logtab <- rbind(logtab12345, logtab22222, logtab54321)

logplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= Loglikelihood,color = seed)) +
  geom_line(aes(x=K, y= Loglikelihood, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw()
logplot
CVplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= CV,color = seed)) +
  geom_line(aes(x=K, y= CV, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw()
CVplot

test2KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/621-merged-snps.2.Q")
merged2KQ <- cbind(test2KQ, popTable)
merged2KQ <- merged2KQ[order(merged2KQ$pop),]
barplot(t(as.matrix(merged2KQ[,1:2])), col=c("steelblue3","orange"), cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged6KQ$id, las =2,cex.axis = 0.2)
axis(3, at=c(140,295,314,325,341,354,356,357,359,366,370,376,383,388,395,400,408,435,464,490,505,523,547,557,595,619), labels = unique(merged10KQ$pop), las =2,cex.axis = 0.5)


test6KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/run3/621-merged-snps.6.Q")
merged6KQ <- cbind(test6KQ, popTable)
merged6KQ <- merged6KQ[order(merged6KQ$pop),]
barplot(t(as.matrix(merged6KQ[,1:6])), col=c("steelblue3","orange","#e41a1c","#5ddb3e","#acdb3e","royalblue2"), cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged6KQ$id, las =2,cex.axis = 0.2)
axis(3, at=c(140,295,314,325,341,354,356,357,359,366,370,376,383,388,395,400,408,435,464,490,505,523,547,557,595,619), labels = unique(merged10KQ$pop), las =2,cex.axis = 0.5)

test10KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/run3/621-merged-snps.10.Q")
merged10KQ <- cbind(test10KQ, popTable)
merged10KQ <- merged10KQ[order(merged10KQ$pop),]
barplot(t(as.matrix(merged10KQ[,1:10])), col=c("steelblue3","orange","#e41a1c","#5ddb3e","#acdb3e","royalblue2","darkcyan","#762a83","#7db777","gray"), cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged10KQ$id, las =2,cex.axis = 0.2)
axis(3, at=c(140,295,314,325,341,354,356,357,359,366,370,376,383,388,395,400,408,435,464,490,505,523,547,557,595,619), labels = unique(merged10KQ$pop), las =2,cex.axis = 0.5)

c17 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1")


test17KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/run3/621-merged-snps.17.Q")
merged17KQ <- cbind(test17KQ, popTable)
merged17KQ <- merged17KQ[order(merged17KQ$pop),]
barplot(t(as.matrix(merged17KQ[,1:17])), col=c17, cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged10KQ$id, las =2,cex.axis = 0.2)
axis(3, at=c(140,295,314,325,341,354,356,357,359,366,370,376,383,388,395,400,408,435,464,490,505,523,547,557,595,619), labels = unique(merged10KQ$pop), las =2,cex.axis = 0.5)

c20 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise")


test20KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/run1/621-merged-snps.20.Q")
merged20KQ <- cbind(test20KQ, popTable)
merged20KQ <- merged20KQ[order(merged20KQ$pop),]
barplot(t(as.matrix(merged20KQ[,1:20])), col=c20, cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged20KQ$id, las =2,cex.axis = 0.2)
axis(3, at=c(140,295,314,325,341,354,356,357,359,366,370,376,383,388,395,400,408,435,464,490,505,523,547,557,595,619), labels = unique(merged10KQ$pop), las =2,cex.axis = 0.5)

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "steelblue"
)

test25KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/run3/621-merged-snps.25.Q")
merged25KQ <- cbind(test25KQ, popTable)
merged25KQ <- merged25KQ[order(merged25KQ$pop),]
barplot(t(as.matrix(merged25KQ[,1:25])), col=c25, cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged20KQ$id, las = 2,cex.axis = 0.2)
axis(3, at=c(140,295,314,325,341,354,356,357,359,366,370,376,383,388,395,400,408,435,464,490,505,523,547,557,595,619), labels = unique(merged10KQ$pop), las =2,cex.axis = 0.5)

c26 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "steelblue"
)
test26KQ <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/1kb/run3/621-merged-snps.26.Q")
merged26KQ <- cbind(test26KQ, popTable)
merged26KQ <- merged26KQ[order(merged26KQ$pop),]
barplot(t(as.matrix(merged26KQ[,1:26])), col=c26, cex.names = 0.2,sort.by.Q = FALSE, border = NA, space = 0, xlab = "id", ylab = "ancestry",las=2)
axis(1, at=1:621, labels = merged20KQ$id, las = 2,cex.axis = 0.2)
axis(3, at=c(140,295,314,325,341,354,356,357,359,366,370,376,383,388,395,400,408,435,464,490,505,523,547,557,595,619), labels = unique(merged10KQ$pop), las =2,cex.axis = 0.5)


### Jacque's
nineK<-read.table("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/jacque-control/jacque-merged-snps.9.Q")
ids<-read.table("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/jacque-control/IdsxPop.txt")
ids_ordered<-ids[c(15,16,19,14,41,43,23:25,73:78,80:86,26:39,42,47:56,3:9,12,13,40,1,2,10,11,17,18,20:22,44:46,57:72,79,87:89),]
pop_order7<-t(as.matrix(nineK[c(15,16,19,14,41,43,23:25,73:78,80:86,26:39,42,47:56,3:9,12,13,40,1,2,10,11,17,18,20:22,44:46,57:72,79,87:89),]))
barplot(pop_order7,col=rainbow(9),ylab="Ancestry",border=NA,space=0,main="K = 9",xaxt="n")
axis(1,at=1:89,labels=ids_ordered$V1,las=2)

tenK<-read.table("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/jacque-control/jacque-merged-snps.10.Q")
ids<-read.table("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/jacque-control/IdsxPop.txt")
ids_ordered<-ids[c(15,16,19,14,41,43,23:25,73:78,80:86,26:39,42,47:56,3:9,12,13,40,1,2,10,11,17,18,20:22,44:46,57:72,79,87:89),]
pop_order8<-t(as.matrix(tenK[c(15,16,19,14,41,43,23:25,73:78,80:86,26:39,42,47:56,3:9,12,13,40,1,2,10,11,17,18,20:22,44:46,57:72,79,87:89),]))
barplot(pop_order8,col=c("steelblue3","orange","#e41a1c","#5ddb3e","#acdb3e","royalblue2","darkcyan","#762a83","#7db777","gray"),ylab="Ancestry",border=NA,space=0,main="K = 10",xaxt="n")
axis(1,at=1:89,labels=ids_ordered$V1,las=2)


logtab12345 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/jacque-control/s12345.txt")
logtab12345$seed <- "Run1"
logtab22222 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/jacque-control/s22222.txt")
logtab22222$seed <- "Run2"
logtab54321 <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/jacque-control/s54321.txt")
logtab54321$seed <- "Run3"
logtab <- rbind(logtab12345, logtab22222, logtab54321)

logplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= Loglikelihood,color = seed)) +
  geom_line(aes(x=K, y= Loglikelihood, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw()
logplot
CVplot <- ggplot(logtab) + 
  geom_point(aes(x=K, y= CV,color = seed)) +
  geom_line(aes(x=K, y= CV, color = seed)) +
  scale_x_continuous(breaks = c(2:30)) +
  theme_bw()
CVplot
