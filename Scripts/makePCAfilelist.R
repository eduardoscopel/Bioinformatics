library(data.table)
library(plyr)
library(scales)
### Read table and subset dataframe excluding unknown categories, and LOWCOV strains
sctab <- fread("sc_table.txt",header=TRUE)
sctab <- sctab[sctab$ecology_category != "Unknown",]
sctab <- sctab[sctab$ploidy != "Unknown",]
sctab <- sctab[sctab$ploidy != "LOWCOV",]
sctab <- sctab[sctab$aneuploidy != "LOWCOV",]
sctab <- sctab[sctab$aneuploidy_binary != "NA",]
sctab <- sctab[sctab$heterozygosity < 0.01,]
sctab <- sctab[sctab$ecology_category %in% names(table(sctab$ecology_category))[table(sctab$ecology_category)>20],]
sctab <- sctab[sctab$ploidy %in% names(table(sctab$ploidy))[table(sctab$ploidy)>20],]
nonMD <- sctab[sctab$monosporic_derivative == "No",]
diploids <- nonMD[nonMD$ploidy == 2,]
diploids[,1]
eco_ind <- fread("eco.ind",header=FALSE)

V2 <- rep("Unknown",773)
V3 <- rep("a",773)
eco_ind <- cbind(eco_ind, V2)
eco_ind <- cbind(eco_ind, V3)
i=0
for(strain in diploids$ID){
  print(strain)
  eco_ind[grep(paste("^",strain,"$",sep=""), eco_ind$V1),2] <- diploids[grep(paste("^",strain,"$",sep=""),diploids$ID),ecology_category]
}

for(x in 1:nrow(eco_ind)){
  eco_ind[x,3] <- diploids[x,8]
}
write.table(eco_ind,"eco.ind")
