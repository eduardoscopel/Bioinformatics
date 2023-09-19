library(data.table)
library(plyr)
library(dplyr)
library(phylolm)
library(scales)
library(epitools)
library(ggplot2)
library(MASS)
library(ape)
library(treeio)
library(rr2)

### Read table and subset dataframe excluding unknown categories, and LOWCOV strains
# Diploids tree
#tree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/tree/alltrees.suptree")
# Diploids and polyploids tree
tree <- read.tree("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/tree/peter_nonad_new.treefile")
rtree <- root(tree, c("EM14S01_3B","EN14S01","GE14S01_7B"))
ntree <- compute.brlen(rtree, method = "Grafen", power=1)

petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt",header=TRUE)
sctab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /sc_table_most_recent.txt", header = TRUE)
nrow(petertab)
#petertab <- petertab[complete.cases(petertab$Aneuploidies),]                                    # Remove strains with no aneuploidy information
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="Euploid","No","Yes")
petertab$amp_response <- ifelse(petertab$Aneuploidy_type =="Euploid","No","Yes")
petertab$amp_binary <- ifelse(petertab$amp_response == "Yes", 1,0)
petertab$admixture <- ifelse(petertab$Clades == "M. Mosaic" |
                               petertab$Clades == "8. Mixed origin",
                             "Admixed","Non-admixed")
petertab$polyploidy <- ifelse(petertab$Ploidy == 3 | 
                                petertab$Ploidy == 4 |
                                petertab$Ploidy == 5, "poly",
                              ifelse(petertab$Ploidy == 1, "haploid",
                                     "diploid"))
petertab$polyploidy <- ifelse(petertab$Ploidy == 3 | 
                                petertab$Ploidy == 4 |
                                petertab$Ploidy == 5, "poly",
                              ifelse(petertab$Ploidy == 1, "haploid",
                                     "diploid"))
nrow(petertab)

# Get important categories from my table into peter df
for(i in 1:nrow(petertab)){
  for(j in 1:nrow(sctab)){
    if(!is.element(petertab$`Isolate name`[i],sctab$ID)){
      petertab$heterozygosity[i] <- NA
      petertab$new_eco[i] <- NA
      petertab$basename[i] <- NA
    }
    else if(petertab$`Isolate name`[i] == sctab$ID[j]){
      petertab$heterozygosity[i] <- sctab$heterozygosity[j]
      petertab$new_eco[i] <- sctab$ecology_category[j]
      petertab$basename[i] <- sctab$basename[j]
    }
  }
}

# Remove HO deleted strains
#petertab <- petertab[petertab$`HO deletion` == "no",]
nrow(petertab)

newdf <- data.frame(row.names = ntree$tip.label)
rownames(newdf)[536] <- "UWOPS034614"
for(i in 1:nrow(newdf)){
  j <- which(petertab$basename == rownames(newdf)[i])
  print(paste("i is ", i))
  print(paste("j is ", j))
  newdf$amp[i] <- petertab$amp_binary[j]
  newdf$amp_no1[i] <- petertab$amp_binary_no1[j]
  newdf$eco[i] <- petertab$`Ecological origins`[j]
  newdf$ploidy[i] <- petertab$polyploidy[j]
  newdf$aneup[i] <- petertab$aneuploidy_response[j]
  newdf$het[i] <- petertab$heterozygosity[j]
}
rownames(newdf)[536] <-ntree$tip.label[536]

#dfnonhet <- data.frame(amp = nonhet$amp_binary, 
#                       eco = nonhet$`Ecological origins`, 
#                       het = nonhet$heterozygosity, 
#                       clade = nonhet$Clades, 
#                       ploidy = nonhet$Ploidy,
#                       row.names = nonhet$basename)

### phylogenetic logistic regression with ecology as independent variable 
eglm <- glm(amp_binary_no1~`Ecological origins`, family=binomial(link="logit"), newdf)
nglm <- glm(amp_response_no1 ~ 1,  family = binomial(link="logit"), newdf)

ecoplpglm <- phyloglm(amp_no1 ~ eco, phy = ntree, data = newdf, boot = 100, btol = 18)
ecoplpglm
summary(ecoplpglm)
ecoplaepglm <- -log(ecoplpglm$bootstrap[,"alpha"])
t.test(ecoplaepglm,  mu = -4)

ecopglm <- phyloglm(amp ~ eco, phy = ntree, data = newdf, btol=18, boot=100)
ecopglm
summary(ecopglm)
ecoaepglm <- -log(ecopglm$bootstrap[,"alpha"])
t.test(ecoaepglm,  mu = -4)

plpglm <- phyloglm(amp_no1 ~ ploidy, phy = ntree, data = newdf, boot = 100)
plpglm
summary(plpglm)
plaepglm <- -log(plpglm$bootstrap[,"alpha"])
t.test(plaepglm,  mu = -4)

hetpglm <- phyloglm(amp_no1 ~ het, phy = ntree, data = newdf, boot =100)
hetpglm
summary(hetpglm)

newdaepglm <- -log(epglm$bootstrap[,"alpha"])
t.test(aepglm,  mu = -4)


# Karyotype imbalance
newdf$KI <- ifelse(newdf$aneup == "Yes" | 
                     newdf$ploidy == "poly", 1, 0)
newdf[584,"KI"] <- 0
table(newdf$eco,newdf$KI)

kpglm <- phyloglm(KI ~ eco, phy = ntree, data = newdf, btol = 18, boot=100)
kiaepglm <- -log(kpglm$bootstrap[,"alpha"])
t.test(kiaepglm,  mu = -4)


# null model
npglm <- phyloglm(amp~1, phy=ntree, data=dfnonhet)

#clade model
cpglm <- phyloglm(amp~clade, phy = ntree, data = dfnonhet,btol=19, boot=100)
acpglm <- -log(cpglm$bootstrap[,"alpha"])
cglm <- glm(amp_response~Clades, family=binomial(link="logit"), nonhet)
t.test(acpglm,  mu = -4)
