library(data.table)
library(ggplot2)
library(fdrtool)

### Read table and subset dataframe excluding unknown categories, and LOWCOV strains
petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt",header=TRUE)
nrow(petertab)
petertab <- petertab[complete.cases(petertab$Aneuploidies),]                                    # Remove strains with no aneuploidy information
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="Euploid","No","Yes")
petertab$amp_response_no1 <- ifelse(petertab$Aneuploidy_type =="Loss" | 
                                      petertab$Aneuploidy_type =="Euploid" |
                                      petertab$Aneuploidy_type == "1only","No","Yes")
petertab$polyploidy <- ifelse(petertab$Ploidy == 1 |
                                petertab$Ploidy == 2,"1-2n","3-5n")

#########################################################################################################################################################
############################################################# CLADE ANALYSIS ###########################################################################
#########################################################################################################################################################

cladematrix <- as.matrix(table(petertab$Clades,petertab$polyploidy), nrow = 27, header = TRUE)

odds <- matrix(c(rep(0,54)), nrow = 27, dimnames = list(rownames(cladematrix),c("1-2n","3-5n")))
pv <- matrix(c(rep(0,54)), nrow = 27, dimnames = list(rownames(cladematrix),c("1-2n","3-5n")))
fdr <- matrix(c(rep(0,54)), nrow = 27, dimnames = list(rownames(cladematrix),c("1-2n","3-5n")))

for(i in 1:nrow(cladematrix)){
  odds[i,1] <- as.numeric(fisher.test(
    matrix(c(cladematrix[i,1],
             cladematrix[i,2],
             sum(cladematrix[,1]),
             sum(cladematrix[,2])),
           nrow=2,
           dimnames = list(c("1-2n","3-5n"),
                           c(rownames(cladematrix)[i],"all")))
  ,alternative = "less")$estimate)
  pv[i,1] <- fisher.test(
    matrix(c(cladematrix[i,1],
             cladematrix[i,2],
             sum(cladematrix[,1]),
             sum(cladematrix[,2])),
           nrow=2,
           dimnames = list(c("1-2n","3-5n"),
                           c(rownames(cladematrix)[i],"all")))
  ,alternative = "less")$p.value
  odds[i,2] <- as.numeric(fisher.test(
    matrix(c(cladematrix[i,1],
             cladematrix[i,2],
             sum(cladematrix[,1]),
             sum(cladematrix[,2])),
           nrow=2,
           dimnames = list(c("1-2n","3-5n"),
                           c(rownames(cladematrix)[i],"all")))
    ,alternative = "greater")$estimate)
  pv[i,2] <- fisher.test(
    matrix(c(cladematrix[i,1],
             cladematrix[i,2],
             sum(cladematrix[,1]),
             sum(cladematrix[,2])),
           nrow=2,
           dimnames = list(c("1-2n","3-5n"),
                           c(rownames(cladematrix)[i],"all")))
    ,alternative = "greater")$p.value
}
mt.rawp2adjp(pv[,1], proc = "BH", alpha = 0.05)
#########################################################################################################################################################
############################################################# ECOLOGY ANALYSIS ###########################################################################
#########################################################################################################################################################

ecomatrix <- as.matrix(table(petertab$`Ecological origins`,petertab$polyploidy), nrow = 23, header = TRUE)

odds <- matrix(c(rep(0,46)), nrow = 23, dimnames = list(rownames(ecomatrix),c("1-2n","3-5n")))
pv <- matrix(c(rep(0,46)), nrow = 23, dimnames = list(rownames(ecomatrix),c("1-2n","3-5n")))
fdr <- matrix(c(rep(0,46)), nrow = 23, dimnames = list(rownames(ecomatrix),c("1-2n","3-5n")))

for(i in 1:nrow(ecomatrix)){
  odds[i,1] <- as.numeric(fisher.test(
    matrix(c(ecomatrix[i,1],
             ecomatrix[i,2],
             sum(ecomatrix[,1]),
             sum(ecomatrix[,2])),
           nrow=2,
           dimnames = list(c("1-2n","3-5n"),
                           c(rownames(ecomatrix)[i],"all")))
    ,alternative = "less")$estimate)
  pv[i,1] <- fisher.test(
    matrix(c(ecomatrix[i,1],
             ecomatrix[i,2],
             sum(ecomatrix[,1]),
             sum(ecomatrix[,2])),
           nrow=2,
           dimnames = list(c("1-2n","3-5n"),
                           c(rownames(ecomatrix)[i],"all")))
    ,alternative = "less")$p.value
  odds[i,2] <- as.numeric(fisher.test(
    matrix(c(ecomatrix[i,1],
             ecomatrix[i,2],
             sum(ecomatrix[,1]),
             sum(ecomatrix[,2])),
           nrow=2,
           dimnames = list(c("1-2n","3-5n"),
                           c(rownames(ecomatrix)[i],"all")))
    ,alternative = "greater")$estimate)
  pv[i,2] <- fisher.test(
    matrix(c(ecomatrix[i,1],
             ecomatrix[i,2],
             sum(ecomatrix[,1]),
             sum(ecomatrix[,2])),
           nrow=2,
           dimnames = list(c("1-2n","3-5n"),
                           c(rownames(ecomatrix)[i],"all")))
    ,alternative = "greater")$p.value
}
mt.rawp2adjp(pv[,1], proc = "BH", alpha = 0.05)
