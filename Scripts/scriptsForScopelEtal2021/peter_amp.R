library(data.table)
library(plyr)
library(scales)
library(epitools)
library(ggplot2)
library(MASS)

### Read table and subset dataframe excluding unknown categories, and LOWCOV strains
petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt",header=TRUE)
sctab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /sc_table_most_recent.txt", header = TRUE)
nrow(petertab)
petertab <- petertab[complete.cases(petertab$Aneuploidies),]                                    # Remove strains with no aneuploidy information
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="euploid","No","Yes")
petertab$amp_response <- ifelse(petertab$Aneuploidy_type =="Deletion" | 
                                  petertab$Aneuploidy_type =="None","No","Yes")
petertab$polyploidy <- ifelse(petertab$Ploidy == 1 |
                                petertab$Ploidy == 2,"1-2n","3-5n")



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
petertab <- petertab[petertab$`HO deletion` == "no",]


#########################################################################################################################################################
############################################################# PLOIDY ANALYSIS ###########################################################################
#########################################################################################################################################################
# Create a ploidy matrix with the entire data set
ploidymatrix <- as.matrix(table(petertab$Ploidy,petertab$amp_response), nrow=5, header = TRUE)
ploidymatrix
fisher.test(ploidymatrix)
ploidymatrix[1:2,]
fisher.test(ploidymatrix[1:2,])
ploidymatrix[3:5,]
fisher.test(ploidymatrix[3:5,])

# Create a ploidy matrix without haploids
nonhap <- petertab[petertab$Ploidy != 1,]
nonhap <- nonhap[nonhap$Aneuploidy_type != "Euploid",]
nrow(nonhap)
ploidymatrix <- as.matrix(table(nonhap$Ploidy,nonhap$Aneuploidy_type), nrow=4, header = TRUE)
ploidymatrix <- ploidymatrix[,c(2,3,1)]
ploidymatrix
fisher.test(ploidymatrix)
# Remove pentaploids
ploidymatrix[1:3,]
fisher.test(ploidymatrix[1:3,])
# Remove tetraploids
ploidymatrix[c(1,2,4),]
fisher.test(ploidymatrix[c(1,2,4),])
# Remove triploids
ploidymatrix[c(1,3,4),]
fisher.test(ploidymatrix[c(1,3,4),])
# Remove diploids
ploidymatrix[c(2,3,4),]
fisher.test(ploidymatrix[c(2,3,4),])

ploidymatrix[2:4,]
fisher.test(ploidymatrix[2:4,])
ploidymatrix[1:2,]
fisher.test(ploidymatrix[1:2,])
ploidymatrix[c(1,3,4),]
fisher.test(ploidymatrix[c(1,3,4),])


# Diploids vs. polyploids
pploidymatrix <- as.matrix(table(nonhap$polyploidy, nonhap$Aneuploidy_type),nrow=2, header=TRUE)
rownames(pploidymatrix)[1] <- "2n"
pploidymatrix <- pploidymatrix[,c(2,3,1)]
pploidymatrix
fisher.test(pploidymatrix)
pploidymatrix[,1:3]
fisher.test(pploidymatrix[,1:3])

# add a column of totals/proportions to the matrix
for(i in 1:nrow(ploidymatrix)){
  if(i==1){
    total_ploidy <- c(sum(ploidymatrix[i,]))
    proportion_ploidy <- c((ploidymatrix[i,2])/(total_ploidy[i]))
  }
  else{
    total_ploidy <- c(total_ploidy,sum(ploidymatrix[i,]))
    proportion_ploidy <- c(proportion_ploidy,(ploidymatrix[i,2])/(total_ploidy[i]))
  }
}
ploidymatrix <- cbind(ploidymatrix,total_ploidy)
ploidymatrix <- cbind(ploidymatrix, proportion_ploidy)

# add columns with the confidence intervals for proportions (binomial test) 
for(i in 1:nrow(ploidymatrix)){
  if(i==1){
    p_CILow <- c(binom.test(ploidymatrix[i,2],ploidymatrix[i,3])$conf.int[1])
    p_CIHigh <- c(binom.test(ploidymatrix[i,2],ploidymatrix[i,3])$conf.int[2])
  }
  else{
    p_CILow <- c(p_CILow,binom.test(ploidymatrix[i,2],ploidymatrix[i,3])$conf.int[1])
    p_CIHigh <- c(p_CIHigh,binom.test(ploidymatrix[i,2],ploidymatrix[i,3])$conf.int[2])
  }
}
ploidymatrix <- cbind(ploidymatrix, p_CILow)
ploidymatrix <- cbind(ploidymatrix, p_CIHigh)

# Set up variables for ploidy proportion plot
p_factors = c()
p_colors = c()
gray_scale = c()
pOdds = c()
pCILow = c()
pCIHigh = c()
pLabels = c()
pcat=row.names(ploidymatrix)

# run a Fisher's Exact Test on each row, and prepare lists to create the odds ratio plot
for(row in 1:nrow(ploidymatrix)){
  if(pcat[row]==1){
    pcat[row]<- "Haploid"
    pLabels <- append(pLabels,paste(pcat[row],", N=",ploidymatrix[row,3],sep=""))
  }
  else if(pcat[row]==2){
    pcat[row]<- "Diploid"
    pLabels <- append(pLabels,paste(pcat[row],", N=",ploidymatrix[row,3],sep=""))
  }
  else if(pcat[row]==3){
    pcat[row]<- "Triploid"
    pLabels <- append(pLabels,paste(pcat[row],", N=",ploidymatrix[row,3],sep=""))
  }
  else if(pcat[row]==4){
    pcat[row]<- "Tetraploid"
    pLabels <- append(pLabels,paste(pcat[row],", N=",ploidymatrix[row,3],sep=""))
  }
  else if(pcat[row]==5){
    pcat[row]<- "Pentaploid"
    pLabels <- append(pLabels,paste(pcat[row],", N=",ploidymatrix[row,3],sep=""))
  }
}

p_df <- data.frame(
  yAxis = length(pLabels[1:5]):1,
  ploidymatrix[,4], ploidymatrix[,5], ploidymatrix[,6])

# create the odds ratio plot
p_plot <- ggplot(p_df, aes(x= ploidymatrix[,4], y = yAxis))
p_plot <- p_plot + 
  geom_vline(aes(xintercept = (sum(ploidymatrix[,2]))/(sum(ploidymatrix[,3]))), size = .5, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ploidymatrix[,6], xmin = ploidymatrix[,5]), size = 2, height = 0.5,color = "gray50") +
  geom_point(size = 5.0) +
  theme_bw() +
  scale_y_continuous(breaks = length(pLabels):1, labels = pLabels) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  ylab("") +
  xlab("Proportion of amplifications") +
  theme(axis.title.x = element_text(h=0.5,size=16))+
  theme(axis.text= element_text(face="bold",color="gray50", size=16))+
  ggtitle("Proportion of amplifications per ploidy")+
  theme(plot.title = element_text(face="bold",h=0.5))
p_plot <- p_plot +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
p_plot <- p_plot + theme(legend.position = "bottom")
png("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/amplifications/ploidy_amp_propplot_peter.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
p_plot
dev.off()

#########################################################################################################################################################

#########################################################################################################################################################
############################################################# HETEROZYGOSITY ANALYSIS ###################################################################
#########################################################################################################################################################
# Plot heterozygosity histogram including highly heterozygous and polyploid strains
het_hist <- ggplot(petertab, aes(x=heterozygosity))
het_hist <- het_hist + geom_histogram(aes(fill=amp_response),
                                      color="black", 
                                      bins=30,
                                      alpha=0.7)+
  scale_x_continuous(breaks=seq(0,0.01,0.001), limits=c(0,0.01)) +
  scale_y_continuous(breaks=seq(0,125,25), limits=c(0,125)) +
  theme_bw() +
  ylab("Total strains") +
  xlab('Heterozygosity') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  theme(legend.position = "bottom")
png("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/amplifications//het_amp_propplot_peter.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
het_hist
dev.off()

# Plot heterozygosity histogram including on diploid strains, excluding highly heterozygous strains (>0.5%)
diploids <- petertab[petertab$Ploidy == 2,]
dip_nonhet <- diploids[diploids$heterozygosity <= 0.005,]
het_hist <- ggplot(dip_nonhet, aes(x=heterozygosity))
het_hist <- het_hist + geom_histogram(aes(fill=amp_response),
                                      color="black", 
                                      bins=30,
                                      alpha=0.7)+
  scale_x_continuous(breaks=seq(0,0.01,0.001), limits=c(0,0.01)) +
  scale_y_continuous(breaks=seq(0,125,25), limits=c(0,125)) +
  theme_bw() +
  ylab("Diploid strains") +
  xlab('Heterozygosity') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  theme(legend.position = "bottom")
png("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/amplifications/2n_het_amp_propplot_peter.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
het_hist
dev.off()

#########################################################################################################################################################

#########################################################################################################################################################
############################################################# PETER ECOLOGY ANALYSIS ####################################################################
#########################################################################################################################################################
# Remove categories with less than 10 strains
oldecopeter <- petertab[petertab$`Ecological origins` %in% names(table(petertab$`Ecological origins`))[table(petertab$`Ecological origins`)>=10],]
oldecomatrix <- as.matrix(table(oldecopeter$`Ecological origins`, oldecopeter$amp_response),
                          nrow=length(unique(oldecopeter$`Ecological origins`)), 
                          header = TRUE)
diploids <- petertab[petertab$Ploidy==2,]
dip_oldecopeter <- diploids[diploids$`Ecological origins` %in% names(table(diploids$`Ecological origins`))[table(diploids$`Ecological origins`)>=10],]
dip_oldecomatrix <- as.matrix(table(dip_oldecopeter$`Ecological origins`, dip_oldecopeter$amp_response),
                          nrow=length(unique(dip_oldecopeter$`Ecological origins`)), 
                          header = TRUE)
oldecomatrix
fisher.test(oldecomatrix, workspace = 2e10)
oldecomatrix[c(1:15, 17:21),]
fisher.test(oldecomatrix[c(1:15, 17:21),], simulate.p.value = TRUE)
fisher.test(oldecomatrix[c(2:15, 17:21),], simulate.p.value = TRUE)
fisher.test(oldecomatrix[c(3:15, 17:21),], simulate.p.value = TRUE)
oldecomatrix[c(3:15, 17:21),]
fisher.test(oldecomatrix[c(3:15, 17:21),], simulate.p.value = TRUE)
oldecomatrix[c(1,2,16),]
fisher.test(oldecomatrix[c(1,2,16),])

dip_oldecomatrix
fisher.test(dip_oldecomatrix, simulate.p.value = TRUE)
dip_oldecomatrix[c(1:15, 17:21),]
fisher.test(dip_oldecomatrix[c(1:15, 17:21),], simulate.p.value = TRUE)


# add a totals column to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){
    Total <- c(sum(scmatrix[i,]))
  }
  else{
    Total <- c(Total, sum(scmatrix[i,]))
  }
}
scmatrix<-cbind(scmatrix,Total)

# add a proportions column to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){proportion <- c((scmatrix[i,2])/(scmatrix[i,3]))}
  else{proportion <- c(proportion,(scmatrix[i,2])/(scmatrix[i,3]))}
}  
scmatrix <- cbind(scmatrix, proportion)

# add a column with the confidence intervales (binomial distribution) to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){
    CILow <- c(binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[1])
    CIHigh <- c(binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[2])
  }
  else{
    CILow <- c(CILow,binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[1])
    CIHigh <- c(CIHigh,binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[2])
    
  }
}  
scmatrix <- cbind(scmatrix, CILow)
scmatrix <- cbind(scmatrix, CIHigh)


# Reorgnaize ecological categories in the scmatrix (domesticated, clinical, wild)
### Including all categories
#scmatrix <- rbind(scmatrix[c(1,2,3,4,5,6,7,12,14,16,17,18,23),],
#                  scmatrix[c(10,11),],
#                  scmatrix[c(8,9,13,15,19,20,22),],
#                  scmatrix[21,])
#rownames(scmatrix)[23] <- "Unknown"

### Ecluding categories with less than 10 strains
scmatrix <- rbind(scmatrix[c(1,2,3,4,5,6,7,12,15,16,21),],
                  scmatrix[c(10,11),],
                  scmatrix[c(8,9,13,14,17,18,20),],
                  scmatrix[19,])
rownames(scmatrix)[21] <- "Unknown"


ecocat <- row.names(scmatrix)

# Set up variables to create the proportion plot
my_factors = c()
my_colors = c()
gray_scale = c()
boxOdds = c()
boxCILow = c()
boxCIHigh = c()
boxLabels = c()

# populate variables to add to the proportion plot
for(row in 1:nrow(scmatrix)){
  if(row>nrow(scmatrix)){break}
  print(ecocat[row])
  print((matrix(c(scmatrix[row,2],
                  sum(scmatrix[,2])-scmatrix[row,2],
                  scmatrix[row,1],
                  sum(scmatrix[,1])-scmatrix[row,1]),
                nrow=2)))
  ftest <- fisher.test(matrix(c(scmatrix[row,2],
                                sum(scmatrix[,2])-scmatrix[row,2],
                                scmatrix[row,1],
                                sum(scmatrix[,1])-scmatrix[row,1]),
                              nrow=2))
  boxCILow <- append(boxCILow, ftest$conf.int[1])
  boxCIHigh <- append(boxCIHigh, ftest$conf.int[2])
  boxOdds <- append(boxOdds, ftest$estimate)
  if(ftest[1]<0.01){
    boxLabels <- append(boxLabels,paste(ecocat[row],"*, N=",scmatrix[row,3],sep=""))
  }
  else{boxLabels <- append(boxLabels,paste(ecocat[row],", N=",scmatrix[row,3],sep=""))}
  if(is.element(ecocat[row], c("Human, clinical", "Human"))){
    my_factors <- append(my_factors, "Clinical")
    my_colors <- append(my_colors, "#FC724F")
  }
  else if(is.element(ecocat[row],c("Flower", "Fruit", "Insect", "Nature", "Soil", "Tree", "Water"))){
    my_factors <- append(my_factors, "Wild")
    my_colors <- append(my_colors, "#316857")
  }
  else if(ecocat[row] == "Unknown"){
    my_factors <- append(my_factors, "Unknown")
    my_colors <- append(my_colors, "black")
  }
  else{
    my_factors <- append(my_factors, "Domesticated")
    my_colors <- append(my_colors, "#225EA8")
  }
  print(ftest)
}

# Create a data frame with proportions and confidence interval values
df <- data.frame(
  yAxis = length(boxLabels):1,
  scmatrix[,4], scmatrix[,5], scmatrix[,6])

# create the odds ratio plot
orplot <- ggplot(df, aes(x= scmatrix[,4], y = yAxis))
orplot <- orplot + 
  geom_vline(aes(xintercept = (sum(scmatrix[,2]))/(sum(scmatrix[,3]))), size = .5, linetype = "dashed") +
  geom_errorbarh(aes(xmax = scmatrix[,6], xmin = scmatrix[,5]), size = 2, height = 0.5,color = "gray50") +
  geom_point(color = my_colors, size = 5.0) +
  scale_color_manual(values=my_colors) +
  theme_bw() +
  scale_y_continuous(breaks = length(boxLabels):1, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) +
  ylab("") +
  xlab('Proportion of amplifications') +
  theme(axis.title.x = element_text(h=0.5))+
  theme(axis.text= element_text(face="bold",color="gray50", size=12))+
  theme(axis.text.y= element_text(color=my_colors))+
  ggtitle("Proportion of amplifications per environment \n (>10 strains) in diploids")+
  theme(plot.title = element_text(face="bold",h=0.5))
orplot <- orplot +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
orplot <- orplot + theme(legend.position = "bottom")
png("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/amplifications/2n_old_eco_amp_pplot_peter.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
orplot
dev.off()
#########################################################################################################################################################

#########################################################################################################################################################
###################################################### NEW ECOLOGY ANALYSIS #############################################################################
#########################################################################################################################################################
newecopeter <- petertab[petertab$new_eco %in% names(table(petertab$new_eco))[table(petertab$new_eco)>=10],]
newecomatrix <- as.matrix(table(newecopeter$new_eco, newecopeter$amp_response), 
                          nrow=length(unique(newecopeter$new_eco)), 
                          header = TRUE)

dip_newecopeter <- diploids[diploids$new_eco %in% names(table(diploids$new_eco))[table(diploids$new_eco)>=10],]
dip_newecomatrix <- as.matrix(table(dip_newecopeter$new_eco, dip_newecopeter$amp_response), 
                          nrow=length(unique(dip_newecopeter$new_eco)), 
                          header = TRUE)
# convert df to matrix
scmatrix <- dip_newecomatrix

# add a totals column to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){
    Total <- c(sum(scmatrix[i,]))
  }
  else{
    Total <- c(Total, sum(scmatrix[i,]))
  }
}
scmatrix<-cbind(scmatrix,Total)

# add a proportions column to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){proportion <- c((scmatrix[i,2])/(scmatrix[i,3]))}
  else{proportion <- c(proportion,(scmatrix[i,2])/(scmatrix[i,3]))}
}  
scmatrix <- cbind(scmatrix, proportion)

# add a column with the confidence intervales (binomial distribution) to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){
    CILow <- c(binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[1])
    CIHigh <- c(binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[2])
  }
  else{
    CILow <- c(CILow,binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[1])
    CIHigh <- c(CIHigh,binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[2])
    
  }
}  
scmatrix <- cbind(scmatrix, CILow)
scmatrix <- cbind(scmatrix, CIHigh)


# Reorgnaize ecological categories in the scmatrix (domesticated, clinical, wild)
### Including all categories
#scmatrix <- rbind(scmatrix[c(1,2,3,6,7,8,12,14,15,16,17,20),],
#                  scmatrix[4,],
#                  scmatrix[c(5,9,10,11,13,18,19),])
#rownames(scmatrix)[13] <- "Clinical"

### Ecluding categories with less than 10 strains
scmatrix <- rbind(scmatrix[c(1,2,3,6,7,11,12,15),],
                  scmatrix[4,],
                  scmatrix[c(5,8,9,10,13,14),])
rownames(scmatrix)[9] <- "Clinical"


ecocat <- row.names(scmatrix)
# Set up variables to create the proportion plot
my_factors = c()
my_colors = c()
gray_scale = c()
boxOdds = c()
boxCILow = c()
boxCIHigh = c()
boxLabels = c()

# populate variables to add to the proportion plot
for(row in 1:nrow(scmatrix)){
  if(row>nrow(scmatrix)){break}
  print(ecocat[row])
  print((matrix(c(scmatrix[row,2],
                  sum(scmatrix[,2])-scmatrix[row,2],
                  scmatrix[row,1],
                  sum(scmatrix[,1])-scmatrix[row,1]),
                nrow=2)))
  ftest <- fisher.test(matrix(c(scmatrix[row,2],
                                sum(scmatrix[,2])-scmatrix[row,2],
                                scmatrix[row,1],
                                sum(scmatrix[,1])-scmatrix[row,1]),
                              nrow=2))
  boxCILow <- append(boxCILow, ftest$conf.int[1])
  boxCIHigh <- append(boxCIHigh, ftest$conf.int[2])
  boxOdds <- append(boxOdds, ftest$estimate)
  if(ftest[1]<0.01){
    boxLabels <- append(boxLabels,paste(ecocat[row],"*, N=",scmatrix[row,3],sep=""))
  }
  else{boxLabels <- append(boxLabels,paste(ecocat[row],", N=",scmatrix[row,3],sep=""))}
  if(is.element(ecocat[row], c("Clinical"))){
    my_factors <- append(my_factors, "Clinical")
    my_colors <- append(my_colors, "#FC724F")
  }
  else if(is.element(ecocat[row],c("Cocoa","Flower", "Fruit", "Insect", "Tree", "Water"))){
    my_factors <- append(my_factors, "Wild")
    my_colors <- append(my_colors, "#316857")
  }
  else{
    my_factors <- append(my_factors, "Domesticated")
    my_colors <- append(my_colors, "#225EA8")
  }
  print(ftest)
}

# Create a data frame with proportions and confidence interval values
df <- data.frame(
  yAxis = length(boxLabels):1,
  scmatrix[,4], scmatrix[,5], scmatrix[,6])

# create the odds ratio plot
orplot <- ggplot(df, aes(x= scmatrix[,4], y = yAxis))
orplot <- orplot + 
  geom_vline(aes(xintercept = (sum(scmatrix[,2]))/(sum(scmatrix[,3]))), size = .5, linetype = "dashed") +
  geom_errorbarh(aes(xmax = scmatrix[,6], xmin = scmatrix[,5]), size = 2, height = 0.5,color = "gray50") +
  geom_point(color = my_colors, size = 5.0) +
  scale_color_manual(values=my_colors) +
  theme_bw() +
  scale_y_continuous(breaks = length(boxLabels):1, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) +
  ylab("") +
  xlab('Proportion of amplifications') +
  theme(axis.title.x = element_text(h=0.5))+
  theme(axis.text= element_text(face="bold",color="gray50", size=10))+
  theme(axis.text.y= element_text(color=my_colors))+
  ggtitle("Proportion of amplifications per \n environment (>10 strains) in diploids")+
  theme(plot.title = element_text(size=11, face="bold", hjust=0.2))
orplot <- orplot +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
orplot <- orplot + theme(legend.position = "bottom")
png("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/amplifications/2n_new_eco_amp_pplot_peter.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
orplot
dev.off()

#########################################################################################################################################################

#########################################################################################################################################################
############################################################ CLADE ANALYSIS #############################################################################
#########################################################################################################################################################
poppeter <- petertab[complete.cases(petertab$Clades),]                                    # Remove strains with unassigned clades
poppeter <- poppeter[poppeter$Clades %in% names(table(poppeter$Clades))[table(poppeter$Clades)>=10],]
popmatrix <- as.matrix(table(poppeter$Clades, poppeter$amp_response),
                       nrow=length(unique(poppeter$Clades)), 
                       header = TRUE)
dip_poppeter <- poppeter[poppeter$Ploidy == 2,]
dip_poppeter <- dip_poppeter[dip_poppeter$Clades %in% names(table(dip_poppeter$Clades))[table(dip_poppeter$Clades)>=10],]
dip_popmatrix <- as.matrix(table(dip_poppeter$Clades, dip_poppeter$amp_response),
                           nrow=length(unique(dip_poppeter$Clades)), 
                           header = TRUE)


# convert df to matrix
scmatrix <- popmatrix
scmatrix
fisher.test(scmatrix, workspace = 2e09)
scmatrix[c(2,4:7,9:13),]
fisher.test(scmatrix[c(2,4:7,9:12),], simulate.p.value = TRUE)

dip_popmatrix
fisher.test(dip_popmatrix, simulate.p.value = TRUE)




# add a totals column to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){
    Total <- c(sum(scmatrix[i,]))
  }
  else{
    Total <- c(Total, sum(scmatrix[i,]))
  }
}
scmatrix<-cbind(scmatrix,Total)

# add a proportions column to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){proportion <- c((scmatrix[i,2])/(scmatrix[i,3]))}
  else{proportion <- c(proportion,(scmatrix[i,2])/(scmatrix[i,3]))}
}  
scmatrix <- cbind(scmatrix, proportion)

# add a column with the confidence intervales (binomial distribution) to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){
    CILow <- c(binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[1])
    CIHigh <- c(binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[2])
  }
  else{
    CILow <- c(CILow,binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[1])
    CIHigh <- c(CIHigh,binom.test(scmatrix[i,2],scmatrix[i,3])$conf.int[2])
    
  }
}  
scmatrix <- cbind(scmatrix, CILow)
scmatrix <- cbind(scmatrix, CIHigh)


# Reorgnaize ecological categories in the scmatrix (domesticated, clinical, wild)
### Including all categories
#scmatrix <- rbind(scmatrix[c(1,3,5,12,18,19,20,22,23,24,26),],
#                  scmatrix[2,],
#                  scmatrix[c(4,6,7,8,9,10,11,13,14,15,16,17,21),],
#                  scmatrix[c(25,27),])
#rownames(scmatrix)[12] <- "10._French_Guiana_human_"

### Excluding categories with less than 10 strains
scmatrix <- rbind(scmatrix[c(1,4,5,7,8,9,10,11),],
                  scmatrix[2,],
                  scmatrix[c(3,6),],
                  scmatrix[c(12,13),])
rownames(scmatrix)[9] <- "10._French_Guiana_human_"

ecocat=row.names(scmatrix)

# Set up variables to create the proportion plot
my_factors = c()
my_colors = c()
gray_scale = c()
boxOdds = c()
boxCILow = c()
boxCIHigh = c()
boxLabels = c()

# populate variables to add to the proportion plot
for(row in 1:nrow(scmatrix)){
  if(row>nrow(scmatrix)){break}
  print(ecocat[row])
  print((matrix(c(scmatrix[row,2],
                  sum(scmatrix[,2])-scmatrix[row,2],
                  scmatrix[row,1],
                  sum(scmatrix[,1])-scmatrix[row,1]),
                nrow=2)))
  ftest <- fisher.test(matrix(c(scmatrix[row,2],
                                sum(scmatrix[,2])-scmatrix[row,2],
                                scmatrix[row,1],
                                sum(scmatrix[,1])-scmatrix[row,1]),
                              nrow=2))
  boxCILow <- append(boxCILow, ftest$conf.int[1])
  boxCIHigh <- append(boxCIHigh, ftest$conf.int[2])
  boxOdds <- append(boxOdds, ftest$estimate)
  if(ftest[1]<0.01){
    boxLabels <- append(boxLabels,paste(ecocat[row],"*, N=",scmatrix[row,3],sep=""))
  }
  else{boxLabels <- append(boxLabels,paste(ecocat[row],", N=",scmatrix[row,3],sep=""))}
  if(row == 9){
    my_factors <- append(my_factors, "Clinical")
    my_colors <- append(my_colors, "#FC724F")
  }
  else if(row >= 12){
    my_factors <- append (my_factors, "Mosaic")
    my_colors <- append(my_colors, "Pink")
  }
  else if(row < 9){
    my_factors <- append(my_factors, "Domesticated")
    my_colors <- append(my_colors, "#225EA8")
  }
  else{
    my_factors <- append(my_factors, "Wild")
    my_colors <- append(my_colors, "#316857")
  }
  print(ftest)
}

# Create a data frame with proportions and confidence interval values
df <- data.frame(
  yAxis = length(boxLabels):1,
  scmatrix[,4], scmatrix[,5], scmatrix[,6])

# create the odds ratio plot
orplot <- ggplot(df, aes(x= scmatrix[,4], y = yAxis))
orplot <- orplot + 
  geom_vline(aes(xintercept = (sum(scmatrix[,2]))/(sum(scmatrix[,3]))), size = .5, linetype = "dashed") +
  geom_errorbarh(aes(xmax = scmatrix[,6], xmin = scmatrix[,5]), size = 2, height = 0.5,color = "gray50") +
  geom_point(color = my_colors, size = 5.0) +
  scale_color_manual(values=my_colors) +
  theme_bw() +
  scale_y_continuous(breaks = length(boxLabels):1, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) +
  ylab("") +
  xlab('Proportion of amplifications') +
  theme(axis.title.x = element_text(hjust = 0.5))+
  theme(axis.text= element_text(face="bold",color="gray50", size=10))+
  theme(axis.text.y= element_text(color=my_colors))+
  ggtitle("Proportion of amplifications per \n clade (>10 strains) in diploids")+
  theme(plot.title = element_text(size=11, face="bold",h=2))
orplot <- orplot +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
orplot <- orplot + theme(legend.position = "bottom")
png("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/amplifications/2n_clade_amp_propplot_peter.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
orplot
dev.off()

