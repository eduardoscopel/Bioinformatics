library(epitools)
library(ggplot2)
library(data.table)
library(MASS)
# Read table and select Peter strains only
sctab <- fread("sc_table_pop_peter.txt",header=TRUE) 
peter <- sctab[sctab$reference == "Peter_et_al_2018",]
peter <- peter[peter$monosporic_derivative == "No",]
peter <- peter[peter$ecology_category %in% names(table(peter$ecology_category))[table(peter$ecology_category)>=10],]

attach(peter)

### Ecology analysis
# convert df to matrix
scmatrix <- as.matrix(table(ecology_category, aneuploidy_binary),nrow=21, header = TRUE)


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

# Set up variables to create the proportion plot
my_factors = c()
my_colors = c()
gray_scale = c()
boxOdds = c()
boxCILow = c()
boxCIHigh = c()
boxLabels = c()

# Reorgnaize ecological categories in the scmatrix (domesticated, clinical, wild)
scmatrix <- rbind(scmatrix[c(1,2,3,6,7,8,12,13,16),],
                  scmatrix[4,],
                  scmatrix[c(5,9,10,11,14,15),])
rownames(scmatrix)[10] <- "Clinical"
ecocat <- row.names(scmatrix)

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
 if(ecocat[row]=="Clinical"){
   my_factors <- append(my_factors, "Clinical")
   my_colors <- append(my_colors, "#FC724F")
 }
 else if(is.element(ecocat[row],c("Cocoa","Flower", "Fruit", "Insect", "Other_plants", "Tree", "Water"))){
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
  xlab('Proportion of aneuploids') +
  theme(axis.title.x = element_text(h=2.4))+
  theme(axis.text= element_text(face="bold",color="gray50", size=12))+
  theme(axis.text.y= element_text(color=my_colors))+
  ggtitle("Proportion of aneuploids per environment")+
  theme(plot.title = element_text(face="bold",h=0.5))
orplot <- orplot +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
orplot <- orplot + theme(legend.position = "bottom")
png("../aneuploidy/plots/propplot_peter_nonMD.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
orplot
dev.off()

# Create a ploidy matrix
peter <- sctab[sctab$reference == "Peter_et_al_2018",]
peter <- peter[peter$monosporic_derivative == "No",]
attach(peter)

ploidymatrix <- as.matrix(table(ploidy,aneuploidy_binary), nrow=7, header = TRUE)
ploidymatrix <- ploidymatrix[1:5,]

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
  else if(pcat[row]==7){
    pcat[row]<- "Heptaploid"
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
  xlab("Proportion of aneuploids") +
  theme(axis.title.x = element_text(h=0.5,size=16))+
  theme(axis.text= element_text(face="bold",color="gray50", size=16))+
  ggtitle("Proportion of aneuploids per ploidy")+
  theme(plot.title = element_text(face="bold",h=0.5))
p_plot <- p_plot +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
p_plot <- p_plot + theme(legend.position = "bottom")
png("../aneuploidy/plots/ploidy_propplot_peter_nonMD.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
p_plot
dev.off()


# Create a monosporic derivatives matrix
diploids <- sctab[sctab$ploidy == 2 | sctab$ploidy == "Unknown",]
attach(diploids)
MDmatrix <- as.matrix(table(monosporic_derivative,aneuploidy_binary), nrow=2, header = TRUE)

# add a column of totals/proportions to the matrix
for(i in 1:nrow(MDmatrix)){
  if(i==1){
    total_MD <- c(sum(MDmatrix[i,]))
    proportion_MD <- c((MDmatrix[i,2])/(total_MD[i]))
  }
  else{
    total_MD <- c(total_MD,sum(MDmatrix[i,]))
    proportion_MD <- c(proportion_MD,(MDmatrix[i,2])/(total_MD[i]))
  }
}
MDmatrix <- cbind(MDmatrix,total_MD)
MDmatrix <- cbind(MDmatrix, proportion_MD)

# add columns with the confidence intervals for proportions (binomial test) 
for(i in 1:nrow(MDmatrix)){
  if(i==1){
    MD_CILow <- c(binom.test(MDmatrix[i,2],MDmatrix[i,3])$conf.int[1])
    MD_CIHigh <- c(binom.test(MDmatrix[i,2],MDmatrix[i,3])$conf.int[2])
  }
  else{
    MD_CILow <- c(MD_CILow,binom.test(MDmatrix[i,2],MDmatrix[i,3])$conf.int[1])
    MD_CIHigh <- c(MD_CIHigh,binom.test(MDmatrix[i,2],MDmatrix[i,3])$conf.int[2])
  }
}
MDmatrix <- cbind(MDmatrix, MD_CILow)
MDmatrix <- cbind(MDmatrix, MD_CIHigh)

# Set up variables for ploidy proportion plot
MD_factors = c()
MD_colors = c()
gray_scale = c()
MDOdds = c()
MDCILow = c()
MDCIHigh = c()
MDLabels = c()
MDcat=row.names(MDmatrix)

# run a Fisher's Exact Test on each row, and prepare lists to create the odds ratio plot
for(row in 1:nrow(MDmatrix)){
  if(MDcat[row]=="No"){
    MDcat[row] <- "Natural strains"
    MDLabels <- append(MDLabels,paste(MDcat[row],", N=",MDmatrix[row,3],sep=""))
  }
  else if(MDcat[row]=="Yes"){
    MDcat[row] <- "Monosporic Derivatives"
    MDLabels <- append(MDLabels,paste(MDcat[row],", N=",MDmatrix[row,3],sep=""))
  }
}

MD_df <- data.frame(
  yAxis = length(MDLabels[1:2]):1,
  MDmatrix[,4], MDmatrix[,5], MDmatrix[,6])

# create the odds ratio plot
MD_plot <- ggplot(MD_df, aes(x= MDmatrix[,4], y = yAxis))
MD_plot <- MD_plot + 
  geom_vline(aes(xintercept = (sum(MDmatrix[,2]))/(sum(MDmatrix[,3]))), size = .5, linetype = "dashed") +
  geom_errorbarh(aes(xmax = MDmatrix[,6], xmin = MDmatrix[,5]), size = 2, height = 0.5,color = "gray50") +
  geom_point(size = 5.0) +
  theme_bw() +
  scale_y_continuous(breaks = length(MDLabels):1, labels = MDLabels) +
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) +
  ylab("") +
  xlab("Proportion of aneuploids") +
  theme(axis.title.x = element_text(h=0.5,size=16))+
  theme(axis.text= element_text(face="bold",color="gray50", size=16))+
  ggtitle("Proportion of aneuploids per monosporic status (any ploidy; Peter only)")+
  theme(plot.title = element_text(face="bold",h=0.5))
MD_plot <- MD_plot +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
MD_plot <- MD_plot + theme(legend.position = "bottom")
png("../aneuploidy/plots/MD_propplot_all_inc_unk.png", width = 12, height = 4, units = 'in', res = 300, bg='transparent')
MD_plot
dev.off()


### Plot heterozygosity histograms on the entire Peter data set
het_hist <- ggplot(peter, aes(x=heterozygosity))
het_hist <- het_hist + geom_histogram(aes(fill=aneuploidy_binary),
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
png("../aneuploidy/plots/het_propplot_peter.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
het_hist
dev.off()



diploids <- peter[peter$ploidy == 2,]

het_hist <- ggplot(diploids, aes(x=heterozygosity))
het_hist <- het_hist + geom_histogram(aes(fill=aneuploidy_binary),
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
png("../aneuploidy/plots/2n_het_propplot_peter.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
het_hist
dev.off()

