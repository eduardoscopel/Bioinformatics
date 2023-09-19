library(data.table)
library(ggplot2)

df <- fread("~/Documents/Papers_Scopel/GWAS/sc_table_most_recent.txt")
### Almeida
almeida <- df[df$reference == "Almeida_et_al_2015",]
# Coverage analysis
ggplot(almeida) + 
  geom_histogram(aes(avg_cov),binwidth = 30, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(almeida$avg_cov)
table(almeida$avg_cov >= 50)

# Contamination
table(almeida$Contaminated)
table(almeida$Contaminated >= 0.05 & almeida$Contaminated != "No")

# Monosporics
table(almeida$monosporic_derivative)

# Ploidy
table(almeida$ploidyBAF)

unique(df$reference)


### Barbosa 2016
barb16 <- df[df$reference == "Barbosa_et_al_2016",]
# Coverage analysis
ggplot(barb16) + 
  geom_histogram(aes(avg_cov),binwidth = 20, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(barb16$avg_cov)
sd(barb16$avg_cov)
table(barb16$avg_cov >= 50)

# Contamination
table(barb16$Contaminated)
table(barb16$Contaminated >= 0.05 & barb16$Contaminated != "No")

# Monosporics
table(barb16$monosporic_derivative)

# Ploidy
table(barb16$ploidyBAF)


### Barbosa 2018
barb18 <- df[df$reference == "Barbosa_et_al_2018",]
# Coverage analysis
ggplot(barb18) + 
  geom_histogram(aes(avg_cov),binwidth = 20, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(barb18$avg_cov)
sd(barb18$avg_cov)
table(barb18$avg_cov >= 50)
# Contamination
table(barb18$Contaminated)
table(barb18$Contaminated >= 0.05 & barb18$Contaminated != "No")
# Monosporics
table(barb18$monosporic_derivative)
# Ploidy
table(barb18$ploidyBAF)

### Bergstrom 2014
bergstrom <- df[df$reference == "Bergstrom_et_al_2014",]
# Coverage analysis
ggplot(bergstrom) + 
  geom_histogram(aes(avg_cov),binwidth = 20, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(bergstrom$avg_cov)
sd(bergstrom$avg_cov)
table(bergstrom$avg_cov >= 50)
# Contamination
table(bergstrom$Contaminated)
table(bergstrom$Contaminated >= 0.05 & bergstrom$Contaminated != "No")
# Monosporics
table(bergstrom$monosporic_derivative)
# Ploidy
table(bergstrom$ploidyBAF)

### Borneman 2016
borneman <- df[df$reference == "Borneman_et_al_2016",]
# Coverage analysis
ggplot(borneman) + 
  geom_histogram(aes(avg_cov),binwidth = 5, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(borneman$avg_cov)
table(borneman$avg_cov >= 50)
sd(borneman$avg_cov)
# Contamination
table(borneman$Contaminated)
table(borneman$Contaminated >= 0.05 & borneman$Contaminated != "No")
# Monosporics
table(borneman$monosporic_derivative)
# Ploidy
table(borneman$ploidyBAF)


### Duan 2018
duan <- df[df$reference == "Duan_et_al_2018",]
# Coverage analysis
ggplot(duan) + 
  geom_histogram(aes(avg_cov),binwidth = 30, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(duan$avg_cov)
table(duan$avg_cov >= 50)
sd(duan$avg_cov)
# Contamination
table(duan$Contaminated)
table(duan$Contaminated >= 0.05 & duan$Contaminated != "No")
# Monosporics
table(duan$monosporic_derivative)
# Ploidy
table(duan$ploidyBAF)

### Fay 2018
fay <- df[df$reference == "Fay_et_al_2018",]
# Coverage analysis
ggplot(fay) + 
  geom_histogram(aes(avg_cov),binwidth = 30, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(fay$avg_cov)
table(fay$avg_cov >= 50)
sd(fay$avg_cov)
# Contamination
table(fay$Contaminated)
table(fay$Contaminated >= 0.05 & fay$Contaminated != "No")
# Monosporics
table(fay$monosporic_derivative)
# Ploidy
table(fay$ploidyBAF)

### Gallone 2016
gallone <- df[df$reference == "Gallone_et_al_2016",]
# Coverage analysis
ggplot(gallone) + 
  geom_histogram(aes(avg_cov),binwidth = 30, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(gallone$avg_cov)
table(gallone$avg_cov >= 50)
sd(gallone$avg_cov)
# Contamination
table(gallone$Contaminated)
table(gallone$Contaminated >= 0.05 & gallone$Contaminated != "No")
# Monosporics
table(gallone$monosporic_derivative)
# Ploidy
table(gallone$ploidyBAF)


### Gayeviskiy 2016
gay <- df[df$reference == "Gayevskiy_et_al_2016",]
# Coverage analysis
ggplot(gay) + 
  geom_histogram(aes(avg_cov),binwidth = 5, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(gay$avg_cov)
table(gay$avg_cov >= 50)
sd(gay$avg_cov)
# Contamination
table(gay$Contaminated)
table(gay$Contaminated >= 0.05 & gay$Contaminated != "No")
# Monosporics
table(gay$monosporic_derivative, useNA = "ifany")
# Ploidy
table(gay$ploidyBAF)


### Goncalves 2016
gonc <- df[df$reference == "Goncalves_et_al_2016",]
# Coverage analysis
ggplot(gonc) + 
  geom_histogram(aes(avg_cov),binwidth = 5, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(gonc$avg_cov)
table(gonc$avg_cov >= 50)
sd(gonc$avg_cov)
# Contamination
table(gonc$Contaminated)
table(gonc$Contaminated >= 0.05 & gonc$Contaminated != "No")
# Monosporics
table(gonc$monosporic_derivative, useNA = "ifany")
# Ploidy
table(gonc$ploidyBAF)


### Han 2021
han <- df[df$reference == "Han_et_al_2021",]
# Coverage analysis
ggplot(han) + 
  geom_histogram(aes(avg_cov),binwidth = 30, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(han$avg_cov)
table(han$avg_cov >= 50)
sd(han$avg_cov)
# Contamination
table(han$Contaminated)
table(han$Contaminated >= 0.05 & han$Contaminated != "No")
# Monosporics
table(han$monosporic_derivative, useNA = "ifany")
# Ploidy
table(han$ploidyBAF)


### Our data
od <- df[df$reference == "Our_Data",]
# Coverage analysis
ggplot(od) + 
  geom_histogram(aes(avg_cov),binwidth = 10, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(od$avg_cov)
table(od$avg_cov >= 50)
sd(od$avg_cov)
# Contamination
table(od$Contaminated)
table(od$Contaminated >= 0.05 & od$Contaminated != "No")
# Monosporics
table(od$monosporic_derivative, useNA = "ifany")
# Ploidy
table(od$ploidyBAF)

### Peter
peter <- df[df$reference == "Peter_et_al_2018",]
# Coverage analysis
ggplot(peter) + 
  geom_histogram(aes(avg_cov),binwidth = 50, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,1000,100),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,1000,100),expand = c(0,0)) + 
  theme_classic()
summary(peter$avg_cov)
table(peter$avg_cov >= 50)
sd(peter$avg_cov,na.rm = TRUE)
# Contamination
table(peter$Contaminated)
table(peter$Contaminated >= 0.05 & peter$Contaminated != "No")
# Monosporics
table(peter$monosporic_derivative, useNA = "ifany")
# Ploidy
table(peter$ploidyBAF)

### Pontes
pontes <- df[df$reference == "Pontes_et_al_2020",]
# Coverage analysis
ggplot(pontes) + 
  geom_histogram(aes(avg_cov),binwidth = 10, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,25),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(pontes$avg_cov)
table(pontes$avg_cov >= 50)
sd(pontes$avg_cov,na.rm = TRUE)
# Contamination
table(pontes$Contaminated)
table(pontes$Contaminated >= 0.05 & pontes$Contaminated != "No")
# Monosporics
table(pontes$monosporic_derivative, useNA = "ifany")
# Ploidy
table(pontes$ploidyBAF)

### Skelly 2013
skelly <- df[df$reference == "Skelly_et_al_2013",]
# Coverage analysis
ggplot(skelly) + 
  geom_histogram(aes(avg_cov),binwidth = 10, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,25),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(skelly$avg_cov)
table(skelly$avg_cov >= 50)
sd(skelly$avg_cov,na.rm = TRUE)
# Contamination
table(skelly$Contaminated)
table(skelly$Contaminated >= 0.05 & skelly$Contaminated != "No")
# Monosporics
table(skelly$monosporic_derivative, useNA = "ifany")
# Ploidy
table(skelly$ploidyBAF)


### Song 2015
song <- df[df$reference == "Song_et_al_2015",]
# Coverage analysis
ggplot(song) + 
  geom_histogram(aes(avg_cov),binwidth = 20, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,450,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,2),expand = c(0,0)) + 
  theme_classic()
summary(song$avg_cov)
table(song$avg_cov >= 50)
sd(song$avg_cov,na.rm = TRUE)
# Contamination
table(song$Contaminated)
table(song$Contaminated >= 0.05 & song$Contaminated != "No")
# Monosporics
table(song$monosporic_derivative, useNA = "ifany")
# Ploidy
table(song$ploidyBAF)


### Strope 2015
strope <- df[df$reference == "Strope_et_al_2015",]
# Coverage analysis
ggplot(strope) + 
  geom_histogram(aes(avg_cov),binwidth = 20, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,1000,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,5),expand = c(0,0)) + 
  theme_classic()
summary(strope$avg_cov)
table(strope$avg_cov >= 50)
sd(strope$avg_cov,na.rm = TRUE)
# Contamination
table(strope$Contaminated)
table(strope$Contaminated >= 0.05 & strope$Contaminated != "No")
# Monosporics
table(strope$monosporic_derivative, useNA = "ifany")
# Ploidy
table(strope$ploidyBAF)


### Yue 2015
yue <- df[df$reference == "Yue_et_al_2017",]
# Coverage analysis
ggplot(yue) + 
  geom_histogram(aes(avg_cov),binwidth = 50, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,1000,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,1),expand = c(0,0)) + 
  theme_classic()
summary(yue$avg_cov)
table(yue$avg_cov >= 50)
sd(yue$avg_cov,na.rm = TRUE)
# Contamination
table(yue$Contaminated)
table(yue$Contaminated >= 0.05 & yue$Contaminated != "No")
# Monosporics
table(yue$monosporic_derivative, useNA = "ifany")
# Ploidy
table(yue$ploidyBAF)


### Zhu 2016
zhu <- df[df$reference == "Zhu_et_al_2016",]
# Coverage analysis
ggplot(zhu) + 
  geom_histogram(aes(avg_cov),binwidth = 20, fill = "white", col = "black") + 
  scale_x_continuous(breaks = seq(0,1000,50),expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,100,20),expand = c(0,0)) + 
  theme_classic()
summary(zhu$avg_cov)
table(zhu$avg_cov >= 50)
sd(zhu$avg_cov,na.rm = TRUE)
# Contamination
table(zhu$Contaminated)
table(zhu$Contaminated >= 0.05 & zhu$Contaminated != "No")
# Monosporics
table(zhu$monosporic_derivative, useNA = "ifany")
# Ploidy
table(zhu$ploidyBAF)


unique(df$reference)
table(df$aneuploidy_binary)
table(df$ecology_category)
table(df$pop_summary)
table(df$ecology_category, df$aneuploidy_type)
table(df$pop_summary, df$aneuploidy_type, useNA = "ifany")
table(df$ploidyBAF, df$aneuploidy_type)


table(df$aneuploidy_binary)
table(df$ecology_category)
table(df$pop_summary)
table(df$ploidy)
table(df$ecology_category, df$aneuploidy_binary)
table(df$pop_summary, df$aneuploidy_binary)
table(df$ploidy, df$aneuploidy_binary, useNA = "always")


scmatrix <- as.matrix(table(df$ecology_category, df$aneuploidy_binary, useNA = "always"), header = TRUE)
for(i in 1:nrow(scmatrix)){
  if(i==1){
    Total <- c(sum(scmatrix[i,]))
  }
  else{
    Total <- c(Total, sum(scmatrix[i,]))
  }
}
scmatrix<-cbind(scmatrix,Total)
scmatrix <- scmatrix[,c(1,2,4)]
for(i in 1:nrow(scmatrix)){
  if(i==1){proportion <- c((scmatrix[i,2])/(scmatrix[i,3]))}
  else{proportion <- c(proportion,(scmatrix[i,2])/(scmatrix[i,3]))}
}  
scmatrix <- cbind(scmatrix, proportion)
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
my_factors = c()
my_colors = c()
gray_scale = c()
boxOdds = c()
boxCILow = c()
boxCIHigh = c()
boxLabels = c()
scmatrix <- rbind(scmatrix[c(1,2,3,6,7,8,9, 13,15,17,18,21),],
                  scmatrix[4,],
                  scmatrix[c(5,10,11,12, 14,16,19,20),],
                  scmatrix[22,])
rownames(scmatrix)[13] <- "Clinical"
rownames(scmatrix)[22] <- "Other"
ecocat <- row.names(scmatrix)
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
  else if(ecocat[row] == "Other"){
    my_factors <- append(my_factors, "Other")
    my_colors <- append(my_colors, "black")
  }
  else if(is.element(ecocat[row],c("Cocoa","Fish","Flower", "Fruit", "Insect", "Other_plants", "Tree", "Water"))){
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
newdf <- data.frame(
  yAxis = length(boxLabels):1,
  scmatrix[,4], scmatrix[,5], scmatrix[,6])
orplot <- ggplot(newdf, aes(x= scmatrix[,4], y = yAxis))
orplot <- orplot + 
  geom_vline(aes(xintercept = (sum(scmatrix[,2]))/(sum(scmatrix[,3]))), size = .5, linetype = "dashed") +
  #geom_errorbarh(aes(xmax = scmatrix[,6], xmin = scmatrix[,5]), size = 2, height = 0.5,color = "gray50") +
  #geom_point(color = my_colors, size = 5.0) +
  geom_bar(aes(x = newdf[,1], y = yAxis),stat = "identity")+
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
orplot

ggplot(newdf) + 
  geom_bar(aes(x = newdf[,1], y = newdf[,2]), stat = "identity") + 
  geom_hline(aes(yintercept = (sum(scmatrix[,2]))/(sum(scmatrix[,3]))), size = .5, linetype = "dashed") +
  theme_classic() + 
  scale_x_continuous(breaks = length(boxLabels):1, labels = boxLabels, expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0,1,0.1), expand = c(0,0),limits = c(0,1))+
  ylab("Proportion of aneuploids")+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))

source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
install.packages("devtools")
BiocManager::install("multtest")
BiocManager::install("snpStats")
install.packages("/Users/es47540/Applications/LDheatmap/LDheatmap_0.99-8.tar.gz", repos = NULL, type = "source")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)

myY <- read.table("http://zzlab.net/GAPIT/data/mdp_traits.txt", head = TRUE)
myGD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
myGM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
myG=read.table("http://zzlab.net/GAPIT/data/mdp_genotype_test.hmp.txt",head=T)

myGAPIT_MLM <- GAPIT(
  Y=myY[,c(1,3)], #fist column is individual ID, the third columns is days to pollination
  GD=myGD,
  GM=myGM,
  PCA.total=3,
  model=c("GLM", "MLM", "CMLM", "FarmCPU", "Blink"),
  Multiple_analysis=TRUE)
