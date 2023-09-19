library(epitools)
library(ggplot2)
library(data.table)
# Read table and subset dataframe excluding unknown categories, and LOWCOV strains
sctab <- fread("sc_table.tsv",header=TRUE)
sctab <- sctab[sctab$ecology_category != "Unknown",]
sctab <- sctab[sctab$ploidy != "LOWCOV",]
sctab <- sctab[sctab$aneuploidy != "LOWCOV",]
sctab <- sctab[sctab$aneuploidy_binary != "NA",]
attach(sctab)
# convert df to matrix
scmatrix <- as.matrix(table(ecology_category,aneuploidy_binary), nrow=21, header = TRUE)
# add a column of totals to the matrix
for(i in 1:nrow(scmatrix)){
  if(i==1){Total <- c(sum(scmatrix[i,]))}
  else{Total <- c(Total,sum(scmatrix[i,]))}
}  
scmatrix <- cbind(scmatrix,Total)

for(i in 1:nrow(scmatrix)){
  if(i==1){proportion <- c((scmatrix[i,2])/(scmatrix[i,3]))}
  else{proportion <- c(proportion,(scmatrix[i,2])/(scmatrix[i,3]))}
}  

scmatrix <- cbind(scmatrix, proportion)
scmatrix <- scmatrix[scmatrix[,"Total"] >20,]
scmatrix
# run a chi-square test on the whole matrix to check for association between variables
chisq.test(scmatrix[,1:2])
my_factors = c()
my_colors = c()
gray_scale = c()
boxOdds = c()
boxCILow = c()
boxCIHigh = c()
boxLabels = c()
scmatrix <- rbind(scmatrix[1:3,],
                       scmatrix[5:7,],
                       scmatrix[11:13,],
                       scmatrix[15,],
                       scmatrix[4,],
                       scmatrix[8:10,],
                       scmatrix[15,])
rownames(scmatrix)[10] <- "Wine"
rownames(scmatrix)[11] <- "Clinical"
rownames(scmatrix)[15] <- "Tree"
ecocat=row.names(scmatrix)
# run a Fisher's Exact Test on each row, and prepare lists to create the odds ratio plot
for(row in 1:nrow(scmatrix)){
  if(row==21){break}
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
    else if(is.element(ecocat[row],c("Flower", "Fruit", "Insect", "Tree"))){
      my_factors <- append(my_factors, "Wild")
      my_colors <- append(my_colors, "#316857")
    }
    else{
      my_factors <- append(my_factors, "Domesticated")
      my_colors <- append(my_colors, "#225EA8")
    }
    print(ftest)
}
# prepare data frame with values from fisher test and ecological categories
df <- data.frame(
  yAxis = length(boxLabels):1,
  boxOdds, boxCILow, boxCIHigh)
# create the odds ratio plot
orplot <- ggplot(df, aes(x= boxOdds, y = yAxis))
orplot <- orplot + 
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 2, height = 0.5,color = "gray50") +
  geom_point(color = my_colors, size = 5.0) +
  scale_color_manual(values=my_colors) +
  theme_bw() +
  scale_y_continuous(breaks = length(boxLabels):1, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,6,1)) +
  coord_trans(x = "sqrt") +
  ylab("") +
  xlab('* Fisher’s Exact test, p-value < 0.01                   Odds ratio (sqrt scale)') +
  theme(axis.title.x = element_text(h=2.4))+
  theme(axis.text= element_text(face="bold",color="gray50", size=12))+
  theme(axis.text.y= element_text(color=my_colors))+
  ggtitle("Association between ecology and aneuploidy in yeast")+
  theme(plot.title = element_text(face="bold",h=1.5))
orplot <- orplot +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
orplot <- orplot + theme(legend.position = "bottom")
png("odds_ratio.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
orplot
dev.off()





# subseting for non-monosporic derivatives
mdtest <- as.matrix(table(monosporic_derivative,aneuploidy_binary), nrow=2,header=TRUE)
mdtest
for(i in 1:nrow(mdtest)){
  if(i==1){total <- c(sum(mdtest[i,]))}
  else{total <- c(total,sum(mdtest[i,]))}
}
mdtest <- cbind(mdtest,total)
mdtest
fisher.test(mdtest[,1:2])

nonmd <- sctab[sctab$monosporic_derivative=="No",]
attach(nonmd)
mdmatrix <- as.matrix(table(ecology_category,aneuploidy_binary), nrow=21, header = TRUE)
mdmatrix
for(i in 1:nrow(mdmatrix)){
  if(i==1){Total <- c(sum(mdmatrix[i,]))}
  else{Total <- c(Total,sum(mdmatrix[i,]))}
}  
mdmatrix <- cbind(mdmatrix,Total)
mdmatrix
mdmatrix <- mdmatrix[mdmatrix[,"Total"] > 20,]
mdmatrix
mdmatrix <- mdmatrix[mdmatrix[,"Total"] > 50,]
mdmatrix
count=0
for(eco in unique(nonmd$ecology_category)){
  eco_cat <- nonmd[nonmd$ecology_category == eco,]
  if(nrow(eco_cat) < 50){next}
  else{
    if(count == 0){
      rand <- eco_cat[sample(nrow(eco_cat),50),]
      count = count+1
    }
    else{
      rand <- rbind(rand, eco_cat[sample(nrow(eco_cat),50),]) 
    }
  }
}
rand
attach(rand)
rand_matrix <- as.matrix(table(ecology_category,aneuploidy_binary),nrow=9,header=TRUE)
rand_matrix
for(i in 1:nrow(rand_matrix)){
  if(i==1){Total <- c(sum(rand_matrix[i,]))}
  else{Total <- c(Total,sum(rand_matrix[i,]))}
}  
rand_matrix <- cbind(rand_matrix,Total)
rand_matrix
chisq.test(rand_matrix[,1:2])
ecocat2=row.names(rand_matrix)
for(row in 1:nrow(rand_matrix)){
  if(row==21){break}
  print(ecocat2[row])
  print(matrix(c(rand_matrix[row,2],
                             sum(rand_matrix[,2])-rand_matrix[row,2],
                             rand_matrix[row,1],
                             sum(rand_matrix[,1])-rand_matrix[row,1]),
                           nrow=2))
  print(fisher.test(matrix(c(rand_matrix[row,2],
                             sum(rand_matrix[,2])-rand_matrix[row,2],
                             rand_matrix[row,1],
                             sum(rand_matrix[,1])-rand_matrix[row,1]),
                           nrow=2)))
}
