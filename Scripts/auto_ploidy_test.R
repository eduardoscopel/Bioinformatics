library(epitools)
sctab <- fread("sc_table.tsv",header=TRUE)
sctab <- sctab[sctab$ecology_category != "Unknown",]
sctab <- sctab[sctab$ploidy != "LOWCOV",]
sctab <- sctab[sctab$aneuploidy != "LOWCOV",]
sctab <- sctab[sctab$aneuploidy_binary != "NA",]
attach(sctab)
scmatrix <- as.matrix(table(ecology_category,aneuploidy_binary), nrow=21, header = TRUE)
for(i in 1:nrow(scmatrix)){
  if(i==1){Total <- c(sum(scmatrix[i,]))}
  else{Total <- c(Total,sum(scmatrix[i,]))}
}  
scmatrix <- cbind(scmatrix,Total)
scmatrix <- scmatrix[scmatrix[,"Total"] >20,]
scmatrix
chisq.test(scmatrix[,1:2])
ecocat=row.names(scmatrix)
for(row in 1:nrow(scmatrix)){
  if(row==21){break}
  print(ecocat[row])
  print((matrix(c(scmatrix[row,2],
                  sum(scmatrix[,2])-scmatrix[row,2],
                  scmatrix[row,1],
                  sum(scmatrix[,1])-scmatrix[row,1]),
                nrow=2)))
  print(fisher.test(matrix(c(scmatrix[row,2],
                             sum(scmatrix[,2])-scmatrix[row,2],
                             scmatrix[row,1],
                             sum(scmatrix[,1])-scmatrix[row,1]),
                           nrow=2)))
}

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
