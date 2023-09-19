library(epitools)
sctab <- fread("sc_table.tsv",header=TRUE)
sctab <- sctab[sctab$ecology_category != "Unknown",]
sctab <- sctab[sctab$ploidy != "LOWCOV",]
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
  print(fisher.test(matrix(c(scmatrix[row,2],sum(scmatrix[,2]),scmatrix[row,1]+scmatrix[row,2],sum(scmatrix[,3])),nrow=2)))
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
fisher.test(mdtest[,2:3])

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

input2 = c("
        Euploids Aneuploids Total
Group_1      52         56   108
Group_2      188        138   326
Group_3      119         87   206
Group_4      248         58   306
Group_5      171         83   254
Group_6       19          14    33
Group_7      175         40   215
Group_8      157         16   173
")
matrix2 = as.matrix(read.table(textConnection(input2), header=TRUE, row.names=1))
matrix2
chisq.test(matrix2)
ecocat2=row.names(matrix2)
for(row in 1:nrow(matrix2)){
  if(row==21){break}
  print(ecocat2[row])
  print(fisher.test(matrix(c(matrix2[row,2],492,matrix2[row,1]+matrix2[row,2],1621),nrow=2)))
}

input3=("
        Euploids Aneuploids Total
Domesticated      607         339   946
Clinical      171         83   254
Wild       351          70    421
")
matrix3 = as.matrix(read.table(textConnection(input3), header=TRUE, row.names=1))
matrix3
chisq.test(matrix3)
ecocat3=row.names(matrix3)
for(row in 1:nrow(matrix3)){
  if(row==21){break}
  print(ecocat3[row])
  print(fisher.test(matrix(c(matrix3[row,2],492,matrix3[row,1]+matrix3[row,2],1621),nrow=2)))
}
