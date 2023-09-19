library(data.table)

dmat <- fread("Documents/GitHub/eduardo/aneuploidy/peter/distmatrix/test", header = FALSE)
strains <- c("strains",dmat$V1)
strains
colnames(dmat) <- strains
rownames(dmat)[1] <- dmat$strains[1]
distmatrix <- as.matrix(dmat, rownames = TRUE)
format(distmatrix["STG_2","STG_3"], scientific = FALSE)
format(distmatrix["STG_1","STG_3"], scientific = FALSE)
format(distmatrix["STG_1","STG_2"], scientific = FALSE)

format(distmatrix["MC6","8_3"], scientific = FALSE)

genthreshold1 <- distmatrix["STG_2","STG_3"]
genthreshold2 <- distmatrix["MC6","8_3"]

clrel <- data.frame(strain = rownames(distmatrix))
a<-c()

for(i in 1:nrow(distmatrix)){
  if(which(distmatrix[i,] <= genthreshold) == i){
    a <- append(a, list(c(rownames(distmatrix)[i],NA)))
  }
  else{
    a <- append(a, list(c(rownames(distmatrix)[i], colnames(distmatrix)[which(distmatrix[i,] <= genthreshold)])))
  }
}  

A <- c()
B <- c()
gendist <- c()


for(i in 1:nrow(distmatrix)){
  for(j in 1:ncol(distmatrix)){
    if(is.null(A[i]) | 
       is.na(A[i])){
      A <- append(A, rownames(distmatrix)[i])
    }
    else if(distmatrix[i,j] <= genthreshold &
       rownames(distmatrix)[i] != colnames(distmatrix)[j] &
       j > i){
      B <- append(B, colnames(distmatrix)[j])
      gendist <- append(gendist, distmatrix[i,j])
      #clrelatives$A[i] <- rownames(distmatrix)[i]
      #clrelatives$B[i] <- colnames(distmatrix)[j]
      #clrelatives$gendist[i] <- distmatrix[i,j]
      #print(c(rownames(distmatrix)[i], colnames(distmatrix)[j], distmatrix[i,j]))
    }
  }
}
clrelatives <- data.frame(A, B, gendist)
clrelatives
