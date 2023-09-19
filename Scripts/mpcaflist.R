library(data.table)
data <- fread("natdip.txt", header=FALSE)
for(i in 1:nrow(data)){
  if(data[i,1] == "peter"){
    data[i,2] <- gsub("-","_",data[i,2])
    data[i,2] <- gsub("\\.","_",data[i,2])
    }
}

df <- fread("filelist.txt",header=FALSE)
df
for(i in 1:nrow(df)){
  if(grepl("peter",df[i])){
    df[i] <- paste(df[i],".sorted.snps.vcf.gz",sep="")
  }
  else{df[i] <- paste(df[i],".snps.vcf.gz",sep="")}
}
write.table(df,"filelist2.txt")
