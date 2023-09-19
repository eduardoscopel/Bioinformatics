library(seqinr)
library(data.table)

### Read tables and match
petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt",header=TRUE)
complete_ssd1 <- read.fasta("Documents/fasta/ssd1/ssd1-1010-prot.fa")

for(seq in complete_ssd1){
  for(i in 1:nrow(petertab)){
    if(attr(seq, "name") == petertab$basename[i]){
      if(seq[1190] == "s" & seq[1196] == "a"){
        petertab$ssd1geno[i] <- "AA"
      }
      else if(seq[1190] == "g" & seq[1196] == "p"){
        petertab$ssd1geno[i] <- "BB"
      }
      else{
        petertab$ssd1geno[i] <- "AB"
      }   
    }
  }
}


for(seq in complete_ssd1){
  for(i in 1:nrow(petertab)){
    if(attr(seq, "name") == petertab$basename[i]){
      if("x" %in% seq){
        petertab$ssd1het[i] <- "het"
      }
      else{
        petertab$ssd1het[i] <- "hom"
      }
    }
  }
}


petertab[petertab$Aneuploidies == "aneu;-1*1;",]
petertab$chr1only <- ifelse(petertab$Aneuploidy_type == "1only", "TRUE", "FALSE")
write.table(petertab,"/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt", row.names = FALSE, sep = "\t")
