library(seqinr)
library(ggplot2)
library(ggseqlogo)


SRS1 <- read.fasta("~/Documents/Papers_Scopel/BrandiWebLogo/fasta/test/SRS1.fasta",
          seqtype = "AA", as.string = TRUE)
cyp51ASRS1list <- c()
for(i in seq(1,length(cyp51ASRS1))){
  cyp51ASRS1list[i] <- cyp51ASRS1[[i]][1]
}


cs1 = make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                      groups=c('gr1', 'gr1', 'gr2', 'gr2'), 
                      cols=c('purple', 'purple', 'blue', 'blue'))

ggseqlogo(cyp51ASRS1list,method = "prob")
function(ggseqlogo)
  
         