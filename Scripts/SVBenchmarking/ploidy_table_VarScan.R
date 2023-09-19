files <- list.files(pattern=".copynumber.called$")
for (i in 1:length(files)){
 strainid<-sub(".copynumber.called.*","",files[i])
 cn<-read.table(files[i],header=TRUE)
 meanratio<-c()
 aneuploidy<-c()
 absmean<-c()
 for(j in c('R',1:9)){
  chr<-cn[cn$chrom == paste0("chr",j),]
  if(j==8){
   meanratio[j]<-min(meanratio[1:8])
   absmean[j]<-min(absmean[1:8])
   aneuploidy[j]<-"N/A"
  }
  else if(j==9){
   meanratio[j]<-max(meanratio[1:8])
   absmean[j]<-max(absmean[1:8])
   aneuploidy[j]<-"N/A"   
  } 
  else{
   meanratio[j] <- mean(chr[,7])
  }
  absmean[j]<-2**meanratio[j]
  if(absmean[j]<1 & abs(absmean[j]-1)>0.35){
   aneuploidy[j]<-"-1"
   }
  else if(absmean[j]>=1.35 & absmean[j]<1.65){
   aneuploidy[j]<-"+1"
   }
  else if(absmean[j]>=1.65 & absmean[j]<2.35){
   aneuploidy[j]<-"+2"
   }
  else if(absmean[j]>=2.35){
   aneuploidy[j]<-"+n"
   }
  else{
  aneuploidy[j]<-"none"
  }
 }
 ratiotable <- matrix(c(c('R',1:7,'Min','Max'),meanratio,absmean,aneuploidy), nrow =10, ncol=4)
 write.table(ratiotable,file=paste0(strainid,".raw.tsv"),sep="\t",row.names=F,col.names=c("Chromosome","MeanLogRatio","MeanAbsRatio","Aneuploidy"))
}

