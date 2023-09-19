files <- list.files(pattern=".copynumber.called$")
for (i in 1:length(files)){
 cn<-read.table(files[i],header=TRUE)
 CNA.object<- CNA(genomdat=cn[,7], chrom=cn[,1], maploc=cn[,2], data.type='logratio', sampleid=sub(".copynumber.called.*","",files[i]))
 CNA.smoothed <- smooth.CNA(CNA.object)
 sdundo.CNA.object <- segment(CNA.smoothed, undo.splits = "sdundo", undo.SD=3, verbose=1)
 a = sub(".copynumber.called.*","",files[i])
 pdf(paste0(a,".pdf"))
 plot(sdundo.CNA.object, plot.type="s", ylim=c(-2,2))
 dev.off()
}

files <- list.files(pattern=".copynumber.called$")
for (i in 1:length(files)){
 strainid<-sub(".copynumber.called.*","",files[i])
 cn<-read.table(files[i],header=TRUE)
 meanratio<-c()
 aneuploidy<-c()
 absmean<-c()
 for(j in 1:18){
  chr<-cn[cn$chrom ==j,]
  if(j==17){
   meanratio[j]<-min(meanratio[1:16])
   absmean[j]<-min(absmean[1:16])
   aneuploidy[j]<-"N/A"
  }
  else if(j==18){
   meanratio[j]<-max(meanratio[1:16])
   absmean[j]<-max(absmean[1:16])
   aneuploidy[j]<-"N/A"   
  } 
  else{
   meanratio[j] <- mean(chr[,7])
#  for(j in 1:length(chr[,4])){
#  normvec[j]<-chr[j,4]*chr[j,7]
#   }        
#  meanratio[i]<-sum(normvec)/sum(chr[,5])
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
 }
 ratiotable <- matrix(c(c(1:16,"Min","Max"),meanratio,absmean,aneuploidy), nrow =18, ncol=4)
 write.table(ratiotable,file=paste0(strainid,".raw.tsv"),sep="\t",row.names=F,col.names=c("Chromosome","MeanLogRatio","MeanAbsRatio","Aneuploidy"))
}

