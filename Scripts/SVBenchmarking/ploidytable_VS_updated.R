files <- list.files(pattern=".raw.tsv$")
for(i in 1:length(files)){
	strainid <- sub(".raw.tsv.*","",files[i])
	cn <- read.table(files[i],header=TRUE)
	cn <- cn[1:8,]
	rawmean <- mean(cn$MeanAbsRatio)
	euploids <- cn[cn$MeanAbsRatio/rawmean < 1.2,]
	euploids <- euploids[euploids$MeanAbsRatio/rawmean > 0.8,]
	adjmean <- mean(euploids$MeanAbsRatio)
	aneuploidy<-c()
	meanratio<-c()
	for(j in c('R',1:7)){
		#meanofothers <- mean(cn[cn$Chromosome != j,]$MeanAbsRatio)
		if(cn[cn$Chromosome == j,]$MeanAbsRatio<adjmean){
			meanratio[j] <- (adjmean/cn[cn$Chromosome == j,]$MeanAbsRatio)-1
			if((adjmean/cn[cn$Chromosome == j,]$MeanAbsRatio)-1 > 0.25){aneuploidy[j]<-"-1"}
			else{aneuploidy[j]<-"NA"}
			}
		else{
			meanratio[j] <- (cn[cn$Chromosome == j,]$MeanAbsRatio/adjmean)-1
			if((cn[cn$Chromosome == j,]$MeanAbsRatio/adjmean)-1 > 0.25 & (cn[cn$Chromosome == j,]$MeanAbsRatio/adjmean)-1 <= 0.75){aneuploidy[j]<-"+1"}
			else if((cn[cn$Chromosome == j,]$MeanAbsRatio/adjmean)-1 > 0.75 & (cn[cn$Chromosome == j,]$MeanAbsRatio/adjmean)-1 <= 1.25){aneuploidy[j]<-"+2"}
			else if((cn[cn$Chromosome == j,]$MeanAbsRatio/adjmean)-1 > 1.25 & (cn[cn$Chromosome == j,]$MeanAbsRatio/adjmean)-1 <= 1.75){aneuploidy[j]<-"+3"}
			else if((cn[cn$Chromosome == j,]$MeanAbsRatio/adjmean)-1 > 1.75){aneuploidy[j]<-"+n"}
			else{aneuploidy[j]<-"NA"}
			}
		}
		ploidytable <- matrix(c(c('R',1:7),meanratio, aneuploidy), nrow =8, ncol=3)
		write.table(ploidytable,file=paste0(strainid,".adjusted.tsv"),sep="\t",row.names=F,col.names=c("Chromosome","meanratio","Aneuploidy"))
	}

