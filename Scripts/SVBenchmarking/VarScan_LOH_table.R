
files <- list.files(pattern=".summary.snp$")
for(i in 1:length(files)){
	strainid <- sub(".summary.snp.*","",files[i])
	df<-read.table(files[i],header=T)
	hetoverall <- c()
	homoverall <- c()
	overall <- c()
	if(mean(df$SamplesHet) <= 0.1){overall[i] <- "homozygous"}
		else{overall[i] <- "heterozygous"}
	hetoverall[i] <- mean(df$SamplesHet)
	homoverall[i] <- mean(df$SamplesHom)
	chromosome <- c()
	hetchrom <- c()
	homchrom <- c()
#	df$Chrom <- as.numeric(df$Chrom)
	for(j in levels(df$Chrom)){
		chr<-df[df$Chrom==j,]
		hetchrom[j] <- mean(chr$SamplesHet)
		homchrom[j] <- mean(chr$SamplesHom)
		if(mean(chr$SamplesHet) <= 0.1){chromosome[j] <- "homozygous"}
		else{chromosome[j]<- "heterozygous"}
		}
	LOHtable <- matrix(c("overall", c(1:7,'R'), hetoverall[i],hetchrom, homoverall[i],homchrom, overall[i], chromosome),nrow = 9, ncol = 4)
	write.table(LOHtable, file=paste0(strainid,".LOH.tsv"),sep="\t", row.names=F, col.names=c("Chromosome", "Het", "Hom", "Zygosity"))
	}

