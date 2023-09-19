files <- list.files(pattern=".bam_ratio.txt$")
for(i in 1:length(files)){
	strainid <- sub(".bam_ratio.txt.*","",files[i])
	cn <- read.table(files[i], header=T)
	ratio <- c()
	medianratio <- c()
	copynumber <- c()
	for(j in c('R',1:7)){
		chr <- cn[cn$Chromosome ==j,]
		ratio[j] <- median(chr$Ratio)
		medianratio[j] <- median(chr$MedianRatio)
		copynumber[j] <- median(chr$CopyNumber)
		}
	ploidytable <- matrix(c(c('R',1:7),ratio, medianratio,copynumber), nrow=8,ncol=4)
	write.table(ploidytable, file=paste0(strainid, ".ploidy.tsv"), sep ="\t", row.names=F, col.names=c("Chromosome", "Ratio", "MedianRatio", "CopyNumber"))
	}
