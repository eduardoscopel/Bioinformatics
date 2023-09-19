df <- data.frame()
files <- list.files(pattern="_ratio.txt$")
for(i in 1:length(files)){
	cn <- read.table(files[i], header = T)
	strain = sub("_CF.bam_ratio.txt", "", files[i]) 
	for(j in 1:16){
		chr <- cn[cn$Chromosome == j,]
		df[j,1]<- j
		df[j,2]<- median(chr$Ratio)
		df[j,3]<- median(chr$MedianRatio)
	}
	write.table(df, file=paste0(strain, "_summary.txt"),sep = " ", row.names = FALSE, col.names = c("Chromosome", "Mean_Ratio", "Mean_Median_Ratio"))
}
