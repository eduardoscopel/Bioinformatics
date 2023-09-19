nquireplot <- function(files){
	files <- list.files(pattern=".txt")
#	header <- c("Chromosome", "FreeModel", "Diploid", "Triploid", "Tetraploid", "deltaDiploid", "deltaTriploid", "deltaTetraploid")
	for(i in 1:length(files)){
		cn <- read.table(files[i],header=T)
		strain = sub(".txt.*","",files[i])
		cn$Chromosome <- as.numeric(cn$Chromosome)
		xrange <- range(cn$Chromosome)
		yrange <- range(cn$deltaDiploid,cn$deltaTriploid, cn$deltaTetraploid)
		pdf(paste0(strain, ".pdf"))
		plot(xrange, yrange, type = "n", xlab = "Chromosome", ylab = "delta Log-Likelihood",xaxp=c(1,16,15),bty="n")
		lines(cn$Chromosome, cn$deltaDiploid, type="o",pch=16, cex=0.6, lwd = 1, col= "darkorchid")
		lines(cn$Chromosome, cn$deltaTriploid, type="o",pch=16,cex=0.6, lwd = 1, col="darkorange3")
		lines(cn$Chromosome, cn$deltaTetraploid, type="o",pch=16,cex=0.6, lwd = 1, col="darkgreen")
		title(main=strain)
		legend("topright", xrange[1], pch=16, legend= c("diploid","triploid","tetraploid"),cex=0.6,col=c("darkorchid","darkorange3","darkgreen"), lty=1, bty="n")
	}
	dev.off()
}
nquireplot()
