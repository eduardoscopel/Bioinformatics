for file in *chr*.bin
do  strain=${file%%-chr*.bin}
/scratch/es47540/apps/nQuire/nQuire lrdmodel ${strain}.bin > ${strain}.txt
done

for file in *chr*.bin
do
strain=${file%%-chr*.bin}
/scratch/es47540/apps/nQuire/nQuire lrdmodel $file | awk 'NR==2 {print$0}' >> ${strain}.txt
done






nquireplot <- function(files){
	files <- list.files(pattern=".txt")
#	header <- c("Chromosome", "FreeModel", "Diploid", "Triploid", "Tetraploid", "deltaDiploid", "deltaTriploid", "deltaTetraploid")
	for(i in 1:length(files)){
		cn <- read.table(files[i],header=T)
		strain = sub(".txt.*","",files[i])
		cn$file <- sort(as.numeric(cn$file))
		xrange <- range(cn$file)
		yrange <- range(cn$d_dip,cn$d_tri, cn$d_tet)
		pdf(paste0(strain, ".pdf"))
		plot(xrange, yrange, type = "n", xlab = "Chromosome", ylab = "delta Log-Likelihood",xaxt="n",bty="n")
		xtick<-c("Overall",1:7,"R")
		axis(side=1,at=c(1:9),labels=xtick)
#		text(x=c(1:9),par("usr")[3.5],labels=xtick,pos=1,xpd=T)
		lines(cn$file, cn$d_dip, type="o",pch=16, cex=0.6, lwd = 1, col= "darkorchid")
		lines(cn$file, cn$d_tri, type="o",pch=16,cex=0.6, lwd = 1, col="darkorange3")
		lines(cn$file, cn$d_tet, type="o",pch=16,cex=0.6, lwd = 1, col="darkgreen")
		title(main=strain)
		legend("topright", xrange[1], pch=16, legend= c("diploid","triploid","tetraploid"),cex=0.6,col=c("darkorchid","darkorange3","darkgreen"), lty=1, bty="n")
		dev.off()
	}
}
nquireplot()
