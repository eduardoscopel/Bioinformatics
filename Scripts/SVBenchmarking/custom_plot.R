files <- list.files()
for(file in 1:length(files)){
files[file]<-sub(".pdf.*","",files[file])
}
strains<-files
chromosomes <- c(1:8,"R")
strains_range <- c(1:length(strains))
xrange <- range(c(1:length(chromosomes))
yrange <- range(strains_range)
xtick <- chromosomes
plot(xrange,yrange, type ="n", xlab="Chromosome", ylab ="Strains", xaxt ="n",yaxt="n")
axis(side=1,at=c(1:9),labels=xtick)
axis(side=2,at=c(1:24),labels=F)
text(y=c(1:24), par("usr")[3]+0.5,labels=strains,cex=0.5,srt=0,pos=2,xpd=TRUE)

