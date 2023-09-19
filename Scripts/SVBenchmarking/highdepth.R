highdepth <- function(file){

	depthtab <- read.table(file, col.names = c("Chromosome", "Position", "Depth"))

	for(i in levels(depthtab$Chromosome)){

		chrm <- depthtab[depthtab$Chromosome == i,2:3]
		attach(chrm)

		W <- 1:as.integer(max(Position)/1000)*1000
		W <- append(W,max(W)+1000)
		meanDepth <- 0
		for(j in 1:length(W)) {meanDepth[j] <- mean(Depth[Position>(W[j]-1000)&Position<=W[j]],na.rm=TRUE,nan.rm=TRUE)} 
		a <- cbind(W,meanDepth)
		write.table(a[meanDepth>200,,drop=FALSE], sep=",",paste0(file,".csv"), append=TRUE,col.names=c(paste(i,"W"),"meanDepth"),row.names=FALSE)
		detach(chrm)

		
	}
}
