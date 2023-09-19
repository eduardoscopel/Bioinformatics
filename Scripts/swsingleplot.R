slidingwindow <- function(files){
	files <- list.files(pattern = "coverage$")
	c = 0
	k = 0
	for (i in 1:length(files)){
		depthtab <- read.table(files[i], col.names = c("Chromosome", "Position", "Depth"))
		cstart <- 0
                cend <- 0
		a = sub(".coverage.*","",files[i])
		b = sub(".bam.*","",a)
		d = sub(".sorted.*","",b)
		f = sub(".merged.*","",d)
		g = sub(".*Depth_","",f)
                pdf(paste0(g,"_SW_SP2.pdf"))
		med = 2*median(depthtab$Depth, na.rm=TRUE)
		a <- levels(depthtab$Chromosome)
		a[5:8] <- levels(depthtab$Chromosome)[7:10]
		a[9] <- levels(depthtab$Chromosome)[5]
		a[10:16] <- levels(depthtab$Chromosome)[11:17]
		a[17] <- levels(depthtab$Chromosome)[6]
		op <- par(mar=c(1,1,4,1)+0.1, oma = c(5,4,4,1)+0.1, mfcol=c(1,16))
                for(i in a){
                        if(grepl("chrM", i)){next}
			chrm <- depthtab[depthtab$Chromosome == i,2:3]
                        attach(chrm)
                        W <- 1:as.integer(max(Position)/1000)*1000
                        W <- append(W,max(W)+1000)
                        meanDepth <- 0
                        for(j in 1:length(W)) {meanDepth[j] <- mean(Depth[Position>(W[j]-1000)&Position<=W[j]])}
			par(mar=c(0,0,2,0))
			if(c%%4==0){color="black"}else if(c%%4==1){color="orange"}else if(c%%4==2){color="purple"}else if(c%%4==3){color="gray"}
                        plot(W[1:max(W)], meanDepth[1:max(W)],pch=19,cex=0.2, ylim=c(0,med), axes = FALSE, col= color)
			title(main = paste0(sub(".*chr", "", i), "\n", median(Depth, na.rm = TRUE), "x"), cex.main = 0.7)
			k = k + max(W)
			axis(side = 1, at = max(W), labels = k/1000, cex.axis = 0.6, pos=0)
			if(c==0){axis(side = 2, at = seq(0,med,(med/5)),pos=0)}
                        #abline(h=mean(Depth),col="red")
                        abline(h=median(Depth, na.rm = TRUE), col="blue")
			abline(h=0, col = "black")
			#if(grepl("chrIII",i)) {cent <- c(114385, 114501)}
                        #else if(grepl("chrII",i)){cent <- c(238207, 238323)}
                        #else if(grepl("chrIV",i)){cent <- c(449711, 449821)}
                        #else if(grepl("chrVIII",i)){cent <- c(105586, 105703)}
                        #else if(grepl("VII",i)){cent <- c(496920, 497038)}
                        #else if(grepl("chrVI",i)){cent <- c(148510,148627)}
                        #else if(grepl("chrIX",i)){cent <- c(355629, 355745)}
                        #else if(grepl("chrXIII",i)){cent <- c(268031, 268149)}
                        #else if(grepl("chrXII",i)){cent <- c(150828, 150947)}
                        #else if(grepl("chrXI",i)){cent <- c(440129, 440246)}
                        #else if(grepl("chrXVI",i)){cent <- c(555957,556073)}
                        #else if(grepl("chrXIV",i)){cent <- c(628758, 628875)}
                        #else if(grepl("chrXV",i)){cent <- c(326584, 326702)}
                        #else if(grepl("chrI",i)){cent <- c(151465, 151582)}
                        #else if(grepl("chrV",i)){cent <- c(151987, 152104)}
                        #else if(grepl("chrX",i)){cent <- c(436307, 436425)}
                        #abline(v=cent[1],col="blue")
                        #abline(v=cent[2],col="blue")
                        detach(chrm)
			c=c+1
		}
		title(main = g, xlab="genomic position (kbp)", ylab="mean depth", outer = TRUE, line = 2)
		par(op)
                dev.off()
		c=0
		k=0
	}
	rm(list=ls())
}
