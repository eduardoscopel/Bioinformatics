makehist <- function(files){

	files <- list.files(pattern = "depth$")

	for (file in 1:length(files)){

		depthtab <- fread(files[file], col.names = c("Chromosome", "Position", "Depth"))
		f = sub(".coverage.*","",files[file])
		g = sub(".*Depth_","",f)
		med <- median(depthtab$Depth)
		print(g)
		print(med)
		pdf(paste0(g,"hist.pdf"))
		par(mfrow = c(1,1))
		hist(depthtab$Depth, 10000, xlim = c(0, med), main = paste0(g, "\n", med/3, "x"), xlab= "Read Depth", ylab = "Frequency")
		abline(v=med/3, col="red",lwd=2)

		title(main=g, outer = TRUE, line = 2)
		dev.off()
	}

			op <- par(mar = c(1,1,4,1) + 0.1, oma = c(5,4,4,1)+0.1, mfrow=c(4,4))

			a <- levels(depthtab$Chromosome)
			a[5:8] <- levels(depthtab$Chromosome)[7:10]
			a[9] <- levels(depthtab$Chromosome)[5]
			a[10:16] <- levels(depthtab$Chromosome)[11:17]
			a[17] <- levels(depthtab$Chromosome)[6]
			for(i in a){
				if(grepl("chrM", i)){next}
				chrm = depthtab[depthtab$Chromosome ==i, 2:3]
				attach(chrm)

				hist(chrm$Depth, 10000, xlim=c(0,med), main = paste0(i, "\n", median(chrm$Depth), "x"),cex.main=0.7)
				abline(v=median(chrm$Depth), col="red", lwd=2)
				#abline(v=mean(chrm$Depth), col="blue", lwd=2)
				#legend(x="topright",c(paste0("Mean = ",as.integer(mean(chrm$Depth))),paste0("Median = ",median(chrm$Depth))),col = c("blue","red"),lwd=c(2,2))
				detach(chrm)
			}
			title(main=g, outer = TRUE, line = 2)
			par(op)
			dev.off()
		}
		rm(list=ls())
	}
