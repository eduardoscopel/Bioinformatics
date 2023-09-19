slidingwindow <- function(files){

files <- list.files(pattern = "coverage$")

for (i in 1:length(files)){
        depthtab <- read.table(files[i], col.names = c("Chromosome", "Position", "Depth"))

                cstart <- 0
                cend <- 0

                pdf(paste0(files[i],"SW.pdf"))

                for(i in levels(depthtab$Chromosome)){

                        chrm <- depthtab[depthtab$Chromosome == i,2:3]
                        attach(chrm)

                        W <- 1:as.integer(max(Position)/1000)*1000
                        W <- append(W,max(W)+1000)
                        meanDepth <- 0
                        for(j in 1:length(W)) {meanDepth[j] <- mean(Depth[Position>(W[j]-1000)&Position<=W[j]])}
                        par(xpd = F, mar = par()$mar + c(0,0,0,7))
                        plot(W[1:max(W)], meanDepth[1:max(W)],pch=19,main= i,sub="1 kb non overlapping", xlab="position",ylab="mean depth", xaxt="n",ylim=c(0,1000))
                        axis(1,xaxp=c(0,max(W),10))
                        abline(h=mean(Depth),col="red")
                        abline(h=median(Depth), col="yellow")
                        if(grepl("chrIII",i)) {
                                cstart <- 114385
                                cend <- 114501
                        }else if(grepl("chrII",i)){
                                cstart <- 238207
                                cend <- 238323
                        }else if(grepl("chrIV",i)){
                                cstart <- 449711
                                cend <- 449821
                        }else if(grepl("chrVIII",i)){
                                cstart <- 105586
                                cend <- 105703
                        }else if(grepl("VII",i)){
                                cstart <- 496920
                                cend <- 497038
                        }else if(grepl("chrVI",i)){
                                cstart <- 148510
                                cend <- 148627
                        }else if(grepl("chrIX",i)){
                                cstart <- 355629
                                cend <- 355745
                        }else if(grepl("chrXIII",i)){
                                cstart <- 268031
                                cend <- 268149
                        }else if(grepl("chrXII",i)){
                                cstart <- 150828
                                cend <- 150947
                        }else if(grepl("chrXI",i)){
                                cstart <- 440129
                                cend <- 440246
                        }else if(grepl("chrXVI",i)){
                                cstart <- 555957
                                cend <- 556073
                        }else if(grepl("chrXIV",i)){
                                cstart <- 628758
                                cend <- 628875
                        }else if(grepl("chrXV",i)){
                                cstart <- 326584
                                cend <- 326702
                        }else if(grepl("chrI",i)){
                                cstart <- 151465
                                cend <- 151582
                        }else if(grepl("chrV",i)){
                                cstart <- 151987
                                cend <- 152104
                        }else if(grepl("chrX",i)){
                                cstart <- 436307
                                cend <- 436425
                        }
                        abline(v=cstart,col="blue")
                        abline(v=cend,col="blue")
                        legend(x="topright", inset=c(-0.4,0),legend= c(paste0("Mean = ", round(mean(Depth))),paste0("Median = ",median(Depth))),lty=c(1,1),col=c("red","yellow"),xpd=T)
                        par(mar=c(5,4,4,2)+0.1, xpd=FALSE)

                        detach(chrm)

                }
                dev.off()
        }
}
