myGD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
myGM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)

X<-myGD[,-1]
taxa<-as.character(myGD[,1])
NQTN = 10
h2 = 0.5
##simulation phyenotype
##-------------------------##
n=nrow(X)
m=ncol(X)
rep=10
rep.power.GLM<-data.frame(matrix(0,rep,6))
rep.FDR.GLM<-data.frame(matrix(0,rep,6))
rep.Power.Alpha.GLM<-data.frame(matrix(0,12,6))

rep.power.MLM<-data.frame(matrix(0,rep,6))
rep.FDR.MLM<-data.frame(matrix(0,rep,6))
rep.Power.Alpha.MLM<-data.frame(matrix(0,12,6))

rep.power.SUPER<-data.frame(matrix(0,rep,6))
rep.FDR.SUPER<-data.frame(matrix(0,rep,6))
rep.Power.Alpha.SUPER<-data.frame(matrix(0,12,6))

rep.power.CMLM<-data.frame(matrix(0,rep,6))
rep.FDR.CMLM<-data.frame(matrix(0,rep,6))
rep.Power.Alpha.CMLM<-data.frame(matrix(0,12,6))

rep.power.MLMM<-data.frame(matrix(0,100,6))
rep.FDR.MLMM<-data.frame(matrix(0,100,6))
rep.Power.Alpha.MLMM<-data.frame(matrix(0,12,6))

rep.power.FarmCPU<-data.frame(matrix(0,100,6))
rep.FDR.FarmCPU<-data.frame(matrix(0,100,6))
rep.Power.Alpha.FarmCPU<-data.frame(matrix(0,12,6))

rep.power.Blink<-data.frame(matrix(0,100,6))
rep.FDR.Blink<-data.frame(matrix(0,100,6))
rep.Power.Alpha.Blink<-data.frame(matrix(0,12,6))

##PCA
##---------------------##

#PCA<-prcomp(X)
#PCVar<-PCA$sdev^2
#myPC<-PCA$x[,1:3]
#m1<-as.data.frame(myPC)

#myCV<-cbind(taxa,m1)
#myCV<-as.data.frame(myCV)

##-----end step 2  for tfam---###
#kcv1<-matrix(1,nrow(myCV),1)
#kcv<-cbind(data.frame(kcv1),myCV)
#write.table(kcv,"pca.txt",row.names = FALSE,col.names = FALSE,sep="\t",quote=FALSE)

for(i in 1:rep){
  mySim=GAPIT.Phenotype.Simulation(GD=myGD,
                                   GM=myGM,
                                   h2=.5,
                                   NQTN=10, 
                                   effectunit =.95,
                                   QTNDist="normal",
                                   category = 0)
  
  myY <- mySim$Y
  max.groups=nrow(myY)
  print(paste("*****************","GWAS by GAPIT...GLM model",i," totle:",rep,sep=""))
  #--------------------------
  myGAPIT_GLM=GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    PCA.total=0,
    file.output=FALSE,
    model="GLM",
    memo="GLM",
    QTN.position=QTN.position,
    threshold.output=0.001,
    iteration.output=TRUE,
  ) 
  
  print(paste("*****************","GWAS by GAPIT...MLM model",i," totle:",rep,sep=""))
  #--------------------------------#
  myGAPIT_MLM=GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    KI=myGAPIT_GLM$KI,
    file.output=FALSE,
    model="MLM",
    memo="MLM",
    QTN.position=QTN.position,
    threshold.output=0.001,
    iteration.output=TRUE,
  ) 
  
  print(paste("*****************","GWAS by GAPIT...SUPER model",i," totle:",rep,sep=""))
  ##--------------------------------#
  myGAPIT_SUPER <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    KI=myGAPIT_GLM$KI,
    #PCA.total=3,
    model="SUPER",
    QTN.position=QTN.position,
    threshold.output=0.001,
    iteration.output=TRUE,
    file.output=FALSE,
  )
  
  print(paste("$$$$$$$$$$$$$$$","GWAS by GAPIT...CMLM model",i," totle:",rep,sep=""))
  #--------------------------------#
  myGAPIT_CMLM=GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    KI=myGAPIT_GLM$KI,
    file.output=FALSE,
    model="CMLM",
    memo="CMLM",
    QTN.position=QTN.position,
    threshold.output=0.001,
    iteration.output=TRUE,
  ) 
  
  print(paste("$$$$$$$$$$$$$$$","GWAS by GAPIT...MLMM model",i," totle:",rep,sep=""))
  #--------------------------------#
  myGAPIT_MLMM=GAPIT(
    Y=myY,
    GD=myGAPIT_GLM$GD,
    GM=myGAPIT_GLM$GM,
    KI=myGAPIT_GLM$KI,
    file.output=FALSE,
    model="MLMM",
    memo="MLMM",
    QTN.position=QTN.position,
    threshold.output=0.001,
    iteration.output=TRUE,
  ) 
  
  print(paste("$$$$$$$$$$$$$$$","GWAS by GAPIT...FarmCPU model",i," totle:",rep,sep=""))
  #--------------------------------#
  myGAPIT_FarmCPU=GAPIT(
    Y=myY,
    GD=myGAPIT_GLM$GD,
    GM=myGAPIT_GLM$GM,
    KI=myGAPIT_GLM$KI,
    file.output=FALSE,
    model="FarmCPU",
    memo="FarmCPU",
    QTN.position=QTN.position,
    threshold.output=0.001,
    iteration.output=TRUE,
  ) 
  print(paste("$$$$$$$$$$$$$$$","GWAS by GAPIT...Blink model",i," totle:",rep,sep=""))
  #--------------------------------#
  myGAPIT_Blink=GAPIT(
    Y=myY,
    GD=myGAPIT_GLM$GD,
    GM=myGAPIT_GLM$GM,
    KI=myGAPIT_GLM$KI,
    file.output=FALSE,
    model="Blink",
    memo="Blink",
    QTN.position=QTN.position,
    threshold.output=0.001,
    iteration.output=TRUE,
  ) 
  
  
  #power #FDR #Power.Alpha
  rep.power.GLM<-rep.power.GLM+myGAPIT_GLM$Power
  rep.FDR.GLM<-rep.FDR.GLM+myGAPIT_GLM$FDR
  rep.Power.Alpha.GLM<-rep.Power.Alpha.GLM+myGAPIT_GLM$Power.Alpha
  
  rep.power.MLM<-rep.power.MLM+myGAPIT_MLM$Power
  rep.FDR.MLM<-rep.FDR.MLM+myGAPIT_MLM$FDR
  rep.Power.Alpha.MLM<-rep.Power.Alpha.MLM+myGAPIT_MLM$Power.Alpha
  
  rep.power.SUPER<-rep.power.SUPER+myGAPIT_SUPER$Power
  rep.FDR.SUPER<-rep.FDR.SUPER+myGAPIT_SUPER$FDR
  rep.Power.Alpha.SUPER<-rep.Power.Alpha.SUPER+myGAPIT_SUPER$Power.Alpha
  
  rep.power.CMLM<-rep.power.CMLM+myGAPIT_CMLM$Power
  rep.FDR.CMLM<-rep.FDR.CMLM+myGAPIT_CMLM$FDR
  rep.Power.Alpha.CMLM<-rep.Power.Alpha.CMLM+myGAPIT_CMLM$Power.Alpha
  
  rep.power.MLMM<-rep.power.MLMM+myGAPIT_MLMM$Power
  rep.FDR.MLMM<-rep.FDR.MLMM+myGAPIT_MLMM$FDR
  rep.Power.Alpha.MLMM<-rep.Power.Alpha.MLMM+myGAPIT_MLMM$Power.Alpha
  
  rep.power.FarmCPU<-rep.power.FarmCPU+myGAPIT_FarmCPU$Power
  rep.FDR.FarmCPU<-rep.FDR.FarmCPU+myGAPIT_FarmCPU$FDR
  rep.Power.Alpha.FarmCPU<-rep.Power.Alpha.FarmCPU+myGAPIT_FarmCPU$Power.Alpha
  
  rep.power.Blink<-rep.power.Blink+myGAPIT_Blink$Power
  rep.FDR.Blink<-rep.FDR.Blink+myGAPIT_Blink$FDR
  rep.Power.Alpha.Blink<-rep.Power.Alpha.Blink+myGAPIT_Blink$Power.Alpha
  
  gc()
}
#mean
rep.power.GLM<-rep.power.GLM/rep
rep.FDR.GLM<-rep.FDR.GLM/rep
rep.Power.Alpha.GLM<-rep.Power.Alpha.GLM/rep

rep.power.MLM<-rep.power.MLM/rep
rep.FDR.MLM<-rep.FDR.MLM/rep
rep.Power.Alpha.MLM<-rep.Power.Alpha.MLM/rep

rep.power.SUPER<-rep.power.SUPER/rep
rep.FDR.SUPER<-rep.FDR.SUPER/rep
rep.Power.Alpha.SUPER<-rep.Power.Alpha.SUPER/rep

rep.power.CMLM<-rep.power.CMLM/rep
rep.FDR.CMLM<-rep.FDR.CMLM/rep
rep.Power.Alpha.CMLM<-rep.Power.Alpha.CMLM/rep

rep.power.MLMM<-rep.power.MLMM/rep
rep.FDR.MLMM<-rep.FDR.MLMM/rep
rep.Power.Alpha.MLMM<-rep.Power.Alpha.MLMM/rep

rep.power.FarmCPU<-rep.power.FarmCPU/rep
rep.FDR.FarmCPU<-rep.FDR.FarmCPU/rep
rep.Power.Alpha.FarmCPU<-rep.Power.Alpha.FarmCPU/rep

rep.power.Blink<-rep.power.Blink/rep
rep.FDR.Blink<-rep.FDR.Blink/rep
rep.Power.Alpha.Blink<-rep.Power.Alpha.Blink/rep

#ouput files power FDR for GLM,MLM,SUPER

myWS=c(1e0,1e3,1e4,1e5,1e6,1e7)
myalpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)

colnames(rep.FDR.GLM)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.GLM)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.GLM)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.MLM)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.MLM)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.MLM)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.SUPER)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.SUPER)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.SUPER)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.CMLM)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.CMLM)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.CMLM)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.MLMM)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.MLMM)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.MLMM)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.FarmCPU)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.FarmCPU)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.FarmCPU)=paste("Power(",myWS,")",sep="")


colnames(rep.FDR.Blink)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.Blink)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.Blink)=paste("Power(",myWS,")",sep="")

write.csv(cbind(rep.FDR.GLM,rep.power.GLM),paste(h2,"_",NQTN,".Power.by.FDR.GLM",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.GLM),paste(h2,"_",NQTN,".Power.by.TypeI.GLM",".csv",sep=""))

write.csv(cbind(rep.FDR.MLM,rep.power.MLM),paste(h2,"_",NQTN,".Power.by.FDR.MLM",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.MLM),paste(h2,"_",NQTN,".Power.by.TypeI.MLM",".csv",sep=""))

write.csv(cbind(rep.FDR.SUPER,rep.power.SUPER),paste(h2,"_",NQTN,".Power.by.FDR.SUPER",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.SUPER),paste(h2,"_",NQTN,".Power.by.TypeI.SUPER",".csv",sep=""))

write.csv(cbind(rep.FDR.CMLM,rep.power.CMLM),paste(h2,"_",NQTN,".Power.by.FDR.CMLM",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.CMLM),paste(h2,"_",NQTN,".Power.by.TypeI.CMLM",".csv",sep=""))

write.csv(cbind(rep.FDR.MLMM,rep.power.MLMM),paste(h2,"_",NQTN,".Power.by.FDR.MLMM",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.MLMM),paste(h2,"_",NQTN,".Power.by.TypeI.MLMM",".csv",sep=""))

write.csv(cbind(rep.FDR.FarmCPU,rep.power.FarmCPU),paste(h2,"_",NQTN,".Power.by.FDR.FarmCPU",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.FarmCPU),paste(h2,"_",NQTN,".Power.by.TypeI.FarmCPU",".csv",sep=""))

write.csv(cbind(rep.FDR.Blink,rep.power.Blink),paste(h2,"_",NQTN,".Power.by.FDR.Blink",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.Blink),paste(h2,"_",NQTN,".Power.by.TypeI.Blink",".csv",sep=""))


write.csv(cbind(rep.FDR.GLM[,6],rep.power.GLM[,6],
                rep.FDR.MLM[,6],rep.power.MLM[,6],
                rep.FDR.CMLM[,6],rep.power.CMLM[,6],
                rep.FDR.MLMM[,6],rep.power.MLMM[,6],
                rep.FDR.SUPER[,6],rep.power.SUPER[,6],
                rep.FDR.FarmCPU[,6],rep.power.FarmCPU[,6],
                rep.FDR.Blink[,6],rep.power.Blink[,6]),
          paste(h2,"_",NQTN,".Power.by.FDR.GLM.MLM.SUPER.MLMM.FarmCPU.Blink",rep,".csv",sep="")
          )
name.of.trait=noquote(names(myY)[2])


pdf(paste("GAPIT.Power ", name.of.trait,".compare to GLM,MLM,CMLM,MLMM,SUPER,FarmCPU,Blink.", ".pdf", sep = ""), width = 4.5, height = 4.5,pointsize=9)
par(mar = c(5,6,5,3))
#win.graph(width=6, height=4, pointsize=9)
#palette(c("blue","red","green4","brown4","orange",rainbow(5)))
palette(c("green4","red","blue","brown4","orange",rainbow(7)))
plot(rep.FDR.SUPER[,6],rep.power.SUPER[,6],bg="lightgray",xlab="FDR",ylab="Power",ylim=c(0,1),xlim=c(0,1),main="Power against FDR",type="o",pch=20,col=1,cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1)
lines(rep.power.Blink[,6]~rep.FDR.Blink[,6], lwd=2,type="o",pch=20,col=2)
lines(rep.power.FarmCPU[,6]~rep.FDR.FarmCPU[,6], lwd=2,type="o",pch=20,col=3)
lines(rep.power.MLMM[,6]~rep.FDR.MLMM[,6], lwd=2,type="o",pch=20,col=4)
lines(rep.power.CMLM[,6]~rep.FDR.CMLM[,6], lwd=2,type="o",pch=20,col=5)
lines(rep.power.MLM[,6]~rep.FDR.MLM[,6], lwd=2,type="o",pch=20,col=6)
lines(rep.power.GLM[,6]~rep.FDR.GLM[,6], lwd=2,type="o",pch=20,col=7)
legend("bottomright",c("SUPER","Blink","FarmCPU","MLMM","CMLM","MLM","GLM"), pch = 20, lty =1,col=c(1:5),lwd=2,cex=1.0,bty="n")
#

dev.off()

###add type I error and power###

kkt<-cbind(rep.Power.Alpha.Blink[,1],rep.Power.Alpha.FarmCPU[,1],rep.Power.Alpha.SUPER[,1],rep.Power.Alpha.MLMM[,1],rep.Power.Alpha.CMLM[,1],rep.Power.Alpha.MLM[,1],rep.Power.Alpha.GLM[,1])
write.csv(cbind(myalpha,rep.Power.Alpha.Blink[,1],rep.Power.Alpha.FarmCPU[,1],rep.Power.Alpha.SUPER[,1],rep.Power.Alpha.MLMM[,1],rep.Power.Alpha.CMLM[,1],rep.Power.Alpha.MLM[,1],rep.Power.Alpha.GLM[,1]),paste(h2,"_",NQTN,".Type I error.Power.by.FDR.GLM.MLM.SUPER.MLMM.FarmCPU.Blink",rep,".csv",sep=""))

myalpha1<-myalpha/10

pdf(paste("GAPIT.Type I error_Power ", name.of.trait,".compare to GLM,MLM,CMLM,MLMM,SUPER,FarmCPU,Blink.", ".pdf", sep = ""), width = 6, height = 4.5,pointsize=9)
par(mar = c(5,6,5,3))
#win.graph(width=6, height=4, pointsize=9)
#palette(c("blue","red","green4","brown4","orange",rainbow(5)))
palette(c("green4","red","blue","brown4","orange",rainbow(7)))
plot(myalpha1,rep.Power.Alpha.SUPER[,1],log="x",bg="lightgray",xlab="Type I error",ylab="Power",main="Power against FDR",type="o",pch=20,col=1,cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1,ylim=c(min(kkt),max(kkt)))
#plot(myalpha1,rep.Power.Alpha.SUPER[,1],bg="lightgray",xlab="Type I error",ylab="Power",ylim=c(0,1),xlim=c(0,1),main="Power against FDR",type="o",pch=20,col=1,cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1)
lines(rep.Power.Alpha.Blink[,1]~myalpha1, lwd=2,type="o",pch=20,col=2)
lines(rep.Power.Alpha.FarmCPU[,1]~myalpha1, lwd=2,type="o",pch=20,col=3)
lines(rep.Power.Alpha.MLMM[,1]~myalpha1, lwd=2,type="o",pch=20,col=4)
lines(rep.Power.Alpha.CMLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=5)
lines(rep.Power.Alpha.MLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=6)
lines(rep.Power.Alpha.GLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=7)
legend("bottomright",c("SUPER","Blink","FarmCPU","MLMM","CMLM","MLM","GLM"), pch = 20, lty =1,col=c(1:7),lwd=2,cex=1.0,bty="n")
#

dev.off()


print(paste("GAPIT.Power ", name.of.trait,".compare to GLM,MLM,CMLM,MLMM,FarmCPU,Blink,SUPER.","successfully!" ,sep = ""))
