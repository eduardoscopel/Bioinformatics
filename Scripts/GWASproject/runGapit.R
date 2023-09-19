library(GAPIT3)
library(data.table)

df <- fread("~/Documents/Papers_Scopel/GWAS/sc_table_most_recent.txt")
nrow(df)
df <- df[df$reference == "Fay_et_al_2018" |
           df$reference == "Gallone_et_al_2016" |
           df$reference == "Han_et_al_2021" |
           df$reference == "Our_Data" |
           df$reference == "Peter_et_al_2018" |
           df$reference == "Pontes_et_al_2020" |
           df$reference == "Zhu_et_al_2016",]
nrow(df)

df$ploidyCons <- ifelse(is.na(df$ploidyFACS),df$ploidyBAF,
                        ifelse(!is.na(df$ploidyFACS) & df$ploidyBAF == "Unknown", df$ploidyFACS,
                               ifelse(df$ploidyFACS == df$ploidyBAF, df$ploidyBAF,
                                      ifelse(df$ploidyFACS != df$ploidyBAF, df$ploidyBAF,"Unknown"))))
df <- df[complete.cases(df$aneuploidy_binary),]
nrow(df)
# Remove low coverage sequences
df <- df[df$aneuploidy_binary != "LOWCOV" | 
           df$avg_cov >= 50,]
nrow(df)
# Remove strains derived from a single spore
df <- df[df$monosporic_derivative != "Yes",]
nrow(df)
# Remove highly contaminated (>10%) strains
df <- df[df$Contaminated < 0.1 |
           df$Contaminated == "No",]
nrow(df)
df$aneuploidy_num <- ifelse(df$aneuploidy_binary == "Yes", 1, 0)
diploids <- df[df$ploidyCons == 2,]
nrow(diploids)
table(diploids$aneuploidy_type, useNA = "ifany")
table(diploids$ecology_category, diploids$aneuploidy_type)
dipGainOnly <- diploids[diploids$aneuploidy_type == "Gain" | diploids$aneuploidy_type == "No",]
nrow(dipGainOnly)


myY <- data.frame(taxa = dipGainOnly$basename, ChrG = dipGainOnly$aneuploidy_num)
myY <- myY[-c(255,257),]
myY$taxa <- gsub('_','',myY$taxa)
nrow(myY)
myG <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.8/gain2n_gSNPs_QC4.hapmap.hmp.txt", header = FALSE,data.table = FALSE)
header <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.2/newIDs.txt")
myG[1,12:914] <- as.character(header$V1)
#myCV <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/PCAgain2n_gSNPs_QC4.evec", col.names = c("taxa","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","pop"),data.table = FALSE)
#myCV <- myCV[-c(255,257),]
#myCV <- myCV[,1:11]
#myCV$taxa <- gsub('_','',myCV$taxa)

#myEval <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/PCAgain2n_gSNPs_QC4.eval")
#sumEval <- sum(myEval$V1)
#currSum <- 0
#for(i in myEval$V1){
#  currSum <- currSum + i
#  if(currSum/sumEval > 0.995){
#    index <- grep(i, myEval$V1)
#    break
#  }
#}
myGAPIT0=GAPIT(G=myG,file.output=FALSE)

p0.05 <- 0.05/index
p0.01 <- 0.01/index
p0.005 <- 0.005/index

myGAPIT_MLM <- GAPIT(
  Y=myY,
  G=myG,
  CV=myCV,
  PCA.total=10,
  model=c("GLM", "MLM", "CMLM", "FarmCPU", "Blink"),
  Multiple_analysis=TRUE,
  cutOff = pThreshold)

simpleMeff <- 70297
p0.05Threshold <- 0.05/simpleMeff
### MLM
CMLM <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.8/GAPIT.MLM.ChrG.GWAS.Results.csv")
head(CMLM)
nrow(CMLM[CMLM$P.value <= p0.05Threshold,])
nrow(CMLM[CMLM$P.value <= p0.01,])
nrow(CMLM[CMLM$P.value <= p0.005,])
nrow(CMLM[CMLM$P.value <= 1e-04,])



### CMLM
CMLM <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/GAPIT.CMLM.ChrG.GWAS.Results.csv")
head(CMLM)
nrow(CMLM[CMLM$P.value <= 1.452513e-07,])
nrow(CMLM[CMLM$P.value <= p0.01,])
nrow(CMLM[CMLM$P.value <= p0.005,])
nrow(CMLM[CMLM$P.value <= 1e-04,])
nrow(CMLM[CMLM$P.value <= p0.05Threshold,])

### SUPER
SUPER <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.8/GAPIT.SUPER.ChrG.GWAS.Results.csv")
head(SUPER)
nrow(SUPER[SUPER$P.value <= p0.05Threshold,])
SUPER[SUPER$P.value <= p0.05Threshold,]

### MLMM
MLMM <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF5.0noRecomb/GAPIT.MLMM.ChrG.GWAS.Results.csv")
head(MLMM)
nrow(MLMM[MLMM$P.value <= p0.05Threshold,])
sigMLMM <- MLMM[MLMM$P.value <= p0.05Threshold,]
write.table(sigMLMM, 
            file = "~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF5.0noRecomb/GAPIT.MLMM.sigSNPs.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
myHapMap <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.8/gain2n_gSNPs_QC4.hapmap.hmp.txt", header = TRUE,data.table = FALSE)
MLMMSNPlist <- as.character(sigMLMM[,1])
# 7:85369
unique(as.character(myHapMap[myHapMap$`rs`== MLMMSNPlist[5],]))
Gind <- grep("G",as.character(myHapMap[myHapMap$`rs` == MLMMSNPlist[5],]))
Sind <- grep("S",as.character(myHapMap[myHapMap$`rs` == MLMMSNPlist[5],]))
myHapMap[myHapMap$`rs#` == MLMMSNPlist[5], Gind]
strainsGind <- colnames(myHapMap[myHapMap$`rs#` == MLMMSNPlist[5], Gind])[2:10]
dipGainOnly[dipGainOnly$basename == "YMD1834",]
myHapMap[myHapMap$`rs#` == MLMMSNPlist[5], Sind]
strainsSind <- colnames(myHapMap[myHapMap$`rs#` == MLMMSNPlist[5], Sind])[2:16]
dipGainOnly[dipGainOnly$basename == "BE021",]
dipGainOnly[which((dipGainOnly$basename %in% strainsSind)==TRUE),]
# 13:211388
unique(as.character(myHapMap[myHapMap$`rs`== MLMMSNPlist[9],]))
Tind <- grep("T",as.character(myHapMap[myHapMap$`rs` == MLMMSNPlist[9],]))
Yind <- grep("Y",as.character(myHapMap[myHapMap$`rs` == MLMMSNPlist[9],]))
myHapMap[myHapMap$`rs#` == MLMMSNPlist[9], Tind]
strainsTind <- colnames(myHapMap[myHapMap$`rs#` == MLMMSNPlist[9], Tind])
dipGainOnly[dipGainOnly$basename == "DBVPG6696",]
myHapMap[myHapMap$`rs#` == MLMMSNPlist[9], Yind]
strainsYind <- colnames(myHapMap[myHapMap$`rs#` == MLMMSNPlist[9], Sind])[2:16]
dipGainOnly[dipGainOnly$basename == "BE021",]

dipGainOnly[which((dipGainOnly$basename %in% c("DBVPG6696","EXF-7197","IMB_53","269521J","YJM464"))==TRUE),]




dipGainOnly[dipGainOnly$basename == "ATCC_52922_1C",]
dipGainOnly[dipGainOnly$basename == "CBS1196",]
dipGainOnly[dipGainOnly$basename == "CBS1592",]
dipGainOnly[dipGainOnly$basename == "CBS1593",]



### FarmCPU
FarmCPU <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.8/GAPIT.FarmCPU.ChrG.GWAS.Results.csv")
head(FarmCPU)
nrow(FarmCPU[FarmCPU$P.value <= p0.05Threshold,])
FarmCPU[FarmCPU$P.value <= p0.05Threshold,]

### Blink
Blink <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF5.0noRecomb/GAPIT.Blink.ChrG.GWAS.Results.csv")
head(Blink)
nrow(Blink[Blink$P.value <= p0.05Threshold,])
Blink[Blink$P.value <= p0.05Threshold,]
head(Blink)
nrow(Blink[Blink$P.value <= p0.05Threshold,])
sigBlink <- Blink[Blink$P.value <= p0.05Threshold,]
write.table(sigBlink, 
            file = "~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF5.0noRecomb/GAPIT.Blink.sigSNPs.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)


### FarmCPU
FarmCPU <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/PC0/GAPIT.FarmCPU.ChrG.GWAS.Results.csv")
head(FarmCPU)
nrow(FarmCPU[FarmCPU$P.value <= 1.452513e-07,])
sigFarmCPU <- FarmCPU[FarmCPU$P.value <= 1.452513e-07,]
unique(as.character(myHapMap[myHapMap$`rs#`== "7:853531",]))
sigFarmCPU$SNP
for(i in sigFarmCPU$SNP){
  print(unique(as.character(myHapMap[myHapMap$`rs#`== i,c(1,2,3,4,5)])))
}
unique(as.character(myHapMap[myHapMap$`rs#`== "4:354897",]))
Aind <- grep("A",as.character(myHapMap[myHapMap$`rs#`== "4:354897",]))
Rind <- grep("R",as.character(myHapMap[myHapMap$`rs#`== "4:354897",]))
myHapMap[myHapMap$`rs#`== "4:354897", Aind]
dipGainOnly[dipGainOnly$basename == "ATCC_52922_1C",]
dipGainOnly[dipGainOnly$basename == "CBS1196",]
dipGainOnly[dipGainOnly$basename == "CBS1592",]
dipGainOnly[dipGainOnly$basename == "CBS1593",]

myHapMap[myHapMap$`rs#`== "4:354897", Rind]
dipGainOnly[dipGainOnly$basename == "BMQ557",]
dipGainOnly[dipGainOnly$basename == "GZJ2",]
dipGainOnly[dipGainOnly$basename == "JSN3",]
dipGainOnly[dipGainOnly$basename == "K12",]
dipGainOnly[dipGainOnly$basename == "SXQ7",]


unique(as.character(myHapMap[myHapMap$`rs#`== "4:1184819",]))
Aind <- grep("A",as.character(myHapMap[myHapMap$`rs#`== "4:1184819",]))
Mind <- grep("M",as.character(myHapMap[myHapMap$`rs#`== "4:1184819",]))
myHapMap[myHapMap$`rs#`== "4:1184819", Aind]
dipGainOnly[dipGainOnly$basename == "CBS1489",]

myHapMap[myHapMap$`rs#`== "4:1184819", Mind]
dipGainOnly[dipGainOnly$basename == "BE011",]
dipGainOnly[dipGainOnly$basename == "BE039",]
dipGainOnly[dipGainOnly$basename == "BE040",]
dipGainOnly[dipGainOnly$basename == "YMD1834",]

unique(as.character(myHapMap[myHapMap$`rs#`== "16:121968",]))
Cind <- grep("C",as.character(myHapMap[myHapMap$`rs#`== "16:121968",]))
Yind <- grep("Y",as.character(myHapMap[myHapMap$`rs#`== "16:121968",]))
myHapMap[myHapMap$`rs#`== "16:121968", Cind]
dipGainOnly[dipGainOnly$basename == "DBVPG1841",]
dipGainOnly[dipGainOnly$basename == "DBVPG1843",]
dipGainOnly[dipGainOnly$basename == "DBVPG6696",]

myHapMap[myHapMap$`rs#`== "16:121968", Yind]
dipGainOnly[dipGainOnly$basename == "SAF26",]
dipGainOnly[dipGainOnly$basename == "SAN10",]

### BLINK
Blink <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/PC0/GAPIT.Blink.ChrG.GWAS.Results.csv")
head(Blink)
nrow(Blink[Blink$P.value <= 1.452513e-07,])
sigBlink <- Blink[Blink$P.value <= 1.452513e-07,]
unique(as.character(myHapMap[myHapMap$`rs#`== "16:205051",]))
sigFarmCPU$SNP
for(i in sigBlink$SNP){
  print(unique(as.character(myHapMap[myHapMap$`rs#`== i,c(1,2,3,4,5)])))
}

unique(as.character(myHapMap[myHapMap$`rs#`== "2:282655",]))
Tind <- grep("T",as.character(myHapMap[myHapMap$`rs#`== "2:282655",]))
Wind <- grep("W",as.character(myHapMap[myHapMap$`rs#`== "2:282655",]))
colnames(myHapMap[myHapMap$`rs#`== "2:282655", Wind])[1]
dipGainOnly[dipGainOnly$basename == "1663",]
dipGainOnly[dipGainOnly$basename == "319_5C",]
dipGainOnly[dipGainOnly$basename == "646_3B",]
dipGainOnly[dipGainOnly$basename == "79",]
dipGainOnly[dipGainOnly$basename == "CBS6216",]
dipGainOnly[dipGainOnly$basename == "CBS8292",]
dipGainOnly[dipGainOnly$basename == "CCY_21_4_97",]
dipGainOnly[dipGainOnly$basename == "CCY_21_4_98",]
dipGainOnly[dipGainOnly$basename == "CECT10131",]
dipGainOnly[dipGainOnly$basename == "CLIB318_1",]
dipGainOnly[dipGainOnly$basename == "CLIB326_1",]
dipGainOnly[dipGainOnly$basename == "CLIB340",]
dipGainOnly[dipGainOnly$basename == "CLQCA_10_386",]
dipGainOnly[dipGainOnly$basename == "SBE_1C",]
dipGainOnly[dipGainOnly$basename == "UWOPS83_883_2",]
dipGainOnly[dipGainOnly$basename == "WI008",]
dipGainOnly[dipGainOnly$basename == "Win_8B",]
dipGainOnly[dipGainOnly$basename == "Y6_b",]
dipGainOnly[dipGainOnly$basename == "YJM1111",]
dipGainOnly[dipGainOnly$basename == "YJM223",]
dipGainOnly[dipGainOnly$basename == "YJM634",]
dipGainOnly[dipGainOnly$basename == "YS15",]
dipGainOnly[dipGainOnly$basename == "YS20",]
dipGainOnly[dipGainOnly$basename == "YS22_E_",]


unique(as.character(myHapMap[myHapMap$`rs#`== "12:732685",]))
Cind <- grep("C",as.character(myHapMap[myHapMap$`rs#`== "12:732685",]))
Yind <- grep("Y",as.character(myHapMap[myHapMap$`rs#`== "12:732685",]))
myHapMap[myHapMap$`rs#`== "12:732685", Cind]
dipGainOnly[dipGainOnly$basename == "1235",]
dipGainOnly[dipGainOnly$basename == "AGME_5I",]
dipGainOnly[dipGainOnly$basename == "DBVPG1107",]
dipGainOnly[dipGainOnly$basename == "NCYC_2743",]

unique(as.character(myHapMap[myHapMap$`rs#`== "13:229233",]))
Cind <- grep("C",as.character(myHapMap[myHapMap$`rs#`== "13:229233",]))
Yind <- grep("Y",as.character(myHapMap[myHapMap$`rs#`== "13:229233",]))
myHapMap[myHapMap$`rs#`== "13:229233", Cind]
dipGainOnly[dipGainOnly$basename == "N37_1A",]
dipGainOnly[dipGainOnly$basename == "YJM428_1b",]
myHapMap[myHapMap$`rs#`== "13:229233", Yind]
dipGainOnly[dipGainOnly$basename == "CBS7836",]
dipGainOnly[dipGainOnly$basename == "YJM308",]


myHapMap <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/gain2n_gSNPs_QC4.hapmap.hmp.txt", header = TRUE,data.table = FALSE)
unique(as.character(myHapMap[myHapMap$`rs#`== "7:157678",]))
Cind <- grep("C",as.character(myHapMap[myHapMap$`rs#`== "7:157678",]))
Mind <- grep("M",as.character(myHapMap[myHapMap$`rs#`== "7:157678",]))
myHapMap[myHapMap$`rs#`== "7:157678", Cind]
myHapMap[myHapMap$`rs#`== "7:157678", Mind]

PCModel <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/PCModelSelection/GAPIT.MLM.ChrG.BIC.Model.Selection.Results.csv", 
                    col.names = c("PC","BIC","logLike"))
