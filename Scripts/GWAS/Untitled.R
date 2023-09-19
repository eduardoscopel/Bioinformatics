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
dipGainOnly <- dipGainOnly[!duplicated(dipGainOnly$basename),]
nrow(dipGainOnly)

MLMM <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/Filter_MLMM.ChrG_GWAS_result.txt")
FarmCPU <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/Filter_FarmCPU.ChrG_GWAS_result.txt")
Blink <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/Filter_Blink.ChrG_GWAS_result.txt")


mergedSNPs <- c(as.character(MLMM$SNP), as.character(FarmCPU$SNP), as.character(Blink$SNP))
mergedSNPs <- mergedSNPs[!duplicated(mergedSNPs)]

myHapMap <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF1.0Recomb0.8/gain2n_gSNPs_QC4.hmp.txt", header = TRUE,data.table = FALSE)
header <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.2/newIDs.txt")
colnames(myHapMap)[12:914] <- as.character(header$V1)

temp <- data.frame(row.names = colnames(myHapMap[,-c(1:11)]))

ref <- strsplit(myHapMap[myHapMap$`rs` ==  mergedSNPs[1],]$alleles,"/")[[1]][1]
mut <- strsplit(myHapMap[myHapMap$`rs` ==  mergedSNPs[1],]$alleles,"/")[[1]][2]
temp <- as.data.frame(t(myHapMap[myHapMap$`rs` ==  mergedSNPs[1],-c(1:11)]))
a <- mergedSNPs[1]
temp$`125262` <- ifelse(temp$`125262` == ref,"ref",
                             ifelse(temp$`125262` == mut, "mut","het"))


StrainSNP <- data.frame(Strain = dipGainOnly$basename)

rbind(StrainSNP,mergedSNPs)
