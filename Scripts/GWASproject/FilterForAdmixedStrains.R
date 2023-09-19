library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

df <- fread("~/Documents/Papers_Scopel/GWAS/sc_table_most_recent.txt", data.table = FALSE,blank.lines.skip = TRUE)
nrow(df)
newdf <- data.frame(matrix(NA, nrow = nrow(df),ncol = 16))
chr <- seq(1,16)
chr <- paste("c", chr, sep = "")
colnames(newdf) <- chr
df <- cbind(df, newdf)
head(df)

first <- strsplit(df$aneuploidy_RD[388],split = ";")
second <- strsplit(first[[1]], split = "@")
third <- strsplit(second[[3]], split = ",")

for(i in as.numeric(rownames(df))){
  first <- strsplit(df$aneuploidy_RD[i], split = ";")
  print(first[[1]])
  if(!(first %in% c(NA, "No", "LOWCOV"))){
    for(j in first){
      second <- strsplit(j, split = "@")
      for(k in second){
        third <- strsplit(k, split = ",")
        type <- third[[1]]
        print(third)
        for(w in third[[2]]){
          column <- as.numeric(w) + 26
          df[i,column]
          df[i,column] <- type
          df[i,column]
        }
      }
    }
  }
}

df

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


######### delete when new tree is finished!!!!!!!!
diploids[diploids$basename == "NIGF1","aneuploidy_type"] <- "No" ######### delete when new tree is finished!!!!!!!!
######### delete when new tree is finished!!!!!!!!


dipGainOnly <- diploids[diploids$aneuploidy_type == "Gain" | diploids$aneuploidy_type == "No",]
nrow(dipGainOnly)
dipGainOnly$lineage <- ifelse(dipGainOnly$population_old %in% c("Europe/wine", "Wine","wine/European","1._Wine/European_(subclade_4)","1._Wine/European_" ,"1._Wine/European_(subclade_1)","1._Wine/European_(subclade_3)","1._Wine/European_(subclade_2)","Wine - Main"),"Wine",
                              ifelse(dipGainOnly$population_old %in% c("Beer2","11._Ale_beer_","Beer 2"), "Ale Beer",
                                     ifelse(dipGainOnly$population_old %in% c("Asian","24._Asian_islands_","Philippines"),"Asian Islands",
                                            ifelse(dipGainOnly$population_old %in% c("Mosaic","Mosaic Honey Wine","Mosaic 3 (Baijiu)","Mosaic 4 (South Africa)","Mosaic 5 (South Africa)","M2._Mosaic_region_2","M3._Mosaic_region_3","M3._Mosaic_region_3_","M1._Mosaic_region_1","M1._Mosaic_region_1_","M2._Mosaic_region_2_"), "Mosaic",
                                                   ifelse(dipGainOnly$population_old %in% c("Mixed 1", "8._Mixed_origin_"), "Mixed Origin",
                                                          ifelse(dipGainOnly$population_old %in% c("Huangjiu","Baijiu","Qingkejiu","Mantou 7","Daqu/Baijiu","26._Asian_fermentation_"), "Asian Ferm.",
                                                                 ifelse(dipGainOnly$population_old %in% c("African Honey Wine"), "Ethiopia",
                                                                        ifelse(dipGainOnly$population_old %in% c("African_Palm_Wine","13._African_palm_wine_"), "African Wine",
                                                                               ifelse(dipGainOnly$population_old %in% c("Mauritius/South Africa"),"South Africa",
                                                                                      ifelse(dipGainOnly$population_old %in% c("West African Beer","South African Beer","6._African_beer_"), "African Beer",
                                                                                             ifelse(dipGainOnly$population_old %in% c("7._Mosaic_beer_"), "Mosaic Beer",
                                                                                                    ifelse(dipGainOnly$population_old %in% c("2._Alpechin_"), "Alpechin",
                                                                                                           ifelse(dipGainOnly$population_old %in% c("4._Mediterranean_oak_"), "Med Oak",
                                                                                                                  ifelse(dipGainOnly$population_old %in% c("12._West_African_cocoa_"), "West African Cocoa",
                                                                                                                         ifelse(dipGainOnly$population_old %in% c("18._Far_East_Asia_"), "Far East Asian",
                                                                                                                                ifelse(dipGainOnly$population_old %in% c("25._Sake_","Sake"), "Sake",
                                                                                                                                       ifelse(dipGainOnly$population_old %in% c("5._French_dairy_","Dairy"), "French Dairy",
                                                                                                                                              ifelse(dipGainOnly$population_old %in% c("10._French_Guiana_human_"), "French Guiana Human",
                                                                                                                                                     ifelse(dipGainOnly$population_old %in% c("21._Ecuadorean_"),"Ecuadorean",
                                                                                                                                                            ifelse(dipGainOnly$population_old %in% c("17._Taiwanese_"), "Taiwan",
                                                                                                                                                                   ifelse(dipGainOnly$population_old %in% c("14._CHNIII_","20._CHN_V_","16._CHNI","15._CHNII_"), "China",
                                                                                                                                                                          ifelse(dipGainOnly$population_old %in% "3._Brazilian_bioethanol_", "Brazilian Bioethanol",
                                                                                                                                                                                 ifelse(dipGainOnly$population_old %in% "9._Mexican_agave_", "Mexican Agave",
                                                                                                                                                                                        ifelse(dipGainOnly$population_old %in% c("22._Far_East_Russian_"), "Far East Russian",
                                                                                                                                                                                               ifelse(dipGainOnly$population_old %in% c("19._Malaysian_"), "Malaysian",
                                                                                                                                                                                                      ifelse(dipGainOnly$population_old %in% c("23._North_American_oak_"),"North American Oak","NA"))))))))))))))))))))))))))

### Reassigning lineages based on NJ clustering 
# Convert lineage to Wine
dipGainOnly[dipGainOnly$basename %in% c("GSY725","CBS9564","YJM1005","GSY723","YJM1095","YJM963","YJM964","YJM965","YJM956","YJM957","YJM967","YJM947",
                                        "YJM955","YJM434","B68019c","GSY1033","CLIB1060","PLU19b.1", "UCD_40_346", "SJ5L14", "995", "YJM332", "DAVAUXa.1", 
                                        "PLU28a.1", "PLU29a.1", "PLU11a.1","PLU22a.1","PLU28b.1","PLU17a.1","PLU12a.1","PLU15b.1","DAVAUXb.1","SON4c.1",
                                        "CLIB1085","CLIB1083","CLIB1059","M9"), "lineage"] <- "Wine"

# Convert lineage to Alpechin
dipGainOnly[dipGainOnly$basename == "CBS2910","lineage"] <- "Alpechin"


# Convert lineage to Brazilian Bioethanol
dipGainOnly[dipGainOnly$basename == "SA_1_5_","basename"] <- "SA_1_5"
dipGainOnly[dipGainOnly$basename %in% c("BI002", "SA_1_5", "BI005", "SP004"), "lineage"] <- "Brazilian Bioethanol"

# Convert lineage to Israel
dipGainOnly[dipGainOnly$basename %in% c("35","33","60","59","13","34"),"lineage"] <- "Israel"

# Convert lineage to Ecuador1
dipGainOnly[dipGainOnly$basename %in% c("CLQCA_10_027","CLQCA_02_003","CLQCA_17_060"),"lineage"] <- "Ecuador1"

# Convert lineage to Ale Beer
dipGainOnly[dipGainOnly$basename == "TUMPI-BA-105","basename"] <- "TUMPI-BB-105"
dipGainOnly[dipGainOnly$basename %in% c("YJM1100","CBS6505","YJM439","1175","NCYC_88","YJM1098","TUMPI-BB-105","CBS1398","YMD1834","WLP570","CBS7539"), "lineage"] <- "Ale Beer"

# Convert lineage to Mosaic 
dipGainOnly[dipGainOnly$basename %in% c("IFO_0877_6_1","YJM1289","CBS4255","YJM223","YJM1143","YJM1101","PYCC4654","B66044","YJM1099","YJM1259","CBS7836","YJM464",
                                        "YJM1094","YJM1096","YJM1097","YJM1135","YJM670","GSY1034","YJM436","AN1f.2.1","PYCC8032", "YJM455","YJM1119",
                                        "CBS7840","YJM1124","YJM1178","EXF_6780", "EXF_5284","UCD_61_190_6A"), "lineage"] <- "Mosaic"

# Convert lineage to French Dairy
dipGainOnly[dipGainOnly$basename %in% c("YJM1141","YJM949"), "lineage"] <- "French Dairy"

# Convert lineage to African Beer
dipGainOnly[dipGainOnly$basename == "N134_7_1_a_","basename"] <- "N134_7_1_a"
dipGainOnly[dipGainOnly$basename == "N134_7_1_a_","lineage"] <- "African Beer"

# Convert lineage to Ethiopia
dipGainOnly[dipGainOnly$basename %in% c("DBVPG1848", "DBVPG1895", "ETPF6","ETPF2"), "lineage"] <- "Ethiopia"

# Convert lineage to FGH
dipGainOnly[dipGainOnly$basename %in% c("CLQCA_20_060", "CEY647"), "lineage"] <- "French Guiana Human"

# Convert lineage to Mexican Agave
dipGainOnly[dipGainOnly$basename %in% c("906","908"),"lineage"] <- "Mexican Agave"

# Convert lineage to Mixed Origin
dipGainOnly[dipGainOnly$basename %in% c("YJM946", "WI008","YJM223", "YJM634", "YJM671", "CBS6308","CBS2165a"),"lineage"] <- "Mixed Origin"

# Convert lineage to Mosaic
dipGainOnly[dipGainOnly$basename %in% c("B68549","YJM669", "YJM1102","AN3e.1.1"), "lineage"] <- "Mosaic"

# Convert lineage to Clinical
dipGainOnly[dipGainOnly$basename %in% c("YJM339","CBS7839","YJM560","YJM467","YJM1125","YJM1115","YJM1112","YJM1116","YJM678","YJM676","CBS7833","YJM560",
                                        "YJM467","YJM1125","YJM1115","YJM1112","YJM1116","YJM678","YJM676","CBS7833","3B4103A","YJM1108","92-123","YJM440",
                                        "YJM1121", "CBS7835","YJM1117","YJM1114","DBVPG6874","CBS7837","YJM521","YJM1111","YJM677","CBS7838","YJM1122",
                                        "384103A", "ATCC_38618_2_2"),"lineage"] <- "Clinical"

# Convert lineage to African West Cocoa
dipGainOnly[dipGainOnly$basename %in% c("WL001","WL003"),"lineage"] <- "West African Cocoa"

# Convert lineages to Siberia
dipGainOnly[dipGainOnly$basename %in% c("N39_7A","N37_1A","N38_4A"), "lineage"] <- "Siberia"

# Convert lineages to Asian Ferm.
dipGainOnly[dipGainOnly$basename %in% c("ATCC_52922_1C", "BI001","BI003","BI004","JSN3", "SAN33","BJQ2","SAN31","HLJU1","AHN1","SXQ1","SCN9","SAN30","SXQ7","SXQ2",
                                        "SAF29", "MAUN4","SAF28","SAF27","SAN32"), "lineage"] <- "Asian Ferm."

# Convert lineages to Sake 
dipGainOnly[dipGainOnly$basename %in% c("SA003","SA007", "SA004","SA001","SA005","SA006"), "lineage"] <- "Sake"

# Convert lineages to African Wine
dipGainOnly[dipGainOnly$basename == "Dji2_2A_a_","basename"] <- "Dji2_2A_a"
dipGainOnly[dipGainOnly$basename %in% c("Dji2_2A_a"), "lineage"] <- "African Wine"

# Convert lineages to North American Oak
dipGainOnly[dipGainOnly$basename %in% c("CLIB414","ATCC_66348_1D","YPS670","N95_5_1A","2163","IY_03_5_30_1_1_1_1","IY_03_5_26_5_1_1_1","N95_6_1C"), "lineage"] <- "North American Oak"

# Convert lineages to Taiwan2
dipGainOnly[dipGainOnly$basename %in% c("ES2M03_7A","ES2M03","ES4M07", "SJ5L12"), "lineage"] <- "Taiwan2"

nrow(dipGainOnly)
dipGainOnly[duplicated(dipGainOnly$basename),]
dipGainOnly <- dipGainOnly[!duplicated(dipGainOnly$basename),]
nrow(dipGainOnly)

myY <- dipGainOnly[, c(20, 45)]
#myY <- data.frame(taxa = dipGainOnly$basename, ChrG = dipGainOnly$aneuploidy_num, clade = dipGainOnly$lineage, specAneuploidy = dipGainOnly$aneuploidy_RD)
#myY <- myY[c(255,257),]
#colnames(myY)[1:4] <- c("taxa", "ChrG", "clade", "specAneuploidy")
#myY$taxa <- gsub('_','',myY$taxa)
nrow(myY)
myY <- myY[!duplicated(myY$basename),]
#aneuploidyDF <- reshape2::melt(myY,id.vars = c("taxa","ChrG","clade", "specAneuploidy"), variable.name = "chr",value.name = "an")
#aneuploidyDF[is.na(aneuploidyDF$an),]$an <- 0
#ancestry <- data.frame(pop =ids$ancestry,strain = ids$strain)
strains <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/RunQC40Miss10MAF0.2Recomb0.2/newIDs.txt")
colnames(strains) <- "basename"
strains$id <- 1:nrow(strains)
setdiff(myY$basename, strains$basename)
setdiff(strains$basename,myY$basename)
myY <- merge(myY, strains,by = "basename",all = FALSE,no.dups = TRUE,sort = TRUE)
nrow(myY)
myY <- myY[order(myY$id),]

Q23 <- read.table("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/admixture/t500bp/s1/AlignedMatrices/1645028616/aligned.files/ready/gain2n_gSNPs_q40_500bp.23.Q.converted.ready")
Q23 <- Q23[!(row.names(Q23) == 602),]
k <- length(Q23)
for(j in 1:nrow(Q23)){
  if(rowSums(Q23[j,1:k] > 0.8) == 0){
    Q23$Panc[j] <- "admix"
  }
  else{
    Q23$Panc[j] <- which.max(Q23[j,1:k])
  }
  Q23$maxV[j] <- max(Q23[j,1:k])
}

Q23 <- cbind(Q23,myY)  
rownames(Q23) <- Q23$basename
nrow(Q23[Q23$Panc == "admix",])
admixedStrains <- Q23[Q23$Panc == "admix",]
admixedDF <- dipGainOnly[dipGainOnly$basename %in% rownames(admixedStrains),]
admixedDF <- admixedDF[!duplicated(admixedDF$basename),]

write.table(x = admixedDF$basename,file = "~/Documents/Papers_Scopel/GWAS/diploids/admixture/AdmixedStrains.txt",quote = FALSE,sep = "\t", row.names = FALSE,col.names = FALSE)

myG <- fread("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/AdmixRunQC40Miss10MAF2.0noRecomb/admixed_gain2n_gSNPs_QC3.hmp.txt", header = TRUE,data.table = FALSE)
head(myG)
p0.05Threshold <- 0.05/114651
Blink <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/AdmixRunQC40Miss10MAF2.0noRecomb/GAPIT.Blink.ChrG.GWAS.Results.csv")
MLM <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/AdmixRunQC40Miss10MAF2.0noRecomb/GAPIT.MLM.ChrG.GWAS.Results.csv")
CMLM <- read.csv("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/AdmixRunQC40Miss10MAF2.0noRecomb/GAPIT.CMLM.ChrG.GWAS.Results.csv")

head(Blink)
nrow(Blink[Blink$P.value <= p0.05Threshold,])
sigBlink <- Blink[Blink$P.value <= p0.05Threshold,]
write.table(sigBlink, 
      file = "~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/AdmixRunQC40Miss10MAF2.0noRecomb/GAPIT.Blink.sigSNPs.tsv",
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE)

head(MLM)
nrow(MLM[MLM$P.value <= p0.05Threshold,])
sigMLM <- MLM[MLM$P.value <= p0.05Threshold,]

head(CMLM)
nrow(CMLM[CMLM$P.value <= p0.05Threshold,])
sigCMLM <- CMLM[CMLM$P.value <= p0.05Threshold,]
