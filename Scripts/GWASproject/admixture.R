library(ggplot2)
library(data.table)
library(cowplot)


logtab1 <- fread("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/S1.txt")
logtab1$seed <- "Run1"
logtab12345 <- fread("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/S12345.txt")
logtab12345$seed <- "Run2"
logtab22222 <- fread("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/S22222.txt")
logtab22222$seed <- "Run3"
logtab54321 <- fread("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/S54321.txt")
logtab54321$seed <- "Run4"
logtab9999 <- fread("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/S9999.txt")
logtab9999$seed <- "Run5"

logtab <- rbind(logtab1, logtab12345, logtab22222, logtab54321, logtab9999)
colnames(logtab) <- c("K","CV","Loglikelihood","seed")
logtab <- logtab[!logtab$K ==1,]
logtab$K <- as.factor(logtab$K)
logplot <- ggplot(logtab, aes(x=K, y= Loglikelihood)) + 
  geom_boxplot(outlier.shape =  NA) + 
  geom_jitter(size = 1.5, alpha = 0.4,position = position_jitter(0.2)) + 
  scale_y_continuous(n.breaks = 8,limits = c(-1900000,-800000), expand = c(0,0))+
  scale_x_discrete(breaks = seq(2,40,2))+
  xlab("Clusters")+
  geom_vline(xintercept = c(21,32), linetype = "dotted")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())
logplot
ggsave("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/logLikePlot.png",
       width=10, 
       height=5, 
       units="in",
       limitsize=FALSE, 
       dpi = 300) 

CVplot <- ggplot(logtab, aes(x=K, y= CV)) + 
  geom_boxplot() + 
  geom_jitter(size = 1.5, alpha = 0.4,position = position_jitter(0.2)) + 
  scale_y_continuous(n.breaks = 9,limits = c(0.045,0.09), expand = c(0,0))+
  xlab("Clusters")+
  ylab("Cross-Validation errors") + 
  geom_vline(xintercept = c(9,19,29,39), linetype = "dotted")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())
CVplot
       
#plot_grid(logplot, CVplot)

logtab <- logtab[order(logtab$seed,logtab$K),]
scaledDeltaList <- c()
for(i in 1:40){
  Pop1Mean <- mean(logtab[logtab$K==i,Loglikelihood])
  Pop2Mean <- mean(logtab[logtab$K==i+1,Loglikelihood])
  Pop2SD <- sd(logtab[logtab$K==i+1,Loglikelihood])
  scaledDelta <- (Pop2Mean - Pop1Mean)/Pop2SD
  scaledDeltaList <- append(scaledDeltaList, scaledDelta)
}
scaledDeltaList <- scaledDeltaList[2:39]
SDLdf <- data.frame(K=seq(3,40,1), scaledDelta = scaledDeltaList)
scaledDeltaPlot <- ggplot(SDLdf, aes(x=K,y=scaledDelta))+
  geom_point() + 
  geom_line() +
  scale_x_continuous(expand = c(0.01,0.01),breaks = seq(4,40,2)) + 
  scale_y_continuous(expand = c(0.01,0.01),breaks = seq(0,12,1)) + 
  theme_bw() +
  xlab("Clusters") +
  ylab("Scaled change in LogLikelihood") +
  geom_vline(xintercept = c(22,33), linetype = "dotted")+
  theme(panel.grid.minor = element_blank())
scaledDeltaPlot

topRowAdmix <- plot_grid(logplot, scaledDeltaPlot, labels = c("(A)","(B)"), label_size = 12, nrow = 1)


ggsave("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/scaledDelta.png",
       width=8, 
       height=5, 
       units="in",
       limitsize=FALSE, 
       dpi = 300) 


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

indOrder <- read.table("../../../../Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/orderIDs.txt")
colnames(indOrder) <- "ID"
indTab <- data.frame(ID = dipGainOnly$basename, Sex = "Unknown", Lineage = dipGainOnly$lineage)
#indTab$ID <- gsub("SA_1_5","SA_1_5_", indTab$ID)
#indTab$ID <- gsub("Dji2_2A_a","Dji2_2A_a_", indTab$ID)
#indTab$ID <- gsub("N134_7_1_a", "N134_7_1_a_", indTab$ID)
indTab2 <- merge(indOrder,indTab, by = "ID",sort = FALSE)
lineageTab <- indTab2[,c(1,3)]
head(lineageTab)

Qlist <- list.files("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/s12345/ready/")
popIndices <- as.numeric(tapply(seq_along(sort(lineageTab$Lineage)),sort(lineageTab$Lineage), max))+0.5
popMean <- as.numeric(tapply(seq_along(sort(lineageTab$Lineage)),sort(lineageTab$Lineage), mean))+0.5
popMean[25] <- popMean[25]-1
popMean[26] <- popMean[26]
popMean[27] <- popMean[27]+1
popMean[28] <- popMean[28]+1
popMean[29]
first_af1 <- "slateblue"
second_af2 <- "slateblue4"
third_fgh <- "darkorange2"
fourth_we1 <- "purple1"
fifth_fd <- "cadetblue2"
sixth_mo <- "deeppink"
seventh_we2 <- "purple4"
eigth_aleb <- "orange"
ninth_afw1 <- "mediumvioletred"
tenth_alpe <- "olivedrab"
eleventh_we3 <- "mediumpurple1"
twelveth_sa <- "red"
thirteenth_other1 <- "black"
fourteenth_we4 <- "mediumpurple4"
fifteenth_NA1 <- "gray30"
sixteenth_afw2 <- "palevioletred4"
seventeenth_afb1 <- "darkkhaki"
eighteenth_ecu <- "blue"
nineteenth_other2 <- "ivory4"
twenteth_NA2 <- "gray50"
twent1_we5 <- "mediumorchid1"
twent2_tai <- "aquamarine"
twent3_NA3 <- "gray70"
twent4_medoak <- "forestgreen"
twent5_M1 <- "hotpink"
twent6_mo2 <- "lightcoral"
twent7_fea <- "aquamarine4"
twent8_we6 <- "orchid4"
twent9_ch <- "lightskyblue"
thirth_NA <- "yellow"

col_list<-c(first_af1,
            second_af2,
            third_fgh,
            fourth_we1,
            fifth_fd,
            sixth_mo,
            seventh_we2,
            eigth_aleb,
            ninth_afw1,
            tenth_alpe,
            eleventh_we3,
            twelveth_sa,
            thirteenth_other1,
            fourteenth_we4,
            fifteenth_NA1,
            sixteenth_afw2,
            seventeenth_afb1,
            eighteenth_ecu,
            nineteenth_other2,
            twenteth_NA2,
            twent1_we5,
            twent2_tai,
            twent3_NA3,
            twent4_medoak,
            twent5_M1,
            twent6_mo2,
            twent7_fea,
            twent8_we6,
            twent9_ch,
            thirth_NA)

plotList <- c()
counter <- 1
for(i in Qlist[2:40]){
  temp <- read.table(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/s1/ready/",i,sep = ""))
  k <- length(temp)
  for(j in 1:nrow(temp)){
    if(rowSums(temp[j,1:k] > 0.8) == 0){
      temp$Panc[j] <- "admix"
    }
    else{
      temp$Panc[j] <- which.max(temp[j,1:k])
    }
    temp$maxV[j] <- max(temp[j,1:k])
  }
  temp <- cbind(temp,lineageTab)  
  rownames(temp) <- temp$ID
  #temp <- temp[order(temp$clade),] # order by lineage
  temp <- temp[order(temp$Lineage,temp$Panc, temp$maxV),] # order by % ancestry
  lineageOrder <- temp$Lineage
  temp2 <- temp[,1:k]
  keys <- rownames(temp2)
  keyDF <- data.frame(key = keys, weigth = 1:length(keys))
  melteddf <- melt(t(temp2))
  head(melteddf)
  if(k %% 10 == 0){
    plotList[[counter]] <- ggplot(melteddf) + 
      geom_bar(aes(x=Var2,y=value,color=Var1,fill=Var1),stat="identity", show.legend = FALSE) +
      geom_vline(xintercept = popIndices, size=.8, size=.4, linetype = "dotted") +
      scale_color_manual(values = c(cladecolors,col_list))+
      scale_fill_manual(values = c(cladecolors,col_list)) +
      scale_x_discrete(expand = c(0,0), breaks = temp$ID[popMean], labels = temp$Lineage[popMean])+
      scale_y_continuous(expand = c(0,0), breaks = NULL,name = paste("k=",k,sep=""))+
      #ggtitle(paste("k = ",k))+
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = 6),
            axis.title.x = element_blank(), 
            plot.margin = unit(c(0,0.1,0,0.1),"cm")) 
    #plot.title = element_text(size = 8,hjust = 0.5,margin = unit(c(0,0,0,0),"cm")))
    #theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())  
  }
  else{
    plotList[[counter]] <- ggplot(melteddf) + 
      geom_bar(aes(x=Var2,y=value,color=Var1,fill=Var1),stat="identity", show.legend = FALSE) +
      geom_vline(xintercept = popIndices, size=.4, linetype = "dotted")+
      scale_color_manual(values = c(cladecolors,col_list))+
      scale_fill_manual(values = c(cladecolors,col_list)) +
      scale_x_discrete(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0),breaks = NULL,name = paste("k=",k,sep=""))+
      #ggtitle(paste("k = ",k, sep = ""))+
      #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4))
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            plot.margin = unit(c(0,0.1,0,0.1),"cm"),
            plot.title = element_blank())
    #plot.title = element_text(size = 8,hjust = 0.5,margin = unit(c(0,0,0,0),"cm")))  
  }
  counter <- counter + 1 
}


plot_grid(plotlist = plotList[20:29], nrow = 10, rel_heights = c(1,1,1,1,1,1,1,1,1,1.8))

plotList <- c()
counter <- 1
NJOrder <- c("European Wine", "Alpechin", "Bioethanol", "Israel", "Ecuador1","Ale Beer","French Dairy", "Ethiopia", "African Beer", "French Guiana Human", 
             "Mexican Agave", "Med Oak", "Mixed Origin", "Clinical", "West African Cocoa", "Siberia", "Sake", "Asian Ferm.", "South Africa", "African Wine", 
             "China", "North American Oak", "Asian Islands", "Taiwan2", "Far East Russia", "Ecuador2", "Malaysian","Far East Asian", "Taiwan", "Mosaic")
for(i in Qlist[29:40]){
  temp <- read.table(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/s12345/ready/",i,sep = ""))
  k <- length(temp)
  for(j in 1:nrow(temp)){
    if(rowSums(temp[j,1:k] > 0.8) == 0){
      temp$Panc[j] <- "admix"
    }
    else{
      temp$Panc[j] <- which.max(temp[j,1:k])
    }
    temp$maxV[j] <- max(temp[j,1:k])
  }
  temp <- cbind(temp,lineageTab)  
  rownames(temp) <- temp$ID
  #temp <- temp[order(temp$clade),] # order by lineage
  temp <- temp[order(temp$Lineage,temp$Panc, temp$maxV),] # order by % ancestry
  lineageOrder <- temp$Lineage
  temp2 <- temp[,1:k]
  keys <- rownames(temp2)
  keyDF <- data.frame(key = keys, weigth = 1:length(keys))
  melteddf <- melt(t(temp2))
  head(melteddf)
  if(k == 40){
    plotList[[counter]] <- ggplot(melteddf) + 
      geom_bar(aes(x=Var2,y=value,color=Var1,fill=Var1),stat="identity", show.legend = FALSE) +
      geom_vline(xintercept = popIndices, size=.4, linetype = "dotted") +
      scale_color_manual(values = c(cladecolors,col_list))+
      scale_fill_manual(values = c(cladecolors,col_list)) +
      scale_x_discrete(expand = c(0,0), breaks = temp$ID[popMean], labels = temp$Lineage[popMean])+
      scale_y_continuous(expand = c(0,0), breaks = NULL,name = k)+
      #ggtitle(paste("k = ",k))+
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = 12),
            axis.title.x = element_blank(), 
            plot.margin = unit(c(0,0.1,0,0.1),"cm")) 
    #plot.title = element_text(size = 8,hjust = 0.5,margin = unit(c(0,0,0,0),"cm")))
    #theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())  
  }
  else{
    plotList[[counter]] <- ggplot(melteddf) + 
      geom_bar(aes(x=Var2,y=value,color=Var1,fill=Var1),stat="identity", show.legend = FALSE) +
      geom_vline(xintercept = popIndices, size=.4, linetype = "dotted")+
      scale_color_manual(values = c(cladecolors,col_list))+
      scale_fill_manual(values = c(cladecolors,col_list)) +
      scale_x_discrete(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0),breaks = NULL,name = k)+
      #ggtitle(paste("k = ",k, sep = ""))+
      #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4))
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            plot.margin = unit(c(0,0.1,0,0.1),"cm"),
            plot.title = element_blank())
    #plot.title = element_text(size = 8,hjust = 0.5,margin = unit(c(0,0,0,0),"cm")))  
  }
  counter <- counter + 1 
}

topRowAdmix <- plot_grid(logplot, scaledDeltaPlot, labels = c("(A)","(B)"), label_size = 12, nrow = 1,hjust = -0.1)
bottomRowAdmix <- plot_grid(plotlist = plotList, nrow = 12, rel_heights = c(rep(1,11),5.7))
AdmixPlots <- plot_grid(topRowAdmix, bottomRowAdmix, nrow = 2,label_size = 12, labels = c(NA,"(C)"), rel_heights = c(1,2.1),hjust = -0.1,vjust = 0)
ggsave("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/AdmixPlotv2.png",
       plot = AdmixPlots,
       width = 12,
       height = 9,
       units = "in",
       limitsize = FALSE,
       dpi = 600)


Q32 <- read.table(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/admixture/s12345/ready/",Qlist[32],sep = ""))
k <- length(Q32)
for(j in 1:nrow(Q32)){
  if(rowSums(Q32[j,1:k] > 0.8) == 0){
    Q32$Panc[j] <- "admix"
  }
  else{
    Q32$Panc[j] <- which.max(Q32[j,1:k])
  }
  Q32$maxV[j] <- max(Q32[j,1:k])
}
Q32 <- cbind(Q32,lineageTab)  
rownames(Q32) <- Q32$ID
#temp <- temp[order(temp$clade),] # order by lineage
Q32 <- Q32[order(Q32$Lineage,Q32$Panc, Q32$maxV),] # order by % ancestry
lineageOrder <- Q32$Lineage
Q32$Lineage <- gsub(pattern = "Mosaic",replacement = "Admixed",x = Q32$Lineage)
Q322 <- Q32[,1:k]
keys <- rownames(Q322)
keyDF <- data.frame(key = keys, weigth = 1:length(keys))
melteddf <- melt(t(Q322))
head(melteddf)

Q32plot <- ggplot(melteddf) + 
    geom_bar(aes(x=Var2,y=value,color=Var1,fill=Var1),stat="identity", show.legend = FALSE) +
    geom_vline(xintercept = popIndices, size=.4, linetype = "dotted") +
    scale_color_manual(values = c(cladecolors,col_list))+
    scale_fill_manual(values = c(cladecolors,col_list)) +
    scale_x_discrete(expand = c(0,0), breaks = Q32$ID[popMean], labels = Q32$Lineage[popMean])+
    scale_y_continuous(expand = c(0,0))+
    #ggtitle(paste("k = ",k))+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = 12),
          axis.title.x = element_blank(), 
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm")) 
  #plot.title = element_text(size = 8,hjust = 0.5,margin = unit(c(0,0,0,0),"cm")))
  #theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())  
Q32plot

k <- length(Q32)
for(j in 1:nrow(Q32)){
  if(rowSums(Q32[j,1:k] > 0.8) == 0){
    Q32$Panc[j] <- "admix"
  }
  else{
    Q32$Panc[j] <- which.max(Q32[j,1:k])
  }
  Q32$maxV[j] <- max(Q32[j,1:k])
}

Q32 <- cbind(Q32,lineageTab)  
rownames(Q32) <- Q32$ID
nrow(Q32[Q32$Panc == "admix",])
Q32$admix <- ifelse(Q32$Panc == "admix", "Yes", "No")
AdmixInfo <- Q32[,c(35:37)]
admixedStrains <- Q32[Q32$Panc == "admix",]
admixedDF <- dipGainOnly[dipGainOnly$basename %in% rownames(admixedStrains),]
admixedDF <- admixedDF[!duplicated(admixedDF$basename),]

colnames(AdmixInfo)[1] <- "basename"
#AdmixInfo$basename <- gsub("SA_1_5_","SA_1_5", AdmixInfo$basename)
#AdmixInfo$basename <- gsub("Dji2_2A_a_","Dji2_2A_a", AdmixInfo$basename)
#AdmixInfo$basename <- gsub("N134_7_1_a_", "N134_7_1_a", AdmixInfo$basename)

dipGainOnly <- merge(dipGainOnly, AdmixInfo, by = "basename")
AdmixedStrains <- dipGainOnly[dipGainOnly$admix == "Yes",]
write.table(AdmixedStrains, "~/Documents/GitHub/eduardo/Dissertation/ScopelThesis/STs/C4ST3.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,na = "NA",sep = "\t")
AnAdmix <- matrix(c(55,238,111,498), nrow = 2, dimnames = list(c("Aneuploid","Euploid"),c("Admixed","non-admixed")))
fisher.test(AnAdmix)

nuc_div <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/nuc_diversity.txt")
nuc_div$pi <- nuc_div$V2/12100000
