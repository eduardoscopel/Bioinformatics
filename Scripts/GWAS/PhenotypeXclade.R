library(data.table)
library(multtest)
library(ggplot2)

### Load data frame with diploids only
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
  #print(first[[1]])
  if(!(first %in% c(NA, "No", "LOWCOV"))){
    for(j in first){
      second <- strsplit(j, split = "@")
      for(k in second){
        third <- strsplit(k, split = ",")
        type <- third[[1]]
   #     print(third)
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

#df

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
df <- df[df$aneuploidy_type == "Gain" | df$aneuploidy_type == "No",]
nrow(df)
df[is.na(df$ecology_category),"ecology_category"] <- "Unknown" 

### Ecology Analys
# Set up eco Vs aneuploid matrix
ecoMatrix <- as.matrix(table(df$ecology_category, df$aneuploidy_binary), 
                       nrow = length(unique(df$ecology_category), 
                                     header = TRUE))
ecoPMatrix <- as.matrix(table(df$ecology_category, df$ploidyCons), 
                        nrow = length(unique(df$ecology_category), 
                                      header = TRUE))
colnames(ecoMatrix) <- c("euploids", "aneuploids")

ecoMatrix

# Identify proportion of aneuploids per clade
Tot <- c()
AnProp <- c()
for(i in 1:nrow(ecoMatrix)){
  if(i == 1){
    Tot <- c(sum(ecoMatrix[i,]))
    AnProp <- c((ecoMatrix[i,2])/(Tot[i]))
  }
  else{
    Tot <- c(Tot, sum(ecoMatrix[i,]))
    AnProp <- c(AnProp,(ecoMatrix[i,2])/(Tot[i]))
  }
}

# Attach total and proportion of aneuploidy to clade V aneuploidy matrix
ecoMatrix <- cbind(ecoMatrix, Tot)
ecoMatrix <- cbind(ecoMatrix, AnProp)
#cladeMatrix <- cladeMatrix[cladeMatrix[,3]>8,] # Exclude categories with 8 or less strains

# set up matrices for pairwise Fisher's tests
OddsRatio <- matrix(
  c(rep(0,nrow(ecoMatrix))),
  nrow = 1, 
  dimnames = list(c("Aneuploids"),rownames(ecoMatrix)))
pValues <- matrix(
  c(rep(0,nrow(ecoMatrix))),
  nrow = 1, 
  dimnames = list(c("Aneuploids"),rownames(ecoMatrix)))
FDR <- matrix(
  c(rep(0,nrow(ecoMatrix))),
  nrow = 1, 
  dimnames = list(c("Aneuploids"),rownames(ecoMatrix)))

# Run Fisher's tests
for(i in 1:nrow(ecoMatrix)){
  OddsRatio[1,i] <- as.numeric(fisher.test(
    matrix(c(ecoMatrix[i,2], 
             ecoMatrix[i,1],
             sum(ecoMatrix[,2])-ecoMatrix[i,2],
             sum(ecoMatrix[,1])-ecoMatrix[i,1]),
           nrow = 2,
           dimnames = list(c("Aneuploids", "Euploids"),
                           c(rownames(ecoMatrix)[i], "All"))))
    $estimate)
  pValues[1,i] <- as.numeric(fisher.test(
    matrix(c(ecoMatrix[i,2], 
             ecoMatrix[i,1],
             sum(ecoMatrix[,2])-ecoMatrix[i,2],
             sum(ecoMatrix[,1])-ecoMatrix[i,1]),
           nrow = 2,
           dimnames = list(c("Aneuploids", "Euploids"),
                           c(rownames(ecoMatrix)[i], "All"))))
    $p.value)
}

# set up variables for eco proportion plot
pLabels = c()

# Prepare DFs to create odds ratio plot
for(row in 1:nrow(ecoMatrix)){
  pLabels <- append(pLabels, 
                    paste(rownames(ecoMatrix)[row], " (", ecoMatrix[row,3],")",sep = ""))
}

ecoDF <- data.frame(xAxis = seq(1, nrow(ecoMatrix),1),
                    eco = rownames(ecoMatrix),
                    prop = ecoMatrix[,4],
                    odds = as.numeric(OddsRatio[1,]),
                    p = as.numeric(pValues[1,]))

# Adjust p-values for FDR using the BH multiple test correction
FDRraw <- mt.rawp2adjp(ecoDF$p, proc = "BH", alpha = 0.05)
index <- FDRraw$index
FDR <- FDRraw$adjp[,2]
ecoDF$FDR[index] <- FDR

# add significance and color labels
ecoDF$sig <- ifelse(ecoDF$FDR <= 0.01, "**", 
                    ifelse(ecoDF$FDR <= 0.05, "*", 
                           ifelse(ecoDF$p <= 0.05, "+", "")))

ecoDF$sigP <- ifelse(ecoDF$p <= 0.01, "**", 
                     ifelse(ecoDF$p <= 0.05, "*", 
                            ifelse(ecoDF$FDR <= 0.11, "+", "")))

ecoDF$xAxis <- factor(ecoDF$xAxis)
ecoDF$lab <- pLabels
avg <- sum(ecoMatrix[,2])/sum(ecoMatrix[,3])
ecoDF$yPos <- ifelse(ecoDF$prop <= avg, avg,
                     ecoDF$prop)
ecoDF$face <- ifelse(ecoDF$sig == "", "plain", "bold")
ecoDF$color <- ifelse(ecoDF$odds > 1 & ecoDF$sig != "", "darkorange",
                      ifelse(ecoDF$odds < 1 & ecoDF$sig != "", "skyblue", "gray70"))

#ecoPlot <- 
ecoPlot <- ggplot(ecoDF, aes(x = lab, y = prop)) + 
  geom_col(fill = ecoDF$color) +
  geom_text(aes(label = paste(lab, sig),y = yPos, fontface = face), 
            angle = 90, 
            na.rm = TRUE, 
            hjust = -0.05,
            size = 4.32) + 
  geom_abline(slope = 0, intercept = sum(ecoMatrix[,2])/sum(ecoMatrix[,3]), color = "gray20", linetype = "dashed") +
  scale_y_continuous(expand = c(0,0), n.breaks = 6, limits = c(0,1.16)) +
  theme_bw()+
  ylab("Proportion of strains with chr gain")+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

ecoPlot
ggsave("../../../../Papers_Scopel/GWAS/diploids/gainOnly/Figure1A.jpg",
       plot = ecoPlot,
       device = "jpg", 
       dpi = 600, 
       limitsize = FALSE,
       width = 5,
       height = 5,
       units = "in")

### Eco Ploidy analysis
df$ploidySum <- ifelse(df$ploidyCons == 1, "haploid", 
                       ifelse(df$ploidyCons == 2, "diploid", 
                              ifelse(df$ploidyCons == "Unknown", "Unknown", "polyploid")))
ecoPMatrix <- as.matrix(table(df$ecology_category, df$ploidySum), 
                        nrow = length(unique(df$ecology_category), 
                                      header = TRUE))
ecoPMatrix <- ecoPMatrix[,c(2,1,3)] # Exclude strains with unknown ploidy
ecoPMatrix <- ecoPMatrix[-6,]
ecoPMatrix <- ecoPMatrix[-9,]


# Identify proportion of aneuploids per clade
Tot <- c()
PolyProp <- c()
for(i in 1:nrow(ecoPMatrix)){
  if(i == 1){
    Tot <- c(sum(ecoPMatrix[i,]))
    PolyProp <- c((ecoPMatrix[i,3])/(Tot[i]))
  }
  else{
    Tot <- c(Tot, sum(ecoPMatrix[i,]))
    PolyProp <- c(PolyProp,(ecoPMatrix[i,3])/(Tot[i]))
  }
}

# Attach total and proportion of aneuploidy to clade V aneuploidy matrix
ecoPMatrix <- cbind(ecoPMatrix, Tot)
ecoPMatrix <- cbind(ecoPMatrix, PolyProp)
#cladeMatrix <- cladeMatrix[cladeMatrix[,3]>8,] # Exclude categories with 8 or less strains


# set up matrices for pairwise Fisher's tests
OddsRatio <- matrix(
  c(rep(0,nrow(ecoPMatrix))),
  nrow = 1, 
  dimnames = list(c("Polyploids"),rownames(ecoPMatrix)))
pValues <- matrix(
  c(rep(0,nrow(ecoPMatrix))),
  nrow = 1, 
  dimnames = list(c("Polyploids"),rownames(ecoPMatrix)))
FDR <- matrix(
  c(rep(0,nrow(ecoPMatrix))),
  nrow = 1, 
  dimnames = list(c("Polyploids"),rownames(ecoPMatrix)))

# Run Fisher's tests
for(i in 1:nrow(ecoPMatrix)){
  OddsRatio[1,i] <- as.numeric(fisher.test(
    matrix(c(ecoPMatrix[i,3], 
             ecoPMatrix[i,1]+ecoPMatrix[i,2],
             sum(ecoPMatrix[,3])-ecoPMatrix[i,3],
             sum(ecoPMatrix[,1])+sum(ecoPMatrix[,2])-ecoPMatrix[i,1]+ecoPMatrix[i,2]),
           nrow = 2,
           dimnames = list(c("Polyploids", "diploids+haploids"),
                           c(rownames(ecoPMatrix)[i], "All"))))
    $estimate)
  pValues[1,i] <- as.numeric(fisher.test(
    matrix(c(ecoPMatrix[i,3], 
             ecoPMatrix[i,1]+ecoPMatrix[i,2],
             sum(ecoPMatrix[,3])-ecoPMatrix[i,3],
             sum(ecoPMatrix[,1])+sum(ecoPMatrix[,2])-ecoPMatrix[i,1]+ecoPMatrix[i,2]),
           nrow = 2,
           dimnames = list(c("Polyploids", "diploids+haploids"),
                           c(rownames(ecoPMatrix)[i], "All"))))
    $p.value)
}

# set up variables for eco proportion plot
pLabels = c()

# Prepare DFs to create odds ratio plot
for(row in 1:nrow(ecoPMatrix)){
  pLabels <- append(pLabels, 
                    paste(rownames(ecoPMatrix)[row], " (", ecoPMatrix[row,4],")",sep = ""))
}

ecoPDF <- data.frame(xAxis = seq(1, nrow(ecoPMatrix),1),
                     eco = rownames(ecoPMatrix),
                     prop = ecoPMatrix[,5],
                     odds = as.numeric(OddsRatio[1,]),
                     p = as.numeric(pValues[1,]))

# Adjust p-values for FDR using the BH multiple test correction
FDRraw <- mt.rawp2adjp(ecoPDF$p, proc = "BH", alpha = 0.05)
index <- FDRraw$index
FDR <- FDRraw$adjp[,2]
ecoPDF$FDR[index] <- FDR

# add significance and color labels
ecoPDF$sig <- ifelse(ecoPDF$FDR <= 0.01, "**", 
                     ifelse(ecoPDF$FDR <= 0.05, "*", 
                            ifelse(ecoPDF$p <= 0.05, "+", "")))

ecoDF$sigP <- ifelse(ecoDF$p <= 0.01, "**", 
                     ifelse(ecoDF$p <= 0.05, "*", 
                            ifelse(ecoDF$FDR <= 0.11, "+", "")))

ecoPDF$xAxis <- factor(ecoPDF$xAxis)
ecoPDF$lab <- pLabels
avgP <- sum(ecoPMatrix[,3])/sum(ecoPMatrix[,4])
ecoPDF$yPos <- ifelse(ecoPDF$prop <= avg, avg,
                      ecoPDF$prop)
ecoPDF$face <- ifelse(ecoPDF$sig == "", "plain", "bold")
ecoPDF$color <- ifelse(ecoPDF$odds > 1 & ecoPDF$sig != "", "darkorange",
                       ifelse(ecoPDF$odds < 1 & ecoPDF$sig != "", "skyblue", "gray70"))

#ecoPlot <- 
ecoPPlot <- ggplot(ecoPDF, aes(x = lab, y = prop)) + 
  geom_col(fill = ecoPDF$color) +
  geom_text(aes(label = paste(lab, sig),y = yPos, fontface = face), 
            angle = 90, 
            na.rm = TRUE, 
            hjust = -0.05,
            size = 4.32) + 
  geom_abline(slope = 0, intercept = avgP, color = "gray20", linetype = "dashed") +
  scale_y_continuous(expand = c(0,0), n.breaks = 6, limits = c(0,1.16)) +
  theme_bw()+
  ylab("Proportion of polyploids")+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

ecoPPlot
ecoGridPlot <- plot_grid(ecoPlot, ecoPPlot,labels = c("(A)","(B)"),label_size = 12,label_x = -0.03)
ggsave("../../Dissertation/ScopelThesis/figures/C4F1.png",
       plot = ecoGridPlot,
       device = "png", 
       dpi = 600, 
       limitsize = FALSE,
       width = 9,
       height = 6,
       units = "in")

### Select diploids only
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

### Ecology Analys (diploids only)
# Set up eco Vs aneuploid matrix
ecoDipMatrix <- as.matrix(table(dipGainOnly$ecology_category, dipGainOnly$aneuploidy_binary), 
                          nrow = length(unique(dipGainOnly$ecology_category), 
                                        header = TRUE))
colnames(ecoDipMatrix) <- c("euploids", "aneuploids")

ecoDipMatrix

# Identify proportion of aneuploids per clade
Tot <- c()
AnProp <- c()
for(i in 1:nrow(ecoDipMatrix)){
  if(i == 1){
    Tot <- c(sum(ecoDipMatrix[i,]))
    AnProp <- c((ecoDipMatrix[i,2])/(Tot[i]))
  }
  else{
    Tot <- c(Tot, sum(ecoDipMatrix[i,]))
    AnProp <- c(AnProp,(ecoDipMatrix[i,2])/(Tot[i]))
  }
}

# Attach total and proportion of aneuploidy to clade V aneuploidy matrix
ecoDipMatrix <- cbind(ecoDipMatrix, Tot)
ecoDipMatrix <- cbind(ecoDipMatrix, AnProp)
#cladeMatrix <- cladeMatrix[cladeMatrix[,3]>8,] # Exclude categories with 8 or less strains

# set up matrices for pairwise Fisher's tests
OddsRatio <- matrix(
  c(rep(0,nrow(ecoDipMatrix))),
  nrow = 1, 
  dimnames = list(c("Aneuploids"),rownames(ecoDipMatrix)))
pValues <- matrix(
  c(rep(0,nrow(ecoDipMatrix))),
  nrow = 1, 
  dimnames = list(c("Aneuploids"),rownames(ecoDipMatrix)))
FDR <- matrix(
  c(rep(0,nrow(ecoDipMatrix))),
  nrow = 1, 
  dimnames = list(c("Aneuploids"),rownames(ecoDipMatrix)))

# Run Fisher's tests
for(i in 1:nrow(ecoDipMatrix)){
  OddsRatio[1,i] <- as.numeric(fisher.test(
    matrix(c(ecoDipMatrix[i,2], 
             ecoDipMatrix[i,1],
             sum(ecoDipMatrix[,2])-ecoDipMatrix[i,2],
             sum(ecoDipMatrix[,1])-ecoDipMatrix[i,1]),
           nrow = 2,
           dimnames = list(c("Aneuploids", "Euploids"),
                           c(rownames(ecoDipMatrix)[i], "All"))))
    $estimate)
  pValues[1,i] <- as.numeric(fisher.test(
    matrix(c(ecoDipMatrix[i,2], 
             ecoDipMatrix[i,1],
             sum(ecoDipMatrix[,2])-ecoDipMatrix[i,2],
             sum(ecoDipMatrix[,1])-ecoDipMatrix[i,1]),
           nrow = 2,
           dimnames = list(c("Aneuploids", "Euploids"),
                           c(rownames(ecoDipMatrix)[i], "All"))))
    $p.value)
}

# set up variables for eco proportion plot
pLabels = c()

# Prepare DFs to create odds ratio plot
for(row in 1:nrow(ecoDipMatrix)){
  pLabels <- append(pLabels, 
                    paste(rownames(ecoDipMatrix)[row], " (", ecoDipMatrix[row,3],")",sep = ""))
}

ecoDipDF <- data.frame(xAxis = seq(1, nrow(ecoDipMatrix),1),
                       eco = rownames(ecoDipMatrix),
                       prop = ecoDipMatrix[,4],
                       odds = as.numeric(OddsRatio[1,]),
                       p = as.numeric(pValues[1,]))

# Adjust p-values for FDR using the BH multiple test correction
FDRraw <- mt.rawp2adjp(ecoDipDF$p, proc = "BH", alpha = 0.05)
index <- FDRraw$index
FDR <- FDRraw$adjp[,2]
ecoDipDF$FDR[index] <- FDR

# add significance and color labels
ecoDipDF$sig <- ifelse(ecoDipDF$FDR <= 0.01, "**", 
                       ifelse(ecoDipDF$FDR <= 0.05, "*", 
                              ifelse(ecoDipDF$p <= 0.05, "+", "")))

ecoDipDF$sigP <- ifelse(ecoDipDF$p <= 0.01, "**", 
                        ifelse(ecoDipDF$p <= 0.05, "*", 
                               ifelse(ecoDipDF$FDR <= 0.11, "+", "")))

ecoDipDF$xAxis <- factor(ecoDipDF$xAxis)
ecoDipDF$lab <- pLabels
avgD <- sum(ecoDipMatrix[,2])/sum(ecoDipMatrix[,3])
ecoDipDF$yPos <- ifelse(ecoDipDF$prop <= avg, avg,
                        ecoDipDF$prop)
ecoDipDF$face <- ifelse(ecoDipDF$sig == "", "plain", "bold")
ecoDipDF$color <- ifelse(ecoDipDF$odds > 1 & ecoDipDF$sig != "", "darkorange",
                         ifelse(ecoDipDF$odds < 1 & ecoDipDF$sig != "", "skyblue", "gray70"))

#ecoPlot <- 
ecoDipPlot <- ggplot(ecoDipDF, aes(x = lab, y = prop)) + 
  geom_col(fill = ecoDipDF$color) +
  geom_text(aes(label = paste(lab, sig),y = yPos, fontface = face), 
            angle = 90, 
            na.rm = TRUE, 
            hjust = -0.05,
            size = 4.32) + 
  geom_abline(slope = 0, intercept = sum(ecoDipMatrix[,2])/sum(ecoDipMatrix[,3]), color = "gray20", linetype = "dashed") +
  scale_y_continuous(expand = c(0,0), n.breaks = 6, limits = c(0,0.8)) +
  theme_bw()+
  ylab("Proportion of diploids with chr gain")+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

ecoDipPlot
ggsave("../../../../Papers_Scopel/GWAS/diploids/gainOnly/Figure1A.png",
       plot = ecoPlot,
       device = "png", 
       dpi = 600, 
       limitsize = FALSE,
       width = 5,
       height = 5,
       units = "in")

topRow <- plot_grid(ecoPlot, ecoPPlot, labels = c("(A)","(B)"), label_size = 12, nrow = 1, label_x = -0.005, label_y = 1)
ecoGridPlot <- plot_grid(topRow, ecoDipPlot, labels = c(NA, "(C)"),label_size = 12,align = "h",nrow = 2 , label_x = -.01, label_y = 1)
ggsave("../../Dissertation/ScopelThesis/figures/C4F1v2.png",
       plot = ecoGridPlot,
       device = "png", 
       dpi = 600, 
       limitsize = FALSE,
       width = 10,
       height = 9.4,
       units = "in")


### Clade Analys in diploids
# Set up clade V aneuploid matrix
cladeMatrix <- as.matrix(table(dipGainOnly$lineage, dipGainOnly$aneuploidy_binary), 
                         nrow = length(unique(dipGainOnly$lineage), 
                                       header = TRUE))
colnames(cladeMatrix) <- c("euploids", "aneuploids")

cladeMatrix

# Identify proportion of aneuploids per clade
Tot <- c()
AnProp <- c()
for(i in 1:nrow(cladeMatrix)){
  if(i == 1){
    Tot <- c(sum(cladeMatrix[i,]))
    AnProp <- c((cladeMatrix[i,2])/(Tot[i]))
  }
  else{
    Tot <- c(Tot, sum(cladeMatrix[i,]))
    AnProp <- c(AnProp,(cladeMatrix[i,2])/(Tot[i]))
  }
}

# Attach total and proportion of aneuploidy to clade V aneuploidy matrix
cladeMatrix <- cbind(cladeMatrix, Tot)
cladeMatrix <- cbind(cladeMatrix, AnProp)
#cladeMatrix <- cladeMatrix[cladeMatrix[,3]>8,] # Exclude categories with 8 or less strains

# set up matrices for pairwise Fisher's tests
OddsRatio <- matrix(
  c(rep(0,nrow(cladeMatrix))),
                   nrow = 1, 
                   dimnames = list(c("Aneuploids"),rownames(cladeMatrix)))
pValues <- matrix(
  c(rep(0,nrow(cladeMatrix))),
  nrow = 1, 
  dimnames = list(c("Aneuploids"),rownames(cladeMatrix)))
FDR <- matrix(
  c(rep(0,nrow(cladeMatrix))),
  nrow = 1, 
  dimnames = list(c("Aneuploids"),rownames(cladeMatrix)))

# Run Fisher's tests
for(i in 1:nrow(cladeMatrix)){
  OddsRatio[1,i] <- as.numeric(fisher.test(
    matrix(c(cladeMatrix[i,2], 
             cladeMatrix[i,1],
             sum(cladeMatrix[,2])-cladeMatrix[i,2],
             sum(cladeMatrix[,1])-cladeMatrix[i,1]),
           nrow = 2,
           dimnames = list(c("Aneuploids", "Euploids"),
                           c(rownames(cladeMatrix)[i], "All"))))
    $estimate)
  pValues[1,i] <- as.numeric(fisher.test(
    matrix(c(cladeMatrix[i,2], 
             cladeMatrix[i,1],
             sum(cladeMatrix[,2])-cladeMatrix[i,2],
             sum(cladeMatrix[,1])-cladeMatrix[i,1]),
           nrow = 2,
           dimnames = list(c("Aneuploids", "Euploids"),
                           c(rownames(cladeMatrix)[i], "All"))))
    $p.value)
}

# set up variables for clade proportion plot
pLabels = c()

# Prepare DFs to create odds ratio plot
for(row in 1:nrow(cladeMatrix)){
    pLabels <- append(pLabels, 
                      paste(rownames(cladeMatrix)[row], " (", cladeMatrix[row,3],")",sep = ""))
}

cladeDF <- data.frame(xAxis = seq(1, nrow(cladeMatrix),1),
                      clade = rownames(cladeMatrix),
                      prop = cladeMatrix[,4],
                      odds = as.numeric(OddsRatio[1,]),
                      p = as.numeric(pValues[1,]))

# Adjust p-values for FDR using the BH multiple test correction
FDRraw <- mt.rawp2adjp(cladeDF$p, proc = "BH", alpha = 0.05)
index <- FDRraw$index
FDR <- FDRraw$adjp[,2]
cladeDF$FDR[index] <- FDR

# add significance and color labels
cladeDF$sig <- ifelse(cladeDF$FDR <= 0.01, "**", 
                  ifelse(cladeDF$FDR <= 0.05, "*", 
                         ifelse(cladeDF$p <= 0.05, "+", "")))
                  
cladeDF$sigP <- ifelse(cladeDF$p <= 0.01, "**", 
                      ifelse(cladeDF$p <= 0.05, "*", 
                             ifelse(cladeDF$FDR <= 0.11, "+", "")))

cladeDF$xAxis <- factor(cladeDF$xAxis)
cladeDF$lab <- pLabels
avg <- sum(cladeMatrix[,2])/sum(cladeMatrix[,3])
cladeDF$yPos <- ifelse(cladeDF$prop <= avg, avg,
                       cladeDF$prop)
cladeDF$face <- ifelse(cladeDF$sig == "", "plain", "bold")
cladeDF$color <- ifelse(cladeDF$odds > 1 & cladeDF$sig != "", "darkorange",
                        ifelse(cladeDF$odds < 1 & cladeDF$sig != "", "skyblue", "gray70"))

#CladePlot <- 
cladeDF[7,]$lab <- "Bioethanol (30)"
cladeDF[8,]$lab <- "CHN1,2,3 (7)"
cladeDF[11,]$lab <- "Ecuador2 (9)"
cladeDF[19,]$lab <- "EU Oak (8)"
cladeDF[22,]$lab <- "Unassigned (87)"
cladeDF[27,]$lab <- "Taiwan1 (3)"

fisherPlot <- ggplot(cladeDF, aes(x = lab, y = prop)) + 
  geom_col(fill = cladeDF$color) +
  geom_text(aes(label = paste(lab, sig),y = yPos, fontface = face), 
            angle = 90, 
            na.rm = TRUE, 
            hjust = -0.05,
            size = 4.32) + 
  geom_abline(slope = 0, intercept = sum(cladeMatrix[,2])/sum(cladeMatrix[,3]), color = "gray20", linetype = "dashed") +
  scale_y_continuous(expand = c(0,0), n.breaks = 7, limits = c(0,0.8)) +
  theme_bw()+
  ylab("Proportion of strains with chr gain")+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

fisherPlot
ggsave("../../Dissertation/ScopelThesis/figures/C4F4v3.png",
       plot = fisherPlot,
       device = "png", 
       dpi = 600, 
       limitsize = FALSE,
       width = 10,
       height = 6,
       units = "in")




