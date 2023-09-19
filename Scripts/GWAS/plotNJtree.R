library(ggtree)
library(ggplot2)
library(ape)
library(treeio)
library(dplyr)
library(ggrepel)
library(cowplot)
library(ggnewscale)
library(data.table)
library(multtest)

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
#dipGainOnly[dipGainOnly$basename == "SA_1_5_","basename"] <- "SA_1_5"
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
#dipGainOnly[dipGainOnly$basename == "N134_7_1_a_","basename"] <- "N134_7_1_a"
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
#dipGainOnly[dipGainOnly$basename == "Dji2_2A_a_","basename"] <- "Dji2_2A_a"
dipGainOnly[dipGainOnly$basename %in% c("Dji2_2A_a"), "lineage"] <- "African Wine"

# Convert lineages to North American Oak
dipGainOnly[dipGainOnly$basename %in% c("CLIB414","ATCC_66348_1D","YPS670","N95_5_1A","2163","IY_03_5_30_1_1_1_1","IY_03_5_26_5_1_1_1","N95_6_1C"), "lineage"] <- "North American Oak"

# Convert lineages to Taiwan2
dipGainOnly[dipGainOnly$basename %in% c("ES2M03_7A","ES2M03","ES4M07", "SJ5L12"), "lineage"] <- "Taiwan2"

nrow(dipGainOnly)
dipGainOnly[duplicated(dipGainOnly$basename),]
dipGainOnly <- dipGainOnly[!duplicated(dipGainOnly$basename),]
nrow(dipGainOnly)

#write.table(dipGainOnly, "~/Documents/GitHub/eduardo/Dissertation/ScopelThesis/STs/C4ST2.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,na = "NA",sep = "\t")

## Generate ind file for PCA
indOrder <- read.table("../../../../Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/orderIDs.txt")
colnames(indOrder) <- "ID"
indTab <- data.frame(ID = dipGainOnly$basename, Sex = "Unknown", Lineage = dipGainOnly$lineage)
#indTab$ID <- gsub("SA_1_5","SA_1_5_", indTab$ID)
#indTab$ID <- gsub("Dji2_2A_a","Dji2_2A_a_", indTab$ID)
#indTab$ID <- gsub("N134_7_1_a", "N134_7_1_a_", indTab$ID)
indTab2 <- merge(indOrder,indTab, by = "ID",sort = FALSE)
indTab2$Lineage <- gsub(" ","_",indTab2$Lineage)
#write.table(indTab2, "../../../../Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/gain2n_gSNPs.ind", quote = FALSE, row.names = FALSE, col.names = FALSE,sep = " ")



#tree <- read.tree("/Users/es47540/Documents/Papers_Scopel/GWAS/diploids/alignment/newMegaNJtreeJC_BS.nwk")
tree <- read.tree("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/SNPeff/newMegaNJtreeJC_BS_correctIDs.nwk")
tree$node.label <- as.numeric(tree$node.label)*100
rtree <- root(tree, c("EM14S01_3B", "EN14S01", "GE14S01_7B"))
newdf <- data.frame(matrix(NA, nrow = nrow(dipGainOnly),ncol = 16))
rownames(newdf) <- rtree$tip.label
chr <- seq(1,16)
chr <- paste("c", chr, sep = "")
colnames(newdf) <- chr

for(i in 1:nrow(newdf)){
  for(j in 1:nrow(dipGainOnly)){
    if(rownames(newdf)[i] == dipGainOnly$basename[j]){
      newdf$aneuploidy_binary[i] <- dipGainOnly$aneuploidy_binary[j]
      newdf$lineage[i] <- dipGainOnly$lineage[j]
      newdf[i,c(1:16)] <- dipGainOnly[j,c(27:42)]
      newdf$admix[i] <- dipGainOnly$admix[j]
    }
  }
}

anBin <- as.matrix(newdf)[,17]
lineage <- as.matrix(newdf)[,18]
anChr <- as.matrix(newdf)[,c(1:16)]
admix <- as.matrix(newdf)[,19]


traitsdf <- data.frame(node = nodeid(rtree, names(anBin)), aneuploidy_binary = anBin, lineage = lineage, chr = anChr, admix = admix)
y <- full_join(rtree, traitsdf, by = "node")

# Group root nodes
#EM14S01_3Bnode <- which(rownames(newdf) == "EM14S01_3B")
#EN14S01node <- which(rownames(newdf) == "EN14S01")
#GE14S01_7Bnode<- which(rownames(newdf) == "GE14S01_7B")
#m<- MRCA(rtree, EM14S01_3Bnode, EN14S01node, GE14S01_7Bnode)
#y <- groupClade(y, m)
#y <- full_join(y, traitsdf, by ='node')
darkred <- "#8B0000"
red <- "#FF0000"
darksalmon <- "#E9967A"
darkorange <- "#FF8C00"
darkgolden <- "#FF8C00"
darkkhaki <- "#BDB76B"
olive <- "#808000"
yellow <- "#ffed6f"
lightyellow <- "#ffff99"
darkgreen <- "#006400"
forestgreen <- "#228B22"
palegreen <- "#98FB98"
mediumgreen <- "#3CB371"
darkgray <- "#2F4F4F"
teal <- "#008080"
aqua <- "#00FFFF"
darkblue <- "#00008B"
blue <- "#0000FF"
lightblue <- "#00BFFF"
darkviolet <- "#9400D3"
indigo <- "#4B0082"
purple <- "#800080"
lightpink <- "#EE82EE"
mediumpink <- "#C71585"
pink <- "#FF1493"
brown <- "#8B4513"
tan <- "#D2B48C"

we1 <- purple
ol2 <- olive
bb3 <- indigo
mo4 <- forestgreen
fd5 <- lightblue
ab6 <- "#bf812d"
mb7 <- "black"
mo8 <- pink
ma9 <- brown 
fgh10 <- darkorange
ab11 <- darkkhaki
wac12 <- tan
apw13 <- darkviolet
ch14 <- palegreen
ch15 <- mediumgreen
ch16 <- darkgreen
t17 <- teal
fea18 <- aqua
m19 <- lightblue
ch20 <- darkgray
e21 <- blue
fer22 <- darkblue
nao23 <- darkkhaki
ai24 <- darksalmon
s25 <- red
af26 <- darkred
israel <- "darkblue"
ecuador1 <- "grey26"
clinical <- "coral2"
siberia <- "darkslategray4"
ethiopia <- "yellowgreen"
taiwan2 <- "springgreen3"
malaysian <- "gray70"

cladecolors <- c(ab11, apw13, ab6, ol2, af26, ai24, bb3, ch14, clinical, ecuador1, e21, ethiopia, fea18, fer22, fd5, fgh10, 
                 israel, malaysian, mo4, ma9, mo8, "black", nao23, s25, siberia, darksalmon, t17, taiwan2, wac12, we1)

ptre <- ggtree(y,size=0.2, layout="rectangular",branch.length = 'none') + 
  geom_treescale() +
  geom_tiplab(size=0.5,aes(label = paste(label,lineage,sep = " "), color=lineage),align = FALSE,linesize = 0.1, show.legend = FALSE)+
  scale_colour_manual(values = cladecolors) + 
  guides(colour = guide_legend(override.aes = list(size=3, label = "")))+
  geom_point2(aes(subset = as.numeric(label) >= 70, label = label), size=1.5, alpha = 0.5,color = "black")
# Finda MRCA for each clade 
wine1 <- MRCA(ptre, "CLI_6","DBVPG4297")
wine2 <- MRCA(ptre, "A_12","WI010")
alpechin <- MRCA(ptre, "UCD_11_601","CBS2909")
wine3 <- MRCA(ptre, "Lib71","CBS3012")
wine4 <- MRCA(ptre, "CLIB1085","DBVPG1145")
bb <- MRCA(ptre, "SM_9_1_AL1", "SP004")
israelNode <- MRCA(ptre, "35","CLIB1077")
eq1 <- MRCA(ptre, "CLQCA_10_027","CLQCA_17_060")
alebeer <- MRCA(ptre, "YJM1100","BE002")
frenchdairy <- MRCA(ptre, "CLIB553","CBS2421")
ethiopiaNode <- MRCA(ptre, "ETPN7","ETPF2")
africanbeer <- MRCA(ptre, "SAN18","NIGF2")
frenchguianahuman <- MRCA(ptre, "HE005", "CLQCA_20_060")
mexicanagave <- MRCA(ptre, "LCBG_3Y8","LCBG_3D2")
medoak <- MRCA(ptre, "STG_2","2162")
mixedorigin <- MRCA(ptre, "CLIB324_2","SBE_1C")
clinicalNode <- MRCA(ptre, "YJM339", "YJM1122")
westAfrica <- MRCA(ptre, "MTF2556", "MTF2546")
siberiaNode <- MRCA(ptre, "N39_7A","N38_4A")
AsianFerm1 <- MRCA(ptre, "SAF29","ZJH7")
sakeNode <- MRCA(ptre, "RIB1006","CBS2270")
AsianFerm2 <- MRCA(ptre, "HBU1","GZU1")
AsianFerm3 <- MRCA(ptre, "SCN5","JSN1")
AsianFerm4 <- MRCA(ptre, "SAN33","JSN3")
AsianFerm5 <- MRCA(ptre, "SXQ7","SXQ2")
AsianFerm6 <- MRCA(ptre, "SCN6","SCN8")
SouthAfrica <- MRCA(ptre, "SAN8", "SAN9")
AfricanWine <- MRCA(ptre, "NPA03_1","DJ74")
CHN1 <- MRCA(ptre, "HN19", "HN10")
naoak <- MRCA(ptre, "YPS154","N95_6_1C")
AsianIslands <- MRCA(ptre, "S8BM_30_2D", "CBS1592")
taiwan2Node <- MRCA(ptre,  "ES2M03_7A", "SJ5L12")
farEastRussian <- MRCA(ptre, "N22_00_5D","N3_00_7A")
ecuador2 <- MRCA(ptre, "CLQCA_20_246","YPS617")
CHN2 <- MRCA(ptre, "HN16", "HN15")
malaysianNode <- MRCA(ptre, "UWOPS03_459_1","UWOPS03_461_4")
farEastAsian <- MRCA(ptre, "WI016","Ksc2_2B")
CHN3 <- MRCA(ptre, "SX3", "SX1")
Taiwan <- MRCA(ptre, "EN14S01","GE14S01_7B")

ptre <- flip(ptre, sakeNode, AsianFerm1)
tipsOnly <- ptre$data[ptre$data$isTip == TRUE,]
strainOrder <- tipsOnly[order(tipsOnly$y),]$label
NJorder <- data.frame(order = seq(1,length(strainOrder)), taxa = strainOrder)

ggsave("../../../../Papers_Scopel/GWAS/diploids/gainOnly/RecNJtree.pdf",
       plot = ptre,
       device = "pdf", 
       dpi = 600, 
       limitsize = FALSE,
       width = 15,
       height = 15,
       units = "in")


circtre <- ggtree(y,size=0.4, layout="circular") + 
  geom_treescale(x=0.008, y=450, offset = 30, width = 0.01,fontsize = 10) + 
  geom_tiplab(size=1,aes(label = label, color = admix),align = TRUE,linesize = 0.02, show.legend = FALSE)+
  scale_color_manual(values = c("black","deeppink"))+
  #geom_tiplab(size=0.5,aes(label = paste(label,lineage,sep = " "), color=lineage),align = TRUE,linesize = 0.1, show.legend = FALSE)+
  #scale_colour_manual(values = cladecolors) + 
  #guides(colour = guide_legend(override.aes = list(size=3, label = "")))+
  geom_point2(aes(subset = as.numeric(label) >= 60, label = label), size = 3, alpha = 0.5,color = "black")
  #geom_text2(aes(subset = as.numeric(label) >= 70, label = label), size = 3, alpha = 0.5,color = "black")

circtre <- flip(circtre, sakeNode, AsianFerm1)

chrgaindf <- data.frame(chr_gain = traitsdf$aneuploidy_binary)
rownames(chrgaindf) <- row.names(traitsdf)
chrgaindf$chr_gain <- ifelse(chrgaindf$chr_gain == "No", NA,"Yes")
circtre <- gheatmap(circtre, chrgaindf, offset = 0.0023, width = 0.025,colnames = FALSE) + 
  scale_fill_manual(values = c("black"),na.translate = FALSE)

circtre <- circtre + new_scale_fill()

#admixDF <- data.frame(admix = traitsdf$admix)
#rownames(admixDF) <- row.names(traitsdf)
#admixDF$admix <- ifelse(admixDF$admix == "No", NA,"Yes")
#circtre <- gheatmap(circtre, admixDF, offset = 0.0031, width = 0.025, colnames = FALSE) + 
#  scale_fill_manual(values = c("deeppink"), na.translate = FALSE)

#circtre <- circtre + new_scale_fill()

circtre <- circtre + geom_hilight(node = wine1, fill = we1, alpha = 0.4, extendto = 0.032) + 
  geom_hilight(node = wine2, fill = we1, alpha = 0.4, extendto = 0.032) + 
  geom_hilight(node = alpechin, fill = ol2, alpha = 0.4, extendto = 0.032)+
  geom_hilight(node = wine3, fill = we1, alpha = 0.4, extendto = 0.032) + 
  geom_hilight(node = wine4, fill = we1, alpha = 0.4, extendto = 0.032) + 
  geom_hilight(node = bb, fill = bb3, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = israelNode, fill = israel, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = eq1, fill = ecuador1, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = alebeer, fill = ab6, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = frenchdairy, fill = fd5, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = ethiopiaNode, fill = ethiopia, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = africanbeer, fill = ab11, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = frenchguianahuman, fill = fgh10, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = mexicanagave, fill = ma9, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = medoak, fill = mo4, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = mixedorigin, fill = mo8, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = clinicalNode, fill = clinical, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = westAfrica, fill = wac12, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = siberiaNode, fill = siberia, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = AsianFerm1, fill = af26, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = sakeNode, fill = s25, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = AsianFerm2, fill = af26, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = AsianFerm3, fill = af26, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = AsianFerm4, fill = af26, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = AsianFerm5, fill = af26, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = AsianFerm6, fill = af26, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = SouthAfrica, fill = darksalmon, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = AfricanWine, fill = apw13, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = CHN1, fill = ch14, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = naoak, fill = nao23, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = AsianIslands, fill = ai24, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = taiwan2Node, fill = taiwan2, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = farEastRussian, fill = fer22, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = ecuador2, fill = e21, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = CHN2, fill = ch14, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = malaysianNode, fill = malaysian, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = farEastAsian, fill = fea18, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = CHN3, fill = ch14, alpha = 0.4, extendto = 0.032) +
  geom_hilight(node = Taiwan, fill = t17, alpha = 0.4, extendto = 0.032)

circtre <- circtre + theme(legend.position = 'none')

circtre <- circtre + geom_strip("CLI_6", "DBVPG1145",label = "European Wine", hjust = 0.5, barsize = 0, offset = 0.0035, angle = 20, fontsize = 5, color = we1, align = TRUE) +
  geom_strip("UCD_11_601","CBS2909", label = "Alpechin", barsize = 0, offset = 0.008, angle = 52, fontsize = 5, color = ol2, align = TRUE) + 
  geom_strip("SM_9_1_AL1", "SP004", label = "Bioethanol", hjust = 0.5, barsize = 0, offset = 0.0035,  angle = 128, fontsize = 5, color = bb3,  align = TRUE) +
  geom_strip("35","CLIB1077", label = "Israel", align = TRUE, color = israel,  fontsize = 5, barsize = 0,offset = 0.007, angle = 35) + 
  geom_strip("CLQCA_10_027","CLQCA_17_060", label = "Ecuador1", align = TRUE, color = ecuador1,  fontsize = 5, barsize = 0,offset = 0.0082, angle = 30) + 
  geom_strip("YJM1100","BE002", label = "Ale Beer", hjust = 0.5, align = TRUE, color = ab6,  fontsize = 5, barsize = 0,offset = 0.0035, angle = 110) + 
  geom_strip("CLIB553","CBS2421", label = "French Dairy", hjust = 0.5, align = TRUE, color = fd5,  fontsize = 5, barsize = 0,offset = 0.0035, angle = 100) + 
  geom_strip("ETPN7","ETPF2", label = "Ethiopia", hjust = 0.5, align = TRUE, color = ethiopia,  fontsize = 5, barsize = 0,offset = 0.0051, angle = 360) + 
  geom_strip("SAN18","NIGF2", label = paste("African", "Beer", sep = "\n"), hjust = 0.6, align = TRUE, color = ab11,  fontsize = 5, barsize = 0,offset = 0.005, angle = 355,vjust = 0.5) + 
  geom_strip("HE005", "CLQCA_20_060", label = paste("French Guiana", "Human", sep = "\n"), hjust = 0.5, align = TRUE, color = fgh10,  fontsize = 5, barsize = 0,offset = 0.0043, angle = 75) + 
  geom_strip("LCBG_3Y8","LCBG_3D2", label = "MX Agave", align = TRUE,color = ma9,  fontsize = 5, barsize = 0,offset = 0.009, angle =340) + 
  geom_strip("STG_2","2162", label = "EU Oak", align = TRUE,color = mo4,  fontsize = 5, barsize = 0,offset = 0.009, angle = 340) + 
  geom_strip("CLIB324_2","SBE_1C", label = "Mixed Origin", hjust = 0.5, align = TRUE,color = mo8, fontsize = 5, barsize = 0,offset = 0.0035, angle = 50) + 
  geom_strip("YJM339", "YJM1122", label = "Clinical", hjust = 0.5, align = TRUE,color = clinical,  fontsize = 5, barsize = 0,offset = 0.0035, angle = 35) + 
  geom_strip("MTF2556", "MTF2546", label = "W. Africa", hjust = 0.5, align = TRUE,color = wac12,  fontsize = 5, barsize = 0,offset = 0.0033, angle = 25) + 
  geom_strip("N39_7A","N38_4A", label = "Siberia", align = TRUE,color = siberia,  fontsize = 5, barsize = 0,offset = 0.004, angle = 110) + 
  geom_strip("SAF29","SCN8", label = "Asian Fermentation", hjust = 0.5, align = TRUE,color = af26,  fontsize = 5, barsize = 0,offset = 0.0032, angle = 328) + 
  geom_strip("RIB1006","CBS2270", label = "Sake", hjust = 0.5, align = TRUE,color = s25,  fontsize = 5, barsize = 0,offset = 0.0035, angle = 355) + 
  geom_strip("SAN8", "SAN9", label = "S. Africa", align = TRUE,color = darksalmon,  fontsize = 5, barsize = 0,offset = 0.0035) + 
  geom_strip("NPA03_1","DJ74", label = paste("African", "Wine", sep = "\n"), hjust = 0.5, align = TRUE,color = apw13,  fontsize = 5, barsize = 0,offset = 0.004, angle = 299) + 
  geom_strip("HN19", "HN10", label = "CHN1", align = TRUE,color = ch14,  fontsize = 5, barsize = 0,offset = 0.0032) + 
  geom_strip("YPS154","N95_6_1C", label = "NA Oak", hjust = 0.5, align = TRUE,color = nao23,  fontsize = 5, barsize = 0,offset = 0.00315, angle = 295) + 
  geom_strip("S8BM_30_2D", "CBS1592", label = "Asian Islands", align = TRUE,color = ai24,  fontsize = 5, barsize = 0,offset = 0.003) +
  geom_strip("ES2M03_7A", "SJ5L12", label = "Taiwan2", align = TRUE,color = taiwan2, fontsize = 5, barsize = 0,offset = 0.003) + 
  geom_strip("N22_00_5D","N3_00_7A", label = "Far East Russia", align = TRUE,color = fer22,  fontsize = 5, barsize = 0,offset = 0.003) + 
  geom_strip("CLQCA_20_246","YPS617", label = "Ecuador2", align = TRUE,color = e21,  fontsize = 5, barsize = 0, offset = 0.003) + 
  geom_strip("HN16", "HN15", label = "CHN2", align = TRUE,color = ch14,  fontsize = 5, barsize = 0, offset = 0.003) + 
  geom_strip("UWOPS03_459_1","UWOPS03_461_4", label = "Malaysian", align = TRUE,color = malaysian,  fontsize = 5, barsize = 0, offset = 0.003) +
  geom_strip("WI016","Ksc2_2B", label = "Far East Asian", align = TRUE,color = fea18,  fontsize = 5, barsize = 0, offset = 0.003) + 
  geom_strip("SX3", "SX1", label = "CHN3", align = TRUE,color = ch14,  fontsize = 5, barsize = 0, offset = 0.003) + 
  geom_strip("EN14S01","GE14S01_7B", label = "Taiwan1", align = TRUE, color = t17, fontsize = 5, barsize = 0, offset = 0.003)

circtre <- circtre + ggtitle("(A)") + theme(title = element_text(size = 16,hjust = 0), plot.margin=grid::unit(c(0,0,0,0), "mm"))


ggsave("../../../../Papers_Scopel/GWAS/diploids/gainOnly/CircNJtreev5.png",
       plot = circtre,
       device = "png", 
       dpi = 600, 
       limitsize = FALSE,
       width = 15,
       height = 15,
       units = "in")

ggsave("../../../../Papers_Scopel/GWAS/diploids/gainOnly/CircNJtree.pdf",
       plot = circtre,
       device = cairo_pdf, 
       dpi = 600, 
       limitsize = FALSE,
       width = 15,
       height = 15,
       units = "in")

system2(command = "pdfcrop", 
        args    = c("../../../../Papers_Scopel/GWAS/diploids/gainOnly/CircNJtree.pdf", 
                    "../../../../Papers_Scopel/GWAS/diploids/gainOnly/CircNJtree.pdf"))






colltre <- ggtree(y,size=0.2, layout="rectangular") + 
  geom_treescale() +
  geom_tiplab(size=0.5,aes(label = paste(label,lineage,sep = " "), color=lineage),align = FALSE,linesize = 0.1, show.legend = FALSE)+
  scale_colour_manual(values = cladecolors) + 
  guides(colour = guide_legend(override.aes = list(size=3, label = "")))+
  geom_point2(aes(subset = as.numeric(label) >= 70, label = label), size=1.5, alpha = 0.5,color = "black")

colltre <- collapse(colltre, wine1, 'max',fill = we1, alpha = 0.8) %>% 
  collapse(wine2, 'max',fill = we1) %>%
  collapse(alpechin, 'max', fill = ol2) %>%
  collapse(wine3, 'max', fill = we1) %>%
  collapse(wine4, 'max', fill = we1) %>%
  collapse(bb, 'max', fill = bb3) %>%
  collapse(israelNode, 'max', fill = israel) %>%
  collapse(eq1, 'max', fill = ecuador1) %>%
  collapse(alebeer, 'max', fill = ab6) %>%
  collapse(frenchdairy, 'max', fill = fd5) %>%
  collapse(ethiopiaNode, 'max', fill = ethiopia) %>%
  collapse(africanbeer, 'max', fill = ab11) %>%
  collapse(frenchguianahuman, 'max', fill = fgh10) %>%
  collapse(mexicanagave, 'max', fill = ma9) %>%
  collapse(medoak, 'max', fill = mo4) %>%
  collapse(mixedorigin, 'max', fill = mo8) %>%
  collapse(clinicalNode, 'max', fill = clinical) %>%
  collapse(westAfrica, 'max', fill = wac12) %>%
  collapse(siberiaNode, 'max', fill = siberia) %>%
  collapse(AsianFerm1, 'max', fill = af26) %>%
  collapse(sakeNode, 'max', fill = s25) %>%
  collapse(AsianFerm2, 'max', fill = af26) %>%
  collapse(AsianFerm3, 'max', fill = af26) %>%
  collapse(AsianFerm4, 'max', fill = af26) %>%
  collapse(AsianFerm5, 'max', fill = af26) %>%
  collapse(AsianFerm6, 'max', fill = af26) %>%
  collapse(SouthAfrica, 'max', fill = darksalmon) %>%
  collapse(AfricanWine, 'max', fill = apw13) %>%
  collapse(CHN1, 'max', fill = ch14) %>%
  collapse(naoak, 'max', fill = nao23) %>%
  collapse(AsianIslands, 'max', fill = ai24) %>%
  collapse(taiwan2Node, 'max', fill = taiwan2) %>%
  collapse(farEastRussian, 'max', fill = fer22) %>%
  collapse(ecuador2, 'max', fill = e21) %>%
  collapse(CHN2, 'max', fill = ch14) %>%
  collapse(malaysianNode, 'max', fill = malaysian) %>%
  collapse(farEastAsian, 'max', fea18) %>%
  collapse(CHN3, 'max', fill = ch14)%>%
  collapse(Taiwan, 'max', fill =t17)

colltre <- flip(colltre, sakeNode, AsianFerm1)


colltre + geom_cladelab(node = wine1, label = "Wine", align = TRUE)
  
ctre + geom_text2(label= "Wine", y=750, x = 0.027, color = we1)
ggsave("../../../../Papers_Scopel/GWAS/diploids/gainOnly/NJtreeCollapsed.pdf",
       plot = ctre,
       device = "pdf", 
       dpi = 600, 
       limitsize = FALSE,
       width = 15,
       height = 15,
       units = "in")



