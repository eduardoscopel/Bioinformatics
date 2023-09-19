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


# Subset haploids
haploids <- df[df$ploidyCons == 1,]
nrow(haploids)
table(haploids$aneuploidy_type, useNA = "ifany")
table(haploids$ecology_category, haploids$aneuploidy_type)

# Subset diploids
diploids <- df[df$ploidyCons == 2,]
nrow(diploids)
table(diploids$aneuploidy_type, useNA = "ifany")
table(diploids$ecology_category, diploids$aneuploidy_type)
dipGainOnly <- diploids[diploids$aneuploidy_type == "Gain" | diploids$aneuploidy_type == "No",]
nrow(dipGainOnly)
popind <- data.frame(ID= dipGainOnly$basename, sex = rep("Unknown", 905), eco = dipGainOnly$ecology_category)
write.table(popind, "gain2n_QC4.ind")
phenotypeMatrix <- data.frame(Taxa = dipGainOnly$basename, ChrGain = dipGainOnly$aneuploidy_binary)
phenotypeMatrix
write.table(dipGainOnly,
            file = "../gainOnly/dipGainOnly.txt", 
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')

# Subset polyploids
polyploids <- df[df$ploidyCons > 2 &
                   df$ploidyCons != "Unknown",]
nrow(polyploids)
table(polyploids$aneuploidy_type, useNA = "ifany")
table(polyploids$ecology_category, polyploids$aneuploidy_type)

# Subset unknown
unknownPloidy <- df[df$ploidyCons =="Unknown",]
nrow(unknownPloidy)
table(unknownPloidy$aneuploidy_type, useNA = "ifany")
table(unknownPloidy$ecology_category, unknownPloidy$aneuploidy_type)

hist(diploids$heterozygosity,breaks = 100, xlim = c(0,0.006))
abline(v= mean(diploids$heterozygosity))
abline(v= mean(diploids$heterozygosity)-3*sd(diploids$heterozygosity),lty = "dashed")
abline(v= mean(diploids$heterozygosity)+3*sd(diploids$heterozygosity), lty = "dashed")


