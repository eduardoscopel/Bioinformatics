library(data.table)

MATA <- read.table("~/Documents/Papers_Scopel/MATtest/red_amb_MATA.bc", col.names = c("basename","TotalACGT","N","Total"))
head(MATA)
MATA$amb <- MATA$Total - MATA$N - MATA$TotalACGT
peter <- df[df$reference == "Peter_et_al_2018",]
MATA <- merge(MATA, peter, by = "basename")
head(MATA)

MATA <- MATA[,c(1:5,15,16)]
mean(MATA[MATA$ploidyBAF == "Unknown","amb"])
mean(MATA[MATA$ploidyFACS == 1,"amb"])
mean(MATA[MATA$ploidyFACS == 2,"amb"])
mean(MATA[MATA$ploidyFACS == 3,"amb"])
mean(MATA[MATA$ploidyFACS == 4,"amb"])
mean(MATA[MATA$ploidyFACS == 1,"N"])
mean(MATA[MATA$ploidyFACS == 2,"N"])
mean(MATA[MATA$ploidyFACS == 3,"N"])
mean(MATA[MATA$ploidyFACS == 4,"N"])

mean(MATA[MATA$ploidyFACS == 2 & MATA$ploidyBAF == 2,"amb"])
mean(MATA[MATA$ploidyFACS == 2 & MATA$ploidyBAF == "Unknown","amb"])
mean(MATA[MATA$ploidyFACS == 1 & MATA$ploidyBAF == "Unknown","amb"])
