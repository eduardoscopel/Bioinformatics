library(data.table)
library(plyr)
library(dplyr)
library(phylolm)
library(scales)
library(epitools)
library(ggplot2)
library(MASS)
library(sjPlot)
library(sjmisc)


### Read tables and match
petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new.txt",header=TRUE)
sctab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /sc_table_most_recent.txt", header = TRUE)
nrow(petertab)
petertab <- petertab[complete.cases(petertab$Aneuploidies),]                                    # Remove strains with no aneuploidy information
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="Euploid","No","Yes")
petertab$amp_response <- ifelse(petertab$Aneuploidy_type =="Loss" | 
                                  petertab$Aneuploidy_type =="Euploid","No","Yes")
petertab$amp_response_no1 <- ifelse(petertab$Aneuploidy_type =="Loss" | 
                                      petertab$Aneuploidy_type =="Euploid" |
                                      petertab$Aneuploidy_type == "1only","No","Yes")
petertab$amp_binary <- ifelse(petertab$amp_response == "Yes", 1,0)
petertab$amp_binary_no1 <- ifelse(petertab$amp_response_no1 == "Yes", 1,0)

petertab$admixture <- ifelse(petertab$Clades == "M. Mosaic" |
                               petertab$Clades == "8. Mixed origin" |
                               petertab$Clades == "7. Mosaic beer","Admixed","Non-admixed")
petertab$polyploidy <- ifelse(petertab$Ploidy == 3 | 
                                petertab$Ploidy == 4 |
                                petertab$Ploidy == 5, "poly",
                              ifelse(petertab$Ploidy == 1, "haploid",
                                     "diploid"))
nrow(petertab)

# Get important categories from my table into peter df
for(i in 1:nrow(petertab)){
  for(j in 1:nrow(sctab)){
    if(!is.element(petertab$`Isolate name`[i],sctab$ID)){
      petertab$heterozygosity[i] <- NA
      petertab$new_eco[i] <- NA
      petertab$basename[i] <- NA
    }
    else if(petertab$`Isolate name`[i] == sctab$ID[j]){
      petertab$heterozygosity[i] <- sctab$heterozygosity[j]
      petertab$new_eco[i] <- sctab$ecology_category[j]
      petertab$basename[i] <- sctab$basename[j]
    }
  }
}

# Remove HO deleted strains
petertab <- petertab[petertab$`HO deletion` == "no",]
nrow(petertab)


### Look at the frequencies  
table(petertab$Ploidy,petertab$amp_response)
table(petertab$`Ecological origins`,petertab$amp_response)
table(petertab$new_eco, petertab$amp_response)
table(petertab$Clades, petertab$amp_response)
table(petertab$admixture, petertab$amp_response)
tapply(petertab$heterozygosity, petertab$amp_response, mean)
table(petertab$wine_eco,petertab$amp_response)

### Factorize the categorical variables
petertab$amp_response <- as.factor(petertab$amp_response)
petertab$amp_response_no1 <- as.factor(petertab$amp_response_no1)
petertab$aneuploidy_response <- as.factor(petertab$aneuploidy_response)
petertab$`Ecological origins` <- as.factor(petertab$`Ecological origins`)
petertab$Ploidy <- as.factor(petertab$Ploidy)
petertab$new_eco <- as.factor(petertab$new_eco)
petertab$admixture <- as.factor(petertab$admixture)
petertab$wine_eco <- as.factor(petertab$wine_eco)

# filter for haploids
haploids <- petertab[petertab$Ploidy == 1,]
# filter for diploids
diploids <- petertab[petertab$Ploidy == 2,]
# filter for triploids
triploids <- petertab[petertab$Ploidy == 3,]
# filter for tetraploids
tetraploids <- petertab[petertab$Ploidy == 4, ]
# filter for pentaploids
pentaploids <- petertab[petertab$Ploidy == 5, ]


######################################################################################################################
###################################################### ECOLOGY MODELS ################################################
######################################################################################################################
diploids <- diploids[diploids$heterozygosity <= 0.005,]
diploids <- diploids[complete.cases(diploids$heterozygosity),]
diploids <- diploids[complete.cases(diploids$amp_response),]
dip_eco <- diploids[complete.cases(diploids$`Ecological origins`),]
nullmodel <- glm(aneuploidy_response ~ 1,
                 family = binomial(link="logit"),
                 dip_eco)
fullmodel <- glm(aneuploidy_response ~ `Ecological origins`*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 dip_eco)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove 2way int term
nmodel <- update(cmodel, ~. -`Ecological origins`:heterozygosity)
summary(nmodel)
anova(nmodel, test="Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel)
anova(cmodel,nmodel, test="Chisq")
finalmodel <- cmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel,nullmodel,test="Chisq")


### Wine/EU only
wineeu <- dip_eco[dip_eco$Clades == "1. Wine/European",]
table(wineeu$Ploidy, wineeu$amp_response)

nullmodel <- glm(amp_response ~ 1,
                 family = binomial(link = "logit"),
                 wineeu)
fullmodel <- glm(amp_response ~ `Ecological origins`*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 wineeu)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove eco:het
nmodel <- update(cmodel, ~. - ecology:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test = "Chisq")
cmodel <- nmodel
# Remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
finalmodel <- nmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test = "Chisq")

### Mixed Origin only
mixor <- dip_eco[dip_eco$Clades == "8. Mixed origin",]
table(mixor$`Ecological origins`, mixor$amp_response)
table(mixor$Ploidy, mixor$amp_response)

nullmodel <- glm(amp_response ~ 1,
                 family = binomial(link = "logit"),
                 mixor)
fullmodel <- glm(amp_response ~ `Ecological origins`*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 mixor)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove eco:het
nmodel <- update(cmodel, ~. - `Ecological origins`:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test = "Chisq")
cmodel <- nmodel
# Remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
finalmodel <- nmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test = "Chisq")

### Mosaics only
dip_mos <- dip_eco[dip_eco$Clades == "M. Mosaic",]
table(dip_mos$`Ecological origins`, dip_mos$amp_response)
nullmodel <- glm(amp_response ~ 1,
                 family = binomial(link = "logit"),
                 dip_mos)
fullmodel <- glm(amp_response ~ `Ecological origins`*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 dip_mos)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove eco:het
nmodel <- update(cmodel, ~. - `Ecological origins`:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test = "Chisq")
cmodel <- nmodel
# Remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
finalmodel <- nmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test = "Chisq")


######################################################################################################################
####################################################### CLADE MODELS #################################################
######################################################################################################################
diploids <- diploids[diploids$heterozygosity <= 0.005,]
diploids <- diploids[complete.cases(diploids$heterozygosity),]
diploids <- diploids[complete.cases(diploids$amp_response),]
dip_clades <- diploids[complete.cases(diploids$Clades),]
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), dip_clades)
fullmodel <- glm(amp_response ~ Clades*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 dip_clades)
summary(fullmodel)
anova(fullmodel)
cmodel <- fullmodel
# Remove 2way int term
nmodel <- update(cmodel, ~. - Clades:heterozygosity)
summary(nmodel)
anova(nmodel)
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel)
anova(cmodel,nmodel, test="Chisq")
finalmodel <- nmodel
summary(finalmodel)
anova(finalmodel, test="Chisq")
anova(finalmodel,nullmodel, test = "Chisq")


### Testing domestication level
dip_clades$domestication <- ifelse(dip_clades$`Ecological origins` == "Wine" | 
                                     dip_clades$`Ecological origins` == "Beer" | 
                                     dip_clades$`Ecological origins` == "Sake" | 
                                     dip_clades$`Ecological origins` == "Bakery" | 
                                     dip_clades$`Ecological origins` == "Fermentation" | 
                                     dip_clades$`Ecological origins` == "Palm wine" | 
                                     dip_clades$`Ecological origins` == "Industrial" | 
                                     dip_clades$`Ecological origins` == "Distillery" | 
                                     dip_clades$`Ecological origins` == "Dairy" | 
                                     dip_clades$`Ecological origins` == "Bioethanol" | 
                                     dip_clades$`Ecological origins` == "Cider" | 
                                     dip_clades$`Ecological origins` == "Probiotic", "domesticated",
                                   ifelse(dip_clades$`Ecological origins` == "Tree" | 
                                            dip_clades$`Ecological origins` == "Nature" | 
                                            dip_clades$`Ecological origins` == "Fruit" | 
                                            dip_clades$`Ecological origins` == "Soil" | 
                                            dip_clades$`Ecological origins` == "Insect" | 
                                            dip_clades$`Ecological origins` == "Flower" | 
                                            dip_clades$`Ecological origins` == "Water", "wild",
                                          ifelse(dip_clades$`Ecological origins` == "Human" | 
                                                   dip_clades$`Ecological origins` == "Human, clinical", "human",
                                                 ifelse(dip_clades$`Ecological origins` == "Lab strain", "lab", "unknown"))))

### Industrial strains
ind_dip <- dip_clades[dip_clades$domestication == "domesticated",]
nrow(ind_dip)
fullmodel <- glm(amp_response ~ Clades*heterozygosity, family = binomial(link = "logit"), ind_dip)
nullmodel <- glm(amp_response ~ 1, family = binomial(link = "logit"), ind_dip)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
summary(nullmodel)
anova(nullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove clades:het
nmodel <- update(cmodel, ~. - Clades:heterozygosity)
summary(nmodel)
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test= "Chisq")
anova(cmodel, nmodel, test = "Chisq")
finalmodel <- nmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test = "Chisq")

### Wild strains
wild_dip <- dip_clades[dip_clades$domestication == "wild",]
nrow(wild_dip)
fullmodel <- glm(amp_response ~ Clades*heterozygosity, family = binomial(link = "logit"), wild_dip)
nullmodel <- glm(amp_response ~ 1, family = binomial(link = "logit"), wild_dip)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
summary(nullmodel)
anova(nullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove clades:het
nmodel <- update(cmodel, ~. - Clades:heterozygosity)
summary(nmodel)
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test= "Chisq")
anova(cmodel, nmodel, test = "Chisq")
finalmodel <- nmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test = "Chisq")

### Human strains
hum_dip <- dip_clades[dip_clades$domestication == "human",]
nrow(hum_dip)
fullmodel <- glm(amp_response ~ Clades*heterozygosity, family = binomial(link = "logit"), hum_dip)
nullmodel <- glm(amp_response ~ 1, family = binomial(link = "logit"), hum_dip)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
summary(nullmodel)
anova(nullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove clades:het
nmodel <- update(cmodel, ~. - Clades:heterozygosity)
summary(nmodel)
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test= "Chisq")
anova(cmodel, nmodel, test = "Chisq")
finalmodel <- nmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test = "Chisq")


######################################################################################################################
############################## Ecology and polyploidy ################################################################
######################################################################################################################
petertab <- petertab[petertab$heterozygosity < 0.004,]
petertab <- petertab[complete.cases(petertab$heterozygosity),]
petertab <- petertab[complete.cases(petertab$amp_response),]
petertab <- petertab[complete.cases(petertab$aamp_response_no1),]
peter_eco <- petertab[complete.cases(petertab$`Ecological origins`),]
peter_eco <- peter_eco[complete.cases(peter_eco$wine_eco)]
nullmodel <- glm(amp_response_no1~1, family = binomial(link="logit"), peter_eco)
fullmodel <- glm(amp_response_no1 ~ `Ecological origins`*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 peter_eco)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove eco:poly
nmodel <- update(cmodel, ~. - `Ecological origins`:polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove polyploidy
nmodel <- update(cmodel, ~. - polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel<-nmodel
# Remove eco
nmodel <- update(cmodel, ~. - `Ecological origins`:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel

nmodel <- update(cmodel, ~. - `Ecological origins`:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel

nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel

nmodel <- update(cmodel, ~. -polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
ephfinal <- nmodel
summary(ephfinal)
anova(ephfinal, test = "Chisq")
anova(ephfinal, nullmodel, test = "Chisq")

# Final model
finalmodel <- cmodel # Reject new model, no more changes to make
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test ="Chisq")


### Wine/EU clade
wine_peter <- peter_eco[peter_eco$Clades == "1. Wine/European",]
table(wine_peter$Ploidy,wine_peter$amp_response)
table(wine_peter$wine_eco,wine_peter$amp_response)
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), wine_peter)
fullmodel <- glm(amp_response ~ wine_eco*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 wine_peter)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove eco:poly
nmodel <- update(cmodel, ~. -wine_eco:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove polyploidy
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove eco
nmodel <- update(cmodel, ~. - wine_eco)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# Final model
finalmodel <- cmodel # Reject new model, no more changes to make
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test ="Chisq")

mixed_or <- peter_eco[peter_eco$Clades == "8. Mixed origin",]
table(mixed_or$Ploidy,mixed_or$amp_response)
table(mixed_or$`Ecological origins`,mixed_or$amp_response)
table(mixed_or$wine_eco, mixed_or$amp_response)
fullmodel <- glm(amp_response ~ wine_eco*polyploidy, family = binomial(link = "logit"), mixed_or)
summary(fullmodel)
anova(fullmodel, test="Chisq")
                 
### Mosaic only
mosaic <- peter_eco[peter_eco$Clades == "M. Mosaic",]
table(mosaic$Ploidy,mosaic$amp_response)
table(mosaic$`Ecological origins`,mosaic$amp_response)
table(mosaic$`Ecological origins`,mosaic$polyploidy)
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), mosaic)
fullmodel <- glm(amp_response ~`Ecological origins`* 
                   polyploidy*heterozygosity, 
                 family = binomial(link = "logit"),
                 mosaic)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove eco:poly:het
nmodel <- update(cmodel, ~. -`Ecological origins`:polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove ploidy:het
nmodel <- update(cmodel, ~. - polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# Remove eco:het
nmodel <- update(cmodel, ~. - `Ecological origins`:heterozygosity)
summary(nmodel)
anova(nmodel, test= "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# Remove eco:ploidy
nmodel <- update(cmodel, ~. - `Ecological origins`:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# Remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test ="Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# Final model
finalmodel <- cmodel # Reject new model, no more changes to make
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test ="Chisq")

### Mosaic + MO
admixed <- peter_eco[peter_eco$Clades == "M. Mosaic" |
                       peter_eco$Clades == "8. Mixed origin" ,]
table(admixed$Ploidy,admixed$amp_response)
table(admixed$`Ecological origins`,admixed$amp_response)
table(admixed$`Ecological origins`,admixed$polyploidy)
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), admixed)
fullmodel <- glm(amp_response ~ 
                   polyploidy*`Ecological origins`, 
                 family = binomial(link = "logit"),
                 admixed)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove eco:poly
nmodel <- update(cmodel, ~. -`Ecological origins`:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove eco
nmodel <- update(cmodel, ~. -`Ecological origins`)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# Remove ploidy
nmodel <- update(cmodel, ~. -polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# Final model
finalmodel <- cmodel # Reject new model, no more changes to make
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test ="Chisq")


######################################################################################################################
############################## Clades and polyploidy #################################################################
######################################################################################################################
petertab <- petertab[petertab$heterozygosity < 0.004,]
petertab <- petertab[complete.cases(petertab$heterozygosity),]
petertab <- petertab[complete.cases(petertab$amp_response),]
peter_clades <- petertab[complete.cases(petertab$Clades),]
nrow(peter_clades)
# No mosaics
#peter_clades <- peter_clades[peter_clades$Clades != "M. Mosaic",]
# No mixed origin
#peter_clades <- peter_clades[peter_clades$Clades != "8. Mixed origin",]
# No admixed (mosaic or mixed origin)
#peter_clades <- peter_clades[peter_clades$Clades != "8. Mixed origin",]
#peter_clades <- peter_clades[peter_clades$Clades != "M. Mosaic",]
nrow(peter_clades)
table(peter_clades$Ploidy, peter_clades$amp_response)
table(peter_clades$Clades, peter_clades$amp_response)
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), peter_clades)
cphfull <- glm(amp_response ~ Clades*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 peter_clades)
cmodel <- cphfull
# Remove eco:poly
nmodel <- update(cmodel, ~. -Clades:polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove polyploidy
nmodel <- update(cmodel, ~. - polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove clades
nmodel <- update(cmodel, ~. -  Clades:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
nmodel <- update(cmodel, ~. - Clades:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cphfinal <- cmodel
anova(cphfinal, nullmodel, test = "Chisq")
# Final model
summary(cphfinal)
anova(cphfinal, test = "Chisq")
anova(cphfinal, nullmodel, test ="Chisq")

### Changing order of variables
# clades*het*ploidy
chpfull <- glm(amp_response ~ Clades*
                 heterozygosity*
                 polyploidy, 
               family = binomial(link = "logit"),
               peter_clades)
summary(chpfull)
anova(chpfull, test = "Chisq")
anova(chpfull, nullmodel, test = "Chisq")
cmodel <- chpfull
# Remove eco:poly
nmodel <- update(cmodel, ~. -Clades:heterozygosity:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove polyploidy
nmodel <- update(cmodel, ~. - heterozygosity:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove clades
nmodel <- update(cmodel, ~. -  Clades:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
nmodel <- update(cmodel, ~. - Clades:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
chpfinal <- cmodel
anova(chpfinal, nullmodel, test = "Chisq")
# Final model
summary(chpfinal)
anova(chpfinal, test = "Chisq")
anova(chpfinal, nullmodel, test ="Chisq")




# het*clades*ploidy
hcpfull <- glm(amp_response ~ heterozygosity*
                 Clades*
                 polyploidy, 
               family = binomial(link = "logit"),
               peter_clades)
summary(hcpfull)
anova(hcpfull, test = "Chisq")
cmodel <- hcpfull
nmodel <- update(cmodel, ~. -heterozygosity:Clades:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove polyploidy
nmodel <- update(cmodel, ~. - Clades:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove clades
nmodel <- update(cmodel, ~. -  heterozygosity:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel

nmodel <- update(cmodel, ~. - heterozygosity:Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
hcpfinal <- cmodel
summary(hcpfinal)
anova(hcpfinal, test = "Chisq")
anova(hcpfinal, nullmodel, test ="Chisq")

# het*ploidy*clades
hpcfull <- glm(amp_response ~ heterozygosity*
                 polyploidy*
                 Clades, 
               family = binomial(link = "logit"),
               peter_clades)
summary(hpcfull)
anova(hpcfull, test = "Chisq")
cmodel <- hpcfull
nmodel <- update(cmodel, ~. -heterozygosity:polyploidy:Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove polyploidy
nmodel <- update(cmodel, ~. - polyploidy:Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove clades
nmodel <- update(cmodel, ~. -  heterozygosity:Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel

nmodel <- update(cmodel, ~. - heterozygosity:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
nmodel <- update(cmodel, ~. - Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
hpcfinal <- cmodel
summary(hpcfinal)
anova(hpcfinal, test = "Chisq")
anova(hpcfinal, nullmodel, test ="Chisq")

# ploidy*het*clades
phcfull <- glm(amp_response ~ polyploidy*
                 heterozygosity*
                 Clades, 
               family = binomial(link = "logit"),
               peter_clades)
summary(phcfull)
anova(phcfull, test = "Chisq")
cmodel <- phcfull
nmodel <- update(cmodel, ~. -polyploidy:heterozygosity:Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove polyploidy
nmodel <- update(cmodel, ~. - heterozygosity:Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove clades
nmodel <- update(cmodel, ~. -  polyploidy:Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel

nmodel <- update(cmodel, ~. - polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
nmodel <- update(cmodel, ~. - Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
phcfinal <- cmodel
summary(phcfinal)
anova(phcfinal, test = "Chisq")
anova(phcfinal, nullmodel, test ="Chisq")

# ploidy*clades*het
pchfull <- glm(amp_response ~ polyploidy*
                 Clades*
                 heterozygosity, 
               family = binomial(link = "logit"),
               peter_clades)
summary(pchfull)
anova(pchfull, test = "Chisq")
cmodel <- pchfull
nmodel <- update(cmodel, ~. -polyploidy:Clades:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove polyploidy
nmodel <- update(cmodel, ~. - Clades:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove clades
nmodel <- update(cmodel, ~. -  polyploidy:heterozygosity )
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel

nmodel <- update(cmodel, ~. - polyploidy:Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
nmodel <- update(cmodel, ~. - Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
pchfinal <- cmodel
summary(pchfinal)
anova(pchfinal, test = "Chisq")
anova(pchfinal, nullmodel, test ="Chisq")


# clades*het*ploidy
chpfull <- glm(amp_response ~ Clades*
                 heterozygosity*
                 polyploidy, 
               family = binomial(link = "logit"),
               peter_clades)
summary(chpfull)
anova(chpfull, test = "Chisq")
cmodel <- chpfull
# remove c:p:h
nmodel <- update(cmodel, ~. -Clades:heterozygosity:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove ph
nmodel <- update(cmodel, ~. - heterozygosity:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove ch
nmodel <- update(cmodel, ~. -  Clades:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# keep ch and remove cp
nmodel <- update(cmodel, ~. - Clades:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# remove p
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# keep p
cphfinal <- cmodel
summary(cphfinal)
anova(cphfinal, test = "Chisq")
anova(cphfinal, nullmodel, test ="Chisq")

# clades*ploidy*het
cphfull <- glm(amp_binary ~ Clades*
                 polyploidy*
                 heterozygosity, 
               family = binomial(link = "logit"),
               peter_clades)
summary(cphfull)
anova(cphfull, test = "Chisq")
cmodel <- cphfull
# remove c:p:h
nmodel <- update(cmodel, ~. -Clades:polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove ph
nmodel <- update(cmodel, ~. - polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove ch
nmodel <- update(cmodel, ~. -  Clades:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# keep ch and remove cp
nmodel <- update(cmodel, ~. - Clades:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove p
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# keep p
cphfinal <- cmodel
summary(cphfinal)
anova(cphfinal, test = "Chisq")
anova(cphfinal, nullmodel, test ="Chisq")

for(i in levels(peter_clades$Clades)){
  auto_p_m <- plot_model(cphfinal, 
                         type = "pred",
                         terms = c("heterozygosity [all]", "polyploidy", paste("Clades [", i ,"]", sep = "")), 
                         ci.lvl = 0.95,
                         title = paste("Strains from the ",i, " clade (", nrow(peter_clades[peter_clades$Clades == i,]),")", sep = ""),
                         axis.title = c("Heterozygosity","Predicted probability of \nchromosome amplification")) + 
    geom_point(data = peter_clades[peter_clades$Clades == i,], 
               aes(x = heterozygosity, y = amp_binary, color = polyploidy), 
               position = position_jitter(0,0.05),
               inherit.aes = FALSE,
               shape = 1) +
    scale_x_continuous(expand = c(0,0.0001), limits = c(0,NA))+
    scale_y_continuous(sec.axis = sec_axis(~., name = "Amplification (yes = 1; no = 0)", breaks = seq(0,1,1)))+
    theme_classic()+
    theme(panel.grid.major.y = element_line(size=.5, linetype = "dashed"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line.y.left = element_line(size=0.51, linetype = "solid"),
          axis.line.x = element_line(size=0.5, linetype = "solid"),
          panel.border = element_blank())
  png(paste("Documents/GitHub/eduardo/aneuploidy/plots/peter/glm/", substr(i,1,2), ".png",sep=""), 
      width = 1800, 
      height = 1100,
      res = 300)
  print(auto_p_m)
  dev.off()
}

auto_p_m <- plot_model(cphfinal, 
           type = "pred",
           terms = c("heterozygosity [all]", "polyploidy", "Clades [8. Mixed origin]"), 
           ci.lvl = 0.95,
           title = paste("Strains from the 8. Mixed origin clade (", nrow(testdf),")", sep = ""),
           axis.title = c("Heterozygosity","Predicted probability of chromosome amplification")) + 
  geom_point(data = peter_clades[peter_clades$Clades == "8. Mixed origin",], 
                      aes(x = heterozygosity, y = amp_binary, color = polyploidy), 
                      position = position_jitter(0,0.05),
                      inherit.aes = FALSE)+
  scale_x_continuous(expand = c(0,0.0001), limits = c(0,NA))+
  scale_y_continuous(sec.axis = sec_axis(~., name = "Amplification (yes = 1; no = 0)", breaks = seq(0,1,1)))+
  theme_classic()+
  theme(panel.grid.major.y = element_line(size=.5, linetype = "dashed"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.y.left = element_line(size=0.51, linetype = "solid"),
        axis.line.y.right = element_line(size=0.51, linetype = "solid"),
        axis.line.x = element_line(size=0.5, linetype = "solid"),
        panel.border = element_blank())
auto_p_m


### Testing domestication level
peter_clades$domestication <- ifelse(peter_clades$`Ecological origins` == "Wine" | 
                                     peter_clades$`Ecological origins` == "Beer" | 
                                     peter_clades$`Ecological origins` == "Sake" | 
                                     peter_clades$`Ecological origins` == "Bakery" | 
                                     peter_clades$`Ecological origins` == "Fermentation" | 
                                     peter_clades$`Ecological origins` == "Palm wine" | 
                                     peter_clades$`Ecological origins` == "Industrial" | 
                                     peter_clades$`Ecological origins` == "Distillery" | 
                                     peter_clades$`Ecological origins` == "Dairy" | 
                                     peter_clades$`Ecological origins` == "Bioethanol" | 
                                     peter_clades$`Ecological origins` == "Cider" | 
                                     peter_clades$`Ecological origins` == "Probiotic", "domesticated",
                                   ifelse(peter_clades$`Ecological origins` == "Tree" | 
                                            peter_clades$`Ecological origins` == "Nature" | 
                                            peter_clades$`Ecological origins` == "Fruit" | 
                                            peter_clades$`Ecological origins` == "Soil" | 
                                            peter_clades$`Ecological origins` == "Insect" | 
                                            peter_clades$`Ecological origins` == "Flower" | 
                                            peter_clades$`Ecological origins` == "Water", "wild",
                                          ifelse(peter_clades$`Ecological origins` == "Human" | 
                                                   peter_clades$`Ecological origins` == "Human, clinical", "human",
                                                 ifelse(peter_clades$`Ecological origins` == "Lab strain", "lab", "unknown"))))

### Industrial strains
ind_peter <- peter_clades[peter_clades$domestication == "domesticated",]
table(ind_peter$Ploidy, ind_peter$amp_response)
table(ind_peter$Clades, ind_peter$amp_response)
nrow(ind_peter)
fullmodel <- glm(amp_response ~ Clades*polyploidy, family = binomial(link = "logit"), ind_peter)
nullmodel <- glm(amp_response ~ 1, family = binomial(link = "logit"), ind_peter)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
summary(nullmodel)
anova(nullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove clades:polyploidy
nmodel <- update(cmodel, ~. - Clades:polyploidy)
summary(nmodel)
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove polyploidy
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test= "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove clades
nmodel <- update(cmodel, ~. -Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# final model
finalmodel <- cmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test = "Chisq")

### Wild strains
wild_peter <- peter_clades[peter_clades$domestication == "wild",]
table(wild_peter$Clades, wild_peter$amp_response)
table(wild_peter$Ploidy, wild_peter$amp_response)

nrow(wild_peter)
fullmodel <- glm(amp_response ~ Clades*polyploidy, family = binomial(link = "logit"), wild_peter)
nullmodel <- glm(amp_response ~ 1, family = binomial(link = "logit"), wild_peter)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
summary(nullmodel)
anova(nullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove clades:polyploidy
nmodel <- update(cmodel, ~. - Clades:polyploidy)
summary(nmodel)
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove polyploidy
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test= "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove clade
nmodel <- update(cmodel, ~. - Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# final model
finalmodel <- cmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test = "Chisq")

### Human strains
hum_peter <- peter_clades[peter_clades$domestication == "human",]
nrow(hum_peter)
table(hum_peter$Clades, hum_peter$amp_response)
table(hum_peter$Ploidy, hum_peter$amp_response)
fullmodel <- glm(amp_response ~ Clades*polyploidy, family = binomial(link = "logit"), hum_peter)
nullmodel <- glm(amp_response ~ 1, family = binomial(link = "logit"), hum_peter)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
summary(nullmodel)
anova(nullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove clades:het
nmodel <- update(cmodel, ~. - Clades:polyploidy)
summary(nmodel)
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove het
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test= "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove clades
nmodel <- update(cmodel, ~. - Clades)
summary(nmodel)
anova(nmodel, test ="Chisq")
anova(cmodel, nmodel, test = "Chisq")
# final model
finalmodel <- cmodel
summary(finalmodel)
anova(finalmodel, test = "Chisq")
anova(finalmodel, nullmodel, test = "Chisq")
