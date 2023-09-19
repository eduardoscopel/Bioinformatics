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
#sctab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /sc_table_most_recent.txt", header = TRUE)
nrow(petertab)
petertab <- petertab[complete.cases(petertab$Aneuploidies),]                                    # Remove strains with no aneuploidy information
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="Euploid","No","Yes")
petertab$amp_response <- ifelse(petertab$Aneuploidy_type =="Euploid","No","Yes")
petertab$amp_binary <- ifelse(petertab$amp_response == "Yes", 1,0)
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
#for(i in 1:nrow(petertab)){
#  for(j in 1:nrow(sctab)){
#    if(!is.element(petertab$`Isolate name`[i],sctab$ID)){
#      petertab$heterozygosity[i] <- NA
#      petertab$new_eco[i] <- NA
#      petertab$basename[i] <- NA
#    }
#    else if(petertab$`Isolate name`[i] == sctab$ID[j]){
#      petertab$heterozygosity[i] <- sctab$heterozygosity[j]
#      petertab$new_eco[i] <- sctab$ecology_category[j]
#      petertab$basename[i] <- sctab$basename[j]
#    }
#  }
#}

# Remove HO deleted strains
petertab <- petertab[petertab$`HO deletion` == "no",]
petertab <- petertab[complete.cases(petertab$heterozygosity),]
petertab <- petertab[complete.cases(petertab$`Ecological origins`),]
petertab <- petertab[complete.cases(petertab$amp_response),]
nrow(petertab)


### Look at the frequencies  
table(petertab$Ploidy,petertab$amp_response)
table(petertab$`Ecological origins`,petertab$amp_response)
table(petertab$Clades, petertab$amp_response)
table(petertab$admixture, petertab$amp_response)
tapply(petertab$heterozygosity, petertab$amp_response, mean)
table(petertab$wine_eco,petertab$amp_response)

### Factorize the categorical variables
petertab$amp_response <- as.factor(petertab$amp_response)
petertab$aneuploidy_response <- as.factor(petertab$aneuploidy_response)
petertab$`Ecological origins` <- as.factor(petertab$`Ecological origins`)
petertab$Ploidy <- as.factor(petertab$Ploidy)
petertab$new_eco <- as.factor(petertab$new_eco)
petertab$admixture <- as.factor(petertab$admixture)
petertab$wine_eco <- as.factor(petertab$wine_eco)

######################################################################################################################
############################## Ecology, Heterozygosity and polyploidy (all strains) ##################################
######################################################################################################################
# Remove chr loss and amplification of chr 1
peter_eco <- petertab
peter_eco <- peter_eco[peter_eco$Aneuploidy_type == "Euploid" | 
                         peter_eco$Aneuploidy_type == "Gain" | 
                         peter_eco$Aneuploidy_type == "Both",]
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), peter_eco)
fullmodel <- glm(amp_response ~ `Ecological origins`*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 peter_eco)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove eco:poly:het
nmodel <- update(cmodel, ~. - `Ecological origins`:polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
# Accept current model as final
eph_no1_noL_final <- cmodel
summary(eph_no1_noL_final)
anova(eph_no1_noL_final, test = "Chisq")
anova(eph_no1_noL_final, nullmodel, test = "Chisq")

### Zooming into the Wine/EU clade
peter_eco <- peter_eco[complete.cases(peter_eco$wine_eco)]
wine_peter <- peter_eco[peter_eco$Clades == "1. Wine/European",]
wine_peter <- wine_peter[wine_peter$`Ecological origins` != "Unknown",]
table(wine_peter$Ploidy,wine_peter$amp_response)
table(wine_peter$wine_eco,wine_peter$amp_response)
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), wine_peter)
fullmodel <- glm(amp_response ~ `Ecological origins`*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 wine_peter)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
anova(fullmodel, nullmodel,test = "Chisq")
# Accept null model as final
wine_ephfinal <- fullmodel
summary(wine_ephfinal)
anova(wine_ephfinal, test = "Chisq")
anova(wine_ephfinal, nullmodel, test = "Chisq")


### Using recategorized wine environments as explanatory variable
wine_peter <- wine_peter[wine_peter$wine_eco != "Unknown",]
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), wine_peter)
fullmodel <- glm(amp_response ~ wine_eco*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 wine_peter)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
anova(fullmodel, nullmodel,test = "Chisq")
# Accept null model as final
recat_wine_ephfinal <- fullmodel

### Collapsing Wine ecologies into Industrial, Nature, Other
wine_peter$collapsed_eco <- ifelse(wine_peter$wine_eco == "Wine" |
                                     wine_peter$wine_eco == "Palm wine" |
                                     wine_peter$wine_eco == "Fermentation" |
                                     wine_peter$wine_eco == "Industrial" |
                                     wine_peter$wine_eco == "Probiotic" |
                                     wine_peter$wine_eco == "Distillery" |
                                     wine_peter$wine_eco == "Beer" |
                                     wine_peter$wine_eco == "Cider" |
                                     wine_peter$wine_eco == "Bakery", "Industrial",
                                   ifelse(wine_peter$wine_eco == "Nature" |
                                            wine_peter$wine_eco == "Fruit" |
                                            wine_peter$wine_eco == "Insect" |
                                            wine_peter$wine_eco == "Tree" |
                                            wine_peter$wine_eco == "Flower" |
                                            wine_peter$wine_eco == "Soil","Nature","Other"))
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), wine_peter)
fullmodel <- glm(amp_response ~ collapsed_eco*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 wine_peter)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
anova(fullmodel, nullmodel,test = "Chisq")
cmodel <- fullmodel
# remove eco:ploidy:het
nmodel <- update(cmodel, ~. - collapsed_eco:polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
#remove ploidy:het
nmodel <- update(cmodel, ~. - polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove eco:het
nmodel <- update(cmodel, ~. - collapsed_eco:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove eco:ploidy
nmodel <- update(cmodel, ~. - collapsed_eco:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
#remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove ploidy
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# remove eco
nmodel <- update(cmodel, ~. - collapsed_eco)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
collaps_wine_ephfinal <- cmodel
summary(collaps_wine_ephfinal)
anova(collaps_wine_ephfinal, test = "Chisq")
anova(collaps_wine_ephfinal, nullmodel,test = "Chisq")

### Mixed origin
mixed_or <- peter_clades[peter_clades$Clades == "8. Mixed origin",]
mixed_or$`Ecological origins` <- as.factor(mixed_or$`Ecological origins`)
mixed_or$amp_response <- as.factor(mixed_or$amp_response)

table(mixed_or$Ploidy,mixed_or$amp_response)
table(mixed_or$`Ecological origins`,mixed_or$amp_response)
nullmodel <- glm(amp_response~1, family = binomial(link = "logit"), mixed_or)
fullmodel <- glm(amp_response ~ `Ecological origins`*polyploidy*heterozygosity, family = binomial(link = "logit"), mixed_or)
summary(fullmodel)
anova(fullmodel, test="Chisq")
mo_ephfinal <- fullmodel
summary(mo_ephfinal)
anova(mo_ephfinal, test="Chisq")
anova(mo_ephfinal, nullmodel,test="Chisq")

### Mosaic only
mosaic <- peter_eco[peter_eco$Clades == "M. Mosaic",]
table(mosaic$Ploidy,mosaic$amp_response)
table(mosaic$`Ecological origins`,mosaic$amp_response)
table(mosaic$`Ecological origins`,mosaic$polyploidy)
nullmodel <- glm(amp_response_no1~1, family = binomial(link="logit"), mosaic)
fullmodel <- glm(amp_response_no1 ~`Ecological origins`* 
                   polyploidy*heterozygosity, 
                 family = binomial(link = "logit"),
                 mosaic)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
anova(fullmodel, nullmodel, test = "Chisq")
mosaic_ephfinal <- fullmodel
summary(mosaic_ephfinal)
anova(mosaic_ephfinal, test = "Chisq")
anova(mosaic_ephfinal, nullmodel, test = "Chisq")

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



#### Removing possible hybrids and repeating the analyses
peter_eco <- peter_eco[peter_eco$heterozygosity < 0.004,]
peter_eco <- peter_eco[complete.cases(peter_eco$heterozygosity),]
peter_eco <- peter_eco[complete.cases(peter_eco$amp_response),]
peter_eco <- peter_eco[complete.cases(peter_eco$`Ecological origins`),]
peter_eco <- peter_eco[complete.cases(peter_eco$wine_eco)]
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), peter_eco)
fullmodel <- glm(amp_response ~ `Ecological origins`*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 peter_eco)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
cmodel <- fullmodel
# Remove eco:poly:het
nmodel <- update(cmodel, ~. - `Ecological origins`:polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove poly:het
nmodel <- update(cmodel, ~. - polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel<-nmodel
# Remove eco:het
nmodel <- update(cmodel, ~. - `Ecological origins`:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel
# remove eco:poly
nmodel <- update(cmodel, ~. - `Ecological origins`:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel
# remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel
# remove poly
nmodel <- update(cmodel, ~. -polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# remove eco
eph_nohyb_final <- cmodel
summary(eph_nohyb_final)
anova(eph_nohyb_final, test = "Chisq")
anova(eph_nohyb_final, nullmodel, test = "Chisq")


### Zooming into the Wine/EU clade
peter_eco <- peter_eco[complete.cases(peter_eco$wine_eco)]
wine_peter <- peter_eco[peter_eco$Clades == "1. Wine/European",]
table(wine_peter$Ploidy,wine_peter$amp_response)
table(wine_peter$wine_eco,wine_peter$amp_response)
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), wine_peter)
fullmodel <- glm(amp_response ~ `Ecological origins`*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 wine_peter)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
anova(fullmodel, nullmodel,test = "Chisq")
# Accept null model as final
nohyb_wine_ephfinal <- fullmodel
summary(nohyb_wine_ephfinal)
anova(nohyb_wine_ephfinal, test = "Chisq")
anova(nohyb_wine_ephfinal, nullmodel, test = "Chisq")



### Mixed origin
mixed_or <- peter_eco[peter_eco$Clades == "8. Mixed origin",]
table(mixed_or$Ploidy,mixed_or$amp_response)
table(mixed_or$`Ecological origins`,mixed_or$amp_response)
table(mixed_or$wine_eco, mixed_or$amp_response)
nullmodel <- glm(amp_response_no1~1, family = binomial(link = "logit"), mixed_or)
fullmodel <- glm(amp_response_no1 ~ `Ecological origins`*polyploidy*heterozygosity, family = binomial(link = "logit"), mixed_or)
summary(fullmodel)
anova(fullmodel, test="Chisq")
nohet_mo_ephfinal <- fullmodel
summary(nohet_mo_ephfinal)
anova(nohet_mo_ephfinal, test="Chisq")
anova(nohet_mo_ephfinal, nullmodel,test="Chisq")

### Mosaic only
mosaic <- peter_eco[peter_eco$Clades == "M. Mosaic",]
table(mosaic$Ploidy,mosaic$amp_response)
table(mosaic$`Ecological origins`,mosaic$amp_response)
table(mosaic$`Ecological origins`,mosaic$polyploidy)
nullmodel <- glm(amp_response_no1~1, family = binomial(link="logit"), mosaic)
fullmodel <- glm(amp_response_no1 ~`Ecological origins`* 
                   polyploidy*heterozygosity, 
                 family = binomial(link = "logit"),
                 mosaic)
summary(fullmodel)
anova(fullmodel, test = "Chisq")
anova(fullmodel, nullmodel, test = "Chisq")
mosaic_ephfinal <- fullmodel
summary(mosaic_ephfinal)
anova(mosaic_ephfinal, test = "Chisq")
anova(mosaic_ephfinal, nullmodel, test = "Chisq")


#### Random 300 samples
petertab <- petertab[complete.cases(petertab$heterozygosity),]
petertab <- petertab[complete.cases(petertab$amp_response),]
peter_eco <- petertab[complete.cases(petertab$`Ecological origins`),]
# Remove chr loss and amplification of chr 1
peter_eco <- peter_eco[peter_eco$Aneuploidy_type == "Euploid" | 
                         peter_eco$Aneuploidy_type == "Gain" | 
                         peter_eco$Aneuploidy_type == "Both",]
peter_eco <- peter_eco[complete.cases(peter_eco$wine_eco)]
set.seed(12)
random_eco <- peter_eco[sample(nrow(peter_eco), 300),]
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), random_eco)
randmodel <- glm(amp_response ~ `Ecological origins`*
                   polyploidy*
                   heterozygosity, 
                 family = binomial(link = "logit"),
                 random_eco)
summary(randmodel)
anova(randmodel, test = "Chisq")
anova(randmodel, nullmodel,test = "Chisq")
cmodel <- randmodel
# remove eco:ploidy:het
# Remove eco:poly:het
nmodel <- update(cmodel, ~. - `Ecological origins`:polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove poly:het
nmodel <- update(cmodel, ~. - polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel<-nmodel
# Remove eco:het
nmodel <- update(cmodel, ~. - `Ecological origins`:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel
# remove eco:poly
nmodel <- update(cmodel, ~. - `Ecological origins`:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel
# remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel
# remove poly
nmodel <- update(cmodel, ~. -polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# remove eco
nmodel <- update(cmodel, ~. -`Ecological origins`)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# final
rand_eph_final <- cmodel
summary(rand_eph_final)
anova(rand_eph_final, test = "Chisq")
anova(rand_eph_final, nullmodel, test = "Chisq")
