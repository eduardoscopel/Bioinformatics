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
nrow(petertab)
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="Euploid","No","Yes")
petertab$amp_response <- ifelse(petertab$Aneuploidy_type == "Euploid","No","Yes")
petertab$amp_binary <- ifelse(petertab$amp_response == "Yes", 1,0)
petertab$admixture <- ifelse(petertab$Clades == "M. Mosaic" |
                               petertab$Clades == "8. Mixed origin",
                             "Admixed","Non-admixed")
petertab$polyploidy <- ifelse(petertab$Ploidy == 3 | 
                                petertab$Ploidy == 4 |
                                petertab$Ploidy == 5, "poly",
                              ifelse(petertab$Ploidy == 1, "haploid",
                                     "diploid"))
petertab$MD <- ifelse(petertab$Reference == 2, "Yes","No")
nrow(petertab)

### Factorize the categorical variables
petertab$amp_response <- as.factor(petertab$amp_response)
petertab$aneuploidy_response <- as.factor(petertab$aneuploidy_response)
petertab$`Ecological origins` <- as.factor(petertab$`Ecological origins`)
petertab$Ploidy <- as.factor(petertab$Ploidy)
petertab$new_eco <- as.factor(petertab$new_eco)
petertab$admixture <- as.factor(petertab$admixture)
petertab$wine_eco <- as.factor(petertab$wine_eco)

test$amp_response <- as.factor(test$amp_response)
petertab$aneuploidy_response <- as.factor(petertab$aneuploidy_response)
petertab$`Ecological origins` <- as.factor(petertab$`Ecological origins`)
petertab$Ploidy <- as.factor(petertab$Ploidy)
petertab$new_eco <- as.factor(petertab$new_eco)
petertab$admixture <- as.factor(petertab$admixture)
petertab$wine_eco <- as.factor(petertab$wine_eco)

# Remove HO deleted strains and strains with no clade, ploidy or heterozygosity information
nrow(petertab)
petertab1 <- petertab[petertab$`HO deletion` == "no",]
nrow(petertab1)
petertab2 <- petertab1[petertab1$MD == "No",]
nrow(petertab2)
petertab3 <- petertab2[complete.cases(petertab2$heterozygosity),]
nrow(petertab3)
petertab4 <- petertab3[complete.cases(petertab3$Clades),]
nrow(petertab4)
petertab5 <- petertab4[complete.cases(petertab4$Aneuploidies),]                                    # Remove strains with no aneuploidy information
nrow(petertab5)
petertab6 <- petertab5[petertab5$Aneuploidy_type == "Euploid" | 
                         petertab5$Aneuploidy_type == "Gain",]
nrow(petertab6)
# Remove strains from the Mosaic clade
petertab7 <- petertab6[petertab6$Clades != "M. Mosaic",]
nrow(petertab7)
#petertab8 <- petertab7[petertab7$Clades != "7. Mosaic beer",]
#nrow(petertab8)
#petertab7[petertab7$Ploidy == 1,"ID"]
petertab8 <- petertab7[petertab7$Ploidy !=1,]
nrow(petertab8)
petertab9 <- petertab8[petertab8$basename != "CBS1593",]
nrow(petertab9)
petertab10 <- petertab9[petertab9$basename != "CBS382",]
nrow(petertab10)
peter_clades <- petertab10
nrow(peter_clades)

# Remove strains that are heterozygous for SSD1
#petertab9 <- petertab8[petertab8$ssd1het == "hom",]
#nrow(petertab9)
#nonhetssd1strains <- petertab9$basename
#nonhetssd1strains
#write(nonhetssd1strains, "/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/ssd1/ssd1_nonhet_481_strains.txt",ncolumns = 1)


# Remove strains that are similar (pairwise genetic distance below 0.000007) 
tokeep <- fread("Documents/GitHub/eduardo/aneuploidy/peter/distmatrix/621strains/outfile.txt", header = FALSE)
petertab11 <- petertab10[petertab10$basename %in% tokeep$V1]
nrow(petertab11)

peter_clades <- petertab11
write.csv(peter_clades,"Documents/GitHub/eduardo/aneuploidy/peter/453strains.csv", row.names = FALSE)
test <- fread("Documents/GitHub/eduardo/aneuploidy/peter/453strains.csv", header = TRUE)



files_list <- fread("Documents/GitHub/eduardo/Strains /strain_list.txt", header = FALSE)
files_list <- files_list$V1
strains455 <- tokeep$V1
strain_list <- peter_clades$basename
new_list <- c()
for(i in 1:length(strain_list)){
  j <- which(strains455 == strain_list[i])
  print(paste("i is ", i))
  print(paste("j is ", j))
  new_list[i] <- paste(strains455[j], ".fa", sep = "")
}
write(new_list, "/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/scopel_etal_tree_455_strains.txt",ncolumns = 1)


### Look at the frequencies  
table(peter_clades$Ploidy,peter_clades$amp_response)
table(peter_clades$Clades, peter_clades$amp_response)
tapply(peter_clades$heterozygosity, peter_clades$amp_response, mean)
table(peter_clades$wine_eco,peter_clades$amp_response)

het_hist <- ggplot(peter_clades, aes(x=heterozygosity))
het_hist <- het_hist + geom_histogram(aes(fill=amp_response),
                                      color="black", 
                                      bins=100,
                                      alpha=0.7)+
  scale_x_continuous(breaks=seq(0,0.007,0.001), limits=c(0,0.007),expand=c(0.0001,0.0001)) +
  scale_y_continuous(breaks=seq(0,45,5), limits=c(0,45),expand=c(0.002,0)) +
  theme_bw() +
  ylab("Total strains") +
  xlab('Heterozygosity') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  theme(legend.position = "bottom")
png("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/455/het_hist.png", width = 6, height = 6, units = 'in', res = 300, bg='transparent')
het_hist
dev.off()


######################################################################################################################
############################## Clades, polyploidy, and het ###########################################################
######################################################################################################################
nrow(peter_clades)
table(peter_clades$Ploidy, peter_clades$amp_response)
table(peter_clades$Clades, peter_clades$amp_response)
cph_619s_null <- glm(amp_response~1, family = binomial(link="logit"), peter_clades)
cph_619s_full <- glm(amp_response ~ Clades*
                 polyploidy*
                 heterozygosity, 
               family = binomial(link = "logit"),
               peter_clades)
model1 <- cph_619s_full
summary(model1)
anova(model1, test = "Chisq")
anova(model1,cphnull, test  ="Chisq")
# Remove eco:poly:het
model2 <- update(model1, ~. -Clades:polyploidy:heterozygosity)
summary(model2)
anova(model2, test = "Chisq")
anova(model1,model2, test="Chisq")
# Keep model 2
# Remove polyploidy:het
model3 <- update(model2, ~. - polyploidy:heterozygosity)
summary(model3)
anova(model3, test = "Chisq")
anova(model2,model3, test= "Chisq")
# Keep model 3 
# Remove clades:poly
model4 <- update(model3, ~. - Clades:polyploidy)
summary(model4)
anova(model4, test = "Chisq")
anova(model3, model4, test = "Chisq")
# Keep model 4
# Remove clades:het
model5 <- update(model4, ~. -  Clades:heterozygosity)
summary(model5)
anova(model5, test = "Chisq")
anova(model4, model5, test = "Chisq")
# Keep model 5
# Remove het
model6 <- update(model5, ~. - heterozygosity)
summary(model6)
anova(model6, test = "Chisq")
anova(model5, model6, test = "Chisq")
# Keep model 6
# remove ploidy
model7 <- update(model6, ~. - polyploidy)
summary(model7)
anova(model7, test = "Chisq")
anova(model6, model7, test = "Chisq")
# Keep model 7
# Accept final model
cph_619s_final <- model7
# Final model
summary(cph_619s_final)
anova(cph_619s_final, test = "Chisq")
anova(cph_619s_final, cph_619s_null, test ="Chisq")

ssd1geno_model <- glm(amp_response ~ ssd1geno, 
                     family = binomial(link = "logit"),
                     peter_clades)
summary(ssd1geno_model)
anova(ssd1geno_model, test = "Chisq")
anova(ssd1geno_model, cph_619s_final, test = "Chisq")

cph_step <- step(cphfull)
step(cph_step)
forw <- step(cph_619s_full)
summary


### Remove possible hybrids (het > 0.4%)
peter_clades <- peter_clades[peter_clades$heterozygosity < 0.004,]
nullmodel <- glm(amp_response~1, family = binomial(link="logit"), peter_clades)
nohet_cphfull <- glm(amp_response ~ Clades*
                 polyploidy*
                 heterozygosity, 
               family = binomial(link = "logit"),
               peter_clades)
cmodel <- nohet_cphfull
summary(cmodel)
anova(cphfull, test = "Chisq")
anova(cmodel,nullmodel, test  ="Chisq")
# Remove eco:poly:het
nmodel <- update(cmodel, ~. -Clades:polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test="Chisq")
cmodel <- nmodel
# Remove polyploidy:het
nmodel <- update(cmodel, ~. - polyploidy:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel,nmodel, test= "Chisq")
cmodel <- nmodel
# Remove clades:het
nmodel <- update(cmodel, ~. -  Clades:heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# Remove clades:poly
nmodel <- update(cmodel, ~. - Clades:polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel
# Remove het
nmodel <- update(cmodel, ~. - heterozygosity)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel<-nmodel
# remove ploidy
nmodel <- update(cmodel, ~. - polyploidy)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
cmodel <- nmodel
# Remove clades
nmodel <- update(cmodel, ~. - Clades)
summary(nmodel)
anova(nmodel, test = "Chisq")
anova(cmodel, nmodel, test = "Chisq")
# Accept final model
no1_noL_nohet_cphfinal <- cmodel
# Final model
summary(no1_noL_nohet_cphfinal)
anova(no1_noL_nohet_cphfinal, test = "Chisq")
anova(no1_noL_nohet_cphfinal, nullmodel, test ="Chisq")



### Run analysis on 455 strains (excluding similar strains (gendist < 0.000007))
nrow(peter_clades)
table(peter_clades$Ploidy, peter_clades$amp_response)
prop.table(table(peter_clades$Clades, peter_clades$amp_response),1)
cph_455s_null <- glm(amp_response~1, family = binomial(link="logit"), test)
cph_455s_full <- glm(amp_response ~ Clades*
                       polyploidy*
                       heterozygosity, 
                     family = binomial(link = "logit"),
                     test)
model8 <- cph_455s_full
summary(model8)
anova(model8, test = "Chisq")
anova(model8,cph_455s_null, test  ="Chisq")
# Remove eco:poly:het
model9 <- update(model8, ~. -Clades:polyploidy:heterozygosity)
summary(model9)
anova(model9, test = "Chisq")
anova(model8,model9, test="Chisq")
# Keep model 9
# Remove polyploidy:het
model10 <- update(model9, ~. - polyploidy:heterozygosity)
summary(model10)
anova(model10, test = "Chisq")
anova(model9,model10, test= "Chisq")
# Keep model 10 
# Remove clades:poly
model11 <- update(model10, ~. - Clades:polyploidy)
summary(model11)
anova(model11, test = "Chisq")
anova(model10, model11, test = "Chisq")
# Keep model 11
# Remove clades:het
model12 <- update(model11, ~. -  Clades:heterozygosity)
summary(model12)
anova(model12, test = "Chisq")
anova(model11, model12, test = "Chisq")
# Keep model 12
# Remove het
model13 <- update(model12, ~. - heterozygosity)
summary(model13)
anova(model13, test = "Chisq")
anova(model12, model13, test = "Chisq")
# Keep model 13
# remove ploidy
model14 <- update(model13, ~. - polyploidy)
summary(model14)
anova(model14, test = "Chisq")
anova(model13, model14, test = "Chisq")
# Keep model 14
# Accept final model
cph_455s_final <- model14
# Final model
summary(cph_455s_final)
anova(cph_455s_final, test = "Chisq")
anova(cph_455s_final, cph_455s_null, test ="Chisq")

cph_455_step <- step(cph_455s_full)
step(cph_455_step)

ce_455s_full <- glm(amp_response ~ Clades*
                       `Ecological origins`, 
                     family = binomial(link = "logit"),
                     peter_clades)

ec_455s_full <- glm(amp_response ~ `Ecological origins`*
                      Clades, 
                    family = binomial(link = "logit"),
                    peter_clades)

cAe_455s_full <- glm(amp_response ~ Clades+
                      `Ecological origins`, 
                    family = binomial(link = "logit"),
                    peter_clades)
summary(cAe_455s_full)
anova(cAe_455s_full, test = "Chisq")
model15<- update(cAe_455s_full, ~. - `Ecological origins`)
summary(model15)
anova(model15, test = "Chisq")
anova(cAe_455s_full, model15, test = "Chisq")

eAc_455s_full <- glm(amp_response ~ `Ecological origins`+
                       Clades, 
                     family = binomial(link = "logit"),
                     peter_clades)
summary(eAc_455s_full)
anova(eAc_455s_full, test = "Chisq")
model16<- update(eAc_455s_full, ~. - `Ecological origins`)
summary(model16)
anova(model16, test = "Chisq")
anova(eAc_455s_full, model16, test = "Chisq")

model17 <- update(eAc_455s_full, ~. - Clades)
summary(model17)
anova(model17, test = "Chisq")
anova(eAc_455s_full, model17, test = "Chisq")


### Categorizing low (<20%), medium (20<=x<50%), high (>=50%) aneuploidy
peter_clades$an_category <- ifelse(peter_clades$Clades == "25. Sake", "High",
                                   ifelse(peter_clades$Clades == "8. Mixed origin" | 
                                            peter_clades$Clades == "9. Mexican agave" | 
                                            peter_clades$Clades == "11. Ale beer" |
                                            peter_clades$Clades == "2. Alpechin" |
                                            peter_clades$Clades == "5. French dairy" | 
                                            peter_clades$Clades == "6. African beer", "Medium","Low"))
cat_cph_455s_null <- glm(amp_response~1, family = binomial(link="logit"), peter_clades)
cat_cph_455s_full <- glm(amp_response ~ an_category*
                       polyploidy*
                       heterozygosity, 
                     family = binomial(link = "logit"),
                     peter_clades)
model18 <- cat_cph_455s_full
summary(model18)
anova(model18, test = "Chisq")
anova(model18,cat_cph_455s_null, test  ="Chisq")
# Remove eco:poly:het
model19 <- update(model18, ~. -an_category:polyploidy:heterozygosity)
summary(model19)
anova(model19, test = "Chisq")
anova(model18,model19, test="Chisq")
# Keep model 19
# Remove polyploidy:het
model20 <- update(model19, ~. - polyploidy:heterozygosity)
summary(model20)
anova(model20, test = "Chisq")
anova(model19,model20, test= "Chisq")
# Keep model 20 
# Remove clades:poly
model21 <- update(model20, ~. - an_category:polyploidy)
summary(model21)
anova(model21, test = "Chisq")
anova(model20, model21, test = "Chisq")
# Keep model 21
# Remove clades:het
model22 <- update(model21, ~. -  an_category:heterozygosity)
summary(model22)
anova(model22, test = "Chisq")
anova(model21, model22, test = "Chisq")
# Keep model 22
# Remove het
model23 <- update(model22, ~. - heterozygosity)
summary(model23)
anova(model23, test = "Chisq")
anova(model22, model23, test = "Chisq")
# Keep model 23
# remove ploidy
model24 <- update(model23, ~. - polyploidy)
summary(model24)
anova(model24, test = "Chisq")
anova(model23, model24, test = "Chisq")
# Keep model 14
# Accept final model
cat_cph_455s_final <- model24
# Final model
summary(cat_cph_455s_final)
anova(cat_cph_455s_final, test = "Chisq")
anova(cat_cph_455s_final, cat_cph_455s_null, test ="Chisq")
anova(cat_cph_455s_final, cph_455s_final, test = "Chisq")
