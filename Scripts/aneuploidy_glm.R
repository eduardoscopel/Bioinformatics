library(data.table)
library(plyr)
library(scales)
### Read table and subset dataframe excluding unknown categories, and LOWCOV strains
sctab <- fread("sc_table.txt",header=TRUE)
sctab <- sctab[sctab$ecology_category != "Unknown",]
sctab <- sctab[sctab$ploidy != "Unknown",]
sctab <- sctab[sctab$ploidy != "LOWCOV",]
sctab <- sctab[sctab$aneuploidy != "LOWCOV",]
sctab <- sctab[sctab$aneuploidy_binary != "NA",]
sctab$aneuploidy_response <- ifelse(sctab$aneuploidy_binary =="Yes",1,0)
sctab$meiosis <- ifelse(sctab$ecology_category == "Beer" | 
                              sctab$ecology_category == "Bread" |
                              sctab$ecology_category == "Commercial" | 
                              sctab$ecology_category == "Sake","No",
                        ifelse(sctab$ecology_category == "Flower" |
                                     sctab$ecology_category == "Fruit" |
                                     sctab$ecology_category == "Insect" |
                                     sctab$ecology_category == "Tree" |
                                 sctab$ecology_category == "Wine" | 
                                 sctab$ecology_category == "Tree","Yes","Unknown"))
sctab$human_assoc <- ifelse(sctab$ecology_category == "Clinical","Clinical",
                        ifelse(sctab$ecology_category == "Flower" |
                                 sctab$ecology_category == "Fruit" |
                                 sctab$ecology_category == "Insect" |
                                 sctab$ecology_category == "Tree","Wild","Domesticated"))
sctab$eco_2 <- as.character(sctab$ecology_category)
sctab$eco_2[sctab$eco_2 == "Beer"] <- "Brewing"
sctab$eco_2[sctab$eco_2 == "Sake"] <- "Brewing"

sctab$eco_3 <- as.character(sctab$ecology_category)
sctab$eco_3[sctab$ecology_category == "Flower" | 
              sctab$ecology_category == "Fruit" |
              sctab$ecology_category == "Insect" |
              sctab$ecology_category == "Tree"] <- "Wild"

sctab$eco_4 <- as.character(sctab$ecology_category)
sctab$eco_4[sctab$ecology_category == "Beer"] <- "Brewing"
sctab$eco_4[sctab$ecology_category == "Sake"] <- "Brewing"
sctab$eco_4[sctab$ecology_category == "Flower" | 
              sctab$ecology_category == "Fruit" |
              sctab$ecology_category == "Insect" |
              sctab$ecology_category == "Tree"] <- "Wild"

sctab$eco_5 <- as.character(sctab$ecology_category)
sctab$eco_5[sctab$ecology_category == "Beer"] <- "Brewing"
sctab$eco_5[sctab$ecology_category == "Sake"] <- "Brewing"
sctab$eco_5[sctab$ecology_category == "Flower" | 
              sctab$ecology_category == "Fruit" |
              sctab$ecology_category == "Insect" |
              sctab$ecology_category == "Tree" | 
              sctab$ecology_category == "Wine"] <- "Seasonal"

sctab$eco_6 <- ifelse(sctab$ecology_category == "Flower" | 
              sctab$ecology_category == "Fruit" |
              sctab$ecology_category == "Insect" |
              sctab$ecology_category == "Tree" | 
              sctab$ecology_category == "Wine","Seasonal", ifelse(sctab$ecology_category == "Beer" | 
                                                                    sctab$ecology_category == "Sake","Brewing", "Other"))

sctab$eco_7 <- as.character(sctab$ecology_category)
sctab$eco_7[sctab$ecology_category == "Beer" | 
              sctab$ecology_category == "Sake" | 
              sctab$ecology_category == "Commercial"] <- "Non-seasonal"
sctab$eco_7[sctab$ecology_category == "Flower" | 
              sctab$ecology_category == "Fruit" |
              sctab$ecology_category == "Insect" |
              sctab$ecology_category == "Tree" | 
              sctab$ecology_category == "Wine"] <- "Seasonal"


sctab$eco_null <- as.character(sctab$ecology_category)
sctab$eco_null[sctab$ecology_category == "Beer" | 
              sctab$ecology_category == "Sake" |
              sctab$ecology_category == "Flower" | 
              sctab$ecology_category == "Fruit" |
              sctab$ecology_category == "Insect" |
              sctab$ecology_category == "Tree" | 
              sctab$ecology_category == "Wine"] <- "Null"


sctab <- sctab[sctab$heterozygosity < 0.01,]
sctab <- sctab[sctab$ecology_category %in% names(table(sctab$ecology_category))[table(sctab$ecology_category)>20],]
sctab <- sctab[sctab$ploidy %in% names(table(sctab$ploidy))[table(sctab$ploidy)>20],]

### Look at the frequencies  
tapply(sctab$aneuploidy_binary, sctab$ploidy,count)
tapply(sctab$aneuploidy_binary, sctab$ecology_category,count)
tapply(sctab$aneuploidy_binary, sctab$monosporic_derivative,count)
tapply(sctab$aneuploidy_binary, sctab$heterozygosity,count)
tapply(sctab$aneuploidy_binary, sctab$human_assoc,count)
tapply(sctab$aneuploidy_binary, sctab$meiosis,count)

### Factorize the categorical variables
sctab$aneuploidy_binary <- as.factor(sctab$aneuploidy_binary)
sctab$ecology_category <- as.factor(sctab$ecology_category)
sctab$monosporic_derivative <- as.factor(sctab$monosporic_derivative)
sctab$ploidy <- as.factor(sctab$ploidy)
sctab$human_assoc <- as.factor(sctab$human_assoc)
sctab$meiosis <- as.factor(sctab$meiosis)
sctab$eco_2 <- as.factor(sctab$eco_2)
sctab$eco_3 <- as.factor(sctab$eco_3)
sctab$eco_4 <- as.factor(sctab$eco_4)
sctab$eco_5 <- as.factor(sctab$eco_5)
sctab$eco_6 <- as.factor(sctab$eco_6)
sctab$eco_null <- as.character.factor(sctab$eco_null)
sctab$eco_7 <- as.factor(sctab$eco_7)

### Add a column with frequency of aneuploids
#sctab <- dcast(setDT(sctab), ecology_category+ploidy+monosporic_derivative+heterozygosity~aneuploidy_binary, length)

# filter for non-monosporic-derivatives
nonMD <- sctab[sctab$monosporic_derivative == "No",]

# filter for monospotic derivatives
MD <- sctab[sctab$monosporic_derivative == "Yes",]

# filter for triploids
haploids <- nonMD[nonMD$ploidy == 1,]

# filter for diploids
diploids <- nonMD[nonMD$ploidy == 2,]

# filter for triploids
triploids <- nonMD[nonMD$ploidy == 3,]

# filter for tetraploids
tetraploids <- nonMD[nonMD$ploidy == 4, ]

### Run glm on diploids (non-MD) this for binary response variable with ecology and heterozygosity as explanatory
model_dip1 <- glm(aneuploidy_response ~ heterozygosity*ecology_category-1, family = binomial(link='logit'), diploids)
summary(model_dip1)
model_dip2 <- update(model_dip1, ~. - heterozygosity:ecology_category)
summary(model_dip2)
anova(model_dip1,model_dip2, test="Chisq")
model_dip3 <- glm(aneuploidy_response ~ heterozygosity + ecology_category -1, family = binomial(),diploids)
summary(model_dip3)
anova(model_dip1,model_dip3, test="Chisq")
model_dip4 <- glm(aneuploidy_response ~ heterozygosity + ecology_category +I(heterozygosity^2),family = binomial(), diploids)
summary(model_dip4)
anova(model_dip3,model_dip4,test="Chisq")
model_dip5 <- glm(aneuploidy_response ~ ecology_category, family = binomial(),diploids)
summary(model_dip5)
anova(model_dip3,model_dip5,test="Chisq")
model_dip6 <- glm(aneuploidy_response ~ heterozygosity, family = binomial(),diploids)
summary(model_dip6)
anova(model_dip3,model_dip6,test="Chisq")

### Convert estimates to odds ratio with confidence intervals
coef_odds_ratio <-exp(cbind(OR = coef(model_dip3), confint(model_dip3)))

### Calculate predicted probability of aneuploidy at each ecological category
# Create new data set keeping mean heterozygosity as a constant
newdf1 <- with(diploids, data.frame(heterozygosity = mean(heterozygosity), ecology_category = levels(ecology_category)))
# Add new column to newdf1 with predicted probabilities using the glm model_dip3
newdf1$ecoP <- predict(model_dip3, newdata = newdf1, type = "response",se=TRUE)
# Create a new data frame with predicted probabilities and confidence intervals varying both het and ecology
newdf2 <- with(diploids, data.frame(heterozygosity = rep(seq(from = 0, to = 0.00674974, length.out = 100),14),
                                    ecology_category = factor(rep(levels(ecology_category), each = 100))))
newdf3 <- cbind(newdf2, predict(model_dip3, newdata = newdf2, type = "link",se=TRUE))
newdf3 <- within(newdf3,{
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
# Plot predicted probabilities of aneuploidy according to ecology maintaining heterozygosity constant
predicted_plot <- ggplot(newdf3, aes(x=heterozygosity, y=PredictedProb)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=ecology_category), alpha=0.2) + 
  scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1),expand=c(0.005,0)) +
  scale_x_continuous(limits = c(0,0.00674974), breaks=seq(0,0.00674974,0.001),expand=c(0.005,0))+
  theme_bw()+
  theme(axis.text= element_text(face="bold",color="gray50", size=14))+
  ylab("Predicted probability of aneuploidy") +
  xlab('Heterozygosity') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  geom_line(aes(colour=ecology_category), size=1)
predicted_plot


### Run glm on diploids (non-MD) this for binary response variable with meiosis and heterozygosity as explanatory
meiosis_m1 <- glm(aneuploidy_response ~ heterozygosity*meiosis-1, family = binomial(link='logit'), diploids)
summary(meiosis_m1)
meiosis_m2 <- update(meiosis_m1, ~. - heterozygosity:meiosis)
summary(meiosis_m2)
anova(meiosis_m1,meiosis_m2, test="Chisq")
meiosis_m3 <- glm(aneuploidy_response ~ heterozygosity + meiosis, family = binomial(),diploids)
summary(meiosis_m3)
anova(meiosis_m1,meiosis_m3, test="Chisq")
meiosis_m4 <- glm(aneuploidy_response ~ heterozygosity + meiosis +I(heterozygosity^2),family = binomial(), diploids)
summary(meiosis_m4)
anova(meiosis_m3,meiosis_m4,test="Chisq")
meiosis_m5 <- update(meiosis_m2, ~. - heterozygosity)
summary(meiosis_m5)
anova(meiosis_m2,meiosis_m5,test="Chisq")
meiosis_m6 <- update(meiosis_m2, ~. - meiosis)
summary(meiosis_m6)
anova(meiosis_m2,meiosis_m6,test="Chisq")

### Convert estimates to odds ratio with confidence intervals
coef_odds_ratio_meiosis <-exp(cbind(OR = coef(meiosis_m2), confint(meiosis_m2)))

### Calculate predicted probability of aneuploidy at each ecological category
# Create new data set keeping mean heterozygosity as a constant
newdf1_meio <- with(diploids, data.frame(heterozygosity = mean(heterozygosity), meiosis = levels(meiosis)))
# Add new column to newdf1 with predicted probabilities using the glm model_dip3
newdf1_meio$meioP <- predict(meiosis_m2, newdata = newdf1_meio, type = "response",se=TRUE)
# Create a new data frame with predicted probabilities and confidence intervals varying both het and ecology
newdf2_meio <- with(diploids, data.frame(heterozygosity = rep(seq(from = 0, to = 0.00674974, length.out = 100),3),
                                    meiosis = factor(rep(levels(meiosis), each = 100))))
newdf3_meio <- cbind(newdf2_meio, predict(meiosis_m2, newdata = newdf2_meio, type = "link",se=TRUE))
newdf3_meio <- within(newdf3_meio,{
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
# Plot predicted probabilities of aneuploidy according to ecology maintaining heterozygosity constant
predicted_plot_meio <- ggplot(newdf3_meio, aes(x=heterozygosity, y=PredictedProb)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=meiosis), alpha=0.2) + 
  scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1),expand=c(0.005,0)) +
  scale_x_continuous(limits = c(0,0.00674974), breaks=seq(0,0.00674974,0.001),expand=c(0.005,0))+
  theme_bw()+
  theme(axis.text= element_text(face="bold",color="gray50", size=14))+
  ylab("Predicted probability of aneuploidy") +
  xlab('Heterozygosity') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  geom_line(aes(colour=meiosis), size=1)
predicted_plot_meio
  
### Run glm on diploids (non-MD) this for binary response variable with domestication and heterozygosity as explanatory
human_assoc_m1 <- glm(aneuploidy_response ~ heterozygosity*human_assoc-1, family = binomial(link='logit'), diploids)
summary(human_assoc_m1)
human_assoc_m2 <- update(human_assoc_m1, ~. - heterozygosity:human_assoc)
summary(human_assoc_m2)
anova(human_assoc_m1,human_assoc_m2, test="Chisq")
human_assoc_m3 <- glm(aneuploidy_response ~ heterozygosity + human_assoc, family = binomial(),diploids)
summary(human_assoc_m3)
anova(human_assoc_m1,human_assoc_m3, test="Chisq")
human_assoc_m4 <- glm(aneuploidy_response ~ heterozygosity + human_assoc +I(heterozygosity^2),family = binomial(), diploids)
summary(human_assoc_m4)
anova(human_assoc_m3,human_assoc_m4,test="Chisq")
human_assoc_m5 <- update(human_assoc_m2, ~. - heterozygosity)
summary(human_assoc_m5)
anova(human_assoc_m2,human_assoc_m5,test="Chisq")
human_assoc_m6 <- update(human_assoc_m2, ~. - human_assoc)
summary(human_assoc_m6)
anova(human_assoc_m2,human_assoc_m6,test="Chisq")

### Convert estimates to odds ratio with confidence intervals
coef_odds_ratio_human_assoc <-exp(cbind(OR = coef(human_assoc_m2), confint(human_assoc_m2)))

### Calculate predicted probability of aneuploidy at each ecological category
# Create new data set keeping mean heterozygosity as a constant
newdf1_human_assoc <- with(diploids, data.frame(heterozygosity = mean(heterozygosity), human_assoc = levels(human_assoc)))
# Add new column to newdf1 with predicted probabilities using the glm model_dip3
newdf1_human_assoc$humanP <- predict(human_assoc_m2, newdata = newdf1_meio, type = "response",se=TRUE)
# Create a new data frame with predicted probabilities and confidence intervals varying both het and ecology
newdf2_human_assoc <- with(diploids, data.frame(heterozygosity = rep(seq(from = 0, to = 0.00674974, length.out = 100),3),
                                         human_assoc = factor(rep(levels(human_assoc), each = 100))))
newdf3_huamn_assoc <- cbind(newdf2_human_assoc, predict(human_assoc_m2, newdata = newdf2_human_assoc, type = "link",se=TRUE))
newdf3_huamn_assoc <- within(newdf3_huamn_assoc,{
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
# Plot predicted probabilities of aneuploidy according to ecology maintaining heterozygosity constant
predicted_plot_human <- ggplot(newdf3_huamn_assoc, aes(x=heterozygosity, y=PredictedProb)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=human_assoc), alpha=0.2) + 
  scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1),expand=c(0.005,0)) +
  scale_x_continuous(limits = c(0,0.00674974), breaks=seq(0,0.00674974,0.001),expand=c(0.005,0))+
  theme_bw()+
  theme(axis.text= element_text(face="bold",color="gray50", size=14))+
  ylab("Predicted probability of aneuploidy") +
  xlab('Heterozygosity') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  geom_line(aes(colour=human_assoc), size=1)
predicted_plot_human


### Brewing test (13 categories)
attempt1 <- glm(aneuploidy_response ~ heterozygosity*eco_2-1, family = binomial(link='logit'), diploids)
summary(attempt1)
attempt2 <- update(attempt1, ~. - heterozygosity:eco_2)
summary(attempt2)
anova(attempt1,attempt2, test="Chisq")
attempt3 <- glm(aneuploidy_response ~ heterozygosity + eco_2 -1, family = binomial(),diploids)
summary(attempt3)
anova(attempt1,attempt3, test="Chisq")
attempt4 <- glm(aneuploidy_response ~ heterozygosity + eco_2 +I(heterozygosity^2),family = binomial(), diploids)
summary(attempt4)
anova(attempt3,attempt4,test="Chisq")
attempt5 <- glm(aneuploidy_response ~ eco_2, family = binomial(),diploids)
summary(attempt5)
anova(attempt1,attempt5,test="Chisq")
attempt6 <- glm(aneuploidy_response ~ heterozygosity, family = binomial(),diploids)
summary(attempt6)
anova(attempt1,attempt6,test="Chisq")

### Wild test (11 categories)
attemptb1 <- glm(aneuploidy_response ~ heterozygosity*eco_3-1, family = binomial(link='logit'), diploids)
summary(attemptb1)
attemptb2 <- update(attemptb1, ~. - heterozygosity:eco_3)
summary(attemptb2)
anova(attemptb1,attemptb2, test="Chisq")
attemptb3 <- glm(aneuploidy_response ~ heterozygosity + eco_3 -1, family = binomial(),diploids)
summary(attemptb3)
anova(attemptb1,attemptb3, test="Chisq")
attemptb4 <- glm(aneuploidy_response ~ heterozygosity + eco_3 +I(heterozygosity^2),family = binomial(), diploids)
summary(attemptb4)
anova(attemptb3,attemptb4,test="Chisq")
attemptb5 <- glm(aneuploidy_response ~ eco_3, family = binomial(),diploids)
summary(attemptb5)
anova(attemptb1,attemptb5,test="Chisq")
attemptb6 <- glm(aneuploidy_response ~ heterozygosity, family = binomial(),diploids)
summary(attemptb6)
anova(attemptb1,attemptb6,test="Chisq")

### Wild + Brewing test (10 categories)
attemptc1 <- glm(aneuploidy_response ~ heterozygosity*eco_4-1, family = binomial(link='logit'), diploids)
summary(attemptc1)
attemptc2 <- update(attemptc1, ~. - heterozygosity:eco_4)
summary(attemptc2)
anova(attemptc1,attemptc2, test="Chisq")
attemptc3 <- glm(aneuploidy_response ~ heterozygosity + eco_4 -1, family = binomial(),diploids)
summary(attemptc3)
anova(attemptc1,attemptc3, test="Chisq")
attemptc4 <- glm(aneuploidy_response ~ heterozygosity + eco_4 +I(heterozygosity^2),family = binomial(), diploids)
summary(attemptc4)
anova(attemptc3,attemptc4,test="Chisq")
attemptc5 <- glm(aneuploidy_response ~ eco_4, family = binomial(),diploids)
summary(attemptc5)
anova(attemptc1,attemptc5,test="Chisq")
attemptc6 <- glm(aneuploidy_response ~ heterozygosity, family = binomial(),diploids)
summary(attemptc6)
anova(attemptc1,attemptc6,test="Chisq")

### Seasonal/brewing test (9 categories)
attemptd1 <- glm(aneuploidy_response ~ heterozygosity*eco_5-1, family = binomial(link='logit'), diploids)
summary(attemptd1)
attemptd2 <- update(attemptd1, ~. - heterozygosity:eco_5)
summary(attemptd2)
anova(attemptd1,attemptd2, test="Chisq")
attemptd3 <- glm(aneuploidy_response ~ heterozygosity + eco_5 -1, family = binomial(),diploids)
summary(attemptd3)
anova(attemptd1,attemptd3, test="Chisq")
attemptd4 <- glm(aneuploidy_response ~ heterozygosity + eco_5 +I(heterozygosity^2),family = binomial(), diploids)
summary(attemptd4)
anova(attemptd3,attemptd4,test="Chisq")
attemptd5 <- glm(aneuploidy_response ~ eco_5, family = binomial(),diploids)
summary(attemptd5)
anova(attemptd1,attemptd5,test="Chisq")
attemptd6 <- glm(aneuploidy_response ~ heterozygosity, family = binomial(),diploids)
summary(attemptd6)
anova(attemptd1,attemptd6,test="Chisq")

### Convert estimates to odds ratio with confidence intervals
coef_odds_ratio_attempt4 <-exp(cbind(OR = coef(attemptd3), confint(attemptd3)))

### Calculate predicted probability of aneuploidy at each ecological category
# Create new data set keeping mean heterozygosity as a constant
newdf1_attemptd3 <- with(diploids, data.frame(heterozygosity = mean(heterozygosity), eco_5 = levels(eco_5)))
# Add new column to newdf1 with predicted probabilities using the glm model_dip3
newdf1_attemptd3$ecoP <- predict(attemptd3, newdata = newdf1_attemptd3, type = "response",se=TRUE)
# Create a new data frame with predicted probabilities and confidence intervals varying both het and ecology
newdf2_attemptd3 <- with(diploids, data.frame(heterozygosity = rep(seq(from = 0, to = 0.00674974, length.out = 100),3),
                                              eco_5 = factor(rep(levels(eco_5), each = 100))))
newdf3_attemptd3 <- cbind(newdf2_attemptd3, predict(attemptd3, newdata = newdf2_attemptd3, type = "link",se=TRUE))
newdf3_attemptd3 <- within(newdf3_attemptd3,{
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
# Plot predicted probabilities of aneuploidy according to ecology maintaining heterozygosity constant
for(i in seq(1,9,1)){
  predicted_plot_attemptd <- ggplot(newdf3_attemptd3[(1+(100*(i-1))) : (100*(i)),], aes(x=heterozygosity, y=PredictedProb)) +
    geom_ribbon(aes(ymin = LL, ymax = UL, fill=eco_5), alpha=0.2) + 
    scale_fill_manual(values=cbbPalette[i])+
    scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1),expand=c(0.005,0)) +
    scale_x_continuous(limits = c(0,0.00674974), breaks=seq(0,0.00674974,0.001),expand=c(0.005,0))+
    theme_bw()+
    theme(axis.text= element_text(face="bold",color="gray50", size=14))+
    #ylab("Predicted probability of aneuploidy") +
    #xlab('Heterozygosity') +
    #theme(axis.title.x = element_text(face="bold", size=14))+
    #theme(axis.title.y = element_text(face="bold", size=14))+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_blank())+
    labs(title = paste0(newdf3_attemptd3[2+100*(i-1),2]," (n= ",length(diploids$eco_5[diploids$eco_5 == newdf3_attemptd3[2+100*(i-1),2]]),")"))+
    theme(plot.title = element_text(face="bold",h=.5,size=24))+
    geom_line(aes(colour=eco_5), size=1)+
    scale_colour_manual(values=cbbPalette[i])+
    theme(legend.position = "none")
  tiff(paste0("plot",i,".tiff"))
  print(predicted_plot_attemptd)
  dev.off()
}
predicted_plot_attemptd <- ggplot(newdf3_attemptd3[newdf3_attemptd3$eco_5 == "Brewing" | newdf3_attemptd3$eco_5 == "Seasonal",], aes(x=heterozygosity, y=PredictedProb))+
  geom_point(data=diploids[diploids$eco_5 == "Brewing" | diploids$eco_5 == "Seasonal",], 
           aes(x=heterozygosity,y=aneuploidy_response,color=eco_5, alpha=0.5,size=1),
             position= position_jitter(width = 0.0002,height = 0))+
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=eco_5), alpha=0.2) + 
  scale_fill_manual(values=c(cbbPalette[3],cbbPalette[8]))+
  scale_y_continuous(name="Predicted probability of aneuploidy", sec.axis = sec_axis(~.,name="Aneuploidy (yes = 1; no = 0)",breaks = seq(0,1,1)),limits = c(0,1), breaks=seq(0,1,0.1),expand=c(0.005,0)) +
  scale_x_continuous(limits = c(0,0.00674974), breaks=seq(0,0.00674974,0.0005),expand=c(0.005,0))+
  theme_bw()+
  theme(axis.text= element_text(face="bold",color="gray50", size=16))+
  #ylab("Predicted probability of aneuploidy") +
  xlab('Heterozygosity') +
  theme(axis.title.x = element_text(face="bold", size=16))+
  theme(axis.title.y = element_text(face="bold", size=16))+
  labs(title = "Brewing (n = 71) and Seasonal (n = 343)")+
  theme(plot.title = element_text(face="bold",h=.5,size=20))+
  geom_line(aes(colour=eco_5), size=1)+
  scale_colour_manual(values=c(cbbPalette[3], cbbPalette[8]))+
  theme(legend.position = "none")
#png("glm_seasonal_brewing.png")
predicted_plot_attemptd
predicted_plot_attemptd + 
  geom_point(data=diploids[diploids$eco_5 == "Brewing" | diploids$eco_5 == "Seasonal",], 
             aes(x=heterozygosity,y=aneuploidy_response,color=eco_5,size=1))+
  geom_jitter()
#dev.off()
ggplot(diploids, aes(x=heterozygosity, y=aneuploidy_response)) + geom_point()

### Seasonal/non-seasonal test (7 categories)
sns1 <- glm(aneuploidy_response ~ heterozygosity*eco_7-1, family = binomial(link='logit'), diploids)
summary(sns1)
sns2 <- update(sns1, ~. - heterozygosity:eco_7)
summary(sns2)
anova(sns1,sns2, test="Chisq")
sns3 <- update(sns2, ~. - heterozygosity)
summary(sns3)
anova(sns2,sns3, test="Chisq")
sns4 <- update(sns2, ~. - eco_7)
summary(sns4)
anova(sns2,sns4, test="Chisq")


### Seasonal, brewing, other test (3 categories)
attempte1 <- glm(aneuploidy_response ~ heterozygosity*eco_6-1, family = binomial(link='logit'), diploids)
summary(attempte1)
attempte2 <- update(attempte1, ~. - heterozygosity:eco_6)
summary(attempte2)
anova(attempte1,attempte2, test="Chisq")
attempte3 <- glm(aneuploidy_response ~ heterozygosity + eco_6 -1, family = binomial(),diploids)
summary(attempte3)
anova(attempte1,attempte3, test="Chisq")
attempte4 <- glm(aneuploidy_response ~ heterozygosity + eco_6 +I(heterozygosity^2),family = binomial(), diploids)
summary(attempte4)
anova(attempte3,attempte4,test="Chisq")
attempte5 <- glm(aneuploidy_response ~ eco_6, family = binomial(),diploids)
summary(attempte5)
anova(attempte1,attempte5,test="Chisq")
attempte6 <- glm(aneuploidy_response ~ heterozygosity, family = binomial(),diploids)
summary(attempte6)
anova(attempte1,attempte6,test="Chisq")


fligner.test(diploids$aneuploidy_response, 
             diploids$eco_5)
fligner.test(diploids$aneuploidy_response[diploids$eco_5 != "Bioethanol" 
                                          & diploids$eco_5 != "Bread"], 
             diploids$eco_5[diploids$eco_5 != "Bioethanol"
                            & diploids$eco_5 != "Bread"])
fligner.test(diploids$aneuploidy_response[diploids$eco_5 != "Bioethanol" 
                                          & diploids$eco_5 != "Bread"
                                          & diploids$eco_5 != "Brewing"], 
             diploids$eco_5[diploids$eco_5 != "Bioethanol"
                            & diploids$eco_5 != "Bread" 
                            & diploids$eco_5 != "Brewing"])
fligner.test(diploids$aneuploidy_response[diploids$eco_5 != "Bioethanol" 
                                          & diploids$eco_5 != "Bread"
                                          & diploids$eco_5 != "Brewing"
                                          & diploids$eco_5 != "Clinical"], 
             diploids$eco_5[diploids$eco_5 != "Bioethanol"
                            & diploids$eco_5 != "Bread" 
                            & diploids$eco_5 != "Brewing"
                            & diploids$eco_5 != "Clinical"])
fligner.test(diploids$aneuploidy_response[diploids$eco_5 != "Bioethanol" 
                                          & diploids$eco_5 != "Bread"
                                          & diploids$eco_5 != "Brewing"
                                          & diploids$eco_5 != "Clinical"
                                          & diploids$eco_5 != "Commercial"], 
             diploids$eco_5[diploids$eco_5 != "Bioethanol"
                            & diploids$eco_5 != "Bread" 
                            & diploids$eco_5 != "Brewing"
                            & diploids$eco_5 != "Clinical"
                            & diploids$eco_5 != "Commercial"])
fligner.test(diploids$aneuploidy_response[diploids$eco_5 != "Bioethanol" 
                                          & diploids$eco_5 != "Bread"
                                          & diploids$eco_5 != "Clinical"
                                          & diploids$eco_5 != "Commercial"
                                          & diploids$eco_5 != "Dairy"
                                          & diploids$eco_5 != "Fermentation"], 
             diploids$eco_5[diploids$eco_5 != "Bioethanol"
                            & diploids$eco_5 != "Bread" 
                            & diploids$eco_5 != "Clinical"
                            & diploids$eco_5 != "Commercial"
                            & diploids$eco_5 != "Dairy"
                            & diploids$eco_5 != "Fermentation"])
fligner.test(diploids$aneuploidy_response[diploids$eco_5 != "Bioethanol" 
                                          & diploids$eco_5 != "Bread"
                                          & diploids$eco_5 != "Clinical"
                                          & diploids$eco_5 != "Commercial"
                                          & diploids$eco_5 != "Dairy"
                                          & diploids$eco_5 != "Fermentation"
                                          & diploids$eco_5 != "SSF"], 
             diploids$eco_5[diploids$eco_5 != "Bioethanol"
                            & diploids$eco_5 != "Bread" 
                            & diploids$eco_5 != "Clinical"
                            & diploids$eco_5 != "Commercial"
                            & diploids$eco_5 != "Dairy"
                            & diploids$eco_5 != "Fermentation"
                            & diploids$eco_5 != "SSF"])
fligner.test(diploids$aneuploidy_response[diploids$eco_5 != "Brewing" 
                                          & diploids$eco_5 != "Seasonal"
                                          & diploids$eco_5 != "SSF"], 
             diploids$eco_5[diploids$eco_5 != "Brewing"
                            & diploids$eco_5 != "Seasonal"
                            & diploids$eco_5 != "SSF"])
fligner.test(diploids$aneuploidy_response[diploids$eco_5 == "Brewing" 
                                          | diploids$eco_5 == "Seasonal"], 
             diploids$eco_5[diploids$eco_5 == "Brewing"
                            | diploids$eco_5 == "Seasonal"])


### Seasonal test with ploidy (9 categories)
ploidy_test1 <- glm(aneuploidy_response ~ heterozygosity*ploidy*eco_5, family = quasibinomial(link='logit'), nonMD)
summary(ploidy_test1)
ploidy_test2 <- update(ploidy_test1, ~. - heterozygosity:ploidy:eco_5)
summary(ploidy_test2)
anova(ploidy_test1,ploidy_test2, test="Chisq")
ploidy_test3 <- update(ploidy_test2, ~. - heterozygosity:ploidy)
summary(ploidy_test3)
anova(ploidy_test2,ploidy_test3, test="Chisq")
ploidy_test4 <- update(ploidy_test3, ~. - ploidy:eco_5)
summary(ploidy_test4)
anova(ploidy_test3,ploidy_test4, test="Chisq")
ploidy_test5 <- update(ploidy_test4, ~. - heterozygosity:eco_5)
summary(ploidy_test5)
anova(ploidy_test4,ploidy_test5, test="Chisq")
ploidy_test6 <- update(ploidy_test4, ~. - ploidy)
summary(ploidy_test6)
anova(ploidy_test4,ploidy_test6, test="Chisq")
ploidy_test7 <- update(ploidy_test4, ~. - heterozygosity)
summary(ploidy_test7)
anova(ploidy_test4,ploidy_test7, test="Chisq")
ploidy_test8 <- update(ploidy_test4, ~. - eco_5)
summary(ploidy_test8)
anova(ploidy_test4,ploidy_test8, test="Chisq")

### Seasonal non-seasonal as null test (8 categories)
null1 <- glm(aneuploidy_response ~ heterozygosity*eco_null-1, family = binomial(link='logit'), diploids)
summary(null1)
null2 <- update(null1, ~. - heterozygosity:eco_null)
summary(null2)
anova(null1,null2, test="Chisq")
null3 <- glm(aneuploidy_response ~ heterozygosity + eco_null +I(heterozygosity^2),family = binomial(), diploids)
summary(null3)
anova(null2,null3,test="Chisq")
null4 <- update(null2, ~. - heterozygosity)
summary(null4)
anova(null3,null4,test="Chisq")
null5 <- update(null2, ~. - eco_null)
summary(null5)
anova(null2,null5,test="Chisq")

anova(null2, attemptd3, test = "Chisq")
