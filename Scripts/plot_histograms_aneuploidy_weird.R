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



### Plot frequencies 
ploidy_hist <- ggplot(sctab, aes(x=ploidy,fill=aneuploidy_binary))+ 
  scale_y_continuous(breaks=seq(0,1000,100))+
  geom_histogram(aes(),
                 stat="count",
                 color="black",
                 binwidth=2)+
  theme_bw()+
  theme(axis.text= element_text(face="bold",color="gray50", size=14))+
  ylab("Total strains") +
  xlab('Ploidy') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))
ploidy_hist

sctab <- within(sctab, ecology_category <- factor(ecology_category, levels = names(sort(table(ecology_category),decreasing = TRUE))))
eco_ploidy_hist <- ggplot(sctab, aes(x=ecology_category,fill=ploidy))+ 
  scale_y_continuous(breaks=seq(0,1000,100))+
  geom_histogram(aes(),
                 stat="count",
                 color="black",
                 binwidth=2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjus=0.5,hjust=1))+
  theme(axis.text= element_text(face="bold",color="gray50", size=14))+
  ylab("Total strains") +
  xlab('Ecology') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))
eco_ploidy_hist


md_hist <- ggplot(subset, aes(monosporic_derivative))+ 
  geom_histogram(aes(fill=aneuploidy_binary),
                 stat="count",
                 color="black",
                 binwidth=2) +
  scale_y_continuous(breaks=seq(0,1200,100))
md_hist

sctab <- within(sctab, ecology_category <- factor(ecology_category, levels = names(sort(table(ecology_category),decreasing = TRUE))))
ecology_hist <- ggplot(sctab, aes(x=ecology_category,fill=aneuploidy_binary))+ 
  geom_histogram(aes(),
                 stat="count",
                 color="black",
                 binwidth=2) + 
  scale_y_continuous(breaks=seq(0,300,50),expand=c(0.005,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,vjus=0.5,hjust=1))+
  theme(axis.text= element_text(face="bold",color="gray50", size=14))+
  ylab("Total strains") +
  xlab('Ecology') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))
ecology_hist

het_hist <- ggplot(sctab, aes(x=heterozygosity))
het_hist + geom_histogram(aes(fill = aneuploidy_binary), alpha = 1.0)+
  scale_x_continuous(breaks=seq(0,0.01,0.001), limits=c(0,0.01)) +
  scale_y_continuous(breaks=seq(0,275,25), limits=c(0,275)) +
  theme_bw() +
  ylab("Total strains") +
  xlab('Heterozygosity') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  theme(legend.position = "bottom")
het_hist+geom_density(aes(fill=aneuploidy_binary),alpha=0.4)+
  scale_x_continuous(breaks=seq(0,0.01,0.001), limits=c(0,0.01))
  

het_hist <- ggplot(sctab, aes(x=heterozygosity))
het_hist <- het_hist + geom_histogram(aes(fill=aneuploidy_binary),
                 color="black", 
                 bins=30,
                 alpha=0.7)+
  stat_bin(aes(y=..count../s))
  scale_x_continuous(breaks=seq(0,0.01,0.001), limits=c(0,0.01)) +
  scale_y_continuous(breaks=seq(0,125,25), limits=c(0,125)) +
  theme_bw() +
  ylab("Total strains") +
  xlab('Heterozygosity') +
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  theme(legend.position = "bottom")
het_hist

het_eco <- ggplot(sctab, aes(x=heterozygosity, fill=ecology_category))+ 
  geom_histogram(aes(),
                 color="black", 
                 bins=10)
het_eco
