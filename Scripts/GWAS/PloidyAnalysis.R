library(data.table)
library(ggplot2)
library(fdrtool)
library(multtest)

### Read table and subset dataframe excluding unknown categories, and LOWCOV strains
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





#########################################################################################################################################################
############################################################# PLOIDY ANALYSIS ###########################################################################
#########################################################################################################################################################

# Create a ploidy matrix without haploids
#nonhap <- petertab[petertab$Ploidy != 1,]
#nonhap <- nonhap[nonhap$Aneuploidy_type != "Euploid",]
#nrow(nonhap)
ploidymatrix <- as.matrix(table(df$ploidyCons,df$aneuploidy_type), header = TRUE)
ploidymatrix <- ploidymatrix[,c(2,3,1,4)]
ploidymatrix
# add a column of totals/proportions to the matrix
for(i in 1:nrow(ploidymatrix)){
  if(i==1){
    total_ploidy <- c(sum(ploidymatrix[i,]))
    proportion_gain <- c((ploidymatrix[i,1])/(total_ploidy[i]))
    proportion_loss <- c((ploidymatrix[i,2])/(total_ploidy[i]))
    proportion_both <- c((ploidymatrix[i,3])/(total_ploidy[i]))
  }
  else{
    total_ploidy <- c(total_ploidy,sum(ploidymatrix[i,]))
    proportion_gain <- c(proportion_gain,(ploidymatrix[i,1])/(total_ploidy[i]))
    proportion_loss <- c(proportion_loss,(ploidymatrix[i,2])/(total_ploidy[i]))
    proportion_both <- c(proportion_both,(ploidymatrix[i,3])/(total_ploidy[i]))
  }
}
ploidymatrix <- cbind(ploidymatrix,total_ploidy)
ploidymatrix <- cbind(ploidymatrix, proportion_gain)
ploidymatrix <- cbind(ploidymatrix, proportion_loss)
ploidymatrix <- cbind(ploidymatrix, proportion_both)


# set up pairwise Fisher's tests
odds_ratio <- matrix(
  c(rep(0,2*nrow(ploidymatrix))),
  nrow = 2,
  dimnames = list(c("Gain", "Loss"),
                  c("Hap","Dip", "Trip","Tet","Pent","Unk")))
p_values <- matrix(
  c(rep(0,2*nrow(ploidymatrix))),
  nrow = 2,
  dimnames = list(c("Gain", "Loss"),
                  c("Hap","Dip", "Trip","Tet","Pent","Unk")))
fdr <- p_values <- matrix(
  c(rep(0,2*nrow(ploidymatrix))),
  nrow = 2,
  dimnames = list(c("Gain", "Loss"),
                  c("Hap","Dip", "Trip","Tet","Pent","Unk")))

for(i in 1:nrow(ploidymatrix)){
  odds_ratio[1,i] <- as.numeric(fisher.test(
    matrix(c(ploidymatrix[i,1]+ploidymatrix[i,3],
             ploidymatrix[i,2]+ploidymatrix[i,4],
             sum(ploidymatrix[,1])+sum(ploidymatrix[,3])-ploidymatrix[i,1]-ploidymatrix[i,3],
             sum(ploidymatrix[,2])+sum(ploidymatrix[,4])-ploidymatrix[i,2]-ploidymatrix[i,4]),
           nrow=2,
           dimnames = list(c("Gain","no gain"),
                           c("Haploids","All")))
  )$estimate)
  p_values[1,i] <- fisher.test(
    matrix(c(ploidymatrix[i,1]+ploidymatrix[i,3],
             ploidymatrix[i,2]+ploidymatrix[i,4],
             sum(ploidymatrix[,1])+sum(ploidymatrix[,3])-ploidymatrix[i,1]-ploidymatrix[i,3],
             sum(ploidymatrix[,2])+sum(ploidymatrix[,4])-ploidymatrix[i,2]-ploidymatrix[i,4]),
           nrow=2,
           dimnames = list(c("Gain","no gain"),
                           c("Haploids","All")))
  )$p.value
  odds_ratio[2,i] <- as.numeric(fisher.test(
    matrix(c(ploidymatrix[i,2]+ploidymatrix[i,3],
             ploidymatrix[i,1]+ploidymatrix[i,4],
             sum(ploidymatrix[,2]+ploidymatrix[,3])-ploidymatrix[i,2]-ploidymatrix[i,3],
             sum(ploidymatrix[,1])+sum(ploidymatrix[,4])-ploidymatrix[i,1]-ploidymatrix[i,4]),
           nrow=2,
           dimnames = list(c("Loss","no loss"),
                           c(i+1,"All")))
  )$estimate)
  p_values[2,i] <- fisher.test(
    matrix(c(ploidymatrix[i,2]+ploidymatrix[i,3],
             ploidymatrix[i,1]+ploidymatrix[i,4],
             sum(ploidymatrix[,2]+ploidymatrix[,3])-ploidymatrix[i,2]+ploidymatrix[i,3],
             sum(ploidymatrix[,1])+sum(ploidymatrix[,4])-ploidymatrix[i,1]+ploidymatrix[i,4]),
           nrow=2,
           dimnames = list(c("Loss","no loss"),
                           c(i+1,"All")))
  )$p.value
  
}

# Set up variables for ploidy proportion plot
pLabels = c()
pcat=row.names(ploidymatrix)

# run a Fisher's Exact Test on each row, and prepare lists to create the odds ratio plot
for(row in 1:nrow(ploidymatrix)){
  if(pcat[row]==1){
    pcat[row]<- "Haploids"
    pLabels <- append(pLabels,paste(pcat[row]," (",ploidymatrix[row,5],")",sep=""))
  }
  if(pcat[row]==2){
    pcat[row]<- "Diploids"
    pLabels <- append(pLabels,paste(pcat[row]," (",ploidymatrix[row,5],")",sep=""))
  }
  else if(pcat[row]==3){
    pcat[row]<- "Triploids"
    pLabels <- append(pLabels,paste(pcat[row]," (",ploidymatrix[row,5],")",sep=""))
  }
  else if(pcat[row]==4){
    pcat[row]<- "Tetraploids"
    pLabels <- append(pLabels,paste(pcat[row]," (",ploidymatrix[row,5],")",sep=""))
  }
  else if(pcat[row]==5){
    pcat[row]<- "Pentaploids"
    pLabels <- append(pLabels,paste(pcat[row]," (",ploidymatrix[row,5],")",sep=""))
  }
  else if(pcat[row]==7){
    pcat[row]<- "Heptaploids"
    pLabels <- append(pLabels,paste(pcat[row]," (",ploidymatrix[row,5],")",sep=""))
  }
  else if(pcat[row]=="Unknown"){
    pcat[row]<- "Unknown"
    pLabels <- append(pLabels,paste(pcat[row]," (",ploidymatrix[row,5],")",sep=""))
  }
}

#p_df <- data.frame(
#  yAxis = pLabels,
#  pgain = ploidymatrix[,6], ploss = ploidymatrix[,7])

ploidydf <- data.frame(xAxis = seq(1,6,1),
                       ploidy = c(rownames(ploidymatrix),rownames(ploidymatrix),rownames(ploidymatrix),rownames(ploidymatrix)),
                       type = c(rep("Gain",6),rep("Loss",6),rep("Both",6),rep("Both",6)), 
                       prop = c(ploidymatrix[,6],-ploidymatrix[,7],ploidymatrix[,8],-ploidymatrix[,8]),
                       odds = c(odds_ratio[1,],odds_ratio[2,],odds_ratio[1,],odds_ratio[2,]),
                       p = c(p_values[1,],p_values[2,],p_values[1,],p_values[2,]))
ploidydf$fdr[mt.rawp2adjp(ploidydf$p, proc = "BH", alpha = 0.05)$index] <- mt.rawp2adjp(ploidydf$p, proc = "BH", alpha = 0.05)$adjp[,2]
ploidydf$sig <- ifelse(ploidydf$fdr <=  0.05, "*",
                       ifelse(ploidydf$p <= 0.05 &
                                ploidydf$fdr <=  0.13, "+",NA))
ploidydf$colour <- ifelse(ploidydf$odds > 1 &
                            !is.na(ploidydf$sig) &
                            ploidydf$type != "Both","#e88446",
                          ifelse(ploidydf$odds > 1 &
                                   !is.na(ploidydf$sig) &
                                   ploidydf$type == "Both", "#F6C5A4",
                                 ifelse(ploidydf$odds < 1 &
                                          !is.na(ploidydf$sig) &
                                          ploidydf$type != "Both", "#0398CE",
                                        ifelse(is.na(ploidydf$sig) &
                                                 ploidydf$type != "Both", "#8EB6D4","#E1EFF9"))))
ploidydf$xAxis <- factor(ploidydf$xAxis)
ploidydf$lab <- pLabels
ploidydf$lab[7:24] <- NA
ploidydf$sig[13:24] <- NA
ploidydf$labjust <- rep(0,24)
ploidydf$labjust[1:6] <- c(-0.2,-0.7,-0.6,-2,-0.2,-1)
ploidydf$labjust[7:24] <- NA
ploidydf$sigjust <- c(0,0,-1.6,-9.5,0,0,0,2,0,11.5,1.5,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)

# create ploidy fisher plot
fisherplot<- ggplot(ploidydf, aes(x= xAxis, y = prop)) + 
  geom_rect(aes(ymin= (-sum(ploidymatrix[,2]) - sum(ploidymatrix[,3]))/sum(ploidymatrix[,5]), ymax=(sum(ploidymatrix[,1]) + sum(ploidymatrix[,3]))/sum(ploidymatrix[,5]),xmin=-Inf,xmax=Inf),
            color="white", 
            alpha=0.01)+
  geom_col(fill=ploidydf$colour, width = 0.5, colour = "darkgray") + 
  geom_text(aes(label = sig, vjust = sigjust),size=7)+
  geom_text(aes(label = lab, vjust = 0.5, angle=90, hjust=labjust),size=3.5)+
  scale_y_continuous(breaks = seq(-0.5,0.6,0.1), 
                     labels = c("50%","40%","30%","20%","10%","0%","10%","20%","30%","40%","50%","60%"), 
                     limits = c(-0.5,1)) + 
  ylab("% Strains with Chr\n\nLosses                             Gains")+
  geom_hline(yintercept = 0, size=.2)+
  geom_segment(x = 0.4, y = -0.5, xend = 0.4, yend = 0.6,size=.2)+
  theme_bw()+
  theme(aspect.ratio = 5/2,
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.title.y = element_text(hjust=0.3, vjust=2),
        axis.line.y.left = element_line(size=0, linetype = "solid"),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(size=.5, linetype = "dashed"),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank())

fisherplot

png("Documents/GitHub/eduardo/aneuploidy/plots/peter/Ploidy/Fig1C_ploidy_chr1.png",width = 3, height = 5, units = "in",bg = "white", res=300)
fisherplot
dev.off()
postscript("Documents/GitHub/eduardo/aneuploidy/plots/peter/Ploidy/Fig1C_ploidy.ps")
fisherplot
dev.off()
