library(data.table)
library(ggplot2)
library(fdrtool)
library(multtest)

### Read table and subset dataframe excluding unknown categories, and LOWCOV strains
petertab <- fread("/Users/es47540/Documents/GitHub/eduardo/Strains /peter_new(11122020).txt",header=TRUE)
nrow(petertab)
petertab <- petertab[complete.cases(petertab$Aneuploidies),]     
#petertab <- petertab[complete.cases(petertab$heterozygosity),]# Remove strains with no aneuploidy information
petertab$aneuploidy_response <- ifelse(petertab$Aneuploidies =="euploid","No","Yes")
petertab$amp_response_no1 <- ifelse(petertab$Aneuploidy_type =="Loss" | 
                                  petertab$Aneuploidy_type =="Euploid" |
                                    petertab$Aneuploidy_type == "1only","No","Yes")
petertab$polyploidy <- ifelse(petertab$Ploidy == 1 |
                                petertab$Ploidy == 2,"1-2n","3-5n")





#########################################################################################################################################################
############################################################# PLOIDY ANALYSIS ###########################################################################
#########################################################################################################################################################

# Create a ploidy matrix without haploids
#nonhap <- petertab[petertab$Ploidy != 1,]
#nonhap <- nonhap[nonhap$Aneuploidy_type != "Euploid",]
#nrow(nonhap)
ploidymatrix <- as.matrix(table(petertab$Ploidy,petertab$Aneuploidy_type), nrow=4, header = TRUE)
ploidymatrix <- ploidymatrix[,c(5,6,3,4,1,2)]
ploidymatrix[,1] <- ploidymatrix[,1] + ploidymatrix[,5]
ploidymatrix[,2] <- ploidymatrix[,2] + ploidymatrix[,6]
ploidymatrix <- ploidymatrix[,c(1,2,3,4)]
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
  c(rep(0,10)),
  nrow = 2,
  dimnames = list(c("Gain", "Loss"),
                  c("Hap","Dip", "Trip","Tet","Pent")))
p_values <- matrix(
  c(rep(0,10)),
  nrow = 2,
  dimnames = list(c("Gain", "Loss"),
                  c("Hap","Dip", "Trip","Tet","Pent")))
fdr <- p_values <- matrix(
  c(rep(0,10)),
  nrow = 2,
  dimnames = list(c("Gain", "Loss"),
                  c("Hap","Dip", "Trip","Tet","Pent")))

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
}

#p_df <- data.frame(
#  yAxis = pLabels,
#  pgain = ploidymatrix[,6], ploss = ploidymatrix[,7])

ploidydf <- data.frame(xAxis = seq(1,5,1),
                       ploidy = c(rownames(ploidymatrix),rownames(ploidymatrix),rownames(ploidymatrix),rownames(ploidymatrix)),
                       type = c(rep("Gain",5),rep("Loss",5),rep("Both",5),rep("Both",5)), 
                       prop = c(ploidymatrix[,6],-ploidymatrix[,7],ploidymatrix[,8],-ploidymatrix[,8]),
                       odds = c(odds_ratio[1,],odds_ratio[2,],odds_ratio[1,],odds_ratio[2,]),
                       p = c(p_values[1,],p_values[2,],p_values[1,],p_values[2,]))
ploidydf$fdr[c(9, 19 , 4, 14 , 3, 13,  2 ,12, 10 ,20,  6, 16,  5 ,15  ,8, 18 , 1, 11 , 7, 17)] <- mt.rawp2adjp(ploidydf$p, proc = "BH", alpha = 0.05)$adjp[,2]
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
ploidydf$lab[6:20] <- NA
ploidydf$sig[11:20] <- NA
ploidydf$labjust <- c(-0.2,-0.2,-0.4,-1,-0.2)
ploidydf$labjust[6:20] <- NA
ploidydf$sigjust <- c(0,0,-1.8,-8,0, 2,2,0,9,2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)

# create ploidy fisher plot
fisherplot<- ggplot(ploidydf, aes(x= xAxis, y = prop)) + 
  geom_rect(aes(ymin= (-sum(ploidymatrix[,2]) - sum(ploidymatrix[,3]))/sum(ploidymatrix[,5]), ymax=(sum(ploidymatrix[,1]) + sum(ploidymatrix[,3]))/sum(ploidymatrix[,5]),xmin=-Inf,xmax=Inf),
            color="#e0e1e2", 
            alpha=0.03)+
  geom_col(fill=ploidydf$colour, width = 0.25, colour = "gray") + 
  geom_text(aes(label = sig, vjust = sigjust),size=3.5)+
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
