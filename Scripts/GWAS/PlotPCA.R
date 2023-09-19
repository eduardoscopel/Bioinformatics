library(ggplot2)
library(gg3D)

### No outliers (>6 sd removed with smartpca)
data <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_noOutliers/PCAgain2n_gSNPs.evec",
                   col.names = c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","lineage"))
head(data)
write(levels(data$lineage),file = "~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/poplist.txt",ncolumns = 1)
eval <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_noOutliers/PCAgain2n_gSNPs.eval")

cladecolors_reduced <- c(ab11, af26, bb3, clinical, ecuador1, e21, ethiopia, fd5, fgh10, israel, mo4, mo8, "black", nao23, s25, darksalmon, wac12, we1)  


ggplot(data,aes(x=PC1, y=PC2, z=PC3, color=lineage, label = lineage)) + 
  #geom_text(aes(color=lineage, label = lineage)) + 
  theme_void() +
  axes_3D(theta=90) +
  stat_3D(theta=90,geom = "text") +
  scale_color_manual(values = cladecolors_reduced) 


for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(data) +
      geom_text(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors_reduced) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_noOutliers/",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(data) +
      geom_point(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors_reduced) +
      xlab(PCx)+
      ylab(PCy)+
      theme(panel.background = element_blank(), 
            legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_noOutliers/point",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

### No outliers (>10 sd removed with smartpca)
data10sd <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_10sd/PCAgain2n_gSNPs_out10sd.evec",
                   col.names = c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","lineage"))
head(data10sd)

cladecolors_10sd <- c(ab11, apw13, ab6, ol2, af26, ai24, bb3, clinical, ecuador1, e21, ethiopia, fd5, fgh10, 
                 israel, mo4, mo8, "black", nao23, s25, darksalmon, taiwan2, wac12, we1)


ggplot(data10sd,aes(x=PC1, y=PC2, z=PC3, color=lineage, label = lineage)) + 
  #geom_text(aes(color=lineage, label = lineage)) + 
  theme_void() +
  axes_3D(theta=90) +
  stat_3D(theta=90,geom = "text") +
  scale_color_manual(values = cladecolors_10sd) 


for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(data10sd) +
      geom_text(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors_10sd) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_10sd/",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(data) +
      geom_point(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors_reduced) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_noOutliers/point",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

### No outliers (>15 sd removed with smartpca)
data15sd <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_15sd/PCAgain2n_gSNPs_out15sd.evec",
                       col.names = c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","lineage"))
head(data15sd)
levels(data15sd$lineage)
cladecolors_15sd <- c(ab11, apw13, ab6, ol2, af26, ai24, bb3, ch14, clinical, ecuador1, e21, ethiopia, fea18, fer22, fd5, fgh10, 
                      israel, malaysian, mo4, ma9, mo8, "black", nao23, s25, siberia, darksalmon, taiwan2, wac12, we1)


ggplot(data15sd,aes(x=PC1, y=PC2, z=PC3, color=lineage, label = lineage)) + 
  #geom_text(aes(color=lineage, label = lineage)) + 
  theme_void() +
  axes_3D(theta=90) +
  stat_3D(theta=90,geom = "text") +
  scale_color_manual(values = cladecolors_15sd) 


for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(data15sd) +
      geom_text(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors_15sd) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_15sd/",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(data15sd) +
      geom_point(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors_15sd) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_15sd/point",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

ggplot(data15sd) +
  geom_point(aes(x = PC1, y = PC2, color=lineage, label = lineage)) + 
  scale_color_manual(values = cladecolors_15sd) +
  xlab(PCx)+
  ylab(PCy)+
  xlim(c(-0.03,0.06))+
  ylim(c(-0.01,0.01))+
  theme(legend.position = 'none')


### With outliers (not removed by smartpca)

dataWO <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_withOutliers/PCAgain2n_gSNPs_noOutRemov.evec",
                   col.names = c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","lineage"))
head(dataWO)
levels(dataWO$lineage)


cladecolors <- c(ab11, apw13, ab6, ol2, af26, ai24, bb3, ch14, clinical, ecuador1, e21, ethiopia, fea18, fer22, fd5, fgh10, 
                 israel, malaysian, mo4, ma9, mo8, "black", nao23, s25, siberia, darksalmon, t17, taiwan2, wac12, we1)

ggplot(dataWO,aes(x=PC1, y=PC2, z=PC3, color=lineage, label = lineage)) + 
  #geom_text(aes(color=lineage, label = lineage)) + 
  theme_void() +
  axes_3D(theta=90) +
  stat_3D(theta=90,geom = "text") +
  scale_color_manual(values = cladecolors) 


for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(dataWO) +
      geom_text(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_withOutliers/",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(dataWO) +
      geom_point(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_withOutliers/point",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}


### With outliers (pops projected)

dataProj <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_projected_withOutliers/PCAgain2n_gSNPs_projected.pca.evec",
                     col.names = c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","lineage"))
head(dataProj)
levels(dataProj$lineage)
evalProj <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_projected_withOutliers/PCAgain2n_gSNPs_projected.eval")
evalProj <- (evalProj/sum(evalProj))*100
head(evalProj,n = 10)

cladecolors <- c(ab11, apw13, ab6, ol2, af26, ai24, bb3, ch14, clinical, ecuador1, e21, ethiopia, fea18, fer22, fd5, fgh10, 
                 israel, malaysian, mo4, ma9, mo8, "black", nao23, s25, siberia, darksalmon, t17, taiwan2, wac12, we1)

ggplot(dataProj,aes(x=PC1, y=PC2, z=PC3, color=lineage, label = lineage)) + 
  #geom_text(aes(color=lineage, label = lineage)) + 
  theme_void() +
  axes_3D(theta=90) +
  stat_3D(theta=90,geom = "text") +
  scale_color_manual(values = cladecolors) 


for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(dataProj) +
      geom_text(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors) +
      xlab(paste(PCx, " (", round(evalProj[i,],digits = 2), "%)", sep = ""))+
      ylab(paste(PCy, " (", round(evalProj[j+1,],digits = 2), "%)", sep = ""))+
      theme_bw()+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_projected_withOutliers/",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(dataProj) +
      geom_point(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage), size = 4.32, shape = 1) + 
      scale_color_manual(values = cladecolors) +
      xlab(paste(PCx, " (", round(evalProj[i,],digits = 2), "%)", sep = ""))+
      ylab(paste(PCy, " (", round(evalProj[j+1,],digits = 2), "%)", sep = ""))+
      theme_bw()+
      theme(panel.grid = element_blank(), 
            legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_projected_withOutliers/point",PCx,PCy,"v2.png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}


AnovaTests <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_projected_withOutliers/AnovaPCA.txt")
tot <- 154
ev1 <- 139
ev2 <- 128
ev3 <- 86
ev4 <- 136
ev5 <- 131
ev6 <- 119
ev7 <- 116
ev8 <- 107
ev9 <- 102
ev10 <- 133

PC1xPC2 <- ggplot(dataProj) +
  geom_point(aes(x = PC1, y = PC2, color=lineage, label = lineage), size = 0.5) + 
  scale_color_manual(values = cladecolors) +
  xlab(paste("PC1", " (", round(evalProj[1,],digits = 2), "%)", sep = ""))+
  ylab(paste("PC2", " (", round(evalProj[2,],digits = 2), "%)", sep = ""))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8),
        legend.position = 'none')

PC1xPC10 <- ggplot(dataProj) +
  geom_point(aes(x = PC1, y = PC10, color=lineage, label = lineage), size = 0.5) + 
  scale_color_manual(values = cladecolors) +
  xlab(paste("PC1", " (", round(evalProj[1,],digits = 2), "%)", sep = ""))+
  ylab(paste("PC10", " (", round(evalProj[10,],digits = 2), "%)", sep = ""))+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8),
        legend.position = 'none')

bottomRow <- plot_grid(PC1xPC2, PC1xPC10, labels = c("(b)","(c)"), label_size = 12, nrow = 1)

ggsave("../../../../Papers_Scopel/GWAS/diploids/gainOnly/PCA.png",
       plot = bottomRow,
       device = "png", 
       dpi = 600, 
       limitsize = FALSE,
       width = 15,
       height = 7,
       units = "in")

### No outliers (QC3, >6 sd removed with smartpca)
dataQC3_6sd <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/QC3_noOutliers/PCAgain2n_gSNPs_QC3.evec",
                   col.names = c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","lineage"))
head(dataQC3_6sd)
levels(dataQC3_6sd$lineage)
cladecolors_QC3_6sd <- c(ab11, ab6, af26, ai24, bb3, ch14, clinical, ecuador1, ethiopia, fea18, fer22, fgh10, 
                         israel, malaysian, mo4, ma9, mo8, "black", nao23, s25, siberia, darksalmon, taiwan2, wac12, we1)

ggplot(dataQC3_6sd,aes(x=PC1, y=PC2, z=PC3, color=lineage, label = lineage)) + 
  #geom_text(aes(color=lineage, label = lineage)) + 
  theme_void() +
  axes_3D(theta=90) +
  stat_3D(theta=90,geom = "text") +
  scale_color_manual(values = cladecolors_QC3_6sd) 


for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(dataQC3_6sd) +
      geom_text(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors_QC3_6sd) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/QC3_noOutliers/",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(data) +
      geom_point(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors_reduced) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_noOutliers/point",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

### With outliers (QC3, no outliers removed with smartpca)
dataQC3 <- read.table("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/QC3_withOutliers/PCAgain2n_gSNPs_QC3_noOut.evec",
                          col.names = c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","lineage"))
head(dataQC3)
levels(dataQC3$lineage)

ggplot(dataQC3,aes(x=PC1, y=PC2, z=PC3, color=lineage, label = lineage)) + 
  #geom_text(aes(color=lineage, label = lineage)) + 
  theme_void() +
  axes_3D(theta=90) +
  stat_3D(theta=90,geom = "text") +
  scale_color_manual(values = cladecolors) 


for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(dataQC3) +
      geom_text(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/QC3_withOutliers/",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}

for(i in seq(1,9,1)){
  PCx <- paste("PC" , i, sep = "")
  for(j in seq(i,9,1)){
    PCy<- paste("PC", j + 1, sep = "")
    PCplot <- ggplot(data) +
      geom_point(aes(x = eval(parse(text = PCx)), y = eval(parse(text = PCy)), color=lineage, label = lineage)) + 
      scale_color_manual(values = cladecolors_reduced) +
      xlab(PCx)+
      ylab(PCy)+
      theme(legend.position = 'none')
    ggsave(paste("~/Documents/Papers_Scopel/GWAS/diploids/gainOnly/newRunQC40Miss10MAF1.0noRecomb/PCA/all_noOutliers/point",PCx,PCy,".png",sep = ""),
           plot = PCplot,
           device = "png", 
           dpi = 600, 
           limitsize = FALSE,
           width = 6,
           height = 6,
           units = "in")
  }
}
