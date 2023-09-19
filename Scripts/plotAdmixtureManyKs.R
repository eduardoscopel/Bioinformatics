library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

col_list<-c("#310041",
            "#6f0286",
            "#ff2b5e",
            "#ffa804",
            "#fcff92",
            "#22ff5b",
            "#5f86d5",
            "#523b73",
            "#183553",
            "#22578d",
            "#b0d7ff",
            "#22cccc",
            "#feda68",
            "#e94963",
            "#a83f5d",
            "#532a48",
            "#351232",
            "#772d46",
            "#cc4452",
            "#e4b79a",
            "#f7e4ac",
            "#dbd7a8",
            "#fca19e",
            "#ff3c7e",
            "#5a2949",
            "#18183a",
            "#9f1793",
            "#f43c7c",
            "#ffb39b",
            "#ffe0a9",
            "#eef663",
            "#9dec5f" 
)

col_list<-c("#9dec5f","#7570B3","#CAB2D6","#FFF2AE","#E5C494","#6A3D9A","#1F78B4","#BF5B17","#FDC086","#E7298A","#BC80BD","#33A02C","#E78AC3",
            "#F2F2F2","#A6761D","#D95F02","#B2DF8A","#984EA3","red","#8DA0CB","#B3E2CD","#FB9A99","gray","orange","forestgreen","steelblue2","darkcyan","goldenrod3","#310041",
            "#ffa804","#ff2b5e")


ids<-read.csv("strain_pop.csv",header=T)
dup <- rownames(ids[duplicated(ids$strain) == TRUE,])
ids[ids$strain == "ETPN8",]
ids[ids$strain == "SJ1L04",]
ids[ids$strain == "N163_01_5A",]
ids <- ids[!duplicated(ids$strain) == TRUE,]
ancestry <- data.frame(pop =ids$ancestry,strain = ids$strain)

Qlist <- list.files("/Users/es47540/Documents/OtherPeople/Jacque/AdmixtureData/Admixture/Q/ready/")
plotList <- c()
counter <- 1
popIndices <- as.numeric(tapply(seq_along(sort(ancestry$pop)),sort(ancestry$pop), max))+0.5
for(i in Qlist[1:20]){
  temp <- read.table(paste("/Users/es47540/Documents/OtherPeople/Jacque/AdmixtureData/Admixture/Q/ready/",i,sep = ""))
  temp <- temp[!row.names(temp) %in% dup,]
  k <- length(temp)
  #for(j in 1:nrow(temp)){
  #  if(rowSums(temp[j,1:k] > 0.5) == 0){
  #    temp$Panc[j] <- "admix"
  #  }
  #  else{
  #    temp$Panc[j] <- which.max(temp[j,1:k])
  #  }
  #  temp$maxV[j] <- max(temp[j,1:k])
  #}
  temp <- cbind(temp,ancestry)  
  rownames(temp) <- temp$strain
  temp <- temp[order(temp$pop),] # order by lineage
  #temp <- temp[order(temp$pop, temp$Panc,temp$maxV),] # order by % ancestry
  temp <- temp[,1:k]
  keys <- rownames(temp)
  keyDF <- data.frame(key = keys, weigth = 1:length(keys))
  melteddf <- melt(t(temp))
  head(melteddf)
  if(counter %% 4 == 0){
    plotList[[counter]] <- ggplot(melteddf) + 
      geom_bar(aes(x=Var2,y=value,color=Var1,fill=Var1),stat="identity", show.legend = FALSE) +
      geom_vline(xintercept = popIndices, size=.8) +
      scale_color_manual(values = col_list)+
      scale_fill_manual(values = col_list) +
      scale_x_discrete(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      ggtitle(paste("k = ",k))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),axis.title = element_blank(), plot.margin = unit(c(0.2,0.1,0,0.1),"cm"), 
            plot.title = element_text(size = 8,hjust = 0.5,margin = unit(c(0,0,0,0),"cm")))
    #theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())  
  }
  else{
    plotList[[counter]] <- ggplot(melteddf) + 
      geom_bar(aes(x=Var2,y=value,color=Var1,fill=Var1),stat="identity", show.legend = FALSE) +
      geom_vline(xintercept = popIndices, size=.8)+
      scale_color_manual(values = col_list)+
      scale_fill_manual(values = col_list) +
      scale_x_discrete(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      ggtitle(paste("k = ",k, sep = ""))+
      #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4))
      theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(0.2,0.1,0.1,0.1),"cm"),
            plot.title = element_text(size = 8,hjust = 0.5,margin = unit(c(0,0,0,0),"cm")))  
  }
  counter <- counter + 1 
}


plot_grid(plotlist = plotList[1:4], nrow = 4, rel_heights = c(1,1,1,1.5))
plot_grid(plotlist = plotList[5:8], nrow = 4, rel_heights = c(1,1,1,1.5))
plot_grid(plotlist = plotList[9:12], nrow = 4, rel_heights = c(1,1,1,1.5))
plot_grid(plotlist = plotList[13:16], nrow = 4, rel_heights = c(1,1,1,1.5))
plot_grid(plotlist = plotList[17:20], nrow = 4, rel_heights = c(1,1,1,1.5))

ggsave(paste("/Users/es47540/Documents/OtherPeople/Jacque/AdmixtureData/Admixture/Q/ready/",k,".png",sep = ""), width=15, height = 15, units="in",limitsize = FALSE, dpi = 600)