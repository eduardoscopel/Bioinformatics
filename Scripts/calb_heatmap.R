library(ggplot2)
library(data.table)

lohtable <- fread("/Users/es47540/Documents/GitHub/eduardo/Calbicans/genes_loh_stop.csv",header=TRUE)
lohtable
lohtable$stop_bi <- ifelse(lohtable$stop == 0, "No", "Yes")
lohtable$stop_tri <- ifelse(lohtable$stop == 0, 0, 
                            ifelse(lohtable$stop == 1,1,">1"))
lohtable <- lohtable[order(lohtable$Phenotype, decreasing = TRUE)]
lohtable$Strain <- factor(lohtable$Strain, levels = unique(lohtable$Strain))

loh_heatmap<- ggplot(lohtable, aes(Gene, Strain, fill=LOH)) + 
  geom_tile(color = c("black")) + 
  scale_fill_manual(values = c("#99FFFF","#FF8000"))+
  ylab("Strain")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1.2, size=4))+
  geom_tile(aes(x=-0.5, y=Strain, color = Phenotype), width = 0, size=2.5, height= 0.9)

loh_heatmap

png("Documents/GitHub/eduardo/Calbicans/LOH_heat_v2.png",width = 10, height = 4, units = "in",bg = "white", res=300)
loh_heatmap
dev.off()

stop_heatmap <- ggplot(lohtable, aes(Gene,Strain,fill=stop_tri)) + 
  geom_tile(color = c("black")) +
  scale_fill_manual(values = c("red","white","gray"))+
  ylab("Strain")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1.2, size=4))+
  geom_tile(aes(x=-0.5, y=Strain, color = Phenotype), width = 0, size=2.5, height= 0.9)
stop_heatmap
png("Documents/GitHub/eduardo/Calbicans/stop_heatmap_v2.png",width = 10, height = 4, units = "in",bg = "white", res=300)
stop_heatmap
dev.off()



library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
