library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)    # Data manipulation
value = c(85,23,112,41,111,10,52,123,44,39,286,20,254,12,64
          ,128,23,173)
df <- data.frame(
  Category = c("Commercial", "Lab", "Beer", "Bioethanol", 
               "Liquid Fermentation", "Industrial","Sake",
               "Bread","Dairy","Solid Fermentation", "Wine", 
               "Coffee/Cocoa", "Clinical", "Other Plants", 
               "Flower", "Fruit", "Insect", "Trees"),
  value = c(85,23,112,41,111,10,52,123,44,39,286,20,254,12,64
            ,128,23,173),
  prop = round(100*(value/1621),digits=1)
)
df <- df %>%
  arrange(desc(Category)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
df
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(df$Category))
bp <- ggplot(df, aes(x="",y=prop,fill=Category))+
  geom_bar(width = 1,stat = "identity",color="white")+
  coord_polar("y",start=0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_void()
bp

par(bg=NA)
value2 = c(421,254,946)
df2 <- data.frame(
  Category = c("Wild","Clinical","Domesticated"),
  value = value2,
  prop = round(100*(value2/1621),digits=1)
)
df2 <- df2 %>%
  arrange(desc(Category)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
df2
mycolors <- c("#FC724F",
              "#225EA8",
              "#316857")
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(df$Category))
bp2 <- ggplot(df2, aes(x="",y=prop,fill=Category))+
  geom_bar(width = 1,stat = "identity",color="white")+
  coord_polar("y",start=0)+
  geom_text(aes(y = lab.ypos, label = value), color = "#fefff0",size=4)+
  scale_fill_manual(values=mycolors) +
  theme_void()+theme(legend.position="bottom",
                     legend.text=element_text(size=6), 
                     legend.title = element_blank(),
                     legend.margin = margin(0,0,0,0),
                     legend.box.margin = margin(-18,0,18,0))
bp2 <- bp2 + theme(
  panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA)
)
tiff("pie_chart.tiff",width = 3, height = 2, units = 'in', res = 300, bg = 'transparent')
bp2
dev.off()




