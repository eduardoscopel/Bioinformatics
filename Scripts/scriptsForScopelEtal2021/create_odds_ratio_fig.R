library(ggplot2)

boxLabels = c("Beer","Bioethanol","Bread","Cocoa/Coffee","Commercial","Dairy",
              "Fermentation","Industrial","Lab","Sake","SSF","Wine",
              "Clinical",
              "Flower","Fruit","Insect","Other plants","Tree","Water")
my_colors = c("Domesticated","Domesticated","Domesticated","Domesticated","Domesticated","Domesticated",
              "Domesticated","Domesticated","Domesticated","Domesticated","Domesticated","Domesticated",
              "Clinical", 
              "Wild","Wild","Wild","Wild","Wild","Wild")
a = c("#225EA8","#225EA8","#225EA8","#225EA8","#225EA8","#225EA8",
              "#225EA8","#225EA8","#225EA8","#225EA8","#225EA8","#225EA8",
              "#FC724F", 
              "#316857","#316857","#316857","#316857","#316857","#316857")

df <- data.frame(
  yAxis = length(boxLabels):1,
  boxOdds = (c(2.029113,0.321546,1.660344,0.3295849,1.976214,1.272808,
            0.9201797,1.317703,0.7163538,1.90005,0.6759632,0.6452242,
            1.076583,
            0.5148992,0.6436411,0.7163538,0.8237529,0.3048032,1.646925)),
  boxCILow = c(1.455527,0.08326178,1.183054,0.0372251,1.347894,0.6754152,
               0.5893619,0.3003143,0.2116077,1.156584,0.2710123,0.4672979,
               0.8132553,
               0.2339044,0.3968007,0.2116077,0.1485614,0.1687242,0.647434
               ),
  boxCIHigh = c(2.813401,0.89435568,2.311772,1.3645965,2.874027,2.2974434,
                1.4010503,4.5933194,1.9400236,3.072945,1.4799832,0.8786156,
                1.4147808,
                1.0194850,1.0075243,1.9400236,3.0693837,0.5155689,3.889906
                )
)

orplot <- ggplot(df, aes(x= boxOdds, y = yAxis))
orplot <- orplot + 
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 2, height = 0.5,color = "gray50") +
  geom_point(aes(color = factor(my_factors)), size = 5.0) +
  scale_color_manual(values=c("#FC724F","#225EA8","#316857")) +
  theme_bw() +
  scale_y_continuous(breaks = length(boxLabels):1, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,5,1)) +
  coord_trans(x = "sqrt") +
  ylab("") +
  xlab("Odds ratio (sqrt scale)\np-value < 0.01") +
  theme(axis.text= element_text(face="bold",color="gray50", size=12))+
  theme(axis.text.y= element_text(color=my_colors))+
  ggtitle("Association between ecology and aneuploidy in yeast")
orplot <- orplot +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
    )
orplot <- orplot + theme(legend.position = "none")
tiff("odds_ratio.tiff", width = 6, height = 8, units = 'in', res = 300, bg='transparent')
orplot
dev.off()
