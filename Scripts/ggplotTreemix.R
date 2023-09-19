library(ggplot2)
library(data.table)
library(cowplot)
library(ape)


setwd("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/k500/s1/")
edges6="new_k500_s1_621-merged-snps.noN.6"
edges7="new_k500_s1_621-merged-snps.noN.7"
edges8="new_k500_s1_621-merged-snps.noN.8"
source("/Users/es47540/Documents/GitHub/eduardo/Scripts/treemix.R")

cladeColors6 <- c(apw13, m19, ol2, nao23, ai24, fgh10, ch15, fea18, wac12, bb3, we1, ma9, af26, e21, t17,ch14, ch16, ch20,fer22, mo4)
aneupCladeColors6 <- c(mo8, s25, ab11,ab6,mb7,fd5)

e6 <- plot.treemix(edges6,plot.nodes = FALSE, branch.colour = "black", clade_color = cladeColors6, aneup_clade_color = aneupCladeColors6)

cladeColors7 <- c(apw13, m19, ol2, nao23, ai24, fgh10, ch15, fea18, wac12,bb3, we1, ma9, af26, e21,t17, ch14, ch16, ch20, fer22, mo4)
aneupCladeColors7 <- c(mo8, s25, ab6, ab11, mb7, fd5)
e7 <- plot.treemix(edges7,plot.nodes = FALSE, branch.colour = "black", clade_color = cladeColors7, aneup_clade_color = aneupCladeColors7)

cladeColors8 <- c(apw13, m19, ol2, nao23, ai24, fgh10, ch15, fea18, wac12,bb3, we1, ma9, af26, e21,t17, ch14, ch16, ch20, fer22, mo4)
aneupCladeColors8 <- c(mo8, s25, ab6, ab11, mb7, fd5)
e8 <- plot.treemix(edges8,plot.nodes = FALSE, branch.colour = "black", clade_color = cladeColors8, aneup_clade_color = aneupCladeColors8, plot.title = "8 migration evens")

lliketab <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/k500/s1/k500_llike.txt")
colnames(lliketab) <- c("M", "llike")
lliketab <- lliketab[lliketab$M < 11,]
llikeplot <- ggplot(lliketab) + 
  geom_point(aes(x=M, y= llike)) +
  geom_line(aes(x=M, y= llike)) +
  scale_x_continuous(breaks = c(0:11),expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(-10000,2000, by = 2000),limits = c(-10000,2000),expand = c(0.01,0.01)) +
  xlab("Migration events")+
  ylab("Log likelihood") + 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.border = element_blank(), axis.line.x = element_line(size=0.5),axis.line.y = element_line(size=0.5))
options(scipen=999)
llikeplot

cladeColors8 <- c(apw13, m19, ol2, nao23, ai24, fgh10, ch15, fea18, wac12,bb3, we1, ma9, af26, e21,t17, ch14, ch16, ch20, fer22, mo4)
aneupCladeColors8 <- c(mo8, s25, ab6, ab11, mb7, fd5)
filtered_e8 <- plot.treemix(edges8,plot.nodes = FALSE, branch.colour = "black", clade_color = cladeColors8, aneup_clade_color = aneupCladeColors8, subset.mig = c(12,25,48,66),
                            plot.title = "Migration events confirmed by f3 tests")
filtered_e8_1 <- plot.treemix(edges8,plot.nodes = FALSE, branch.colour = "black", plot.migration = FALSE,clade_color = cladeColors8, aneup_clade_color = aneupCladeColors8, subset.mig = c(12,25,48,66),
                            plot.title = "Migration events confirmed by f3 tests")



toprow <- plot_grid(llikeplot, e8, nrow = 1, ncol = 2, rel_widths = c(2,3))
toprow
suppFigS4_4 <- plot_grid(toprow, filtered_e8, nrow=2,ncol=1, rel_heights = c(2,3))
suppFigS4_4
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v5.png", width=12, height = 12, units="in",limitsize = FALSE, dpi = 600)
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v5.pdf", width=10, height = 10, units="in",limitsize = FALSE, dpi = 600)
ggsave("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/plots/peter/trees/FigS4.v5.eps",device=cairo_ps, width=12, height = 12, units="in",limitsize = FALSE, dpi = 600)


f3pop <- fread("/Users/es47540/Documents/GitHub/eduardo/aneuploidy/peter/admixture/treemix/k500/k500-f3pop-test.txt")
f3pop[f3pop$f3stat < 0 & f3pop$`Z-score` < -2, ][order(`3pop`)]
f3pop[startsWith(f3pop$`3pop`, "Ale_beer"),][order(`3pop`)]
