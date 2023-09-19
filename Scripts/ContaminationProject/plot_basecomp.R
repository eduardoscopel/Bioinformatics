library(data.table)
library(ggplot2)
df <- fread("bc_summary(transpose).txt",header=TRUE)
hap_hap <- df[1:6,2:32]
dip_hap <- df[7:12,2:32]
tri_hap <- as.matrix(df[13:18,2:32])
tet_hap <- as.matrix(df[19:24,2:32])

attach(df)
ACGT <- ggplot(df,aes(x=Contamination, y=TotalACGT, group=Mix))
ACGT <- ACGT + 
  geom_point(size = 2.0, aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_y_continuous(labels = comma) +
  ylab("High-quality ACGT calls")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
ACGT

DEG <- ggplot(df,aes(x=Contamination, y=TotalDEG, group=Mix))
DEG <- DEG +
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(0,110000,by=10000), labels = comma) +
  ylab("High-quality ambiguous calls")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
DEG

acgt <- ggplot(df,aes(x=Contamination, y=TotalLacgt, group=Mix))
acgt <- acgt + 
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 600000, by=100000), limits = c(0, 600000), labels = comma) +
  ylab("Low-quality acgt calls")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
acgt

deg <- ggplot(df,aes(x=Contamination, y=Totaldeg, group=Mix))
deg <- deg + 
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(0,20000,by=2500), labels = comma) +
  ylab("Low-quality ambiguous calls (not N)")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
deg

n <- ggplot(df,aes(x=Contamination, y=n, group=Mix))
n <- n + 
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix))+
  theme_bw() +
  scale_y_continuous(breaks=seq(0,250000,by=20000),labels = comma) +
  ylab("Low-quality N calls (any base)")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
n

lqdeg <- ggplot(df,aes(x=Contamination, y=`LQ-deg`, group=Mix))
lqdeg <- lqdeg + 
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix))+
  theme_bw() +
  scale_y_continuous(breaks=seq(0,250000,by=10000),labels = comma) +
  ylab("Low-quality ambiguous calls (total)")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
lqdeg

ACGT_het <- ggplot(df,aes(x=het, y=TotalACGT, group=Mix))
ACGT_het <- ACGT_het + 
  geom_point(size = 2.0, aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_y_continuous(labels = comma) +
  ylab("High-quality ACGT calls")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
ACGT_het

DEG_het <- ggplot(df,aes(x=het, y=TotalDEG, group=Mix))
DEG_het <- DEG_het +
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(0,110000,by=10000), labels = comma) +
  ylab("High-quality ambiguous calls")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
DEG_het

acgt_het <- ggplot(df,aes(x=het, y=TotalLacgt, group=Mix))
acgt_het <- acgt_het + 
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 600000, by=100000), limits = c(0, 600000), labels = comma) +
  ylab("Low-quality acgt calls")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
acgt_het

deg_het <- ggplot(df,aes(x=het, y=Totaldeg, group=Mix))
deg_het <- deg_het + 
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(0,20000,by=2500), labels = comma) +
  ylab("Low-quality ambiguous calls (not N)")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
deg_het

n_het <- ggplot(df,aes(x=het, y=n, group=Mix))
n_het <- n_het + 
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix))+
  theme_bw() +
  scale_y_continuous(breaks=seq(0,250000,by=20000),labels = comma) +
  ylab("Low-quality N calls (any base)")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
n_het

lqdeg_het <- ggplot(df,aes(x=het, y=`LQ-deg`, group=Mix))
lqdeg_het <- lqdeg_het + 
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix))+
  theme_bw() +
  scale_y_continuous(breaks=seq(0,250000,by=10000),labels = comma) +
  ylab("Low-quality ambiguous calls (total)")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
lqdeg_het


het_cont <- ggplot(df,aes(x=Contamination, y=het, group=Mix))
het_cont <- het_cont + 
  geom_point(size=2.0,aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix))+
  theme_bw() +
  scale_y_continuous(breaks=seq(0,0.01,by=0.001)) +
  ylab("High-quality heterozygosity")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
het_cont

attach(df)
SNPs <- ggplot(df,aes(x=Contamination, y=SNPs, group=Mix))
SNPs <- SNPs + 
  geom_point(size = 2.0, aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_y_continuous(labels = comma) + 
  ylab("SNPs")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
SNPs

pDiff <- ggplot(df,aes(x=Contamination, y=pDiff, group=Mix))
pDiff <- pDiff + 
  geom_point(size = 2.0, aes(color=Mix)) +
  geom_line(size=0.5, aes(color=Mix)) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,50,5))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  ylab("pDiffs")+
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
pDiff
