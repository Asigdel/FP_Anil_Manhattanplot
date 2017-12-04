# FP_Anil_Manhattanplot
```

# The way to merge 6 plots 
# Remember, make 6 plots and just get the snp_sol_inter_1 for 29 chromosomes and 
# the final way is to combine all into one 
# Bingo 

rm(list=ls())
setwd("D:/Random_Regression_Testday_Model/Reproduction_Traits/")
library(biomaRt)
database = useMart("ensembl")
genome = useDataset("btaurus_gene_ensembl", mart = database)
gene = getBM(c("ensembl_gene_id", "entrezgene", "external_gene_name",
               "start_position", "end_position", "chromosome_name"), mart = genome); dim(gene)
# CREATE SNP FILTER - upload file freqdata.count.after.clean 

filter = read.table("freqdata.count.after.clean", header = F)
colnames(filter) = c("SNP", "Frequency", "Filter")
f = filter$Filter == 0; table(f)

# Manhattan plot 

require(ggplot2)
require(gridExtra)
require(grid)
require(cowplot)

SOL=read.table("snp_sol_inter_1")
colnames(SOL)<-c("trait", "effect", "SNP", "Chromosome", "Position", "SNP solution", "weight", "Variance")
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$within = seq(1:nrow(SNP))

gap = 0; xx = as.numeric(); x0 = 0

for(w in 1:29){
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)

coloR <- c(rep(c("blue2","darkorange2"),14),"blue2")
plot1 <- ggplot(SNP, aes(x = within, y = Var, col = factor(Chr))) +
  scale_color_manual(values = coloR) + 
  geom_point( ) + 
  ylim(0,2.0) +
  ylab("Genetic Variance (%)")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x  = element_text(size = 12)) + 
  theme(axis.text.y  = element_text(size = 12)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=12))


SOL=read.table("snp_sol_slope_1")
colnames(SOL)<-c("trait", "effect", "SNP", "Chromosome", "Position", "SNP solution", "weight", "Variance")
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$within = seq(1:nrow(SNP))

gap = 0; xx = as.numeric(); x0 = 0

for(w in 1:29){
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)

coloR <- c(rep(c("blue2","darkorange2"),14),"blue2")
plot2 <- ggplot(SNP, aes(x = within, y = Var, col = factor(Chr))) +
  scale_color_manual(values = coloR) + 
  geom_point( ) + 
  ylim(0,1.20) +
  ylab("Genetic Variance (%)")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x  = element_text(size = 12)) + 
  theme(axis.text.y  = element_text(size = 12)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=12))


SOL=read.table("snp_sol_inter_2")
colnames(SOL)<-c("trait", "effect", "SNP", "Chromosome", "Position", "SNP solution", "weight", "Variance")
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$within = seq(1:nrow(SNP))

gap = 0; xx = as.numeric(); x0 = 0

for(w in 1:29){
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)

coloR <- c(rep(c("blue2","darkorange2"),14),"blue2")
plot3 <- ggplot(SNP, aes(x = within, y = Var, col = factor(Chr))) +
  scale_color_manual(values = coloR) + 
  geom_point( ) + 
  ylim(0,1.80) +
  ylab("Genetic Variance (%)")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x  = element_text(size = 12)) + 
  theme(axis.text.y  = element_text(size = 12)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=12))


SOL=read.table("snp_sol_slope_2")
colnames(SOL)<-c("trait", "effect", "SNP", "Chromosome", "Position", "SNP solution", "weight", "Variance")
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$within = seq(1:nrow(SNP))

gap = 0; xx = as.numeric(); x0 = 0

for(w in 1:29){
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)

coloR <- c(rep(c("blue2","darkorange2"),14),"blue2")
plot4 <- ggplot(SNP, aes(x = within, y = Var, col = factor(Chr))) +
  scale_color_manual(values = coloR) + 
  geom_point( ) + 
  ylim(0,1.60) +
  ylab("Genetic Variance (%)")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x  = element_text(size = 12)) + 
  theme(axis.text.y  = element_text(size = 12)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=12))


SOL=read.table("snp_sol_inter_3")
colnames(SOL)<-c("trait", "effect", "SNP", "Chromosome", "Position", "SNP solution", "weight", "Variance")
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$within = seq(1:nrow(SNP))

gap = 0; xx = as.numeric(); x0 = 0

for(w in 1:29){
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)

coloR <- c(rep(c("blue2","darkorange2"),14),"blue2")
plot5 <- ggplot(SNP, aes(x = within, y = Var, col = factor(Chr))) +
  scale_color_manual(values = coloR) + 
  geom_point( ) + 
  ylim(0,1.80) +
  ylab("Genetic Variance (%)")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x  = element_text(size = 12)) + 
  theme(axis.text.y  = element_text(size = 12)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=12))

SOL=read.table("snp_sol_slope_3")
colnames(SOL)<-c("trait", "effect", "SNP", "Chromosome", "Position", "SNP solution", "weight", "Variance")
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$within = seq(1:nrow(SNP))

gap = 0; xx = as.numeric(); x0 = 0

for(w in 1:29){
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)

coloR <- c(rep(c("blue2","darkorange2"),14),"blue2")
plot6 <- ggplot(SNP, aes(x = within, y = Var, col = factor(Chr))) +
  scale_color_manual(values = coloR) + 
  geom_point( ) + 
  ylim(0,1.40) +
  ylab("Genetic Variance (%)")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x  = element_text(size = 12)) + 
  theme(axis.text.y  = element_text(size = 12)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=12))


tiff("Figure 1.tiff", width = 14, height = 12, units = 'in', res = 300)
plot_grid(plot1, plot2,plot3,plot4, plot5,plot6, align = c("hv"), nrow = 3,  
          labels = c("A", "B","C","D","E","F"), label_size= 20, label_colour = "darkgreen")
dev.off()
```
