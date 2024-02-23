
#GeoMx R
library(ggplot2)
library(plotly)
library(tidyverse)
library(readxl)
library(ggrepel)
library(RColorBrewer)
library(ggpubr)

#PCA-plot
data <- read.csv("PCA_data_R.csv")
PC1 <- data$PC1
PC2 <- data$PC2
PC3 <- data$PC3
Condition <- data$Condtion

axx <- list(title = "PC1(0.352)") #overskrifter til figur
axy <- list(title = "PC2(0.281)")
axz <- list(title = "PC3(0.146)")

figPCA <- plot_ly(data, x=PC1, y=PC2, z=PC3, 
                  color=Condition,
                  marker = list(size = 5))
figPCA <- figPCA %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
figPCA


#export
library(kaleido)
library(reticulate)

kaleido(figPCA, "PCAfigureX.pdf")
save.image(file=figPCA, "PCAfigureX.pdf")


p <- plot_ly(x = 1:10)
plot_ly(fig, type="scatter3d"[, ...])
add_trace(p, type="scatter3d"[, ...])
reticulate::py_run_string("import sys")
save_image(p, "./pic.png")
save_image(figPCA, "./PCA-plot.pdf")





#########################################

#volcano plots 

#exocrine
exo <- read_excel("Volcano plot, exocrine.xlsx")
ggplot(data=exo, aes(x=`Fold changes`, y=`-log10 adjusted pvalue`)) + 
  geom_point(size = 2) +
  ggtitle("Volcano plot, exocrine") +
  xlim(-6,2) +
  ylim(-0.05, 4.2) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0) + 
  theme_light(base_size = 14) +
  geom_text_repel(data=subset(exo, `-log10 adjusted pvalue` > 1.2),
            aes(`Fold changes`,`-log10 adjusted pvalue`,label=`Target name`),
            nudge_y = 0.2)


#islet distal vs proximal
islet_dp <- read_excel("Volcano plot, islet d vs p.xlsx")
ggplot(data=islet_dp, aes(x=`Fold changes`, y=`-log10 adjusted pvalue`)) +
  geom_point(size = 2, position=position_jitter(w=0.05)) + #jitter added to separate S6 and GZMB to improve clarity
  ggtitle("Volcano plot, islet distal vs proximal") +
  xlim(-2,1.5) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0) + 
  theme_light(base_size = 14) +
  geom_text_repel(data=subset(islet_dp, `-log10 adjusted pvalue` > 1.2),
                  aes(`Fold changes`,`-log10 adjusted pvalue`,label=`Target name`),
                  nudge_y = 0.2)



#islet ctr vs proximal
islet_cp <- read_excel("Volcano plot, islet c vs p.xlsx")
ggplot(data=islet_cp, aes(x=`Fold changes`, y=`-log10 adjusted pvalue`)) +
  geom_point(size = 2) +
  ggtitle("Volcano plot, islet ctr vs proximal") +
  xlim(-2.8,1.5) +
  ylim(-0.05, 2.38) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0) + 
  theme_light(base_size = 14) +
  geom_text_repel(data=subset(islet_cp, `-log10 adjusted pvalue` > 1.3),
                  aes(`Fold changes`,`-log10 adjusted pvalue`,label=`Target name`),
                  nudge_y = 0.2)


#islet ctr vs distal, no differences found after statistical testing
islet_cd <- read_excel("Volcano plot, islet c vs d.xlsx")
ggplot(data=islet_cd, aes(x=`Fold changes`, y=`-log10 adjusted pvalue`)) +
  geom_point(size = 2) +
  ggtitle("Volcano plot, islet ctr vs distal") +
  ylim(-0.05, 0.48) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0) + 
  theme_light(base_size = 14) +
  geom_text_repel(data=subset(islet_cd, `-log10 adjusted pvalue` > 1.3),
                  aes(`Fold changes`,`-log10 adjusted pvalue`,label=`Target name`),
                  nudge_y = 0.2)



########################################################################3

#boxplots

#summary of signif differences
#exocrine: CD45, SMA, CD4, HistonH3 higher in proximal. PD-1 borderline signif, higher in distal
#islet p vs d: GZMB and S6 higher in distal. CD45 and CD8a higher in proximal
#islet c vs p: CD45, CD4, CD11c higher in proximal
#islet c vs d: no differences


#exocrine boxplots
#adj.pvals: CD4 0.00697, CD45 0.00009, HistonH3 0.01248, PD-1 0.05125, SMA 0.02024

#CD4 exocrine
bexCD4 <- read_excel("ex-CD4.xlsx")
bexCD4p <- tibble::tribble(~group1, ~group2, ~p.adj, "distal", "proximal", 0.00697)

tiff("ex-CD4.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(bexCD4, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(bexCD4p, y.position = 5.8, step.increase = 0.1) +
  ylim(0,6) +
  theme_grey() +
  theme(legend.position="none") +
  ggtitle("Exocrine CD4")
dev.off()


#CD45 exocrine
bexCD45 <- read_excel("ex-CD45.xlsx")
bexCD45p <- tibble::tribble(~group1, ~group2, ~p.adj, "distal", "proximal", 0.00009)

tiff("ex-CD45.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(bexCD45, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(bexCD45p, y.position = 32, step.increase = 0.1) +
  ylim(0,33) +
  theme_grey() +
  theme(legend.position="none") +
  ggtitle("Exocrine CD45")
dev.off()


#HistonH3 exocrine (housekeeping)
bexH3 <- read_excel("ex-HistoneH3.xlsx")
bexH3p <- tibble::tribble(~group1, ~group2, ~p.adj, "distal", "proximal", 0.01248)

tiff("ex-H3.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(bexH3, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(bexH3p, y.position = 290, step.increase = 0.1) +
  ylim(40,300) +
  theme_grey() +
  theme(legend.position="none") +
  ggtitle("Exocrine Histone H3")
dev.off()


#PD-1 exocrine (borderline signif)
bexPD1 <- read_excel("ex-PD-1 border signif.xlsx")
bexPD1p <- tibble::tribble(~group1, ~group2, ~p.adj, "distal", "proximal", 0.05125)

tiff("ex-PD-1.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(bexPD1, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(bexPD1p, y.position = 6, step.increase = 0.1) +
  ylim(0,6.3) +
  theme_grey() +
  theme(legend.position="none") +
  ggtitle("Exocrine PD-1")
dev.off()


#SMA exocrine 
bexSMA <- read_excel("ex-SMA.xlsx")
bexSMAp <- tibble::tribble(~group1, ~group2, ~p.adj, "distal", "proximal", 0.02024)

tiff("ex-SMA.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(bexSMA, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(bexSMAp, y.position = 135, step.increase = 0.1) +
  ylim(0,140) +
  theme_grey() +
  theme(legend.position="none") +
  ggtitle("Exocrine SMA")
dev.off()



#islet boxplots


#CD4
#adjp ctr vs distal 0.5456, ctr vs proximal 0.00685, distal vs proximal 0.28367
biCD4 <- read_excel("is-CD4.xlsx")
biCD4p <- tibble::tribble(~group1, ~group2, ~p.adj, 
                          "ctr", "distal", 0.00446,
                          "ctr", "proximal", 0.0001,
                          "distal", "proximal", 0.039566)

tiff("is-CD4.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(biCD4, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(biCD4p, y.position = 2.5, step.increase = 0.14, label = "p.adj", size = 2.7) +
  ylim(0,3.25) +
  theme_grey() +
  theme(legend.position="none") + 
  ggtitle("Islet CD4")
dev.off()


#CD8a
#adjp ctr vs distal, 0.73671, ctr vs proximal 0.19685, distal vs proximal 0.0377
biCD8a <- read_excel("is-CD8a.xlsx")
biCD8ap <- tibble::tribble(~group1, ~group2, ~p.adj, 
                          "ctr", "distal", 0.73671,
                          "ctr", "proximal", 0.19685,
                          "distal", "proximal", 0.0377)

tiff("is-CD8a.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(biCD8a, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(biCD8ap, y.position = 3.5, step.increase = 0.14, label = "p.adj", size = 2.7) +
  ylim(0.5,4.5) +
  theme_grey() +
  theme(legend.position="none") + 
  ggtitle("Islet CD8a")
dev.off()


#CD11c
#adjp pval ctr vs distal 0.33286, ctr vs proximal 0.00685, distal vs proximal 0.28367
biCD11c <- read_excel("is-CD11c.xlsx")
biCD11cp <- tibble::tribble(~group1, ~group2, ~p.adj, 
                           "ctr", "distal", 0.33286,
                           "ctr", "proximal", 0.00685,
                           "distal", "proximal", 0.28367)

tiff("is-CD11c.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(biCD11c, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(biCD11cp, y.position = 7, step.increase = 0.14, label = "p.adj", size = 2.7) +
  ylim(1,8.8) +
  theme_grey() +
  theme(legend.position="none") + 
  ggtitle("Islet CD11c")
dev.off()


#CD45
#adjp pval ctr vs distal 0.5456, ctr vs proximal 0.00685, distal vs proximal 0.00175
biCD45 <- read_excel("is-CD45.xlsx")
biCD45p <- tibble::tribble(~group1, ~group2, ~p.adj, 
                            "ctr", "distal", 0.5456,
                            "ctr", "proximal", 0.00685,
                            "distal", "proximal", 0.00175)

tiff("is-CD45.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(biCD45, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(biCD45p, y.position = 35, step.increase = 0.14, label = "p.adj", size = 2.7) +
  ylim(2.33,45) +
  theme_grey() +
  theme(legend.position="none") + 
  ggtitle("Islet CD45")
dev.off()


#GZMB
#adjp pval ctr vs distal 0.93915, ctr vs proximal 0.37189, distal vs proximal 0.01094
biGZMB <- read_excel("is-GZMB.xlsx")
biGZMBp <- tibble::tribble(~group1, ~group2, ~p.adj, 
                           "ctr", "distal", 0.93915,
                           "ctr", "proximal", 0.37189,
                           "distal", "proximal", 0.01094)

tiff("is-GZMB.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(biGZMB, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(biGZMBp, y.position = 55, step.increase = 0.14, label = "p.adj", size = 2.7) +
  ylim(5,70) +
  theme_grey() +
  theme(legend.position="none") + 
  ggtitle("Islet GZMB")
dev.off()


#S6 (housekeeping protein)
#adjp pval ctr vs distal 0.84139, ctr vs proximal 0.28279, distal vs proximal 0.01094
biS6 <- read_excel("is-S6.xlsx")
biS6p <- tibble::tribble(~group1, ~group2, ~p.adj, 
                           "ctr", "distal", 0.84139,
                           "ctr", "proximal", 0.28279,
                           "distal", "proximal", 0.01094)

tiff("is-S6.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
ggboxplot(biS6, x="Condition", y="Count", fill = "Condition", size=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.25), size = 1) +
  stat_pvalue_manual(biS6p, y.position = 630, step.increase = 0.14, label = "p.adj", size = 2.7) +
  ylim(90,800) +
  theme_grey() +
  theme(legend.position="none") + 
  ggtitle("Islet S6")
dev.off()



###########################################################

#OLD GGPLOT BOXPLOT VERSION
#tiff("ex-CD4.tiff", units="cm", width=4.5, height=8, res=800, compression = "lzw")
#ggplot(bexCD4, aes(x=Condition, y=Count, fill = Condition)) +
#  ylim(0,5.8) + 
#  geom_boxplot() +
#  #theme_minimal(base_size = 14) +
#  theme(legend.position="none") + 
#  geom_jitter(shape=16, position=position_jitter(0.25)) +
#  annotate("text", x=1.5, y=5.8, label= "adj.p = 0.00697", size = 2.7)
#  dev.off() 


##New PCA with saveable figure demo
library(readxl)
library(scatterplot3d)
data2 = read_excel('PCAdata2.xlsx')

data2$Condtion <- as.factor(data2$Condtion)
scatterplot3d(x = data2$PC1, y = data2$PC2, z = data2$PC3,
              cex.symbols = 0.6, pch =16, angle = -165,
              color = as.numeric(data2$Condtion),
              xlab="PC1 (0.352)", ylab="PC2 (0.281)", zlab="PC3 (0.146)")
        #export as pdf 4x4 inches så det passer med størrelse på datapunkter.
        #legend exporteres separat og sættes på i photoshop

#legend was saved separately and then added to the PCA plot in Gimp.
legend("center", legend = levels(data2$Condtion), fill = unique(data2$Condtion), title = "Condition",
       x.intersp = 0.8, y.intersp = 0.8, cex = 0.8)




