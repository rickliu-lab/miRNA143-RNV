library(ggplot2)
library(dplyr)
setwd("/Users/vikrants/Desktop/OIR_New/OIR_mir143_Scrambled/")
oir <- read.delim("Sleuth.DE_OIR_scrambled_vs_OIR_miR143.txt")
#new volcano plot with angiogensis gene highlights
##############################################
mir143_angio <- read.csv("~/Desktop/OIR_New/OIR_mir143_Scrambled/oir_mir143_angio.txt", sep="")
###############################################
###############################################
####miR126#####
##############################################
a<- data.frame(oir$ext_gene, oir$qval,oir$b)
colnames(a)[1]<- "Gene"
colnames(a)[2]<- "Q.Value"
colnames(a)[3]<- "logFC"
res<-a
#################################################
genes<- res
genes<- na.omit(genes)
genes$Significant <- ifelse(genes$Q.Value < 0.05, "Q.Value < 0.05", "Not Sig")
#################################################
aa<- merge(genes,mir143_angio, by="Gene", ignore.case=FALSE)
aa1 <- subset(aa, aa$Significant == "Q.Value < 0.05")
aa1$Significant [aa1$Significant == "Q.Value < 0.05"] <- "Angiogenesis" 

aa2 <- rbind(genes,aa1)
##################################################
ppp <- ggplot(aa2, aes(x = logFC, y = -log10(Q.Value))) 
ppp +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("blue", "black", "red"), labels=c("Angiogenesis","Not Significant","Q.Value <0.05")) +
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom") +
  geom_point( alpha = 1, na.rm = T) +
  geom_point(aes(color=Significant))
ggsave("OIR_miR143_volcanoplot.tiff",
       dpi = 150,
       width = 6,
       height = 6,
       units = "in")
