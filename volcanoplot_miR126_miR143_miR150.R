library(ggplot2)
library(RColorBrewer)
setwd("/Users/vikrants/Desktop/OIR_New/mir143_mir150_mir126/")

###########################################
#new volcano plot with angiogensis gene highlights
##############################################
mir126_angio <- read.csv("~/Desktop/OIR_New/mir143_mir150_mir126/mir126_angio.txt", sep="")
mir143_angio <- read.csv("~/Desktop/OIR_New/mir143_mir150_mir126/mir143_angio.txt", sep="")
mir150_angio <- read.csv("~/Desktop/OIR_New/mir143_mir150_mir126/mir150_angio.txt", sep="")
###############################################
oir <- read.delim("Sleuth.DE_Control_vs_miR126.txt")
oir1 <- read.delim("Sleuth.DE_Control_vs_miR143.txt")
oir2 <- read.delim("Sleuth.DE_Control_vs_miR150.txt")
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
genes$Significant <- ifelse(genes$Q.Value < 0.01, "Q.Value < 0.01", "Not Sig")
#################################################
aa<- merge(genes,mir126_angio, by="Gene")
aa1 <- subset(aa, aa$Significant == "Q.Value < 0.01")
aa1$Significant [aa1$Significant == "Q.Value < 0.01"] <- "Angiogenesis" 

aa2 <- rbind(genes,aa1)
##################################################
ppp <- ggplot(aa2, aes(x = logFC, y = -log10(Q.Value))) 
ppp +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("blue", "black", "red"), labels=c("Angiogenesis","Not Significant","Q.Value <0.01")) +
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom") +
  geom_point( alpha = 1, na.rm = T) +
  geom_point(aes(color=Significant))
ggsave("miR126_volcanoplot.tiff",
       dpi = 300,
       width = 4,
       height = 10,
       units = "in")
  
########################################################################
########################################################################
########################################################################
###
###miR143
########################################################################
########################################################################
b<- data.frame(oir1$ext_gene, oir1$qval,oir1$b)
colnames(b)[1]<- "Gene"
colnames(b)[2]<- "Q.Value"
colnames(b)[3]<- "logFC"
res1<-b
#################################################
genes1<- res1
genes1<- na.omit(genes1)
genes1$Significant <- ifelse(genes1$Q.Value < 0.01, "Q.Value < 0.01", "Not Sig")
#################################################
bb<- merge(genes1,mir143_angio, by="Gene")
bb1 <- subset(bb, bb$Significant == "Q.Value < 0.01")
bb1$Significant [bb1$Significant == "Q.Value < 0.01"] <- "Angiogenesis" 

bb2 <- rbind(genes,bb1)
##################################################

ppp1 <- ggplot(bb2, aes(x = logFC, y = -log10(Q.Value))) 
ppp1 +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("blue", "black", "red"), labels=c("Angiogenesis","Not Significant","Q.Value <0.01")) +
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom") +
  geom_point( alpha = 1, na.rm = T) +
  geom_point(aes(color=Significant))
ggsave("miR143_volcanoplot.tiff",
       dpi = 150,
       width = 6,
       height = 6,
       units = "in")
########################################################################
########################################################################
########################################################################
###
###miR150
########################################################################
########################################################################
c<- data.frame(oir2$ext_gene, oir2$qval,oir2$b)
colnames(c)[1]<- "Gene"
colnames(c)[2]<- "Q.Value"
colnames(c)[3]<- "logFC"
res2<-c
#################################################
genes2<- res2
genes2<- na.omit(genes2)
genes2$Significant <- ifelse(genes2$Q.Value < 0.01, "Q.Value < 0.01", "Not Sig")
#################################################
cc<- merge(genes2,mir150_angio, by="Gene")
cc1 <- subset(cc, cc$Significant == "Q.Value < 0.01")
cc1$Significant [cc1$Significant == "Q.Value < 0.01"] <- "Angiogenesis" 

cc2 <- rbind(genes,cc1)
##################################################

ppp2 <- ggplot(cc2, aes(x = logFC, y = -log10(Q.Value))) 
ppp2 +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("blue", "black", "red"), labels=c("Angiogenesis","Not Significant","Q.Value <0.01")) +
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom") +
  geom_point( alpha = 1, na.rm = T) +
  geom_point(aes(color=Significant))
ggsave("miR150_volcanoplot.tiff",
       dpi = 150,
       width = 6,
       height = 6,
       units = "in")
