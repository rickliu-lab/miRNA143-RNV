library(ggplot2)
library(ggpubr)
library(cowplot)
library(gridExtra)

setwd(dir = "/Users/vikrants/Desktop/bubbleplot_ipa/new_1_06_2019/")
df<- read.delim("mir150_IPA.txt")

a <- ggplot(df, aes(y = df$Diseases.or.Functions.Annotation, x = -log10(df$p.value), label = df$Diseases.or.Functions.Annotation)) +
  geom_jitter(aes(size = df$Molecules, colour = df$Diseases.or.Functions.Annotation),alpha = .8,show.legend=F) +
  #geom_point(aes(size=gene))+)
  geom_text(hjust = 0, size = 11) +
  scale_size(range = c(6,30)) +
  xlab("-log10(p value)")+ ylab("miR150") +
  xlim(c(2,7))+
  #ylim(c(0,df$Name))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(3)))+
  theme(
    axis.title=element_text(size=40,face="bold"),
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
#dimensions
#900X1200
#############################################################################
#############################################################################
#############################################################################
#############################################################################
df1<- read.delim("mir126_IPA.txt")

b <- ggplot(df1, aes(y = df1$Diseases.or.Functions.Annotation, x = -log10(df1$p.value), label = df1$Diseases.or.Functions.Annotation)) +
  geom_jitter(aes(size = df1$Molecules, colour = df1$Diseases.or.Functions.Annotation),alpha = .8,show.legend=F) +
  #geom_point(aes(size=gene))+
  geom_text(hjust = 0, size = 11) +
  scale_size(range = c(10,160)) +
  xlab("-log10(p value)")+ ylab("miR126") +
  xlim(c(3,22))+
  #ylim(c(0,df$Name))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(3)))+
  theme(
    axis.title=element_text(size=40,face="bold"),
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
#dimension
#2450X3725
#############################################################################
#############################################################################
#############################################################################
#############################################################################
df2<- read.delim("mir143_IPA.txt")

c <- ggplot(df2, aes(y = df2$Diseases.or.Functions.Annotation, x = -log10(df2$p.value), label = df2$Diseases.or.Functions.Annotation)) +
  geom_jitter(aes(size = df2$Molecules, colour = df2$Diseases.or.Functions.Annotation),alpha = .8,show.legend=F) +
  #geom_point(aes(size=gene))+
  geom_text(hjust = 0, size = 11) +
  scale_size(range = c(6,60)) +
  xlab("-log10(p value)")+ ylab("miR143") +
  xlim(c(2,10))+
  #ylim(c(0,df$Name))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(3)))+
  theme(
    axis.title=element_text(size=40,face="bold"),
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
#dimensions
#900X1200
#############################################################################
#############################################################################
#############################################################################
#############################################################################
df3<- read.delim("invivo_mir143_IPA.txt")

d <- ggplot(df3, aes(y = df3$Diseases.or.Functions.Annotation, x = -log10(df3$p.value), label = df3$Diseases.or.Functions.Annotation)) +
  geom_jitter(aes(size = df3$Molecules, colour = df3$Diseases.or.Functions.Annotation),alpha = .8,show.legend=F) +
  #geom_point(aes(size=gene))+
  geom_text(hjust = 0, size = 12) +
  scale_size(range = c(5,20)) +
  xlab("-log10(p value)")+ ylab("OIR-miR143") +
  xlim(c(3,6))+
  #ylim(c(0,df$Name))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(2)))+
  theme(
    axis.title=element_text(size=48,face="bold"),
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
#dimensions
#1200X900
#############################################################################
#############################################################################
#############################################################################
#############################################################################
df4<- read.delim("invivo_mir143_neuron_IPA.txt")

e <- ggplot(df4, aes(y = df4$Specific.diseases.or.functions.annotation, x = -log10(df4$p.value), label = df4$Specific.diseases.or.functions.annotation)) +
  geom_jitter(aes(size = df4$Molecules, colour = df4$Specific.diseases.or.functions.annotation),alpha = .8,show.legend=F) +
  #geom_point(aes(size=gene))+
  geom_text(hjust = 0, size = 12) +
  scale_size(range = c(5,30)) +
  xlab("-log10(p value)")+ ylab("OIR-miR143") +
  xlim(c(3,13))+
 
  #ylim(c(0,df$Name))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(2)))+
  theme(
    axis.title=element_text(size=48,face="bold"),
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
#dimension
#900X1200
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#ggarrange(a,b,c, nrow = 2, ncol = 2)

grid.arrange(arrangeGrob(a,c, ncol=1, nrow=2),
             arrangeGrob(b, ncol=1, nrow=1), heights=c(10,1), widths=c(1,1))

#dimension
#3200X2900

