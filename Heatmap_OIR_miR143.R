setwd("~/Desktop/OIR_New/OIR_mir143_Scrambled/kallisto/")
library(reshape) 
library(ggplot2) 
library(gplots)
library(RColorBrewer)
library(viridis) 

a<- read.delim("299_R/kallisto/abundance.tsv")
b<- read.delim("301_R/kallisto/abundance.tsv")
c<- read.delim("309_L/kallisto/abundance.tsv")
d<- read.delim("310_L/kallisto/abundance.tsv")
oir_mir143_angio <- read.csv("~/Desktop/OIR_New/OIR_mir143_Scrambled/Kallisto/oir_mir143_angio.txt", sep="")
All_data <- read.delim("~/Desktop/OIR_New/OIR_mir143_Scrambled/Sleuth.DE_OIR_scrambled_vs_OIR_miR143.txt")


a1<- data.frame(a$target_id,a$tpm)
colnames(a1)[1]<- "target_id"
colnames(a1)[2]<-"299_R"
b1<- data.frame(b$target_id,b$tpm)
colnames(b1)[1]<- "target_id"
colnames(b1)[2]<-"301_R"
c1<- data.frame(c$target_id,c$tpm)
colnames(c1)[1]<- "target_id"
colnames(c1)[2]<-"309_L"
d1<- data.frame(d$target_id,d$tpm)
colnames(d1)[1]<- "target_id"
colnames(d1)[2]<-"310_L"


a2<- merge(a1,b1, by="target_id")
c2<- merge(c1,d1,by="target_id")

oir_tpm <- merge(a2,c2,by="target_id")

ss <- data.frame(All_data$ext_gene,All_data$target_id,All_data$qval)
colnames(ss)[1]<- "Gene"
colnames(ss)[2]<-"target_id"
colnames(ss)[3]<- "Q.value"
oir_tpm1<- merge(oir_tpm,ss,by="target_id")

hhoir <- merge(oir_tpm1,oir_mir143_angio, by="Gene", all=FALSE)
hhoir<- na.omit(hhoir)
hhoir1<-subset(hhoir,hhoir$Q.value <=0.05)
heatoir<-data.frame(hhoir1$Gene,hhoir1$`299_R`,hhoir1$`301_R`,hhoir1$`309_L`,hhoir1$`310_L`)
colnames(heatoir)[1]<- "Gene"
colnames(heatoir)[2]<- "OIR_Scrambled_1"
colnames(heatoir)[3]<- "OIR_Scrambled_2"
colnames(heatoir)[4]<- "OIR_miR143_1"
colnames(heatoir)[5]<- "OIR_miR143_2"

write.csv(heatoir, file = "Supp_OIR_mir143_Angio.csv",quote = FALSE, row.names = FALSE)

data<- heatoir
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames  
heatmap.2(mat_data, trace="none", margins=c(8,8), dendrogram="both", 
          sepwidth=c(0.01,0.01),
          sepcolor="azure3",
          colsep=1:ncol(mat_data),
          rowsep=1:nrow(mat_data),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="complete"),
          cexRow =1,cexCol = 1 ,scale="row", 
          col=brewer.pal(11,"RdBu")
         )
######################################################
#######
Supp_OIR_mir143_Angio <- read.csv("~/Desktop/OIR_New/OIR_mir143_Scrambled/Kallisto/Supp_OIR_mir143_Angio.csv")

mata_data <- Supp_OIR_mir143_Angio
data<- mata_data
rnames <- data[,1]      
mat_data4 <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data4) <- rnames  


#####correcting color for Up and Down regulating genes###
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames  
heatmap.2(mat_data4, trace="none", margins=c(8,8), dendrogram="both", 
          sepwidth=c(0.01,0.01),
          sepcolor="azure3",
          colsep=1:ncol(mat_data4),
          rowsep=1:nrow(mat_data4),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="complete"),
          cexRow =1,cexCol = 1 ,scale="row", 
          col=rev(brewer.pal(11,"RdBu")),
          key= FALSE,
          symbreaks = TRUE
)