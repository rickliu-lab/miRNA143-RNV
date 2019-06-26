setwd("~/Desktop/OIR_New/mir143_mir150_mir126/kallisto/")
library(reshape) 
library(ggplot2) 
library(gplots)
library(RColorBrewer)
#####miR126#####
mir126b<- read.delim("~/Desktop/tak_new_data/results_1/126_2/kallisto/abundance.tsv")
mir126c <- read.delim("~/Desktop/tak_new_data/results_1/126_3/kallisto/abundance.tsv")
mir126d <- read.delim("~/Desktop/tak_new_data/results_1/126_4/kallisto/abundance.tsv")

a<- data.frame(mir126b$target_id,mir126b$tpm)
b<- data.frame(mir126c$target_id,mir126c$tpm)
c<- data.frame(mir126d$target_id,mir126d$tpm)


colnames(a)[1]<- "target_id"
colnames(b)[1]<- "target_id"
colnames(c)[1]<- "target_id"

mir126aa<- merge(a,b, by="target_id")
mir126_main <- merge(mir126aa,c, by="target_id")

colnames(mir126_main)[2] <- "126_2"
colnames(mir126_main)[3] <- "126_3"
colnames(mir126_main)[4] <- "126_4"

######################################################

#####miR150#####
mir150b<- read.delim("~/Desktop/tak_new_data/results_1/150_2/kallisto/abundance.tsv")
mir150c <- read.delim("~/Desktop/tak_new_data/results_1/150_3/kallisto/abundance.tsv")
mir150d <- read.delim("~/Desktop/tak_new_data/results_1/150_4/kallisto/abundance.tsv")

a<- data.frame(mir150b$target_id,mir150b$tpm)
b<- data.frame(mir150c$target_id,mir150c$tpm)
c<- data.frame(mir150d$target_id,mir150d$tpm)


colnames(a)[1]<- "target_id"
colnames(b)[1]<- "target_id"
colnames(c)[1]<- "target_id"

mir150aa<- merge(a,b, by="target_id")
mir150_main <- merge(mir150aa,c, by="target_id")

colnames(mir150_main)[2] <- "150_2"
colnames(mir150_main)[3] <- "150_3"
colnames(mir150_main)[4] <- "150_4"


######################################################

#####miR143#####
mir143b<- read.delim("~/Desktop/tak_new_data/results_1/143_2/kallisto/abundance.tsv")
mir143c <- read.delim("~/Desktop/tak_new_data/results_1/143_3/kallisto/abundance.tsv")
mir143d <- read.delim("~/Desktop/tak_new_data/results_1/143_4/kallisto/abundance.tsv")

a<- data.frame(mir143b$target_id,mir143b$tpm)
b<- data.frame(mir143c$target_id,mir143c$tpm)
c<- data.frame(mir143d$target_id,mir143d$tpm)


colnames(a)[1]<- "target_id"
colnames(b)[1]<- "target_id"
colnames(c)[1]<- "target_id"

mir143aa<- merge(a,b, by="target_id")
mir143_main <- merge(mir143aa,c, by="target_id")

colnames(mir143_main)[2] <- "143_2"
colnames(mir143_main)[3] <- "143_3"
colnames(mir143_main)[4] <- "143_4"


##############################################

#####control#####
controlb<- read.delim("~/Desktop/tak_new_data/results_1/C_2/kallisto/abundance.tsv")
controlc <- read.delim("~/Desktop/tak_new_data/results_1/C_3/kallisto/abundance.tsv")
controld <- read.delim("~/Desktop/tak_new_data/results_1/C_4/kallisto/abundance.tsv")

a<- data.frame(controlb$target_id,controlb$tpm)
b<- data.frame(controlc$target_id,controlc$tpm)
c<- data.frame(controld$target_id,controld$tpm)


colnames(a)[1]<- "target_id"
colnames(b)[1]<- "target_id"
colnames(c)[1]<- "target_id"

controlaa<- merge(a,b, by="target_id")
control_main <- merge(controlaa,c, by="target_id")

colnames(control_main)[2] <- "Control_2"
colnames(control_main)[3] <- "Control_3"
colnames(control_main)[4] <- "Control_4"



#################################################
con126 <- merge(control_main,mir126_main, by="target_id")

con143 <- merge(control_main,mir143_main, by="target_id")

con150 <- merge(control_main,mir150_main, by="target_id")
##########################################################
mir126angio<- read.delim("mir126_angio.txt")
colnames(mir126angio)[1]<- "ext_gene"
mir143angio<- read.delim("mir143_angio.txt")
colnames(mir143angio)[1]<- "ext_gene"
mir150angio<- read.delim("mir150_angio.txt")
colnames(mir150angio)[1]<- "ext_gene"
#########################################################
controlmir126 <- read.delim("~/Desktop/tak_new_data/sleuth.DE_transcripts.qval_0.01_controlvs_miR126.txt")
controlmir143 <- read.delim("~/Desktop/tak_new_data/sleuth.DE_transcripts.qval_0.01_controlvs_miR143.txt")
controlmir150 <- read.delim("~/Desktop/tak_new_data/sleuth.DE_transcripts.qval_0.01_controlvs_miR150.txt")
#########################################################
hhmir126<- merge(controlmir126,mir126angio, by="ext_gene")
hhmir126a <- data.frame(hhmir126$target_id)
colnames(hhmir126a)[1]<- "target_id"
heatmir126 <- merge(con126,hhmir126a, by="target_id")
#########################################################
hhmir143<- merge(controlmir143,mir143angio,by="ext_gene")
hhmir143a <- data.frame(hhmir143$target_id)
colnames(hhmir143a)[1]<- "target_id"
heatmir143 <- merge(con143,hhmir143a, by="target_id")
#########################################################
hhmir150<- merge(controlmir150,mir150angio,by="ext_gene")
hhmir150a<- data.frame(hhmir150$target_id)
colnames(hhmir150a)[1]<- "target_id"
heatmir150 <- merge(con150,hhmir150a, by="target_id")
#####################################################
amir126 <- merge(heatmir126,controlmir126,by="target_id")
amir126_1<- data.frame(amir126$ext_gene,amir126$Control_2,amir126$Control_3,amir126$Control_4,amir126$`126_2`,amir126$`126_3`,amir126$`126_4`)
colnames(amir126_1)[1]<- "Gene"
colnames(amir126_1)[2]<-"Control_2"
colnames(amir126_1)[3]<-"Control_3"
colnames(amir126_1)[4]<-"Control_4"
colnames(amir126_1)[5]<-"miR126_2"
colnames(amir126_1)[6]<-"miR126_3"
colnames(amir126_1)[7]<-"miR126_4"
#####################################################
#####################################################
amir143<- merge(heatmir143,controlmir143,by="target_id")
amir143_1<- data.frame(amir143$ext_gene,amir143$Control_2,amir143$Control_3,amir143$Control_4,amir143$`143_2`,amir143$`143_3`,amir143$`143_4`)
colnames(amir143_1)[1]<- "Gene"
colnames(amir143_1)[2]<-"Control_2"
colnames(amir143_1)[3]<-"Control_3"
colnames(amir143_1)[4]<-"Control_4"
colnames(amir143_1)[5]<-"miR143_2"
colnames(amir143_1)[6]<-"miR143_3"
colnames(amir143_1)[7]<-"miR143_4"

#####################################################
#####################################################
amir150<- merge(heatmir150,controlmir150,by="target_id")
amir150_1<- data.frame(amir150$ext_gene,amir150$Control_2,amir150$Control_3,amir150$Control_4,amir150$`150_2`,amir150$`150_3`,amir150$`150_4`)
colnames(amir150_1)[1]<- "Gene"
colnames(amir150_1)[2]<-"Control_2"
colnames(amir150_1)[3]<-"Control_3"
colnames(amir150_1)[4]<-"Control_4"
colnames(amir150_1)[5]<-"miR150_2"
colnames(amir150_1)[6]<-"miR150_3"
colnames(amir150_1)[7]<-"miR150_4"


write.table(amir126_1, file='Supp_miR126_heatmap.csv', sep="\t",row.names=F, quote=F)
write.table(amir143_1, file='Supp_miR143_heatmap.csv', sep="\t",row.names=F, quote=F)
write.table(amir150_1, file='Supp_miR150_heatmap.csv', sep="\t",row.names=F, quote=F)
#####################################################
#####################################################
##Heatmaps for supplemenatry figures################
######################################################
data1<- amir126_1
rnames <- data1[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data1[,2:ncol(data1)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames  
heatmap.2(mat_data, trace="none", margins=c(8,5), labRow = NA, dendrogram="both", scale="row", col=colorRampPalette(c("red","white","blue")))
######################################################
######################################################
data2<- amir143_1
rnames <- data2[,1]                            # assign labels in column 1 to "rnames"
mat_data1 <- data.matrix(data2[,2:ncol(data2)])  # transform column 2-5 into a matrix
rownames(mat_data1) <- rnames  
heatmap.2(mat_data1, trace="none", margins=c(8,5),  dendrogram="both", scale="row", col=colorRampPalette(c("red","white","blue")))
######################################################
######################################################
data3<- amir150_1
rnames <- data3[,1]                            # assign labels in column 1 to "rnames"
mat_data2 <- data.matrix(data3[,2:ncol(data3)])  # transform column 2-5 into a matrix
rownames(mat_data2) <- rnames  
heatmap.2(mat_data2, trace="none", margins=c(8,5), dendrogram="both", scale="row", col=colorRampPalette(c("red","white","blue")))
######################################################
######################################################
#correcting color order for heatmaps
###############################
heatmap.2(mat_data2, trace="none", margins=c(8,8), dendrogram="both", 
          sepwidth=c(0.01,0.01),
          sepcolor="azure3",
          colsep=1:ncol(mat_data2),
          rowsep=1:nrow(mat_data2),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="complete"),
          cexRow =1,cexCol = 1 ,scale="row", 
          col=rev(brewer.pal(11,"RdBu")),
          key= FALSE,
          symbreaks = TRUE
)
