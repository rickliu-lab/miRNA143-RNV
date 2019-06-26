setwd("~/Desktop/tak_new_data/")
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
###For Correlation plot####
con126$control <- (con126$Control_2+con126$Control_3+con126$Control_4)/3
con126$mir126<- (con126$`126_2`+con126$`126_3`+con126$`126_4`)/3

con143$control <- (con143$Control_2+con143$Control_3+con143$Control_4)/3
con143$mir143<- (con143$`143_2`+con143$`143_3`+con143$`143_4`)/3

con150$control <- (con150$Control_2+con150$Control_3+con150$Control_4)/3
con150$mir150<- (con150$`150_2`+con150$`150_3`+con150$`150_4`)/3
############################################################
##correlation plot###
##########################################################

####correlation plot for scrambled###
#######mir126######
ggplot(con126, aes(log10(con126$control),log10(con126$mir126))) + xlim(-3,6) + ylim(-3,6) + xlab("log10(Control_readcounts)") +
  ylab("log10(miR126_readcounts)") +  
  #ggtitle("Control vs miR126") +
  geom_smooth(method=glm,se=FALSE,fullrange=TRUE, color="darkred") +
  geom_point(data=con126,colour = "cornflowerblue", size =0.005) + 
  #geom_point(data = mlscr,aes(log(mlscr$hi_li_scr_main.Control.Scrambled.RNA),log(mlscr$hi_li_scr_main.Control.miR143),labels=mlscr$hi_li_scr_main.tracking_id), colour="lightpink4", size=4)+
  #geom_text(data = mlscr,aes(log(mlscr$hi_li_scr_main.Control.Scrambled.RNA),log(mlscr$hi_li_scr_main.Control.miR143),label=mlscr$hi_li_scr_main.tracking_id), colour="midnightblue", size=4,hjust=0,vjust=1,angle=-15)+
  
  #geom_smooth(se = FALSE)+
  #geom_bar(stat = "identity", width = .5) +
  #geom_text(hjust=-0,vjust=0, stat = "identity")+
  annotate(x=3, y=6, 
           label=paste("R = ", round(cor(con126$control, con126$mir126),2)), 
           geom="text", size=5, colour= "black")+
  theme_classic()
######################
#######mir143######
ggplot(con143, aes(log10(con143$control),log10(con143$mir143))) + xlim(-3,6) + ylim(-3,6) + xlab("log10(Control_readcounts)") +
  ylab("log10(miR143_readcounts)") +  
  #ggtitle("Control vs miR126") +
  geom_smooth(method=glm,se=FALSE,fullrange=TRUE, color="darkred") +
  geom_point(data=con143,colour = "cornflowerblue", size =0.005) + 
  #geom_point(data = mlscr,aes(log(mlscr$hi_li_scr_main.Control.Scrambled.RNA),log(mlscr$hi_li_scr_main.Control.miR143),labels=mlscr$hi_li_scr_main.tracking_id), colour="lightpink4", size=4)+
  #geom_text(data = mlscr,aes(log(mlscr$hi_li_scr_main.Control.Scrambled.RNA),log(mlscr$hi_li_scr_main.Control.miR143),label=mlscr$hi_li_scr_main.tracking_id), colour="midnightblue", size=4,hjust=0,vjust=1,angle=-15)+
  
  #geom_smooth(se = FALSE)+
  #geom_bar(stat = "identity", width = .5) +
  #geom_text(hjust=-0,vjust=0, stat = "identity")+
  annotate(x=3, y=6, 
           label=paste("R = ", round(cor(con143$control, con143$mir143),2)), 
           geom="text", size=5, colour= "black")+
  theme_classic()
######################
#######mir150######
ggplot(con150, aes(log10(con150$control),log10(con150$mir150))) + xlim(-3,6) + ylim(-3,6) + xlab("log10(Control_readcounts)") +
  ylab("log10(miR150_readcounts)") +  
  #ggtitle("Control vs miR126") +
  geom_smooth(method=glm,se=FALSE,fullrange=TRUE, color="darkred") +
  geom_point(data=con126,colour = "cornflowerblue", size =0.005) + 
  #geom_point(data = mlscr,aes(log(mlscr$hi_li_scr_main.Control.Scrambled.RNA),log(mlscr$hi_li_scr_main.Control.miR143),labels=mlscr$hi_li_scr_main.tracking_id), colour="lightpink4", size=4)+
  #geom_text(data = mlscr,aes(log(mlscr$hi_li_scr_main.Control.Scrambled.RNA),log(mlscr$hi_li_scr_main.Control.miR143),label=mlscr$hi_li_scr_main.tracking_id), colour="midnightblue", size=4,hjust=0,vjust=1,angle=-15)+
  
  #geom_smooth(se = FALSE)+
  #geom_bar(stat = "identity", width = .5) +
  #geom_text(hjust=-0,vjust=0, stat = "identity")+
  annotate(x=3, y=6, 
           label=paste("R = ", round(cor(con150$control, con150$mir150),2)), 
           geom="text", size=5, colour= "black")+
  theme_classic()
######################
#######################################################
controlmir126 <- read.delim("~/Desktop/tak_new_data/sleuth.DE_transcripts.qval_0.01_controlvs_miR126.txt")
controlmir143 <- read.delim("~/Desktop/tak_new_data/sleuth.DE_transcripts.qval_0.01_controlvs_miR143.txt")
controlmir150 <- read.delim("~/Desktop/tak_new_data/sleuth.DE_transcripts.qval_0.01_controlvs_miR150.txt")


##########################################

mir126_heat <- merge(con126,controlmir126, by="target_id")
mir143_heat <- merge(con143,controlmir143, by="target_id")
mir150_heat<- merge(con150,controlmir150, by="target_id")

mir126_heata <- data.frame(mir126_heat$ext_gene,mir126_heat$Control_2,mir126_heat$Control_3,mir126_heat$Control_4,mir126_heat$`126_2`,mir126_heat$`126_3`,mir126_heat$`126_4`)
colnames(mir126_heata)[1]<- "Gene"
colnames(mir126_heata)[2]<- "Control_2"
colnames(mir126_heata)[3]<- "Control_3"
colnames(mir126_heata)[4]<- "Control_4"
colnames(mir126_heata)[5] <- "126_2"
colnames(mir126_heata)[6] <- "126_3"
colnames(mir126_heata)[7] <- "126_4"
mir143_heata <- data.frame(mir143_heat$ext_gene,mir143_heat$Control_2,mir143_heat$Control_3,mir143_heat$Control_4,mir143_heat$`143_2`,mir143_heat$`143_3`,mir143_heat$`143_4`)
colnames(mir143_heata)[1]<- "Gene"
colnames(mir143_heata)[2]<- "Control_2"
colnames(mir143_heata)[3]<- "Control_3"
colnames(mir143_heata)[4]<- "Control_4"
colnames(mir143_heata)[5] <- "143_2"
colnames(mir143_heata)[6] <- "143_3"
colnames(mir143_heata)[7] <- "143_4"

mir150_heata <- data.frame(mir150_heat$ext_gene,mir150_heat$Control_2,mir150_heat$Control_3,mir150_heat$Control_4,mir150_heat$`150_2`,mir150_heat$`150_3`,mir150_heat$`150_4`)
colnames(mir150_heata)[1]<- "Gene"
colnames(mir150_heata)[2]<- "Control_2"
colnames(mir150_heata)[3]<- "Control_3"
colnames(mir150_heata)[4]<- "Control_4"
colnames(mir150_heata)[5] <- "150_2"
colnames(mir150_heata)[6] <- "150_3"
colnames(mir150_heata)[7] <- "150_4"

write.table(mir126_heata, file='controlvs_miR126_heatmap.csv', sep="\t",row.names=F, quote=F)
write.table(mir143_heata, file='controlvs_miR143_heatmap.csv', sep="\t",row.names=F, quote=F)
write.table(mir150_heata, file='controlvs_miR150_heatmap.csv', sep="\t",row.names=F, quote=F)

######################

x<- data.frame(mir126_heata$Gene)
colnames(x)[1]<- "Gene"
x<- unique(x)

y<- data.frame(mir143_heata$Gene)
colnames(y)[1]<- "Gene"
y<- unique(y) 
z<- data.frame(mir150_heata$Gene)
colnames(z)[1]<- "Gene"
z<- unique(z)
x1<- merge(x,y, by="Gene")
x2<- merge(x,z, by="Gene")
x3<- merge(y,z, by="Gene")
#########

xx<- merge(x1,x2, by="Gene", all = FALSE)
xx2<- merge(xx,x3, by= "Gene", all= FALSE)
########################################
mir150 <- read.delim("~/Desktop/tak_new_data/mir150.txt")
mir143 <- read.delim("~/Desktop/tak_new_data/mir143.txt")
mir126 <- read.delim("~/Desktop/tak_new_data/mir126.txt")

##############
mir150_1<- data.frame(mir150$Gene.Symbol)
colnames(mir150_1)[1]<-"Gene"

mir143_1<- data.frame(mir143$Gene.Symbol)
colnames(mir143_1)[1]<-"Gene"

mir126_1<- data.frame(mir126$Gene.Symbol)
colnames(mir126_1)[1]<-"Gene"

mirdata126 <- merge(x,mir126_1, by="Gene", all = FALSE)
write.table(mirdata126, file='mirdatabase_mir126_comparsion_overlap_genes.csv', sep="\t",row.names=F, quote=F)

mirdata143<- merge(y,mir143_1, by="Gene", all = FALSE)
write.table(mirdata143, file='mirdatabase_mir143_comparsion_overlap_genes.csv', sep="\t",row.names=F, quote=F)

mirdata150<- merge(z,mir150_1, by="Gene", all = FALSE)
write.table(mirdata150, file='mirdatabase_mir150_comparsion_overlap_genes.csv', sep="\t",row.names=F, quote=F)

