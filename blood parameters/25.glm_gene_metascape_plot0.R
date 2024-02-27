rm(list = ls())
library(readxl)

library(reshape2)
library(gridExtra)
library(scales)
library(ggsci)
library(ggplot2)
library(stringr)
library(grid)
library(RColorBrewer)

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
show_col(pal)


##for enrichment pathway for four lists of linear genes  
#construct dataframe for all target GO terms 
Decidual_Down <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Metascape_Enrichment_result/glm/decidua_Down_gene_glm/metascape_result.xlsx",sheet = 2))
Decidual_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Metascape_Enrichment_result/glm/decidua_Up_gene_glm/metascape_result.xlsx",sheet = 2))
Villus_Down <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Metascape_Enrichment_result/glm/Vill_Down_gene_glm/metascape_result.xlsx",sheet = 2))
Villus_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Metascape_Enrichment_result/glm/Vill_Up_gene_glm/metascape_result.xlsx",sheet = 2))

Decidual_Down$GeneList<-"Decidual_Down"
Decidual_Up$GeneList<-"Decidual_Up"
Villus_Down$GeneList<-"Villus_Down"
Villus_Up$GeneList<-"Villus_Up"

####merge
Enrich_all_plot<-as.data.frame(rbind(Decidual_Down[grep("Summary",Decidual_Down$GroupID),],Decidual_Up[grep("Summary",Decidual_Up$GroupID),],
                                     Villus_Down[grep("Summary",Villus_Down$GroupID),],Villus_Up[grep("Summary",Villus_Up$GroupID),]))

dim(Enrich_all_plot)#80 10
head(Enrich_all_plot)
Enrich_all_plot$count <- as.numeric(sapply(strsplit(Enrich_all_plot$InTerm_InList, "/"), function(x) x[1]))

Enrich_all_plot2<-Enrich_all_plot[,c("Term","Description","count","LogP","Log(q-value)","Symbols","GeneList")];dim(Enrich_all_plot2)#340  10

colnames(Enrich_all_plot2)<-c("Category","Description","Count","Log_pvalue","Log_qvalue","Hits","group")
table(Enrich_all_plot2$group)
head(Enrich_all_plot2,n=20)
write.table(Enrich_all_plot2, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Metascape_Enrichment_result/glm/enrichment_merge_TOP20_four_lists.txt",row.names=T, col.names=T) 


##plot enrichment
#construct dataframe for all target GO terms 
Decidual_Down <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Decidual_Down/metascape_result.xlsx",sheet = 2))
Decidual_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Decidual_Up/metascape_result.xlsx",sheet = 2))
Villus_Down <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Villus_Down/metascape_result.xlsx",sheet = 2))
Villus_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Villus_Up/metascape_result.xlsx",sheet = 2))
Decidual_Up_Villus_Down <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Decidual_Up_Villus_Down/metascape_result.xlsx",sheet = 2))
Decidual_Down_Villus_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Decidual_Down_Villus_Up/metascape_result.xlsx",sheet = 2))

Decidual_Down$GeneList<-"Decidual_Down"
Decidual_Up$GeneList<-"Decidual_Up"
Villus_Down$GeneList<-"Villus_Down"
Villus_Up$GeneList<-"Villus_Up"
Decidual_Up_Villus_Down$GeneList<-"Decidual_Up_Villus_Down"
Decidual_Down_Villus_Up$GeneList<-"Decidual_Down_Villus_Up"

####Up_down
Enrich_all_plot<-as.data.frame(rbind(Decidual_Up_Villus_Down[grep("Summary",Decidual_Up_Villus_Down$GroupID),],Decidual_Down_Villus_Up[grep("Summary",Decidual_Down_Villus_Up$GroupID),]))

dim(Enrich_all_plot)#22 10
head(Enrich_all_plot)
Enrich_all_plot$count <- as.numeric(sapply(strsplit(Enrich_all_plot$InTerm_InList, "/"), function(x) x[1]))

Enrich_all_plot2<-Enrich_all_plot[,c("Term","Description","count","LogP","Log(q-value)","Symbols","GeneList")];dim(Enrich_all_plot2)#340  10

colnames(Enrich_all_plot2)<-c("Category","Description","Count","Log_pvalue","Log_qvalue","Hits","group")

Enrich_all_plot2$pathway<-paste0(Enrich_all_plot2$Category,":",Enrich_all_plot2$Description)
head(Enrich_all_plot2)
Enrich_all_plot2<-na.omit(Enrich_all_plot2)
Enrich_all_plot2$Log_pvalue<-c(-Enrich_all_plot2$Log_pvalue)
Enrich_all_plot2$Log_qvalue<-c(-Enrich_all_plot2$Log_qvalue)

Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_pvalue)),]
Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_qvalue)),]

range(Enrich_all_plot2$Log_pvalue)# 2.012332 4.294106
range(Enrich_all_plot2$Count)# 3 8
head(Enrich_all_plot2)

#2.论文画法：
table(Enrich_all_plot2$group)
Enrich_all_plot_down_up<-Enrich_all_plot2[which(Enrich_all_plot2$group == "Decidual_Down_Villus_Up"),]
Enrich_all_plot_up_down<-Enrich_all_plot2[which(Enrich_all_plot2$group == "Decidual_Up_Villus_Down"),]

#自定义主题：
Enrich_all_plot_down_up$text_x <- rep(0.03,nrow(Enrich_all_plot_down_up)) 
Enrich_all_plot_down_up<-Enrich_all_plot_down_up[order(Enrich_all_plot_down_up$Log_pvalue,decreasing = T),]
Enrich_all_plot_down_up$pathway<-factor(Enrich_all_plot_down_up$pathway,levels = rev(c(unique(as.character(Enrich_all_plot_down_up$pathway)))))

p1<- ggplot(data = Enrich_all_plot_down_up,aes(x = Log_pvalue, y = pathway)) +
  geom_bar(aes(fill = Log_pvalue), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "Log_pvalue", y = "pathway", title = "Metascape enrichment barplot:Decidual_Down_Villus_Up") +
  geom_text(aes(x = text_x, label= pathway),hjust= 0)+ #hjust=0，左对齐
  scale_x_continuous(limits = c(0,5), breaks = seq(0, 5, by = 1))+  # 设置刻度范围和间隔
  theme_classic()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title = element_text(size = 13),axis.text = element_text(size = 11),plot.title = element_text(size = 14,hjust= 0.5,face= "bold"),legend.title = element_text(size = 13),legend.text = element_text(size = 11))
p1

#自定义主题：
Enrich_all_plot_up_down$text_x <- rep(0.03,nrow(Enrich_all_plot_up_down)) 
Enrich_all_plot_up_down<-Enrich_all_plot_up_down[order(Enrich_all_plot_up_down$Log_pvalue,decreasing = T),]
Enrich_all_plot_up_down$pathway<-factor(Enrich_all_plot_up_down$pathway,levels = rev(c(unique(as.character(Enrich_all_plot_up_down$pathway)))))

p2<- ggplot(data = Enrich_all_plot_up_down,aes(x = Log_pvalue, y = pathway)) +
  geom_bar(aes(fill = Log_pvalue), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Log_pvalue", y = "pathway", title = "Metascape enrichment barplot:Decidual_Down_Villus_Up") +
  geom_text(aes(x = text_x, label= pathway),hjust= 0)+ #hjust=0，左对齐
  scale_x_continuous(limits = c(0,4), breaks = seq(0, 5, by = 1))+  # 设置刻度范围和间隔
  theme_classic()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title = element_text(size = 13),axis.text = element_text(size = 11),plot.title = element_text(size = 14,hjust= 0.5,face= "bold"),legend.title = element_text(size = 13),legend.text = element_text(size = 11))
p2
combined_plot <- grid.arrange(p1, p2, ncol = 1, heights  = c(12,10))

combined_plot <- grid.arrange(p1, p2, ncol = 1)  # 将两个图按2列排列
ggsave(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Down_Up_Enrich_Metascape_result.pdf",combined_plot,width = 6, height =6)

###输出其他图
Enrich_all_number<-data.frame(table(as.character(Enrich_all_plot2$Description)))
#Enrich_all_number2<-Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot3<-merge(Enrich_all_plot2,Enrich_all_number)
dim(Enrich_all_plot3)#  80  8
Enrich_all_plot3<-Enrich_all_plot3[order(Enrich_all_plot3$freq,Enrich_all_plot3$Description,Enrich_all_plot3$group,Enrich_all_plot3$Hits,decreasing = T),]
head(Enrich_all_plot3)
write.table(Enrich_all_plot3, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Contract_Decidua_vill_up_down_lm_gene_enrichment_merge.txt",row.names=T, col.names=T) 

dim(Enrich_all_plot3[which(Enrich_all_plot3$Log_pvalue > -log10(0.05)),])#22 9
dim(Enrich_all_plot3[which(Enrich_all_plot3$Log_qvalue > -log10(0.05)),])#0  9

Enrich_all_plot3$group<-factor(Enrich_all_plot3$group,levels = c("Decidual_Down_Villus_Up","Decidual_Up_Villus_Down"))
Enrich_all_plot3<-Enrich_all_plot3[order(Enrich_all_plot3$group,decreasing = T ),]
Enrich_all_plot3$Description<-factor(Enrich_all_plot3$Description,levels = rev(c(unique(as.character(Enrich_all_plot3$Description)))))

plot_BP<-ggplot(Enrich_all_plot3, aes(x=group,y=Description,size=Count,colour=Log_pvalue))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,2,3,4,5,6,7,8),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,1,2,3,4,5),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="metascape enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Contract_metascape_Down_Up_Enrich_Metascape_result.pdf",plot_BP,width = 6, height =6)
write.table(Enrich_all_plot3, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/Decidua_vill_up_down_lm_gene_GO_BP_merge.txt",row.names=T, col.names=T) 

#commom
Enrich_all_plot<-as.data.frame(rbind(Decidual_Down[grep("Summary",Decidual_Down$GroupID),],Decidual_Up[grep("Summary",Decidual_Up$GroupID),],
                                     Villus_Down[grep("Summary",Villus_Down$GroupID),],Villus_Up[grep("Summary",Villus_Up$GroupID),]))

dim(Enrich_all_plot)#80 10
head(Enrich_all_plot)
Enrich_all_plot$count <- as.numeric(sapply(strsplit(Enrich_all_plot$InTerm_InList, "/"), function(x) x[1]))

Enrich_all_plot2<-Enrich_all_plot[,c("Term","Description","count","LogP","Log(q-value)","Symbols","GeneList")];dim(Enrich_all_plot2)#340  10

colnames(Enrich_all_plot2)<-c("Category","Description","Count","Log_pvalue","Log_qvalue","Hits","group")
head(Enrich_all_plot2)
Enrich_all_plot2<-na.omit(Enrich_all_plot2)
Enrich_all_plot2$Log_pvalue<-c(-Enrich_all_plot2$Log_pvalue)
Enrich_all_plot2$Log_qvalue<-c(-Enrich_all_plot2$Log_qvalue)
Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_pvalue)),]
Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_qvalue)),]

range(Enrich_all_plot2$Log_pvalue)#3.04682 25.06279
range(Enrich_all_plot2$Count)#  4 182

head(Enrich_all_plot2)
Enrich_all_number<-data.frame(table(as.character(Enrich_all_plot2$Description)))
#Enrich_all_number2<-Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot3<-merge(Enrich_all_plot2,Enrich_all_number)
dim(Enrich_all_plot3)#  80  8
Enrich_all_plot3<-Enrich_all_plot3[order(Enrich_all_plot3$freq,Enrich_all_plot3$Description,Enrich_all_plot3$group,Enrich_all_plot3$Hits,decreasing = T),]
head(Enrich_all_plot3)
write.table(Enrich_all_plot3, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/result_output/Decidua_vill_up_down_lm_gene_enrichment_merge.txt",row.names=T, col.names=T) 

dim(Enrich_all_plot3[which(Enrich_all_plot3$Log_pvalue > -log10(0.05)),])#80  8
dim(Enrich_all_plot3[which(Enrich_all_plot3$Log_qvalue > -log10(0.05)),])#80  8

Enrich_all_plot3$group<-factor(Enrich_all_plot3$group,levels = c("Villus_Up","Villus_Down","Decidual_Up","Decidual_Down"))
Enrich_all_plot3<-Enrich_all_plot3[order(Enrich_all_plot3$freq,Enrich_all_plot3$group,decreasing = T ),]
Enrich_all_plot3$Description<-factor(Enrich_all_plot3$Description,levels = rev(c(unique(as.character(Enrich_all_plot3$Description)))))
length(unique(as.character(Enrich_all_plot3[which(Enrich_all_plot3$freq>1),]$Description)))#16 
range(Enrich_all_plot3$Count)#4 182
range(Enrich_all_plot3$Log_pvalue)#3.04682 25.06279

plot_BP<-ggplot(Enrich_all_plot3, aes(x=group,y=Description,size=Count,colour=Log_pvalue))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,5,10,20,40,80,120,160),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,3,6,9,12,15,18,21,24),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="AMA related mom and kids commmon DEGs: BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/result_output/Four_up_down_metascape_Enrich_Metascape_result.pdf",plot_BP,width = 8, height =12)
write.table(Enrich_all_plot3, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/glm_relationship/metascape_result/result_output/Four_up_down_Decidua_vill_lm_gene_Metascape_result_merge.txt",row.names=T, col.names=T) 