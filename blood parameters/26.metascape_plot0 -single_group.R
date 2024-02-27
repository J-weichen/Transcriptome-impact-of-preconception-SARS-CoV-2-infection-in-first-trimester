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

##plot enrichment
#construct dataframe for all target GO terms 
Decidual_D7_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/decidua_result/Day_group_DEGs/final_single_group/metascape/Decidua_D7_Up/metascape_result.xlsx",sheet = 2))
Decidual_D7_Down <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/decidua_result/Day_group_DEGs/final_single_group/metascape/Decidua_D7_Down/metascape_result.xlsx",sheet = 2))
Decidual_D77_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/decidua_result/Day_group_DEGs/final_single_group/metascape/Decidua_D77_Up/metascape_result.xlsx",sheet = 2))

Vill_D7_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/vill_result/Day_group_DEGs/New/Final_result_single_group/metascape_result/Vill_D7_Up/metascape_result.xlsx",sheet = 2))
Vill_D7_Down <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/vill_result/Day_group_DEGs/New/Final_result_single_group/metascape_result/Vill_D7_Down/metascape_result.xlsx",sheet = 2))
Vill_D63_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/vill_result/Day_group_DEGs/New/Final_result_single_group/metascape_result/Vill_D63_Up/metascape_result.xlsx",sheet = 2))

Decidual_D7_Up$GeneList<-"Decidual_D7_Up"
Decidual_D7_Down$GeneList<-"Decidual_D7_Down"
Decidual_D77_Up$GeneList<-"Decidual_D77_Up"
Vill_D7_Up$GeneList<-"Vill_D7_Up"
Vill_D7_Down$GeneList<-"Vill_D7_Down"
Vill_D63_Up$GeneList<-"Vill_D63_Up"

####For target six group
Enrich_all_plot<-as.data.frame(rbind(Decidual_D7_Up[grep("Summary",Decidual_D7_Up$GroupID),],
                                     Decidual_D7_Down[grep("Summary",Decidual_D7_Down$GroupID),],
                                     Decidual_D77_Up[grep("Summary",Decidual_D77_Up$GroupID),],
                                     Vill_D7_Up[grep("Summary",Vill_D7_Up$GroupID),],
                                     Vill_D7_Down[grep("Summary",Vill_D7_Down$GroupID),],
                                     Vill_D63_Up[grep("Summary",Vill_D63_Up$GroupID),]))

dim(Enrich_all_plot)#111  10
head(Enrich_all_plot)
Enrich_all_plot$count <- as.numeric(sapply(strsplit(Enrich_all_plot$InTerm_InList, "/"), function(x) x[1]))
write.table(Enrich_all_plot, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Metascape_Enrichment_result/DEGs/supplementary3_Decidua_vill_DEG_Day7_Day77_Day63_enrichment_merge_single_group.txt",row.names=T, col.names=T) 


Enrich_all_plot2<-Enrich_all_plot[,c("GroupID","Term","Description","count","LogP","Log(q-value)","Symbols","GeneList")];dim(Enrich_all_plot2)#340  10

colnames(Enrich_all_plot2)<-c("GroupID","Category","Description","Count","Log_pvalue","Log_qvalue","Hits","group")

Enrich_all_plot2$pathway<-paste0(Enrich_all_plot2$Category,":",Enrich_all_plot2$Description)
head(Enrich_all_plot2)
Enrich_all_plot2<-na.omit(Enrich_all_plot2)
Enrich_all_plot2$Log_pvalue<-c(-Enrich_all_plot2$Log_pvalue)
Enrich_all_plot2$Log_qvalue<-c(-Enrich_all_plot2$Log_qvalue)

Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_pvalue)),]
Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_qvalue)),]

range(Enrich_all_plot2$Log_pvalue)# 2.432989 25.150316
range(Enrich_all_plot2$Log_qvalue)# 0.00000 20.80385

range(Enrich_all_plot2$Count)# 3 130
head(Enrich_all_plot2)

###输出其他图
Enrich_all_number<-data.frame(table(as.character(Enrich_all_plot2$Description)))
#Enrich_all_number2<-Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot3<-merge(Enrich_all_plot2,Enrich_all_number)
dim(Enrich_all_plot3)#104    10
Enrich_all_plot3<-Enrich_all_plot3[order(Enrich_all_plot3$freq,Enrich_all_plot3$Description,Enrich_all_plot3$group,Enrich_all_plot3$Hits,decreasing = T),]
head(Enrich_all_plot3)
write.table(Enrich_all_plot3, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Metascape_Enrichment_result/DEGs/Decidua_vill_DEG_Day7_Day77_Day63_enrichment_merge_single_group.txt",row.names=T, col.names=T) 

dim(Enrich_all_plot3[which(Enrich_all_plot3$Log_pvalue > -log10(0.05)),])# 104  10
dim(Enrich_all_plot3[which(Enrich_all_plot3$Log_qvalue > -log10(0.05)),])#84   10
table(Enrich_all_plot3$group)
#Decidual_D7_Down   Decidual_D7_Up  Decidual_D77_Up      Vill_D63_Up     Vill_D7_Down       Vill_D7_Up 
#           4               20               20               20               20               20

##提取各组前top6
Enrich_all_plot4<-Enrich_all_plot3[which(Enrich_all_plot3$GroupID %in% c("1_Summary","2_Summary","3_Summary","4_Summary","5_Summary")),]
Enrich_all_plot4$GroupID<-factor(Enrich_all_plot4$GroupID,levels = c("1_Summary","2_Summary","3_Summary","4_Summary","5_Summary"))
Enrich_all_plot4$group<-factor(Enrich_all_plot4$group,levels = c("Decidual_D7_Up","Decidual_D7_Down","Decidual_D77_Up","Vill_D7_Up","Vill_D7_Down","Vill_D63_Up"))
Enrich_all_plot4<-Enrich_all_plot4[order(Enrich_all_plot4$group,Enrich_all_plot4$GroupID,decreasing = F),]
Enrich_all_plot4$pathway<-factor(Enrich_all_plot4$pathway,levels = c(unique(as.character(Enrich_all_plot4$pathway))))
Enrich_all_plot4$group<-factor(Enrich_all_plot4$group,levels = rev(c("Decidual_D7_Up","Decidual_D7_Down","Decidual_D77_Up","Vill_D7_Up","Vill_D7_Down","Vill_D63_Up")))

length(unique(as.character(Enrich_all_plot4[which(Enrich_all_plot4$freq>1),]$pathway)))# 12
range(Enrich_all_plot4$Count)#  4 130
range(Enrich_all_plot4$Log_pvalue)#  2.432989 25.150316

plot_BP<-ggplot(Enrich_all_plot4, aes(x=pathway,y=group,size=Count,colour=Log_pvalue))+geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,20,40,60,80,100,120,140),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,5,10,15,20,25),name='-log10(pvalue)')+ 
  scale_x_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_bw()+labs(x="",y="GO terms",title="metascape enrichment")+
  theme(legend.text = element_text(size = 8, colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 8,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 8,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Metascape_Enrichment_result/DEGs/Decidua_vill_DEG_Day7_Day77_Day63_enrichment_merge_Metascape_result_single_group.pdf",plot_BP,width =10, height =5)
write.table(Enrich_all_plot4, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Metascape_Enrichment_result/DEGs/Decidua_vill_DEG_Day7_Day77_Day63_enrichment_merge_TOP5_single_group.txt",row.names=T, col.names=T) 

########################single group common###########################
#construct dataframe for all target GO terms 
#Decidual_common_Up_D35_D77 <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/decidua_result/Day_group_DEGs/final_single_group/metascape/Decidual_DEGs_common_Up_D35_D77/metascape_result.xlsx",sheet = 2))
Villi_gene_common_Up_D21_D49_D63_D77 <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/vill_result/Day_group_DEGs/New/Final_result_single_group/metascape_result/Villi_gene_common_Up_D21_D49_D63_D77/metascape_result.xlsx",sheet = 2))
Villi_gene_common_Up_D49_D63_D77_sep <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/vill_result/Day_group_DEGs/New/Final_result_single_group/metascape_result/Villi_gene_common_Up_D49_D63_D77_sep/metascape_result.xlsx",sheet = 2))

Villi_gene_common_Up_D21_D49_D63_D77$GeneList<-"Villi_gene_common_Up_D21_D49_D63_D77"
Villi_gene_common_Up_D49_D63_D77_sep$GeneList<-"Villi_gene_common_Up_D49_D63_D77_sep"

####For target six group
Enrich_all_plot<-as.data.frame(rbind(Villi_gene_common_Up_D21_D49_D63_D77[grep("Summary",Villi_gene_common_Up_D21_D49_D63_D77$GroupID),],
                                     Villi_gene_common_Up_D49_D63_D77_sep[grep("Summary",Villi_gene_common_Up_D49_D63_D77_sep$GroupID),]))

dim(Enrich_all_plot)#26 10
head(Enrich_all_plot)
Enrich_all_plot$count <- as.numeric(sapply(strsplit(Enrich_all_plot$InTerm_InList, "/"), function(x) x[1]))

Enrich_all_plot2<-Enrich_all_plot[,c("GroupID","Term","Description","count","LogP","Log(q-value)","Symbols","GeneList")];dim(Enrich_all_plot2)#340  10
colnames(Enrich_all_plot2)<-c("GroupID","Category","Description","Count","Log_pvalue","Log_qvalue","Hits","group")

Enrich_all_plot2$pathway<-paste0(Enrich_all_plot2$Category,":",Enrich_all_plot2$Description)
head(Enrich_all_plot2)
Enrich_all_plot2<-na.omit(Enrich_all_plot2)
Enrich_all_plot2$Log_pvalue<-c(-Enrich_all_plot2$Log_pvalue)
Enrich_all_plot2$Log_qvalue<-c(-Enrich_all_plot2$Log_qvalue)

Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_pvalue)),]
Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_qvalue)),]

range(Enrich_all_plot2$Log_pvalue)#  2.050510 9.342688
range(Enrich_all_plot2$Log_qvalue)#   0.000000 4.996217

range(Enrich_all_plot2$Count)# 3 30
head(Enrich_all_plot2)
table(Enrich_all_plot2$group)
# Villi_gene_common_Up_D21_D49_D63_D77 Villi_gene_common_Up_D49_D63_D77_sep 
#           6                                   20 
write.table(Enrich_all_plot2, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/metascape_result/ClassI_II_Day_common_Enrich_Metascape_result_single_group.txt",row.names=T, col.names=T) 

#2.论文画法：
table(Enrich_all_plot2$group)
Enrich_all_plot_common_one<-Enrich_all_plot2[which(Enrich_all_plot2$group %in% c("Villi_gene_common_Up_D21_D49_D63_D77")),]
Enrich_all_plot_common_two<-Enrich_all_plot2[which(Enrich_all_plot2$group == "Villi_gene_common_Up_D49_D63_D77_sep"),]

#自定义主题：
Enrich_all_plot_common_one$text_x <- rep(0.03,nrow(Enrich_all_plot_common_one)) 
Enrich_all_plot_common_one$group<-factor(Enrich_all_plot_common_one$group,levels =  c("Villi_gene_common_Up_D21_D49_D63_D77","Villi_gene_common_Up_D49_D63_D77_sep"))

Enrich_all_plot_common_one<-Enrich_all_plot_common_one[order(Enrich_all_plot_common_one$group,Enrich_all_plot_common_one$Log_pvalue,decreasing = T),]
Enrich_all_plot_common_one$pathway<-factor(Enrich_all_plot_common_one$pathway,levels = rev(c(unique(as.character(Enrich_all_plot_common_one$pathway)))))
range(Enrich_all_plot_common_one$Log_pvalue)#2.050510 4.813017
p1<- ggplot(data = Enrich_all_plot_common_one,aes(x = Log_pvalue, y = pathway)) +
  geom_bar(aes(fill = Log_pvalue), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(x = "Log_pvalue", y = "pathway", title = "Metascape enrichment barplot:Decidual_Down_Villus_Up") +
  geom_text(aes(x = text_x, label= pathway),hjust= 0)+ #hjust=0，左对齐
  scale_x_continuous(limits = c(0,5), breaks = seq(0, 5, by = 1))+  # 设置刻度范围和间隔
  theme_classic()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title = element_blank(),axis.text = element_text(size = 11),plot.title = element_blank(),legend.title = element_text(size = 13),legend.text = element_text(size = 11))
p1

#自定义主题：
Enrich_all_plot_common_two$text_x <- rep(0.03,nrow(Enrich_all_plot_common_two)) 
Enrich_all_plot_common_two<-Enrich_all_plot_common_two[order(Enrich_all_plot_common_two$Log_pvalue,decreasing = T),]
Enrich_all_plot_common_two$pathway<-factor(Enrich_all_plot_common_two$pathway,levels = rev(c(unique(as.character(Enrich_all_plot_common_two$pathway)))))
range(Enrich_all_plot_common_two$Log_pvalue)# 2.949528 9.342688

p2<- ggplot(data = head(Enrich_all_plot_common_two,n=15),aes(x = Log_pvalue, y = pathway)) +
  geom_bar(aes(fill = Log_pvalue), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "Purples",direction = 1) +
  labs(x = "Log_pvalue", y = "pathway", title = "Metascape enrichment barplot:Decidual_Down_Villus_Up") +
  geom_text(aes(x = text_x, label= pathway),hjust= 0)+ #hjust=0，左对齐
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10, by = 2))+  # 设置刻度范围和间隔
  theme_classic()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title = element_blank(),axis.text = element_text(size = 11),plot.title = element_blank(),legend.title = element_text(size = 13),legend.text = element_text(size = 11))
p2
combined_plot <- grid.arrange(p1, p2, ncol = 1, heights  = c(6.5,15))
ggsave(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/metascape_result/Day_common_Enrich_Metascape_result_single_group.pdf",combined_plot,width = 6, height =6)

#vill and decidual commom day DEGs
#construct dataframe for all target GO terms 
Decidual_Vill_D7_common_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/D7_common_Up/metascape_result.xlsx",sheet = 2))
Decidual_Vill_D77_common_Up <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/D77_common_Up/metascape_result.xlsx",sheet = 2))

Decidual_Vill_D7_common_Up$GeneList<-"Decidual_Vill_D7_common_Up"
Decidual_Vill_D77_common_Up$GeneList<-"Decidual_Vill_D77_common_Up"

####For target six group
Enrich_all_plot<-as.data.frame(rbind(Decidual_Vill_D7_common_Up[grep("Summary",Decidual_Vill_D7_common_Up$GroupID),],
                                     Decidual_Vill_D77_common_Up[grep("Summary",Decidual_Vill_D77_common_Up$GroupID),]))

dim(Enrich_all_plot)# 4 10
head(Enrich_all_plot)
Enrich_all_plot$count <- as.numeric(sapply(strsplit(Enrich_all_plot$InTerm_InList, "/"), function(x) x[1]))

Enrich_all_plot2<-Enrich_all_plot[,c("GroupID","Term","Description","count","LogP","Log(q-value)","Symbols","GeneList")];dim(Enrich_all_plot2)#12  8
colnames(Enrich_all_plot2)<-c("GroupID","Category","Description","Count","Log_pvalue","Log_qvalue","Hits","group")

Enrich_all_plot2$pathway<-paste0(Enrich_all_plot2$Category,":",Enrich_all_plot2$Description)
head(Enrich_all_plot2)
Enrich_all_plot2<-na.omit(Enrich_all_plot2)
Enrich_all_plot2$Log_pvalue<-c(-Enrich_all_plot2$Log_pvalue)
Enrich_all_plot2$Log_qvalue<-c(-Enrich_all_plot2$Log_qvalue)

Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_pvalue)),]
Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_qvalue)),]

range(Enrich_all_plot2$Log_pvalue)# 2.204272 8.790356
range(Enrich_all_plot2$Log_qvalue)# 0.000000 4.444943

range(Enrich_all_plot2$Count)# 3 7
head(Enrich_all_plot2)
table(Enrich_all_plot2$group)
# Decidual_Vill_D7_common_Up Decidual_Vill_D77_common_Up 
#            2                           2 

#2.论文画法：
table(Enrich_all_plot2$group)
Enrich_all_plot_common_one<-Enrich_all_plot2[which(Enrich_all_plot2$group %in% c("Decidual_Vill_D7_common_Up")),]
Enrich_all_plot_common_two<-Enrich_all_plot2[which(Enrich_all_plot2$group == "Decidual_Vill_D77_common_Up"),]

#自定义主题：
Enrich_all_plot_common_one$text_x <- rep(0.03,nrow(Enrich_all_plot_common_one)) 

Enrich_all_plot_common_one<-Enrich_all_plot_common_one[order(Enrich_all_plot_common_one$Log_pvalue,decreasing = T),]
Enrich_all_plot_common_one$pathway<-factor(Enrich_all_plot_common_one$pathway,levels = rev(c(unique(as.character(Enrich_all_plot_common_one$pathway)))))
range(Enrich_all_plot_common_one$Log_pvalue)#5.361353 7.103928
p1<- ggplot(data = Enrich_all_plot_common_one,aes(x = Log_pvalue, y = pathway)) +
  geom_bar(aes(fill = Log_pvalue), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "OrRd", direction = 1) +
  labs(x = "Log_pvalue", y = "pathway", title = "Metascape enrichment barplot:Decidual_Down_Villus_Up") +
  geom_text(aes(x = text_x, label= pathway),hjust= 0)+ #hjust=0，左对齐
  scale_x_continuous(limits = c(0,9), breaks = seq(0,9, by = 1))+  # 设置刻度范围和间隔
  theme_classic()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title = element_blank(),axis.text = element_text(size = 11),plot.title = element_blank(),legend.title = element_text(size = 13),legend.text = element_text(size = 11))
p1

#自定义主题：
Enrich_all_plot_common_two$text_x <- rep(0.03,nrow(Enrich_all_plot_common_two)) 
Enrich_all_plot_common_two<-Enrich_all_plot_common_two[order(Enrich_all_plot_common_two$Log_pvalue,decreasing = T),]
Enrich_all_plot_common_two$pathway<-factor(Enrich_all_plot_common_two$pathway,levels = rev(c(unique(as.character(Enrich_all_plot_common_two$pathway)))))
range(Enrich_all_plot_common_two$Log_pvalue)#  2.157623 7.150843

p2<- ggplot(data = head(Enrich_all_plot_common_two,n=10),aes(x = Log_pvalue, y = pathway)) +
  geom_bar(aes(fill = Log_pvalue), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "OrRd",direction = 1) +
  labs(x = "Log_pvalue", y = "pathway", title = "Metascape enrichment barplot:Decidual_Down_Villus_Up") +
  geom_text(aes(x = text_x, label= pathway),hjust= 0)+ #hjust=0，左对齐
  scale_x_continuous(limits = c(0,9), breaks = seq(0,9, by = 1))+  # 设置刻度范围和间隔
  theme_classic()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title = element_blank(),axis.text = element_text(size = 11),plot.title = element_blank(),legend.title = element_text(size = 13),legend.text = element_text(size = 11))
p2
combined_plot <- grid.arrange(p1, p2, ncol = 1, heights  = c(2,2))
ggsave(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/metascape_result/Day_gene_common_Enrich_Metascape_result_D7_D77_single_group.pdf",combined_plot,width = 6, height =5)
