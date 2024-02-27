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
Decidua_N15 <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/Decidua_N15/metascape_result.xlsx",sheet = 2))
Decidua_N68 <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/Decidua_N68/metascape_result.xlsx",sheet = 2))
Decidua_N74 <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/Decidua_N74/metascape_result.xlsx",sheet = 2))
Decidua_N89 <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/Decidua_N89/metascape_result.xlsx",sheet = 2))

Vill_N37 <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/Vill_N37/metascape_result.xlsx",sheet = 2))
Vill_N72 <-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/Vill_N72/metascape_result.xlsx",sheet = 2))

Decidua_N15$GeneList<-"Decidua_N15"
Decidua_N68$GeneList<-"Decidua_N68"
Decidua_N74$GeneList<-"Decidua_N74"
Decidua_N89$GeneList<-"Decidua_N89"
Vill_N37$GeneList<-"Vill_N37"
Vill_N72$GeneList<-"Vill_N72"

####For target six group
Enrich_all_plot<-as.data.frame(rbind(Decidua_N15[grep("Summary",Decidua_N15$GroupID),],
                                     Decidua_N68[grep("Summary",Decidua_N68$GroupID),],
                                     Decidua_N74[grep("Summary",Decidua_N74$GroupID),],
                                     Decidua_N89[grep("Summary",Decidua_N89$GroupID),],
                                     Vill_N37[grep("Summary",Vill_N37$GroupID),],
                                     Vill_N72[grep("Summary",Vill_N72$GroupID),]))

dim(Enrich_all_plot)#120  10
head(Enrich_all_plot)
Enrich_all_plot$count <- as.numeric(sapply(strsplit(Enrich_all_plot$InTerm_InList, "/"), function(x) x[1]))

Enrich_all_plot2<-Enrich_all_plot[,c("GroupID","Term","Description","count","LogP","Log(q-value)","Symbols","GeneList")];dim(Enrich_all_plot2)# 120   8

colnames(Enrich_all_plot2)<-c("GroupID","Category","Description","Count","Log_pvalue","Log_qvalue","Hits","group")

Enrich_all_plot2$pathway<-paste0(Enrich_all_plot2$Category,":",Enrich_all_plot2$Description)
head(Enrich_all_plot2)
Enrich_all_plot2<-na.omit(Enrich_all_plot2)
Enrich_all_plot2$Log_pvalue<-c(-Enrich_all_plot2$Log_pvalue)
Enrich_all_plot2$Log_qvalue<-c(-Enrich_all_plot2$Log_qvalue)

Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_pvalue)),]
Enrich_all_plot2[which(is.na(Enrich_all_plot2$Log_qvalue)),]

range(Enrich_all_plot2$Log_pvalue)# 2.460629 14.832786
range(Enrich_all_plot2$Log_qvalue)# 0.1353477 10.4863157

range(Enrich_all_plot2$Count)#  3 75
head(Enrich_all_plot2)

###输出其他图
Enrich_all_number<-data.frame(table(as.character(Enrich_all_plot2$Description)))
#Enrich_all_number2<-Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot3<-merge(Enrich_all_plot2,Enrich_all_number)
dim(Enrich_all_plot3)# 120  10
Enrich_all_plot3<-Enrich_all_plot3[order(Enrich_all_plot3$freq,Enrich_all_plot3$Description,Enrich_all_plot3$group,Enrich_all_plot3$Hits,decreasing = T),]
head(Enrich_all_plot3)
write.table(Enrich_all_plot3, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/Decidua_vill_SWAN_gene_N15_N68_N74_N89_N37_N72_enrichment_merge_single_group.txt",row.names=T, col.names=T) 

dim(Enrich_all_plot3[which(Enrich_all_plot3$Log_pvalue > -log10(0.05)),])# 120  10
dim(Enrich_all_plot3[which(Enrich_all_plot3$Log_qvalue > -log10(0.05)),])#94 10
table(Enrich_all_plot3$group)
#Decidua_N15 Decidua_N68 Decidua_N74 Decidua_N89    Vill_N37    Vill_N72 
#        20          20          20          20          20          20 

##提取各组前top5
Enrich_all_plot4<-Enrich_all_plot3[which(Enrich_all_plot3$GroupID %in% c("1_Summary","2_Summary","3_Summary","4_Summary","5_Summary")),]
Enrich_all_plot4$GroupID<-factor(Enrich_all_plot4$GroupID,levels = c("1_Summary","2_Summary","3_Summary","4_Summary","5_Summary"))

Enrich_all_plot4$group<-factor(Enrich_all_plot4$group,levels = rev(c("Decidua_N15","Decidua_N68","Decidua_N74","Decidua_N89","Vill_N37","Vill_N72")))
Enrich_all_plot4<-Enrich_all_plot4[order(Enrich_all_plot4$group,Enrich_all_plot4$GroupID,decreasing = F),]
Enrich_all_plot4$pathway<-factor(Enrich_all_plot4$pathway,levels = c(unique(as.character(Enrich_all_plot4$pathway))))
Enrich_all_plot4$group<-factor(Enrich_all_plot4$group,levels = c("Decidua_N15","Decidua_N68","Decidua_N74","Decidua_N89","Vill_N37","Vill_N72"))

length(unique(as.character(Enrich_all_plot4[which(Enrich_all_plot4$freq>1),]$pathway)))# 10
range(Enrich_all_plot4$Count)#  10 59
range(Enrich_all_plot4$Log_pvalue)#  3.58568 14.83279

plot_BP<-ggplot(Enrich_all_plot4, aes(y=pathway,x=group,size=Count,colour=Log_pvalue))+geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,10,20,30,40,50,60),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,3,6,9,12,15),name='-log10(pvalue)')+ 
  scale_x_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_bw()+labs(x="",y="GO terms",title="metascape enrichment")+
  theme(legend.text = element_text(size = 8, colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 8,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 12),legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/Decidua_vill_SWAN_gene_N15_N68_N74_N89_N37_N72_enrichment_merge_single_group.pdf",plot_BP,width =6, height =6)
write.table(Enrich_all_plot4, file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/DEGs_new_1012_single_group/metascape_result/Decidua_vill_SWAN_gene_N15_N68_N74_N89_N37_N72_enrichment_merge_single_group_wait.txt",row.names=T, col.names=T) 
