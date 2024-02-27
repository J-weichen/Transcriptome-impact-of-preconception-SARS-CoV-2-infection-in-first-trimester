rm(list = ls())

library(reshape2)
library(RColorBrewer)
library(grid)
library(scales)
library(ggsci)
library(stringr)
library(ggpubr) 
library(dplyr)
require(gridExtra)
library(readxl)
library(ggrepel)
library(moRandi)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

##building metadata for samples
clinical_data<-data.frame(sample=c("A60","A61","A65","A76","A96","B80","B91","C38","C39","C44","C63","C64"),
                          Age=c(25,26,31,28,27,27,29,20,25,24,23,28),
                          Pregnant_day=c(54,42,46,47,54,61,52,70,51,44,53,47), 
                          Infect_day=c(29,33,37,39,15,47,64,42,61,69,92,89),
                          patient=c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12"))

##reading FACS index
FACS_index<-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/1.data/4.Facs_data/流式各群占比计算.xlsx", sheet =3, col_names = T, col_types = NULL, na = "", skip = 0))
head(FACS_index)
##B91的流式结果存在异常
FACS_index<-FACS_index[which(!(FACS_index$sample %in% c("B91","B91-2"))),]
data_anlysis<-merge(clinical_data,FACS_index,by="sample")
head(data_anlysis);dim(FACS_index);dim(data_anlysis)#9 37 41

##calculate the ratio for target cells
Cell_infor<-data.frame(Details=c("T_cells_in_CD45_cell","CD8_cells_in_CD3_Ts","CD4_cells_in_CD3_Ts","CD8_cells_in_live_PBMCs","CD4_cells_in_live_PBMCs",
                                 "CD8_cells_in_CD45_cell","CD4_cells_in_CD45_cell",
                                 "CD69_resident_cells_in_CD4_cell","CD69_resident_cells_in_CD8_cell",
                                 "Ki67_proliferation_cells_in_CD4_cell","Ki67_proliferation_cells_in_CD8_cell",
                                 "Naive_Tn_in_CD4_cell","Central_Memory_Tcm_in_CD4_cell","Effector_Memory_Tem_in_CD4_cell","Terminally-differentiated_Effector_Memory_Temra_in_CD4_cell",
                                 "Naive_Tn_in_CD8_cell","Central_Memory_Tcm_in_CD8_cell","Effector_Memory_Tem_in_CD8_cell","Terminally-differentiated_Effector_Memory_Temra_in_CD8_cell",
                                 "CD69_resident_in_CD4_Naive_Tn","CD69_resident_in_CD4_Central_Memory_Tcm","CD69_resident_in_CD4_Effector_Memory_Tem","CD69_resident_in_CD4_Terminally-differentiated_Effector_Memory_Temra",
                                 "CD69_resident_in_CD8_Naive_Tn","CD69_resident_in_CD8_Central_Memory_Tcm","CD69_resident_in_CD8_Effector_Memory_Tem","CD69_resident_in_CD8_Terminally-differentiated_Effector_Memory_Temra",
                                 "Ki67_proliferation_in_CD4_Naive_Tn","Ki67_proliferation_in_CD4_Central_Memory_Tcm","Ki67_proliferation_in_CD4_Effector_Memory_Tem","Ki67_proliferation_in_CD4_Terminally-differentiated_Effector_Memory_Temra",
                                 "Ki67_proliferation_in_CD8_Naive_Tn","Ki67_proliferation_in_CD8_Central_Memory_Tcm","Ki67_proliferation_in_CD8_Effector_Memory_Tem","Ki67_proliferation_in_CD8_Terminally-differentiated_Effector_Memory_Temra"),
                       Facs_index=paste0("Index_",1:35))

data_anlysis$Index_1<-round((data_anlysis$`T`/data_anlysis$P4)*100,3)
data_anlysis$Index_2<-round((data_anlysis$CD8/data_anlysis$`T`)*100,3)
data_anlysis$Index_3<-round((data_anlysis$CD4/data_anlysis$`T`)*100,3)
data_anlysis$Index_4<-round((data_anlysis$CD8/data_anlysis$P3)*100,3)
data_anlysis$Index_5<-round((data_anlysis$CD4/data_anlysis$P3)*100,3)
data_anlysis$Index_6<-round((data_anlysis$CD8/data_anlysis$P4)*100,3)
data_anlysis$Index_7<-round((data_anlysis$CD4/data_anlysis$P4)*100,3)
data_anlysis$Index_8<-round((data_anlysis$P6/data_anlysis$CD4)*100,3)
data_anlysis$Index_9<-round((data_anlysis$P8/data_anlysis$CD8)*100,3)
data_anlysis$Index_10<-round((data_anlysis$P7/data_anlysis$CD4)*100,3)
data_anlysis$Index_11<-round((data_anlysis$P9/data_anlysis$CD8)*100,3)
data_anlysis$Index_12<-round((data_anlysis$'Q2-UR'/data_anlysis$CD4)*100,3)
data_anlysis$Index_13<-round((data_anlysis$'Q2-UL'/data_anlysis$CD4)*100,3)
data_anlysis$Index_14<-round((data_anlysis$'Q2-LL'/data_anlysis$CD4)*100,3)
data_anlysis$Index_15<-round((data_anlysis$'Q2-LR'/data_anlysis$CD4)*100,3)
data_anlysis$Index_16<-round((data_anlysis$'Q1-UR'/data_anlysis$CD8)*100,3)
data_anlysis$Index_17<-round((data_anlysis$'Q1-UL'/data_anlysis$CD8)*100,3)
data_anlysis$Index_18<-round((data_anlysis$'Q1-LL'/data_anlysis$CD8)*100,3)
data_anlysis$Index_19<-round((data_anlysis$'Q1-LR'/data_anlysis$CD8)*100,3)

data_anlysis$Index_20<-round((data_anlysis$P22/data_anlysis$'Q2-UR')*100,3)
data_anlysis$Index_21<-round((data_anlysis$P23/data_anlysis$'Q2-UL')*100,3)
data_anlysis$Index_22<-round((data_anlysis$P24/data_anlysis$'Q2-LL')*100,3)
data_anlysis$Index_23<-round((data_anlysis$P25/data_anlysis$'Q2-LR')*100,3)
data_anlysis$Index_24<-round((data_anlysis$P18/data_anlysis$'Q1-UR')*100,3)
data_anlysis$Index_25<-round((data_anlysis$P19/data_anlysis$'Q1-UL')*100,3)
data_anlysis$Index_26<-round((data_anlysis$P20/data_anlysis$'Q1-LL')*100,3)
data_anlysis$Index_27<-round((data_anlysis$P21/data_anlysis$'Q1-LR')*100,3)

data_anlysis$Index_28<-round((data_anlysis$P14/data_anlysis$'Q2-UR')*100,3)
data_anlysis$Index_29<-round((data_anlysis$P15/data_anlysis$'Q2-UL')*100,3)
data_anlysis$Index_30<-round((data_anlysis$P16/data_anlysis$'Q2-LL')*100,3)
data_anlysis$Index_31<-round((data_anlysis$P17/data_anlysis$'Q2-LR')*100,3)
data_anlysis$Index_32<-round((data_anlysis$P10/data_anlysis$'Q1-UR')*100,3)
data_anlysis$Index_33<-round((data_anlysis$P11/data_anlysis$'Q1-UL')*100,3)
data_anlysis$Index_34<-round((data_anlysis$P12/data_anlysis$'Q1-LL')*100,3)
data_anlysis$Index_35<-round((data_anlysis$P13/data_anlysis$'Q1-LR')*100,3)

write.table(as.data.frame(data_anlysis), file="D:/PROJECT/新冠/manuscript/3.Table/Facs/Facs_PBMCs_T1_index_all_sample_infection_data.txt",sep="\t",quote=F, row.names=F, col.names=T) 
##数据准备
head(data_anlysis)
#data_anlysis2<-data_anlysis[which(data_anlysis$sample != "B91"),]
#analysis_used<-data_anlysis2[,-c(1,6:25)]
analysis_used<-data_anlysis[,-c(1,6:41)]

head(analysis_used)
analysis_used$Infect_group<-ifelse(analysis_used$Infect_day<50,"Group_one","Group_two")
analysis_final <-melt(analysis_used,id=c("patient","Age","Pregnant_day","Infect_day","Infect_group"),variable.name="Facs_index",value.name="ratio")
head(analysis_final)
analysis_final2<-merge(Cell_infor,analysis_final,by="Facs_index")
head(analysis_final2);dim(analysis_final2)
analysis_final2$Infect_group<-factor(analysis_final2$Infect_group,levels = c("Group_one","Group_two"))
analysis_final2$Details<-factor(analysis_final2$Details,levels = Cell_infor$Details)

##绘图
Facs_index_boxplot<-ggboxplot(analysis_final2,x="Infect_group", y="ratio", add = "jitter",size =0.5,fill="Infect_group") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  scale_fill_manual(values=c("navy","red"))+ylab("% ratio")+labs(title = "PBMCs_T1_index")+
  #  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE)+
  stat_compare_means(label="p.format", method="wilcox.test",hide.ns = F,label.x = 1.5)+
  stat_summary(aes(group= Infect_group), fun = "mean", geom = "point",shape=23,size=2,fill="white",position=position_dodge(0.8)) +
  #stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ Details, scales = "free",ncol =6)
Facs_index_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/all_sample/Facs_PBMCs_T1_index_all_sample_compare_stat_box_infection_data_split.pdf",Facs_index_boxplot,width=24, height=26)
##绘制相关性图
head(analysis_final2)
range(analysis_final2$Infect_day)# 15 92
range(analysis_final2$Pregnant_day)# 44 70

index_plot<-ggplot(analysis_final2, aes(x = Infect_day, y = ratio)) + #color = group,
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "PBMCs_T1_index")+
  scale_size(breaks =  c(40,45,50,55,60,65,70),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(10,100,10)) + facet_wrap(~ Details, scales = "free",ncol =6)+
  scale_color_manual(values=my_morandi_colors)+
  # scale_color_brewer(palette = "Dark2") + 
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot1<-index_plot+  geom_smooth(method = "auto",size=2) +  geom_point(alpha = 1,aes(color = patient,size =Pregnant_day)) 
index_plot2<-index_plot+  geom_line(size=1,color="blue")+geom_point(alpha = 0.9,aes(color = patient,size=2)) 
index_plot3<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 1,aes(color = patient,size =Pregnant_day)) +  ggpubr::stat_cor(label.x = 15,method = "spearman")
index_plot4<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(label.x = 15,method = "spearman")
index_plot3+geom_text(stat="identity",aes(label=patient), color="black", size=2)

#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
ggsave(index_plot1,file="D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/all_sample/Facs_PBMCs_T1_index_Dot_smooth_plot_for_all_sample_infection_data_split.pdf",width = 24, height =26)
ggsave(index_plot2,file="D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/all_sample/Facs_PBMCs_T1_index_Dot_line_plot_for_all_sample_infection_data_split.pdf",width = 24, height =26)
ggsave(index_plot3,file="D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/all_sample/Facs_PBMCs_T1_index_Dot_lm_plot_for_all_sample_infection_data_split.pdf",width = 24, height =26)
ggsave(index_plot4,file="D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/all_sample/Facs_PBMCs_T1_index_Dot_lm_plot_for_all_sample_infection_data_split2.pdf",width = 24, height =26)

##去除异常点
##考虑到P5中月经末期距离感染日期仅15天,如果纳入一定的记忆偏差,此时患者可能还处于急性感染期内,因此考虑去除
clinical_data
analysis_final3<-analysis_final2[which(!(analysis_final2$patient %in% c("P5"))),]

Facs_index_boxplot2<-ggboxplot(analysis_final3,x="Infect_group", y="ratio", add = "jitter",size =0.5,fill="Infect_group") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  scale_fill_manual(values=c("navy","red"))+ylab("% ratio")+labs(title = "PBMCs_T1_index")+
  #  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE)+
  stat_compare_means(label="p.format", method="wilcox.test",hide.ns = F,label.x = 1.5)+
  stat_summary(aes(group= Infect_group), fun = "mean", geom = "point",shape=23,size=2,fill="white",position=position_dodge(0.8)) +
  #stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ Details, scales = "free",ncol =6)
Facs_index_boxplot2
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/P5_remove/Facs_PBMCs_T1_index_P5_remove_sample_compare_stat_box_infection_data_split.pdf",Facs_index_boxplot2,width=24, height=26)

##绘制相关性图
range(analysis_final3$Infect_day)#37 92
range(analysis_final3$Pregnant_day)#44 70
##输出相关性统计结果
write.table(as.data.frame(analysis_final3[,c("Details","patient","Infect_day","ratio")]), file="D:/PROJECT/新冠/manuscript/3.Table/Facs/Cell_ratio_PBMCs_T1_index_for_P5_remove_sample_infection_data.txt",sep = "\t",row.names=F) 
correlation_results <- analysis_final3 %>%  group_by(Details) %>% summarize(correlation = cor(Infect_day, ratio,method = "spearman"), p_value = cor.test(Infect_day, ratio,method = "spearman")$p.value)
write.table(as.data.frame(correlation_results), file="D:/PROJECT/新冠/manuscript/3.Table/Facs/lm_spearman_PBMCs_T1_index_for_P5_remove_sample_infection_data.txt",sep = "\t",row.names=F) 
correlation_results[which(correlation_results$p_value<0.05),]

index_plot<-ggplot(analysis_final3, aes(x = Infect_day, y = ratio)) + #color = group,
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "PBMCs_T1_index")+
  scale_size(breaks =  c(40,45,50,55,60,65,70),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(35,100,10)) + facet_wrap(~ Details, scales = "free",ncol =6)+
  scale_color_manual(values=my_morandi_colors)+
  # scale_color_brewer(palette = "Dark2") + 
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot10<-index_plot+  geom_smooth(method = "auto",size=2) +  geom_point(alpha = 1,aes(color = patient,size =Pregnant_day)) 
index_plot20<-index_plot+  geom_line(size=1,color="blue")+geom_point(alpha = 0.9,aes(color = patient,size=2)) 
index_plot30<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 1,aes(color = patient,size =Pregnant_day)) +  ggpubr::stat_cor(label.x = 35,method = "spearman")
index_plot40<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(label.x = 35,method = "spearman")
index_plot30+geom_text(stat="identity",aes(label=patient), color="black", size=2)

#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
ggsave(index_plot10,file="D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/P5_remove/Facs_PBMCs_T1_index_Dot_smooth_plot_for_P5_remove_sample_infection_data_split.pdf",width = 24, height =26)
ggsave(index_plot20,file="D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/P5_remove/Facs_PBMCs_T1_index_Dot_line_plot_for_P5_remove_sample_infection_data_split.pdf",width = 24, height =26)
ggsave(index_plot30,file="D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/P5_remove/Facs_PBMCs_T1_index_Dot_lm_plot_for_P5_remove_sample_infection_data_split.pdf",width = 24, height =26)
ggsave(index_plot40,file="D:/PROJECT/新冠/manuscript/2.Figure/Facs/plot_figure/PBMCs/P5_remove/Facs_PBMCs_T1_index_Dot_lm_plot_for_P5_remove_sample_infection_data_split2.pdf",width = 24, height =26)
