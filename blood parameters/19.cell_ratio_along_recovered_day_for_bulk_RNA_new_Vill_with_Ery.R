rm(list = ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggpmisc)
library(moRandi)
library(scales)
library(Laurae)
library(dplyr)
##set colour
x <- c(30,4,1,2,3,20,26,29,37,41,6,7,8,51,39,42,56,52,60, 43,58,12,50)
my_morandi_colors<-morandi_diy(my_colors = x)
show_col(my_morandi_colors)
#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

###########evaluation using psudo bulk RNA####################
True_sc_Vill_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/True_sc_Vill_cell_ratio.txt",header = T,sep = "\t")
#
flag_tag<-"manu_sig_no_1K"
CIBERSORTx_pesudo_Vill <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/seurat_result/manu_1K_no_adjsut/CIBERSORTx_pesudo_Mixture_Vill.txt",header = T,sep = "\t")

head(True_sc_Vill_cell_ratio)
head(CIBERSORTx_pesudo_Vill)

True_data <-melt(True_sc_Vill_cell_ratio,id=c("sample_code"),variable.name="Cell",value.name="ratio_true")
head(True_data)
CB_data <-melt(CIBERSORTx_pesudo_Vill,id=c("Mixture"),variable.name="Cell",value.name="ratio_CB")
head(CB_data)
CB_data$ratio_CB<-CB_data$ratio_CB*100
colnames(True_data)<-c("Mixture","Cell","ratio_true")
com_data<-merge(True_data,CB_data,by=c("Mixture","Cell"))
head(com_data)
table(com_data$Cell)
##绘制相关性图
##输出相关性统计结果
correlation_results <- com_data %>%  group_by(Cell) %>% summarize(correlation = cor(ratio_true, ratio_CB,method = "spearman"), p_value = cor.test(ratio_true, ratio_CB,method = "spearman")$p.value)
correlation_results[which(correlation_results$p_value<0.05),]

write.table(as.data.frame(com_data),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/seurat_result/manu_1K_no_adjsut/file1_True_SC_Cell_ratio_CsX_Ratio_value_",flag_tag,".txt"),sep = "\t",row.names=F) 
write.table(as.data.frame(correlation_results),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/seurat_result/manu_1K_no_adjsut/file2_True_SC_Cell_ratio_CsX_correlationship_data_",flag_tag,".txt"),sep = "\t",row.names=F) 

index_plot<-ggplot(com_data, aes(x = ratio_true, y =  ratio_CB)) +
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(method = "spearman")+
  xlab("sc_RNA_Cell ratio(%)")+ylab("CsX_Cell ratio(%)") +labs(title = "Vill_correlationship")+
  facet_wrap(~ Cell, scales = "free",ncol =4)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/seurat_result/manu_1K_no_adjsut/file3_True_SC_Cell_ratio_CsX_correlationship_",flag_tag,".pdf"),width = 12, height =10)

##线性相关性分析1）计算相关性 reference:https://zhuanlan.zhihu.com/p/543768987
#ref:https://www.coder.work/article/6652294#google_vignette
#ref:https://www.saoniuhuo.com/question/detail-2441094.html
#https://rdrr.io/cran/ggpmisc/man/stat_poly_eq.html
cor.test(com_data$ratio_true,com_data$ratio_CB,data=com_data)

###根据结果可以知道相关性系数为0.9374051，且呈正相关，P值为2.2e-16 <0.05
###回归系数的95% 置信区间为[ 0.9072704 0.9579627]
df_cor<-lm(ratio_CB~ratio_true,data=com_data)
summary(df_cor)
###回归方程为ratio_CB=0.9934360*ratio_CB+0.0003455 
###拟合优度R^2=0.8787,即拟合度很好 
##模型 p 值，而不是斜率 p 值(通常) p-value: < 2.2e-16

range(com_data$ratio_true)# 0.0000 74.7855
range(com_data$ratio_CB)#  0.00000 83.21935
formula <- y ~ x
index_plot<-ggplot(com_data, aes(x = ratio_true, y = ratio_CB)) + #color = group,
  xlab("True SC Vill Cell ratio")+ylab("CIBERSORTx Vill Cell ratio") +labs(title = "Vill correlationship(spearman)")+
  scale_x_continuous(breaks = seq(0,100,10)) +  scale_color_manual(values=ppCor)+
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 0.5,aes(color = Mixture,size =2)) +  
  ggpubr::stat_cor(method = "spearman",size =4,label.x = 10,label.y = 60)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,..p.value.label..,sep ="~~~")),size =4,colour = "black", formula = formula, label.x = 0.1,label.y = 0.7, parse = TRUE)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/seurat_result/manu_1K_no_adjsut/file4_True_SC_Cell_ratio_CsX_correlationship_no_split_",flag_tag,".pdf"),width = 12, height =10)


##########################################################
##for Vill sample
Vill_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_manu.txt",header = T,sep = "\t")
Vill_cell_ratio2<-Vill_cell_ratio[,-((ncol(Vill_cell_ratio)-3):ncol(Vill_cell_ratio))]
head(Vill_cell_ratio2)
rowSums(Vill_cell_ratio2[,-1])

###remian all cell
Vill_cell_ratio3<-as.data.frame(t(apply(Vill_cell_ratio2[,2:ncol(Vill_cell_ratio2)],1,function(x) round(prop.table(x)*100,8))))
Vill_cell_ratio3$sample_code<-as.character(Vill_cell_ratio2$Mixture)
head(Vill_cell_ratio3)

rowSums(Vill_cell_ratio3[,-ncol(Vill_cell_ratio3)])
analysis_data <-melt(Vill_cell_ratio3,id=c("sample_code"),variable.name="Cell",value.name="ratio")
head(analysis_data)
colnames(analysis_data)<-c("sample","Cell","ratio")
#analysis_data$sample <- factor(analysis_data$sample,levels=paste0("T"))

Ratio_0<-ggplot(analysis_data,aes(sample,ratio,fill=Cell))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors)+guides(fill=guide_legend(title=NULL))
Ratio_0
ggsave(Ratio_0,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_barplot_all_cell.pdf",width = 16, height =8)

##合并康复日期信息
colData_used <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/Final_fielt_Decidua_Vill_analysis_metadata_genderadd.txt",sep="\t",header =T)
rownames(colData_used)<-as.character(colData_used$sample_code)
dim(colData_used)#245  13
head(colData_used)

colData_Vill<-colData_used[which(colData_used$sample_code %in% as.character(Vill_cell_ratio2$Mixture)),]
head(colData_Vill);dim(colData_Vill)#125  13
###############################
colData_Vill$group_day2<-ifelse(colData_Vill$infect_state =="No","CTRL",
                                   ifelse(colData_Vill$LMP_infect <= 7,"D7",
                                          ifelse(colData_Vill$LMP_infect <= 21,"D21",
                                                 ifelse(colData_Vill$LMP_infect <= 35,"D35",
                                                        ifelse(colData_Vill$LMP_infect <= 49,"D49",
                                                               ifelse(colData_Vill$LMP_infect <= 63,"D63",
                                                                      ifelse(colData_Vill$LMP_infect <= 77,"D77","Dmore77")))))))
table(colData_Vill$group_day2)
class_order<-c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77")
colData_Vill$group_day2<-factor(colData_Vill$group_day2,levels = class_order)

analysis_Vill<-merge(colData_Vill,Vill_cell_ratio3,by="sample_code")
analysis_Vill[1:5,1:16]
analysis_Vill$infect_state <- factor(analysis_Vill$infect_state,levels=c("Infect","No"))
analysis_Vill$LMP_infect<-as.numeric(analysis_Vill$LMP_infect)
head(analysis_Vill)
dim(distinct(analysis_Vill[,c("sample","infect_state")]))# 120  2

##感染与未感染的组间比较
#因子在感染与非感染组的组间比较
analysis_collect0<-analysis_Vill[,c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2",colnames(Vill_cell_ratio3[,-c(ncol(Vill_cell_ratio3))]))]
head(analysis_collect0)
data_anlysis0 <- melt(analysis_collect0,variable.name="Cell",value.name = "Cell_ratio",id.vars = c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2"))
data_anlysis0$age_group<-ifelse(data_anlysis0$Age>=35,"AMA","YMA")
head(data_anlysis0)
#data_anlysis0<-data_anlysis0[which(data_anlysis0$Cell_ratio>0),]
#data_anlysis0$Cell_ratio<-data_anlysis0$Cell_ratio*100

stat_data<-compare_means(Cell_ratio~infect_state, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]

stat_boxplot<-ggplot(data_anlysis0, aes(x=infect_state,y=Cell_ratio,fill=infect_state))+
  geom_boxplot(position=position_dodge(),width=0.5)+geom_jitter(width = 0.2,size=1,color="grey",alpha=0.8)+ 
  scale_color_manual(values=ppCor)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_compare_means(label="p.format", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= infect_state), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ Cell, scales = "free",ncol =4)
stat_boxplot
#ggsave(stat_boxplot,file=paste0("/mnt/data/chenwei/covid19/bulk_RNA_result/1.rm_long/Vill_result/Vill_compare_stat_box_infect_noremove_",index_name,".pdf"),width = 16, height =8)
ggsave(stat_boxplot,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_compare_stat_box_infect_state_allcell.pdf",width = 16, height =12)

################# 比较不同区间分组
head(data_anlysis0)
stat_data<-compare_means(Cell_ratio~group_day2, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]
# A tibble:  30 × 9
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/intergroup_compare_pvalue_wilcox_test_each_cell_Vill_data_allcell.txt",sep = "\t",row.names=F) 

head(data_anlysis0)
Daygroup_stat_boxplot1 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =4)
Daygroup_stat_boxplot1
Daygroup_stat_boxplot10<-Daygroup_stat_boxplot1+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

Daygroup_stat_boxplot2 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =4)
Daygroup_stat_boxplot2
Daygroup_stat_boxplot20<-Daygroup_stat_boxplot2+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

ggsave(Daygroup_stat_boxplot1,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare_allcell.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot2,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_D7_compare_allcell.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare2_allcell.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot20,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_D7_compare2_allcell.pdf",width = 25, height =22)

##绘制康复曲线
data_anlysis1<-data_anlysis0[which(data_anlysis0$infect_state == "Infect"),]
range(na.omit(data_anlysis1$LMP_infect))##-2 110
##输出相关性统计结果
correlation_results <- data_anlysis1 %>%  group_by(Cell) %>% summarize(correlation = cor(LMP_infect, Cell_ratio,method = "spearman"), p_value = cor.test(LMP_infect, Cell_ratio,method = "spearman")$p.value)
write.table(as.data.frame(correlation_results), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/Spearman_pvalue_each_cell_Vill_data_allcell.txt",sep = "\t",row.names=F) 
correlation_results[which(correlation_results$p_value<0.05),]
#  Cell  correlation p_value
#1 STBs        0.222 0.0208 
#2 EVTs       -0.247 0.0100 
#3 STCs       -0.244 0.0108 
#4 HCs        -0.280 0.00336

index_plot<-ggplot(data_anlysis1, aes(x = LMP_infect, y = Cell_ratio)) + #color = group,
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "Vill_cell_ratio")+
  geom_rug(position = "jitter", size = 0.1, color = "black")+
  scale_x_continuous(breaks = seq(-5,120,10)) + facet_wrap(~ Cell, scales = "free",ncol =4)+
  scale_color_manual(values=my_morandi_colors)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot10<-index_plot+  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)+ geom_point(alpha = 0.5,aes(color = age_group,size =2))
index_plot30<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 0.5,aes(color = age_group,size =2)) +  ggpubr::stat_cor(label.x = 3,method = "spearman")
index_plot31<-index_plot+ geom_smooth(aes(color = age_group,fill=age_group), method = lm,se = TRUE, linewidth=2, alpha = 0.2, fullrange = TRUE)+ geom_point(alpha = 0.5,aes(color = age_group),size = 2) +   ggpubr::stat_cor(aes(color = age_group), label.x = 0,method = "spearman")

index_plot32<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+geom_smooth(method = lm, se = TRUE, fill = "grey",size=2, alpha = 0.5, fullrange = TRUE) +  ggpubr::stat_cor(label.x = 0,method = "spearman")
index_plot33<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+ geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)

#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
ggsave(index_plot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_sm_trend_plot_allcell.pdf",width = 16, height =12)
ggsave(index_plot30,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_lm_single_plot_allcell.pdf",width = 16, height =12)
ggsave(index_plot31,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_lm_age_split_plot_allcell.pdf",width = 16, height =12)
ggsave(index_plot32,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_lm_no_age_split_plot_allcell.pdf",width = 16, height =12)
ggsave(index_plot33,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_sm_no_age_split_plot_allcell.pdf",width = 16, height =12)

###remove STBs and Ery cell
head(Vill_cell_ratio2)
Vill_cell_ratio4<-Vill_cell_ratio2[,c("Mixture","CTBs_1","CTBs_2","EVTs","STCs","Endos","Epis","MyCs","HCs","TNK")]
Vill_cell_ratio3<-as.data.frame(t(apply(Vill_cell_ratio4[,2:ncol(Vill_cell_ratio4)],1,function(x) round(prop.table(x)*100,8))))
Vill_cell_ratio3$sample_code<-as.character(Vill_cell_ratio4$Mixture)
head(Vill_cell_ratio3)

rowSums(Vill_cell_ratio3[,-ncol(Vill_cell_ratio3)])
analysis_data <-melt(Vill_cell_ratio3,id=c("sample_code"),variable.name="Cell",value.name="ratio")
head(analysis_data)
colnames(analysis_data)<-c("sample","Cell","ratio")
#analysis_data$sample <- factor(analysis_data$sample,levels=paste0("T"))

Ratio_0<-ggplot(analysis_data,aes(sample,ratio,fill=Cell))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors[c(1:2,4:10)])+guides(fill=guide_legend(title=NULL))
Ratio_0
ggsave(Ratio_0,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_barplot_subset_cell.pdf",width = 16, height =8)

##合并康复日期信息
colData_used <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/Final_fielt_Decidua_Vill_analysis_metadata_genderadd.txt",sep="\t",header =T)
rownames(colData_used)<-as.character(colData_used$sample_code)
dim(colData_used)#245  13
head(colData_used)

colData_Vill<-colData_used[which(colData_used$sample_code %in% as.character(Vill_cell_ratio2$Mixture)),]
head(colData_Vill);dim(colData_Vill)#125  13
###############################
colData_Vill$group_day2<-ifelse(colData_Vill$infect_state =="No","CTRL",
                                ifelse(colData_Vill$LMP_infect <= 7,"D7",
                                       ifelse(colData_Vill$LMP_infect <= 21,"D21",
                                              ifelse(colData_Vill$LMP_infect <= 35,"D35",
                                                     ifelse(colData_Vill$LMP_infect <= 49,"D49",
                                                            ifelse(colData_Vill$LMP_infect <= 63,"D63",
                                                                   ifelse(colData_Vill$LMP_infect <= 77,"D77","Dmore77")))))))
table(colData_Vill$group_day2)
class_order<-c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77")
colData_Vill$group_day2<-factor(colData_Vill$group_day2,levels = class_order)

head(Vill_cell_ratio3)
#Vill_cell_ratio3$sample_code <-Vill_cell_ratio3$Mixture
analysis_Vill<-merge(colData_Vill,Vill_cell_ratio3,by="sample_code")
analysis_Vill[1:5,1:16]
analysis_Vill$infect_state <- factor(analysis_Vill$infect_state,levels=c("Infect","No"))
analysis_Vill$LMP_infect<-as.numeric(analysis_Vill$LMP_infect)
head(analysis_Vill)
dim(distinct(analysis_Vill[,c("sample","infect_state")]))# 120  2

##感染与未感染的组间比较
#因子在感染与非感染组的组间比较
analysis_collect0<-analysis_Vill[,c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2",colnames(Vill_cell_ratio3[,-ncol(Vill_cell_ratio3)]))]
head(analysis_collect0)
data_anlysis0 <- melt(analysis_collect0,variable.name="Cell",value.name = "Cell_ratio",id.vars = c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2"))
data_anlysis0$age_group<-ifelse(data_anlysis0$Age>=35,"AMA","YMA")
head(data_anlysis0)
#data_anlysis0<-data_anlysis0[which(data_anlysis0$Cell_ratio>0),]
#data_anlysis0$Cell_ratio<-data_anlysis0$Cell_ratio*100

stat_data<-compare_means(Cell_ratio~infect_state, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]

stat_boxplot<-ggplot(data_anlysis0, aes(x=infect_state,y=Cell_ratio,fill=infect_state))+
  geom_boxplot(position=position_dodge(),width=0.5)+geom_jitter(width = 0.2,size=1,color="grey",alpha=0.8)+ 
  scale_color_manual(values=ppCor)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_compare_means(label="p.format", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= infect_state), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ Cell, scales = "free",ncol =3)
stat_boxplot
#ggsave(stat_boxplot,file=paste0("/mnt/data/chenwei/covid19/bulk_RNA_result/1.rm_long/Vill_result/Vill_compare_stat_box_infect_noremove_",index_name,".pdf"),width = 16, height =8)
ggsave(stat_boxplot,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_compare_stat_box_infect_state_subset_cell.pdf",width = 11, height =11)

################# 比较不同区间分组
head(data_anlysis0)
stat_data<-compare_means(Cell_ratio~group_day2, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]
# A tibble:  11 × 9
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/intergroup_compare_pvalue_wilcox_test_each_cell_Vill_data_subset_cell.txt",sep = "\t",row.names=F) 

head(data_anlysis0)
Daygroup_stat_boxplot1 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =3)
Daygroup_stat_boxplot1
Daygroup_stat_boxplot10<-Daygroup_stat_boxplot1+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

Daygroup_stat_boxplot2 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =3)
Daygroup_stat_boxplot2
Daygroup_stat_boxplot20<-Daygroup_stat_boxplot2+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

ggsave(Daygroup_stat_boxplot1,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare_subset_cell.pdf",width = 15, height =12)
ggsave(Daygroup_stat_boxplot2,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_D7_compare_subset_cell.pdf",width = 15, height =12)
ggsave(Daygroup_stat_boxplot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare2_subset_cell.pdf",width = 15, height =12)
ggsave(Daygroup_stat_boxplot20,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_D7_compare2_subset_cell.pdf",width = 15, height =12)

##绘制康复曲线
data_anlysis1<-data_anlysis0[which(data_anlysis0$infect_state == "Infect"),]
range(na.omit(data_anlysis1$LMP_infect))##-2 110
##输出相关性统计结果
correlation_results <- data_anlysis1 %>%  group_by(Cell) %>% summarize(correlation = cor(LMP_infect, Cell_ratio,method = "spearman"), p_value = cor.test(LMP_infect, Cell_ratio,method = "spearman")$p.value)
write.table(as.data.frame(correlation_results), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/Spearman_pvalue_each_cell_Vill_data_subset_cell.txt",sep = "\t",row.names=F) 
correlation_results[which(correlation_results$p_value<0.05),]
#  Cell  correlation p_value
#1 CTBs_2       0.300 0.00158

index_plot<-ggplot(data_anlysis1, aes(x = LMP_infect, y = Cell_ratio)) + #color = group,
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "Vill_cell_ratio")+
  geom_rug(position = "jitter", size = 0.1, color = "black")+
  scale_x_continuous(breaks = seq(-5,120,10)) + facet_wrap(~ Cell, scales = "free",ncol =3)+
  scale_color_manual(values=my_morandi_colors)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot10<-index_plot+  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)+ geom_point(alpha = 0.5,aes(color = age_group,size =2))
index_plot30<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 0.5,aes(color = age_group,size =2)) +  ggpubr::stat_cor(label.x = 3,method = "spearman")
index_plot31<-index_plot+ geom_smooth(aes(color = age_group,fill=age_group), method = lm,se = TRUE, linewidth=2, alpha = 0.2, fullrange = TRUE)+ geom_point(alpha = 0.5,aes(color = age_group),size = 2) +   ggpubr::stat_cor(aes(color = age_group), label.x = 0,method = "spearman")
index_plot32<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+geom_smooth(method = lm, se = TRUE, fill = "grey",size=2, alpha = 0.5, fullrange = TRUE) +  ggpubr::stat_cor(label.x = 0,method = "spearman")
index_plot33<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+ geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)

#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
ggsave(index_plot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_sm_trend_plot_subset_cell.pdf",width = 12, height =12)
ggsave(index_plot30,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_lm_single_plot_subset_cell.pdf",width = 12, height =12)
ggsave(index_plot31,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_lm_age_split_plot_subset_cell.pdf",width = 12, height =12)
ggsave(index_plot32,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_lm_no_age_split_plot_subset_cell.pdf",width = 12, height =12)
ggsave(index_plot33,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result_vill_with_Ery_final/bulk_result/manu_1k_no_adjsut/CIBERSORTx_Vill_bulk_cell_ratio_sm_no_age_split_plot_subset_cell.pdf",width = 12, height =12)
