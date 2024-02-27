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

##for Decidua sample
Decidua_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk.txt",header = T,sep = "\t")
Decidua_cell_ratio2<-Decidua_cell_ratio[,-((ncol(Decidua_cell_ratio)-2):ncol(Decidua_cell_ratio))]
rowSums(Decidua_cell_ratio2[,-1])
analysis_data <-melt(Decidua_cell_ratio2,id=c("Mixture"),variable.name="Cell",value.name="ratio")
head(analysis_data)
colnames(analysis_data)<-c("sample","Cell","ratio")
Ratio_0<-ggplot(analysis_data,aes(sample,ratio,fill=Cell))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors)+guides(fill=guide_legend(title=NULL))
Ratio_0
ggsave(Ratio_0,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk_cell_ratio_barplot.pdf",width = 16, height =8)

##合并康复日期信息
colData_used <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/Final_fielt_Decidua_Vill_analysis_metadata_genderadd.txt",sep="\t",header =T)
rownames(colData_used)<-as.character(colData_used$sample_code)
dim(colData_used)#245  13
head(colData_used)

colData_Decidua<-colData_used[which(colData_used$sample_code %in% as.character(Decidua_cell_ratio2$Mixture)),]
head(colData_Decidua)
Decidua_cell_ratio2$sample_code <-Decidua_cell_ratio2$Mixture


analysis_Decidua<-merge(colData_Decidua,Decidua_cell_ratio2,by="sample_code")
analysis_Decidua[1:5,1:16]

class_order<-c("CTRL","D7","D21","D35","D49","D63","D77","D91","Dmore91")
analysis_Decidua$group_day1 <- factor(analysis_Decidua$group_day1,levels=class_order)
analysis_Decidua$infect_state <- factor(analysis_Decidua$infect_state,levels=c("Infect","No"))
analysis_Decidua$LMP_infect<-as.numeric(analysis_Decidua$LMP_infect)

table(analysis_Decidua$Year_month)
head(analysis_Decidua)
dim(distinct(analysis_Decidua[,c("sample","infect_state")]))# 125  2


##感染与未感染的组间比较
#因子在感染与非感染组的组间比较
analysis_collect0<-analysis_Decidua[,c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day1",colnames(Decidua_cell_ratio2[,-c(1,ncol(Decidua_cell_ratio2))]))]
head(analysis_collect0)
data_anlysis0 <- melt(analysis_collect0,variable.name="Cell",value.name = "Cell_ratio",id.vars = c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day1"))
data_anlysis0$age_group<-ifelse(data_anlysis0$Age>=35,"AMA","YMA")
head(data_anlysis0)
#data_anlysis0<-data_anlysis0[which(data_anlysis0$Cell_ratio>0),]

stat_data<-compare_means(Cell_ratio~infect_state, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]

stat_boxplot<-ggplot(data_anlysis0, aes(x=infect_state,y=Cell_ratio,fill=infect_state))+
  geom_boxplot(position=position_dodge(),width=0.5)+geom_jitter(width = 0.2,size=1,color="grey",alpha=0.8)+ 
  scale_color_manual(values=ppCor)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_compare_means(label="p.format", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= infect_state), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ Cell, scales = "free",ncol =5)
stat_boxplot
#ggsave(stat_boxplot,file=paste0("/mnt/data/chenwei/covid19/bulk_RNA_result/1.rm_long/Decidua_result/Decidua_compare_stat_box_infect_noremove_",index_name,".pdf"),width = 16, height =8)
ggsave(stat_boxplot,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk_cell_ratio_compare_stat_box_infect_state.pdf",width = 16, height =12)

################# 比较不同区间分组
head(data_anlysis0)
stat_data<-compare_means(Cell_ratio~group_day1, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]
# A tibble:  69 x 9

head(data_anlysis0)
Daygroup_stat_boxplot1 <-ggboxplot(data_anlysis0, x="group_day1", y="Cell_ratio", color="group_day1",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =5)
Daygroup_stat_boxplot1
Daygroup_stat_boxplot10<-Daygroup_stat_boxplot1+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

Daygroup_stat_boxplot2 <-ggboxplot(data_anlysis0, x="group_day1", y="Cell_ratio", color="group_day1",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =5)
Daygroup_stat_boxplot2
Daygroup_stat_boxplot20<-Daygroup_stat_boxplot2+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

ggsave(Daygroup_stat_boxplot1,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot2,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_D7_compare.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare2.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot20,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_D7_compare2.pdf",width = 25, height =22)


##绘制康复曲线
data_anlysis1<-data_anlysis0[which(data_anlysis0$infect_state == "Infect"),]
range(na.omit(data_anlysis1$LMP_infect))##-2 110


index_plot<-ggplot(data_anlysis1, aes(x = LMP_infect, y = Cell_ratio)) + #color = group,
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "Decidua_cell_ratio")+
  geom_rug(position = "jitter", size = 0.1, color = "black")+
  scale_x_continuous(breaks = seq(-5,120,10)) + facet_wrap(~ Cell, scales = "free",ncol =5)+
  scale_color_manual(values=my_morandi_colors)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot10<-index_plot+  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)+ geom_point(alpha = 0.5,aes(color = age_group,size =2))
index_plot30<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 0.5,aes(color = age_group,size =2)) +  ggpubr::stat_cor(label.x = 3,method = "spearman")
index_plot31<-index_plot+ geom_smooth(aes(color = age_group,fill=age_group), method = lm,se = TRUE, linewidth=2, alpha = 0.2, fullrange = TRUE)+ geom_point(alpha = 0.5,aes(color = age_group),size = 2) +   ggpubr::stat_cor(aes(color = age_group), label.x = 0,method = "spearman")

#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
ggsave(index_plot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk_cell_ratio_sm_trend_plot.pdf",width = 25, height =20)
ggsave(index_plot30,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk_cell_ratio_lm_single_plot.pdf",width = 25, height =20)
ggsave(index_plot31,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Decidua_bulk_cell_ratio_lm_age_split_plot.pdf",width = 25, height =20)

#my.formula <- y ~ x
#formula <- y ~ poly(x, 3, raw = TRUE)
#index_plot2<-index_plot+stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(rr.label),sep = '~~~')), formula =formula,parse =T,label.x.npc = "right", label.y.npc = "top",size = 3, vstep = 0.06) 
#ggsave(index_plot2,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/raw_Dot_line_plot_for_45_clinical_index.pdf",width = 40, height =40)

##for Vill sample
Vill_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk.txt",header = T,sep = "\t")
Vill_cell_ratio2<-Vill_cell_ratio[,-((ncol(Vill_cell_ratio)-2):ncol(Vill_cell_ratio))]
rowSums(Vill_cell_ratio2[,-1])
analysis_data <-melt(Vill_cell_ratio2,id=c("Mixture"),variable.name="Cell",value.name="ratio")
head(analysis_data)
colnames(analysis_data)<-c("sample","Cell","ratio")
Ratio_0<-ggplot(analysis_data,aes(sample,ratio,fill=Cell))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors)+guides(fill=guide_legend(title=NULL))
Ratio_0
ggsave(Ratio_0,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk_cell_ratio_barplot.pdf",width = 16, height =8)

##合并康复日期信息
colData_used <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/Final_fielt_Decidua_Vill_analysis_metadata_genderadd.txt",sep="\t",header =T)
rownames(colData_used)<-as.character(colData_used$sample_code)
dim(colData_used)#245  13
head(colData_used)

colData_Vill<-colData_used[which(colData_used$sample_code %in% as.character(Vill_cell_ratio2$Mixture)),]
head(colData_Vill)
Vill_cell_ratio2$sample_code <-Vill_cell_ratio2$Mixture


analysis_Vill<-merge(colData_Vill,Vill_cell_ratio2,by="sample_code")
analysis_Vill[1:5,1:16]

class_order<-c("CTRL","D7","D21","D35","D49","D63","D77","D91","Dmore91")
analysis_Vill$group_day1 <- factor(analysis_Vill$group_day1,levels=class_order)
analysis_Vill$infect_state <- factor(analysis_Vill$infect_state,levels=c("Infect","No"))
analysis_Vill$LMP_infect<-as.numeric(analysis_Vill$LMP_infect)

table(analysis_Vill$Year_month)
head(analysis_Vill)
dim(distinct(analysis_Vill[,c("sample","infect_state")]))# 120   2


##感染与未感染的组间比较
#因子在感染与非感染组的组间比较
analysis_collect0<-analysis_Vill[,c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day1",colnames(Vill_cell_ratio2[,-c(1,ncol(Vill_cell_ratio2))]))]
head(analysis_collect0)
data_anlysis0 <- melt(analysis_collect0,variable.name="Cell",value.name = "Cell_ratio",id.vars = c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day1"))
data_anlysis0$age_group<-ifelse(data_anlysis0$Age>=35,"AMA","YMA")
head(data_anlysis0)
#data_anlysis0<-data_anlysis0[which(data_anlysis0$Cell_ratio>0),]

stat_data<-compare_means(Cell_ratio~infect_state, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]

stat_boxplot<-ggplot(data_anlysis0, aes(x=infect_state,y=Cell_ratio,fill=infect_state))+
  geom_boxplot(position=position_dodge(),width=0.5)+geom_jitter(width = 0.2,size=1,color="grey",alpha=0.8)+ 
  scale_color_manual(values=ppCor)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_compare_means(label="p.format", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= infect_state), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ Cell, scales = "free",ncol =5)
stat_boxplot
#ggsave(stat_boxplot,file=paste0("/mnt/data/chenwei/covid19/bulk_RNA_result/1.rm_long/vill_result/vill_compare_stat_box_infect_noremove_",index_name,".pdf"),width = 16, height =8)
ggsave(stat_boxplot,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk_cell_ratio_compare_stat_box_infect_state.pdf",width = 16, height =9)

################# 比较不同区间分组
head(data_anlysis0)
stat_data<-compare_means(Cell_ratio~group_day1, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]
# A tibble:  69 x 9

head(data_anlysis0)
Daygroup_stat_boxplot1 <-ggboxplot(data_anlysis0, x="group_day1", y="Cell_ratio", color="group_day1",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =5)
Daygroup_stat_boxplot1
Daygroup_stat_boxplot10<-Daygroup_stat_boxplot1+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

Daygroup_stat_boxplot2 <-ggboxplot(data_anlysis0, x="group_day1", y="Cell_ratio", color="group_day1",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =5)
Daygroup_stat_boxplot2
Daygroup_stat_boxplot20<-Daygroup_stat_boxplot2+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

ggsave(Daygroup_stat_boxplot1,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare.pdf",width = 25, height =18)
ggsave(Daygroup_stat_boxplot2,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_D7_compare.pdf",width = 25, height =18)
ggsave(Daygroup_stat_boxplot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare2.pdf",width = 25, height =18)
ggsave(Daygroup_stat_boxplot20,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk_cell_ratio_box_static_compare_group_day_D7_compare2.pdf",width = 25, height =18)


##绘制康复曲线
data_anlysis1<-data_anlysis0[which(data_anlysis0$infect_state == "Infect"),]
range(na.omit(data_anlysis1$LMP_infect))##-2 110


index_plot<-ggplot(data_anlysis1, aes(x = LMP_infect, y = Cell_ratio)) + #color = group,
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "Vill_cell_ratio")+
  geom_rug(position = "jitter", size = 0.1, color = "black")+
  scale_x_continuous(breaks = seq(-5,120,10)) + facet_wrap(~ Cell, scales = "free",ncol =5)+
  scale_color_manual(values=my_morandi_colors)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot10<-index_plot+  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)+ geom_point(alpha = 0.5,aes(color = age_group,size =2))
index_plot30<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 0.5,aes(color = age_group,size =2)) +  ggpubr::stat_cor(label.x = 3,method = "spearman")
index_plot31<-index_plot+ geom_smooth(aes(color = age_group,fill=age_group), method = lm,se = TRUE, linewidth=2, alpha = 0.2, fullrange = TRUE)+ geom_point(alpha = 0.5,aes(color = age_group),size = 2) +   ggpubr::stat_cor(aes(color = age_group), label.x = 0,method = "spearman")

#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
ggsave(index_plot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk_cell_ratio_sm_trend_plot.pdf",width = 25, height =15)
ggsave(index_plot30,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk_cell_ratio_lm_single_plot.pdf",width = 25, height =15)
ggsave(index_plot31,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_Vill_bulk_cell_ratio_lm_age_split_plot.pdf",width = 25, height =15)

#my.formula <- y ~ x
#formula <- y ~ poly(x, 3, raw = TRUE)
#index_plot2<-index_plot+stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(rr.label),sep = '~~~')), formula =formula,parse =T,label.x.npc = "right", label.y.npc = "top",size = 3, vstep = 0.06) 
#ggsave(index_plot2,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/raw_Dot_line_plot_for_45_clinical_index.pdf",width = 40, height =40)
