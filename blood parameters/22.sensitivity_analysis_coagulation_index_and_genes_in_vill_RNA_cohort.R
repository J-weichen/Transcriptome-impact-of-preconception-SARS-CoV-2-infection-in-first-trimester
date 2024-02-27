rm(list = ls())
library(MatchIt)

library(reshape2)
library(RColorBrewer)
library(grid)
library(scales)
library(ggsci)
library(stringr)
library(ggpubr) 
library(ggstatsplot)
library(dplyr)
#library(export)
require(gridExtra)
library(readxl)
library(ggrepel)
library(VennDiagram)

library(circlize)
library(ggcorrplot)
library(scales)
library(ggpmisc)
library(cli)
#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor_all)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]



#################################################################################
##正式的绘制
#################################################################################
analysis_index_wide <- read.table(file= "D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/analysis_index_wide.txt",sep="\t",header = T)
analysis_index_wide$ID_full<-str_pad(analysis_index_wide$ID_full,width =8 ,side = c("left"),pad = "0")
head(analysis_index_wide);dim(analysis_index_wide)# 132  48
###选择凝血相关指标
Coagulation_factor<-c( "PT","APTT","APTT_R","A","INR","Fib","TT","R")
analysis_index_wide<-analysis_index_wide[,c("ID_full", "Vill_code", "Decidua_code",Coagulation_factor)]
head(analysis_index_wide)

###读取最终用于计算的转录组信息
RNA_meta <- read.table(file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/Final_fielt_Decidua_Vill_analysis_metadata.txt",sep="\t",header = T)
head(RNA_meta);dim(RNA_meta)#245   8

code_data <-melt(analysis_index_wide[,1:3],id=c("ID_full"),variable.name="type",value.name="sample_code")
code_data2 <-na.omit(code_data[,c(1,3)])
RNA_meta2<-merge(code_data2,RNA_meta,by="sample_code")
RNA_meta3<-distinct(RNA_meta2[,c(2,4:6)])
head(RNA_meta3);dim(RNA_meta3)#132   4

analysis_index_wide2<-merge(analysis_index_wide,RNA_meta3,by="ID_full")
head(analysis_index_wide2);dim(analysis_index_wide2)#132   51

#####输入表达矩阵
count_table_all<-read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidual_Vill_final_245_sample_gene_expression_matrix_log2_normalized_counts.txt",sep="\t", row.names=1,header =T)
count_table_all[1:4,1:4]

#####输入凝血因子编码基因
##读取目标基因
clinical_index<-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/3.Table/Coagulation/1.Blood_coagulation_gene.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_index)
table(clinical_index$Source)
# GO_term_GO_0007596 HALLMARK_COAGULATION 
#     42                  138 
index_maker0<-unique(clinical_index$Gene_symbols)
index_maker1<-index_maker0[which(index_maker0 %in% colnames(count_table_all))]
length(index_maker0);length(index_maker1)## 159 ##158

##读取分泌蛋白列表
secretome_index1 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_SPOCTOPUS.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))
secretome_index2 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_SignalP.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))
secretome_index3 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_Secreted.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))
secretome_index4 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_Phobius.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))
secretome_index5 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_Predicted.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))

index_maker0<-unique(c(secretome_index1$Gene,secretome_index2$Gene,secretome_index3$Gene,secretome_index4$Gene,secretome_index5$Gene))
index_maker2<-index_maker0[which(index_maker0 %in% colnames(count_table_all))]
length(index_maker0);length(index_maker2)## 4741 ##4143

index_maker3<-Reduce(intersect, list(index_maker1,index_maker2))
length(index_maker3)#107

index_maker1<-c(index_maker1,"IL6",'CXCL8',"CFD")
######for all coagulation_gene detected############
count_table_all[1:4,1:4]
Vill_count <-count_table_all[as.character(na.omit(analysis_index_wide2$Vill_code)),index_maker1]
Vill_count<-Vill_count[,which(colSums(Vill_count)>0)]
dim(Vill_count)#120 158
Vill_count[1:4,1:4]

###############合并凝血相关基因与外周凝血指标数据#############
Vill_index<-analysis_index_wide2[which(!is.na(analysis_index_wide2$Vill_code)),-c(1,3)]
rownames(Vill_index)<-Vill_index$Vill_code;Vill_index<-Vill_index[,-1]
Vill_marix<-merge(Vill_index,Vill_count,by=0)

##线性拟合
##for APTT APTT_R CFD FURIN
index<-c("APTT","APTT_R","CFD","FURIN","LMP_infect")
#index<-c("PT","INR","A","TT","R","Fib", "APTT","APTT_R","CFD","FURIN","LMP_infect")
index<-c("APTT","APTT_R","IL6",'CXCL8',"LMP_infect")

head(Vill_marix)
data_test<-Vill_marix[,index]
data_analysis <-melt(data_test,id=c("LMP_infect"),variable.name="index",value.name="value")
data_analysis$index<-factor(data_analysis$index,levels = index)

#data_analysis<-na.omit(data_analysis)
head(data_analysis)
#range(data_analysis$LMP_infect)#-2 110
index_plot01<-ggplot(data_analysis, aes(x = LMP_infect, y = value)) + #color = group,
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE) + 
  # geom_smooth(aes(color = age_group), method = lm,linewidth=2, alpha = 0.5,se = FALSE, fullrange = TRUE)+
  geom_point(alpha = 0.5,size = 1) + 
  geom_rug(aes(color =age_group),position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_x_continuous(breaks = seq(-5,115,10)) + facet_wrap(~ index , scales = "free",ncol =2)+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  ggpubr::stat_cor(label.x = 0,method = "spearman")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot01
#ggsave(index_plot01,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/raw_lm_line_plot_for_45_clinical_index_single_group.pdf",width = 40, height =40)


index_plot1 <-ggplot(data_test, aes(y = APTT, x = CFD)) + theme_bw()+
  geom_point(alpha = 0.5,size = 2) + ggpubr::stat_cor(method = "spearman")+
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE) 
index_plot1
index_plot2<-ggplot(data_test, aes(y = APTT, x = FURIN)) + theme_bw()+
  geom_point(alpha = 0.5,size = 2) + ggpubr::stat_cor(method = "spearman")+
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE,color="red") 
index_plot2
combined_plot <- grid.arrange(index_plot1, index_plot2, ncol = 2)
#ggsave(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidua_Vill_inter/metascape_result/Day_gene_common_Enrich_Metascape_result_D7_D77_single_group.pdf",combined_plot,width = 6, height =5)

data_analysis$group_day<-ifelse(is.na(data_analysis$LMP_infect),"CTRL",
                                    ifelse(data_analysis$LMP_infect <= 7,"D7",
                                           ifelse(data_analysis$LMP_infect <= 21,"D21",
                                                  ifelse(data_analysis$LMP_infect <= 35,"D35",
                                                         ifelse(data_analysis$LMP_infect <= 49,"D49",
                                                                ifelse(data_analysis$LMP_infect <= 63,"D63",
                                                                       ifelse(data_analysis$LMP_infect <= 77,"D77","Dmore77")))))))
table(data_analysis$group_day)
length(unique(data_analysis$E_name))
data_analysis$group_day<-factor(data_analysis$group_day,levels = c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77"))
head(data_analysis)

stat_data<-compare_means(value~group_day, data=data_analysis,group.by = "index")
stat_data[which(stat_data$p.signif != "ns"),]
#A tibble: 70 × 9
#write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Clinical_index_collected_sample_compare_stat_among_groups_three_months.txt",row.names=T, col.names=T) 

stat_boxplot <-ggboxplot(data_analysis, x="group_day", y="value", color="group_day",add = "jitter") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ index, scales = "free",ncol =2)
stat_boxplot1<-stat_boxplot+  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)
stat_boxplot2<-stat_boxplot+  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)
stat_boxplot1
#ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Clinical_index_collected_sample_compare_stat_box_CTRL_reference_three_months.pdf",Daygroup_stat_boxplot,width=38, height=35,limitsize = F)

data_analysis$group_day2<-ifelse(data_analysis$group_day== "CTRL","CTRL","Infect")
data_analysis$group_day2<-factor(data_analysis$group_day2,levels = c("CTRL","Infect"))
stat_boxplot3 <-ggboxplot(data_analysis, x="group_day2", y="value", color="group_day2",add = "jitter") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ index, scales = "free",ncol =4)
stat_boxplot3

######for CFD output##
data_test<-Vill_marix[,index]

index_plot1 <-ggplot(data_test, aes(y = APTT, x = CFD)) + theme_bw()+
  geom_point(alpha = 0.5,size = 2) + ggpubr::stat_cor(method = "spearman")+
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE)+
  theme(axis.text= element_text(size = 10,colour = 'black'))
index_plot1
index_plot2 <-ggplot(data_test, aes(y = APTT_R, x = CFD)) + theme_bw()+
  geom_point(alpha = 0.5,size = 2) + ggpubr::stat_cor(method = "spearman")+
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE)+
  theme(axis.text= element_text(size = 10,colour = 'black'))
index_plot2

data_final<-data_analysis[which(data_analysis$index == "CFD"),]
##intergroup comparison
CFD_boxplot <-ggboxplot(data_final, x="group_day", y="value", color="group_day",add = "jitter") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)
CFD_boxplot

#ggsave(index_plot01,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/raw_lm_line_plot_for_45_clinical_index_single_group.pdf",width = 40, height =40)

index_plot01<-ggplot(data_final, aes(x = LMP_infect, y = value)) + 
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE) + 
  geom_point(alpha = 0.5,size = 2) + 
  xlab("infect day to last menstrual period") + theme_bw()+
  scale_x_continuous(breaks = seq(-5,110,10)) +
  ggpubr::stat_cor(label.x = 0,method = "spearman")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot01

combined_plot0 <- grid.arrange(index_plot1, index_plot2,CFD_boxplot,index_plot01, ncol = 4)
ggsave(combined_plot0,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/figure3_CFD_plot_final.pdf",width = 40, height =8)
