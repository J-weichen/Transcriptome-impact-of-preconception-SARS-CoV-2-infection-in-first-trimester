library(reshape2)
library(RColorBrewer)
library(grid)
library(scales)
library(ggsci)
library(stringr)
library(ggpubr) 
library(dplyr)
#library(export)
require(gridExtra)
library(readxl)
library(ggrepel)
library(pheatmap)

#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

RNA_meta <- read.table(file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/Final_fielt_Decidua_Vill_analysis_metadata_genderadd.txt",sep="\t",header = T)
head(RNA_meta);dim(RNA_meta)#245  13
range(na.omit(RNA_meta$LMP_infect))# -2 110
RNA_meta$age_group<-ifelse(RNA_meta$Age>=35,"AMA","YMA")
RNA_meta$GS_week <-RNA_meta$LMP_operate/7

###############################
RNA_meta$group_day<-ifelse(RNA_meta$infect_state =="No","CTRL",
                                ifelse(RNA_meta$LMP_infect <= 7,"D7",
                                       ifelse(RNA_meta$LMP_infect <= 21,"D21",
                                              ifelse(RNA_meta$LMP_infect <= 35,"D35",
                                                     ifelse(RNA_meta$LMP_infect <= 49,"D49",
                                                            ifelse(RNA_meta$LMP_infect <= 63,"D63",
                                                                   ifelse(RNA_meta$LMP_infect <= 77,"D77","Dmore77")))))))
class_order<-c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77")
RNA_meta$group_day<-factor(RNA_meta$group_day,levels = class_order)
dim(RNA_meta)
table(RNA_meta$group_day)


RNA_meta_sample<-distinct(RNA_meta[,c("family","Age","LMP_operate","GS_week","group_day","infect_state")])
head(RNA_meta_sample);dim(RNA_meta_sample)#132   6
RNA_meta_sample$group_day<-factor(RNA_meta_sample$group_day,levels = c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77"))

##for age
compare_means(Age~group_day, data=RNA_meta_sample,method = "t.test")
stat_data<-compare_means(Age~group_day, data=RNA_meta_sample,method = "wilcox.test")
stat_data[which(stat_data$p.signif != "ns"),]
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Age_RNA_132_sample_compare_stat_among_groups.txt",row.names=T, col.names=T) 

table(RNA_meta_sample$group_day)
Daygroup_age_boxplot <-ggboxplot(RNA_meta_sample, x="group_day", y="Age", color="group_day",add = "jitter") +  scale_color_manual(values=ppCor)+ #, palette = "jco"
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_compare_means(label.y = 45)+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
Daygroup_age_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Age_RNA_132_sample_diff_group_compare_stat_box_CTRL_ref.pdf",Daygroup_age_boxplot,width=6, height=4,limitsize = F)

##for LMP_operate
compare_means(LMP_operate ~ group_day, data=RNA_meta_sample,method = "t.test")
stat_data<-compare_means(LMP_operate~group_day, data=RNA_meta_sample,method = "wilcox.test")
stat_data[which(stat_data$p.signif != "ns"),]
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/pregancy_days_RNA_132_sample_compare_stat_among_groups.txt",row.names=T, col.names=T) 

table(RNA_meta_sample$group_day)
Daygroup_pregancy_day_boxplot <-ggboxplot(RNA_meta_sample, x="group_day", y="LMP_operate", color="group_day",add = "jitter") +  scale_color_manual(values=ppCor)+ #, palette = "jco"
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_compare_means(label.y = 85)+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
Daygroup_pregancy_day_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/pregancy_days_RNA_132_sample_diff_group_compare_stat_box_CTRL_ref.pdf",Daygroup_pregancy_day_boxplot,width=6, height=4,limitsize = F)

##for LMP_operate
compare_means(GS_week ~ group_day, data=RNA_meta_sample,method = "t.test")
stat_data<-compare_means(GS_week~group_day, data=RNA_meta_sample,method = "wilcox.test")
stat_data[which(stat_data$p.signif != "ns"),]
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/GS_week_RNA_132_sample_compare_stat_among_groups.txt",row.names=T, col.names=T) 

table(RNA_meta_sample$group_day)
Daygroup_GS_week_boxplot <-ggboxplot(RNA_meta_sample, x="group_day", y="GS_week", color="group_day",add = "jitter") +  scale_color_manual(values=ppCor)+ #, palette = "jco"
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_compare_means(label.y = 12)+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
Daygroup_GS_week_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/GS_week_RNA_132_sample_diff_group_compare_stat_box_CTRL_ref.pdf",Daygroup_GS_week_boxplot,width=6, height=4,limitsize = F)


############meta data output#####
report_meta <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.report_file/merge_read_report.txt",sep="\t",header = T)
head(report_meta);dim(report_meta)#274   3
report_meta2<-merge(report_meta,RNA_meta,by="sample_code")
head(report_meta2);dim(report_meta2)#245  15
report_meta3<-report_meta2[,c("sample_code","Raw_read_pair","Mapped_reads","Age","LMP_infect","LMP_operate","infect_state","tissue","gender","age_group","GS_week","group_day" )]
genenumber_meta <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.report_file/Final_fielt_Decidua_Vill_raw_count_data_gene_number.txt",sep="\t",header = T)
head(genenumber_meta);dim(genenumber_meta)
genenumber_meta$sample_code<-rownames(genenumber_meta)
report_meta4<-merge(report_meta3,genenumber_meta,by="sample_code")
head(report_meta4);dim(report_meta4)#245  15
write.table(as.data.frame(report_meta4), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.report_file/Final_fielt_Decidua_Vill_analysis_metadata_output.txt",row.names=T, col.names=T) 
