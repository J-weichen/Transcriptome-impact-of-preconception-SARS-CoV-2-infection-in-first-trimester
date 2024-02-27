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
analysis_used <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/clinical_index_data_for_analysis_used_final.txt",header = T,sep = "\t")
dim(analysis_used )# 78529    15

analysis_used[which(analysis_used$E_name == "\xa6\xc3-GT"),]$E_name <-"γ-GT"
index_order<-c("WBC","RBC","HGB","HCT","MCV","MCH","MCHC","PLT","RDW_CV",
               "PCT","LYMPH_pct","LYMPH_av","Neut_pct","Neut_av","EO_pct","EO_av","BASO_pct","BASO_av","MONO_pct","MONO_av",
               "MPV","PDW","P-LCR","ALT","AST","TBIL","γ-GT","CK","CK-MB","BUN_Urea","Cr","UA","T-CHO","TG",
               "HDL-C","LDL-C","GLU","PT","A","INR","Fib","APTT","APTT_R","TT","R")
length(index_order)#45
analysis_used$E_name <- factor(analysis_used$E_name,levels=index_order)
analysis_used$value<-as.numeric(analysis_used$value)
head(analysis_used)
analysis_used$value<-as.numeric(analysis_used$value)


##读取被记录的康复信息完整的样本内容
collect_sample_info_final <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/collect_sample_infor.txt",header = T,sep = "\t")
collect_sample_info_final$ID_full<-str_pad(collect_sample_info_final$ID,width =8 ,side = c("left"),pad = "0")
head(collect_sample_info_final);dim(collect_sample_info_final)#257   6
length(unique(collect_sample_info_final$ID_full))#257

analysis_collect<-analysis_used[which(analysis_used$ID_full %in% unique(collect_sample_info_final$ID_full)),]
dim(analysis_collect)#11120    15


##去除非2023年2月后样本
distinct(analysis_collect[which(analysis_collect$OPR_Year == 2022),c("sample","ID_full")])
#      sample  ID_full
#51913   1209 08425191
#54047   1254 22141225
#3130      85 13575102

analysis_collect<-analysis_collect[which(analysis_collect$Year_month %in% c("20232","20233","20234")),]
length(unique(analysis_collect$ID_full))##249
dim(analysis_collect)# 10985        14
distinct(analysis_collect[which(analysis_collect$ID_full %in% c("08425191","22141225","13575102")),c("sample","ID_full")])
#sample  ID_full
#  1768 13575102
#  1878 22141225
#  86433   2034 08425191
dim(collect_sample_info_final)#249   6
table(collect_sample_info_final$Infect_state)
# no yes 
# 17 232 

analysis_collect2<-merge(analysis_collect,collect_sample_info_final,by="ID_full")
analysis_collect2<-analysis_collect2[which(analysis_collect2$Year_month %in% c("20232","20233","20234")),]

length(unique(analysis_collect2$E_name))##45
analysis_collect2$E_name <- factor(analysis_collect2$E_name,levels=index_order)
analysis_collect2$Infect_state <- factor(analysis_collect2$Infect_state,levels=c("yes","no"))
analysis_collect2$LMP_infect<-as.numeric(analysis_collect2$LMP_infect)

table(analysis_collect2$Year_month)
head(analysis_collect2)
dim(distinct(analysis_collect2[,c("sample","Infect_state")]))#249   2
##比较感染与非感染组
table(collect_sample_info_final$Infect_state)
# no yes 
# 17 232
sample_info3<-distinct(analysis_collect2[,c("ID_full","Age","Infect_state")])
dim(sample_info3)##249   3
stat_data0<-compare_means(Age~Infect_state, data=sample_info3)
stat_data0[which(stat_data0$p.signif != "ns"),]
#  .y.   group1 group2       p  p.adj p.format p.signif method  
#  Age   yes    no     0.00742 0.0074 0.0074   **       Wilcoxon

month_stat_boxplot12<-ggplot(sample_info3, aes(x=Infect_state,y=Age,fill=Infect_state))+
  geom_boxplot(position=position_dodge(),width=0.5)+geom_jitter(width = 0.2,color="grey",alpha=1)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= Infect_state), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))
month_stat_boxplot12
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Age_collect_sample_compare_stat_box_infect_compare.pdf",month_stat_boxplot12,width=6, height=6)

##因子在感染与非感染组的组间比较
stat_data<-compare_means(value~Infect_state, data=analysis_collect2, group.by = "E_name")
stat_data[which(stat_data$p.signif != "ns"),]
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Index_collect_sample_compare_stat_box_infect_compare.txt",row.names=T, col.names=T) 

month_stat_boxplot13<-ggplot(analysis_collect2, aes(x=Infect_state,y=value,fill=Infect_state))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= Infect_state), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =7)
month_stat_boxplot13
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Index_collect_sample_compare_stat_box_infect_compare.pdf",month_stat_boxplot13,width=30, height=30)

##收集样本中康复样本的康复日期分布以及分组
head(analysis_collect2)
sample_info4<-distinct(analysis_collect2[which(analysis_collect2$Infect_state == "yes"),c("ID_full","LMP_infect")])
head(sample_info4);dim(sample_info4)#232   2
range(na.omit(sample_info4$LMP_infect))##-31 156
sample_info4[which(is.na(sample_info4$LMP_infect)),]
#     ID_full LMP_infect
#126 20265792         NA :纵迪迪 尽管可查到检测信息，但是病历中未记录其具体末次月经信息

sample_info4<-sample_info4[!(is.na(sample_info4$LMP_infect)),]
dim(sample_info4)##231   2 ==> 具备详细信息样本为 231
###
sample_distribution<-ggplot(sample_info4,aes(x=LMP_infect))+
  geom_histogram(binwidth=14,color='black',fill='#03B0AB',cex=1)+
  stat_bin(binwidth=14, geom='text', color='black', size=4,
           aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+
  scale_x_continuous(breaks = seq(-35, 165, 14))
sample_distribution
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/distribution_of_lmp_infection_samples_infect.pdf",sample_distribution,width=6, height=5)

analysis_collect2$group_day<-ifelse(analysis_collect2$Infect_state =="no","CTRL",
                                    ifelse(analysis_collect2$LMP_infect <= 7,"D7",
                                           ifelse(analysis_collect2$LMP_infect <= 21,"D21",
                                                  ifelse(analysis_collect2$LMP_infect <= 35,"D35",
                                                         ifelse(analysis_collect2$LMP_infect <= 49,"D49",
                                                                ifelse(analysis_collect2$LMP_infect <= 63,"D63",
                                                                       ifelse(analysis_collect2$LMP_infect <= 77,"D77","Dmore77")))))))
table(analysis_collect2$group_day)
length(unique(analysis_collect2$E_name))
analysis_collect2$group_day<-factor(analysis_collect2$group_day,levels = c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77"))

analysis_collect2[!(is.na(analysis_collect2$group_day)),]$sample #2085  纵迪迪
analysis_collect3<-analysis_collect2[!(is.na(analysis_collect2$group_day)),]
dim(analysis_collect3)

stat_data<-compare_means(value~group_day, data=analysis_collect3,group.by = "E_name")
stat_data[which(stat_data$p.signif != "ns"),]
#A tibble: 55 × 9
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Clinical_index_collected_sample_compare_stat_among_groups.txt",row.names=T, col.names=T) 

Daygroup_stat_boxplot <-ggboxplot(analysis_collect3, x="group_day", y="value", color="group_day") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ E_name, scales = "free",ncol =7)
Daygroup_stat_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Clinical_index_collected_sample_compare_stat_box_CTRL_reference.pdf",Daygroup_stat_boxplot,width=35, height=35,limitsize = F)
Daygroup_stat_boxplot <-ggboxplot(analysis_collect3, x="group_day", y="value", color="group_day") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ E_name, scales = "free",ncol =7)
Daygroup_stat_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Clinical_index_collected_sample_compare_stat_box_D7_reference.pdf",Daygroup_stat_boxplot,width=35, height=35,limitsize = F)


sample_info5<-distinct(analysis_collect3[,c("ID_full","Age","LMP_infect","group_day","Infect_state")])
head(sample_info5);dim(sample_info5)#246   5
sample_info5$group_day<-factor(sample_info5$group_day,levels = c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77"))
stat_data<-compare_means(Age~group_day, data=sample_info5)
stat_data[which(stat_data$p.signif != "ns"),]
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Age_collected_sample_compare_stat_among_groups.txt",row.names=T, col.names=T) 

table(sample_info5$group_day)
Daygroup_age_boxplot <-ggboxplot(sample_info5, x="group_day", y="Age", color="group_day") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
Daygroup_age_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Age_collected_sample_diff_group_compare_stat_box.pdf",Daygroup_age_boxplot,width=6, height=6,limitsize = F)




##绘制康复曲线
analysis_collect3<-analysis_collect2[which(analysis_collect2$Infect_state == "yes"),]
analysis_collect3<-analysis_collect3[!(is.na(analysis_collect3$group_day)),]

head(analysis_collect3);dim(analysis_collect3)# 10184    21
range(na.omit(analysis_collect3$LMP_infect))##-31 156
range(na.omit(analysis_collect3$day_detected))##29 81
##去除首尾跨度过大样本
analysis_collect4<-analysis_collect3[which(analysis_collect3$LMP_infect >(-31) & analysis_collect3$LMP_infect<110),]
#analysis_collect4<-analysis_collect3[which(!(analysis_collect3$LMP_infect %in% c(-31,156))),]
range(na.omit(analysis_collect4$day_detected))##29 81

index_plot<-ggplot(analysis_collect3, aes(x = LMP_infect, y = value,size =day_detected)) + #color = group,
  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5) + 
  geom_point(alpha = 0.5,aes(color = age_group)) + 
  geom_rug(position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(-35,160,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_brewer(palette = "Dark2") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot

#my.formula <- y ~ x
formula <- y ~ poly(x, 3, raw = TRUE)
index_plot2<-index_plot+stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(rr.label),sep = '~~~')), formula =formula,parse =T,label.x.npc = "right", label.y.npc = "top",size = 3, vstep = 0.06) 
ggsave(index_plot2,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/Dot_line_plot_for_45_clinical_index.pdf",width = 40, height =40)

###去掉首尾跨度过大样本
index_plot2<-ggplot(analysis_collect4, aes(x = LMP_infect, y = value,size =day_detected)) + #color = group,
  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5) + 
  geom_point(alpha = 0.5,aes(color = age_group)) + 
  geom_rug(position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(0,110,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_brewer(palette = "Dark2") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot2

#my.formula <- y ~ x
formula <- y ~ poly(x, 3, raw = TRUE)
index_plot3<-index_plot2+stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(rr.label),sep = '~~~')), formula =formula,parse =T,label.x.npc = "right", label.y.npc = "top",size = 3, vstep = 0.06) 
ggsave(index_plot3,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/Dot_line_plot_for_45_clinical_index_limit110.pdf",width = 40, height =40)

##线性拟合
index_plot01<-ggplot(analysis_collect3, aes(x = LMP_infect, y = value,size =day_detected)) + #color = group,
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE) + 
  # geom_smooth(aes(color = age_group), method = lm,linewidth=2, alpha = 0.5,se = FALSE, fullrange = TRUE)+
  geom_point(alpha = 0.5,aes(color = age_group)) + 
  geom_rug(aes(color =age_group),position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(-35,160,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  ggpubr::stat_cor(label.x = 0,method = "spearman")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot01
ggsave(index_plot01,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/lm_line_plot_for_45_clinical_index_single_group.pdf",width = 40, height =40)

index_plot02<-ggplot(analysis_collect3, aes(x = LMP_infect, y = value,size =day_detected)) + #color = group,
  geom_smooth(aes(color = age_group), method = lm,linewidth=2, alpha = 0.5,se = FALSE, fullrange = TRUE)+
  geom_point(alpha = 0.5,aes(color = age_group)) + 
  geom_rug(aes(color =age_group),position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(-35,160,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  ggpubr::stat_cor(aes(color = age_group), label.x = 0,method = "spearman")+#pearson
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot02
ggsave(index_plot02,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/lm_line_plot_for_45_clinical_index_two_group.pdf",width = 40, height =40)

##for limit
index_plot04<-ggplot(analysis_collect4, aes(x = LMP_infect, y = value,size =day_detected)) + #color = group,
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE) + 
  # geom_smooth(aes(color = age_group), method = lm,linewidth=2, alpha = 0.5,se = FALSE, fullrange = TRUE)+
  geom_point(alpha = 0.5,aes(color = age_group)) + 
  geom_rug(aes(color =age_group),position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(0,110,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  ggpubr::stat_cor(label.x = 3,method = "spearman")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot04
ggsave(index_plot04,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/lm_line_plot_for_45_clinical_index_limit_single_group.pdf",width = 40, height =40)

index_plot4<-ggplot(analysis_collect4, aes(x = LMP_infect, y = value,size =day_detected)) + #color = group,
  geom_smooth(aes(color = age_group), method = lm,linewidth=2, alpha = 0.5,se = FALSE, fullrange = TRUE)+
  geom_point(alpha = 0.5,aes(color = age_group)) + 
  geom_rug(aes(color =age_group),position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(0,110,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  ggpubr::stat_cor(aes(color = age_group), label.x = 3,method = "spearman")+#pearson
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot4
ggsave(index_plot4,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/lm_line_plot_for_45_clinical_index_limit_two_group.pdf",width = 40, height =40)

index_plot41<-ggplot(analysis_collect4, aes(x = LMP_infect, y = value)) +
  geom_smooth(aes(color = age_group,fill=age_group), method = lm,se = TRUE, linewidth=2, alpha = 0.2, fullrange = TRUE)+
  geom_point(alpha = 0.5,aes(color = age_group),size = 3) + 
  geom_rug(aes(color =age_group),position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(0,110,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  ggpubr::stat_cor(aes(color = age_group), label.x = 3,method = "spearman")+#pearson
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot41
ggsave(index_plot41,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/lm_line_plot_for_45_clinical_index_limit_two_group_2.pdf",width = 40, height =40)


#################仅绘制单组的康复轨迹
##绘制康复曲线
##不标记妊娠天数以及年龄分组
index_plot<-ggplot(analysis_collect3, aes(x = LMP_infect, y = value)) + #color = group,
  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5) + 
  geom_point(alpha = 0.5,size = 1) + 
  geom_rug(position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(-35,160,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_brewer(palette = "Dark2") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot

#my.formula <- y ~ x
formula <- y ~ poly(x, 3, raw = TRUE)
index_plot2<-index_plot+stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(rr.label),sep = '~~~')), formula =formula,parse =T,label.x.npc = "right", label.y.npc = "top",size = 3, vstep = 0.06) 
ggsave(index_plot2,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/raw_Dot_line_plot_for_45_clinical_index.pdf",width = 40, height =40)

index_plot2<-ggplot(analysis_collect4, aes(x = LMP_infect, y = value)) + #color = group,
  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5) + 
  geom_point(alpha = 0.5,size = 1)+ 
  geom_rug(position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(0,110,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_brewer(palette = "Dark2") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot2

#my.formula <- y ~ x
formula <- y ~ poly(x, 3, raw = TRUE)
index_plot3<-index_plot2+stat_poly_eq(aes(label = paste(after_stat(eq.label),after_stat(rr.label),sep = '~~~')), formula =formula,parse =T,label.x.npc = "right", label.y.npc = "top",size = 3, vstep = 0.06) 
ggsave(index_plot3,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/raw_Dot_line_plot_for_45_clinical_index_limit110.pdf",width = 40, height =40)

##线性拟合
index_plot01<-ggplot(analysis_collect3, aes(x = LMP_infect, y = value)) + #color = group,
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE) + 
  # geom_smooth(aes(color = age_group), method = lm,linewidth=2, alpha = 0.5,se = FALSE, fullrange = TRUE)+
  geom_point(alpha = 0.5,size = 1) + 
  geom_rug(aes(color =age_group),position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(-35,160,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  ggpubr::stat_cor(label.x = 0,method = "spearman")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot01
ggsave(index_plot01,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/raw_lm_line_plot_for_45_clinical_index_single_group.pdf",width = 40, height =40)

##for limit in selected sample
index_plot04<-ggplot(analysis_collect4, aes(x = LMP_infect, y = value)) + #color = group,
  geom_smooth(method = lm, se = TRUE, fill = "grey",linewidth=2, alpha = 0.5, fullrange = TRUE) + 
  # geom_smooth(aes(color = age_group), method = lm,linewidth=2, alpha = 0.5,se = FALSE, fullrange = TRUE)+
  geom_point(alpha = 0.5,size = 1) + 
  geom_rug(aes(color =age_group),position = "jitter", size = 0.1, color = "black") + 
  xlab("infect day to last menstrual period") +#labs(title = target)+ylab(target)+
  scale_size(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85),range = c(1,6),name='pregnancy days')+
  scale_x_continuous(breaks = seq(0,110,20)) + facet_wrap(~ E_name , scales = "free",ncol =7)+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  ggpubr::stat_cor(label.x = 3,method = "spearman")+
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) 
index_plot04
ggsave(index_plot04,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/recover_Time_group/raw_lm_line_plot_for_45_clinical_index_limit_single_group.pdf",width = 40, height =40)

##输出相关性统计结果
head(analysis_collect4)
write.table(as.data.frame(analysis_collect4[,c("ID_full","sample","LMP_infect","E_name","value")]), file="D:/PROJECT/新冠/manuscript/3.Table/Figure1/data_45_clinical_index_limit_sample.txt",sep = "\t",row.names=F) 
correlation_results <- analysis_collect4 %>%  group_by(E_name) %>% summarize(correlation = cor(LMP_infect, value,method = "spearman"), p_value = cor.test(LMP_infect, value,method = "spearman")$p.value)
write.table(as.data.frame(correlation_results), file="D:/PROJECT/新冠/manuscript/3.Table/Figure1/spearman_lm_line_for_45_clinical_index_limit_single_group.txt",sep = "\t",row.names=F) 
correlation_results[which(correlation_results$p_value<0.05),]
#  E_name   correlation p_value
#1 RDW_CV        -0.153  0.0217
#2 Neut_pct      -0.144  0.0306
#3 Cr            -0.149  0.0272
#4 TT             0.155  0.0197
#5 R              0.152  0.0219

###explore until three months 20231221#######################
head(analysis_collect2)
Col_data<- distinct(analysis_collect2[,c("ID_full","sample","LMP_infect","group_day")])
dim(Col_data)
table(Col_data$group_day)
#CTRL      D7     D21     D35     D49     D63     D77 Dmore77 
# 17       8      18      52      44      40      35      34 

analysis_collect2$group_day2<-ifelse(analysis_collect2$Infect_state =="no","CTRL",
                                     ifelse(analysis_collect2$LMP_infect <= 7,"D7",
                                            ifelse(analysis_collect2$LMP_infect <= 21,"D21",
                                                   ifelse(analysis_collect2$LMP_infect <= 35,"D35",
                                                          ifelse(analysis_collect2$LMP_infect <= 49,"D49",
                                                                 ifelse(analysis_collect2$LMP_infect <= 63,"D63",
                                                                        ifelse(analysis_collect2$LMP_infect <= 77,"D77",
                                                                               ifelse(analysis_collect2$LMP_infect <= 91,"D91","Dmore91"))))))))

table(analysis_collect2$group_day2)
length(unique(analysis_collect2$E_name))
analysis_collect2$group_day2<-factor(analysis_collect2$group_day2,levels = c("CTRL","D7","D21","D35","D49","D63","D77","D91","Dmore91"))
Col_data<- distinct(analysis_collect2[,c("ID_full","sample","LMP_infect","group_day2")])
#Col_data[which(Col_data$ID_full ==  "20265792"),]
#     ID_full sample LMP_infect group_day2
#132 20265792   2085         NA       <NA>

table(Col_data$group_day)
#CTRL      D7     D21     D35     D49     D63     D77     D91 Dmore91 
#  17       8      18      52      44      40      35      21      13 
unique(analysis_collect2[!(is.na(analysis_collect2$group_day2)),]$sample) #2085  纵迪迪
analysis_collect3<-analysis_collect2[!(is.na(analysis_collect2$group_day2)),]
dim(analysis_collect3)#10940    22

stat_data<-compare_means(value~group_day2, data=analysis_collect3,group.by = "E_name")
stat_data[which(stat_data$p.signif != "ns"),]
#A tibble: 70 × 9
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Clinical_index_collected_sample_compare_stat_among_groups_three_months.txt",row.names=T, col.names=T) 

Daygroup_stat_boxplot <-ggboxplot(analysis_collect3, x="group_day2", y="value", color="group_day2") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ E_name, scales = "free",ncol =7)
Daygroup_stat_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Clinical_index_collected_sample_compare_stat_box_CTRL_reference_three_months.pdf",Daygroup_stat_boxplot,width=38, height=35,limitsize = F)
Daygroup_stat_boxplot <-ggboxplot(analysis_collect3, x="group_day2", y="value", color="group_day2") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ E_name, scales = "free",ncol =7)
Daygroup_stat_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Clinical_index_collected_sample_compare_stat_box_D7_reference_three_months.pdf",Daygroup_stat_boxplot,width=38, height=35,limitsize = F)


sample_info5<-distinct(analysis_collect3[,c("ID_full","Age","LMP_infect","group_day2","Infect_state")])
head(sample_info5);dim(sample_info5)#248   5
sample_info5$group_day2<-factor(sample_info5$group_day2,levels =  c("CTRL","D7","D21","D35","D49","D63","D77","D91","Dmore91"))
stat_data<-compare_means(Age~group_day2, data=sample_info5)
stat_data[which(stat_data$p.signif != "ns"),]
#A tibble: 6 × 8
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Age_collected_sample_compare_stat_among_groups_three_months.txt",row.names=T, col.names=T) 

table(sample_info5$group_day2)
Daygroup_age_boxplot <-ggboxplot(sample_info5, x="group_day2", y="Age", color="group_day2") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
Daygroup_age_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Infect_state_group/Age_collected_sample_diff_group_compare_stat_box_three_months.pdf",Daygroup_age_boxplot,width=6, height=6,limitsize = F)
