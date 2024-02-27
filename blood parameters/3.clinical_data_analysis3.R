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

analysis_used <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/clinical_index_data_for_analysis_used_final.txt",header = T,sep = "\t")
dim(analysis_used)# 78529    15
index_order<-c("WBC","RBC","HGB","HCT","MCV","MCH","MCHC","PLT","RDW_CV",
               "PCT","LYMPH_pct","LYMPH_av","Neut_pct","Neut_av","EO_pct","EO_av","BASO_pct","BASO_av","MONO_pct","MONO_av",
               "MPV","PDW","P-LCR","ALT","AST","TBIL","γ-GT","CK","CK-MB","BUN_Urea","Cr","UA","T-CHO","TG",
               "HDL-C","LDL-C","GLU","PT","A","INR","Fib","APTT","APTT_R","TT","R")
length(index_order)#45
analysis_used$E_name <- factor(analysis_used$E_name,levels=index_order)
analysis_used$value<-as.numeric(analysis_used$value)
################################

month_select<-c("20221","20231","20222","20232","20223","20233","20224","20234")
analysis_select<-analysis_used[which(analysis_used$Year_month %in%month_select),]
head(analysis_select)
analysis_select$OPR_month <- factor(analysis_select$OPR_month,levels=1:4)
analysis_select$OPR_Year <- factor(analysis_select$OPR_Year,levels=c(2022,2023))

analysis_select$E_name_mounth<-paste0(analysis_select$E_name,":",analysis_select$OPR_month)
head(analysis_select)

stat_data<-compare_means(value~OPR_Year,data=analysis_select,group.by = "E_name")
stat_data[which(stat_data$p.signif != "ns"),]
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Clinical_target_index_all_used_sample_compare_four_months_year_compare_group_by_year.txt",row.names=T, col.names=T) 

month_stat_boxplot08<-ggplot(analysis_select, aes(x=OPR_Year,y=value,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =7)
month_stat_boxplot08
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Clinical_target_index_all_used_sample_compare_stat_box_four_months_year_compare_group_by_year.pdf",month_stat_boxplot08,width=30, height=30)

stat_data_sig<-stat_data[which(stat_data$p.signif != "ns"),]
stat_data_sig$E_name <- as.character(stat_data_sig$E_name)
analysis_sig<-analysis_select[which(analysis_select$E_name %in% unique(stat_data_sig$E_name)),]
length(unique(stat_data_sig$E_name))##16

month_stat_boxplot18<-ggplot(analysis_sig, aes(x=OPR_Year,y=value,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =4)
month_stat_boxplot18
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Clinical_target_index_significant_used_sample_compare_stat_box_four_months_year_compare_group_by_year.pdf",month_stat_boxplot18,width=20, height=20)

##age comparison
stat_data<-compare_means(Age~OPR_Year,data=analysis_select,group.by = "E_name")
stat_data[which(stat_data$p.signif != "ns"),]
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Age_all_used_sample_compare_stat_box_four_months_year_compare_group_by_year.txt",row.names=T, col.names=T) 

month_stat_boxplot09<-ggplot(analysis_select, aes(x=OPR_Year,y=Age,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =7)
month_stat_boxplot09
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Age_all_used_sample_compare_stat_box_four_months_year_compare_group_by_year.pdf",month_stat_boxplot09,width=30, height=30)

################################
#提取2022年度前四个月信息
month_select1<-c("20221","20222","20223","20224")
analysis_select_2022<-analysis_used[which(analysis_used$Year_month %in%month_select1),]
head(analysis_select_2022)

##读取被记录的康复信息完整的样本内容
collect_sample_infor <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/collect_sample_infor.txt",header = T,sep = "\t")
collect_sample_infor$ID_full<-str_pad(collect_sample_infor$ID,width =8 ,side = c("left"),pad = "0")
head(collect_sample_infor);dim(collect_sample_infor)#249   6
length(unique(collect_sample_infor$ID_full))#249

analysis_collect<-analysis_used[which(analysis_used$ID_full %in% unique(collect_sample_infor$ID_full)),]
table(analysis_collect$Year_month)
analysis_collect<-analysis_collect[which(analysis_collect$Year_month %in% c("20232","20233","20234")),]
dim(analysis_collect)# 10985    15

collect_sample_info2<-collect_sample_infor[which(collect_sample_infor$ID_full %in% unique(analysis_collect$ID_full)),]
head(collect_sample_info2)
head(analysis_collect)
dim(collect_sample_info2)#249   6
table(collect_sample_info2$Infect_state)
# no yes 
# 17 232 

analysis_collect2<-merge(analysis_collect,collect_sample_info2,by="ID_full")
analysis_collect2<-analysis_collect2[which(analysis_collect2$Year_month %in% c("20232","20233","20234")),]

length(unique(analysis_collect2$E_name))##45
analysis_collect2$E_name <- factor(analysis_collect2$E_name,levels=index_order)
analysis_collect2$Infect_state <- factor(analysis_collect2$Infect_state,levels=c("yes","no"))
analysis_collect2$LMP_infect<-as.numeric(analysis_collect2$LMP_infect)

##提取2023年收集样本中明确感染信息的数据
head(analysis_collect2)
analysis_collect3<-analysis_collect2[which(analysis_collect2$Infect_state =="yes"),]
analysis_select_2023<-analysis_used[which(analysis_used$sample %in% unique(analysis_collect3$sample)),]
head(analysis_select_2023)

##合并两个年度信息数据
analysis_select_2022_23<-rbind(analysis_select_2022,analysis_select_2023)
head(analysis_select_2022_23)
analysis_select_2022_23$GS_week<-floor(analysis_select_2022_23$day_detected /7)
analysis_select_2022_23$GS_week2<-floor(analysis_select_2022_23$pregnancy_day_final/7)
##
analysis_select_2022_23$OPR_Year <- factor(analysis_select_2022_23$OPR_Year,levels=c("2022","2023"))
table(analysis_select_2022_23$OPR_Year)
head(analysis_select_2022_23)

write.table(as.data.frame(analysis_select_2022_23), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/Infect_or_not_clinical_index_4month_all.txt",sep = "\t",row.names=F) 

####数据查看
analysis_select_2022_23 <- read.table(file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/Infect_or_not_clinical_index_4month_all.txt",header = T,sep = "\t")

##年龄在感染与非感染组的组间比较
Age_data<-distinct(analysis_select_2022_23[,c("sample","Age","OPR_Year")])
stat_data<-compare_means(Age~OPR_Year, data=Age_data)
stat_data
#.y.   group1 group2     p p.adj p.format p.signif method  
#Age   2022   2023   0.983  0.98 0.98     ns       Wilcoxon

##检测时候的孕龄在感染与非感染组的组间比较
stat_data<-compare_means(day_detected~OPR_Year, data=analysis_select_2022_23, group.by = "E_name")
stat_data[which(stat_data$p.signif != "ns"),]
#37 x 9
##检测时候的孕周在感染与非感染组的组间比较
stat_data<-compare_means(GS_week~OPR_Year, data=analysis_select_2022_23, group.by = "E_name")
stat_data[which(stat_data$p.signif != "ns"),]
##37 x 9
##手术时候的孕天/孕周在感染与非感染组的组间比较

OPR_data<-distinct(analysis_select_2022_23[,c("sample","pregnancy_day_final","GS_week2","OPR_Year")])
dim(OPR_data)# 629   4
stat_data<-compare_means(pregnancy_day_final~OPR_Year, data=OPR_data)
stat_data[which(stat_data$p.signif != "ns"),]
##.y.                 group1 group2             p        p.adj p.format p.signif method 
#1 pregnancy_day_final 2022   2023   0.00000000691 0.0000000069 6.9e-09  ****     Wilcoxon
stat_data<-compare_means(GS_week2~OPR_Year, data=OPR_data)
stat_data[which(stat_data$p.signif != "ns"),]
#1 GS_week2 2022   2023   0.000000201 0.0000002 2e-07    ****     Wilcoxon

##因子在感染与非感染组的组间比较
stat_data<-compare_means(value~OPR_Year, data=analysis_select_2022_23, group.by = "E_name")
stat_data[which(stat_data$p.signif != "ns"),]
##13 x 9

###进行两组间样本匹配分析：重点参考：https://www.cnblogs.com/ayueme/articles/16844697.html
##首先可以看一下原始数据的基线资料表，用的是tableone这个包，它能计算SMD（后面会介绍这个SMD的作用）
library(tableone)

evaluat_data<-distinct(analysis_select_2022_23[,c("sample","Age","pregnancy_day_final","GS_week2","OPR_Year")])## 619   2
head(evaluat_data);dim(evaluat_data)#628   5
table2 <- CreateTableOne(vars = c('Age', 'pregnancy_day_final'),data = evaluat_data, strata = 'OPR_Year', smd=TRUE)
                        # factorVars = c('x.Gender', 'CVD'),
table2 <- print(table2,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table2

##              Stratified by OPR_Year
##                                level 2022           2023           p        test SMD    
## n                               ""    "397"          "232"          ""       ""   ""     
##Age (mean (SD))                 ""    "32.00 (6.11)" "31.93 (6.32)" "0.890"  ""   "0.011"
##pregnancy_day_final (mean (SD)) ""    "54.29 (8.29)" "50.64 (7.99)" "<0.001" ""   "0.448"

table3 <- CreateTableOne(vars = c('Age', 'GS_week2'),data = evaluat_data, strata = 'OPR_Year', smd=TRUE)
table3 <- print(table3,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table3
#write.csv(table2, file = "Table2_before_matching.csv")
#Stratified by OPR_Year
#                      level 2022           2023           p        test SMD    
# n                    ""    "397"          "232"          ""       ""   ""     
#Age (mean (SD))      ""    "32.00 (6.11)" "31.93 (6.32)" "0.890"  ""   "0.011"
#GS_week2 (mean (SD)) ""    "7.30 (1.21)"  "6.81 (1.19)"  "<0.001" ""   "0.413"

##结果中可以看出Age在两组间没有差异而pregnancy_day_final/GS_week变量在两组间是有差异的，其中SMD(standardized mean differences)可以用来衡量协变量在不同组间的差异
evaluat_data0<-na.omit(evaluat_data)
evaluat_data0$treat<-ifelse(evaluat_data0$OPR_Year==2022,0,1)
table(evaluat_data0$treat)
# 0   1 
#394 231 
####有一个没有记录lmp,因此孕天和孕周无法计算

##https://www.cnblogs.com/ayueme/articles/16844697.html
##distance = "logit" # 选择logistic回归
#主要匹配方法选择
#在确定了使用哪种算法计算PS后，匹配方法也是需要注意的一个问题，需要注意以下几个方面，首先是匹配方法的选择（method），然后是采样手段（有无放回），相似度的度量（卡钳值或其他），匹配比例（1:1或1：多）。

##caliper = 0.05,    #卡钳值
##ratio = 1,         #1：N匹配
##replace = F
###1:1匹配策略1：GS_week和Age均采用最近邻匹配
m.out <-matchit(treat ~ Age + GS_week2, method = "nearest",data = evaluat_data0,ratio = 1,replace = F)
summary(m.out,standardize = TRUE)

###1:1匹配策略1：GS_week 和Age均采用精确匹配
m.out <-matchit(treat ~ Age + GS_week2, method = "nearest",data = evaluat_data0,ratio = 1,replace = F,exact = ~ Age + GS_week2)
summary(m.out,standardize = TRUE)

###1:1匹配策略1：最近邻匹配，GS_week采用精确匹配而 控制在Age两岁以内
##https://cloud.tencent.com/developer/ask/sof/1059835/answer/1489437
##使用卡尺，您可以使用完全匹配的卡尺来保留尽可能多的单位。完全匹配通过最小化处理单元和控制单元之间的总层内距离来创建恰好具有一个处理单元或恰好一个控制单元的地层。卡尺确保每一层中的处理单元和控制单元之间的距离得到控制。你可以直接在协变量上设置一个卡尺。例如，为了确保层内的单元年龄不超过2岁，并且年龄和诊断年龄之间的距离不超过2年，您可以使用以下代码：
library(optmatch)
m.out <-matchit(treat ~ Age + GS_week2,method = "nearest",exact= "GS_week2",caliper = c(Age = 2),data = evaluat_data0,ratio = 1,replace = F)
summary(m.out,standardize = TRUE)

mdata <- match.data(m.out)
mdata2<-as.data.frame(mdata[,c("Age","subclass")])
mdata2$subclass<-as.character(mdata2$subclass)
mdata2$group<-rep(c("CT","TR"),length(unique(mdata2$subclass)))
mdata2$group<-factor(mdata2$group,levels=c( "CT","TR"))
mdata3<-as.data.frame(acast(mdata2, subclass ~ group,value.var="Age"))
range(mdata3$CT -mdata3$TR)##-2  2

###1:1匹配策略1：GS_week 精确匹配；Age采用最近邻匹配
##MatchIt包：结合“最近邻”匹配和“精确”匹配：https://it.cha138.com/wen3/show-15901221.html
##https://zhuanlan.zhihu.com/p/596876245
head(evaluat_data0)
m.out <- matchit(treat ~ Age + GS_week2, method = "nearest",distance='logit', exact  = "GS_week2", data = evaluat_data0,ratio = 1,replace = F)
#上面结果说明你匹配方法是1:1无放回最近邻匹配，计算方法是logistic回归，615例总样本中匹配了430例。
summary(m.out,standardize = TRUE)
#####1:1匹配策略1：pregnancy_day_final 和Age采用最近邻匹配
m.out <- matchit(treat ~ Age + pregnancy_day_final, method = "nearest",distance='logit', exact  = "pregnancy_day_final", data = evaluat_data0,ratio = 1,replace = F)
summary(m.out,standardize = TRUE)
###1:1匹配策略1：pregnancy_day_final 和Age均采用精确匹配
m.out <-matchit(treat ~ Age + pregnancy_day_final, method = "nearest",data = evaluat_data0,ratio = 1,replace = F,exact = ~ Age + pregnancy_day_final)
summary(m.out,standardize = TRUE)

###1:1匹配策略1：pregnancy_day_final 采用精确匹配和Age控制在两岁
m.out <-matchit(treat ~ Age + pregnancy_day_final,method = "nearest",exact= "pregnancy_day_final",caliper = c(Age = 2),data = evaluat_data0,ratio = 1,replace = F)
summary(m.out,standardize = TRUE)

###1:1匹配策略1：pregnancy_day_final 采用最近邻匹配和Age控制在两岁
m.out <-matchit(treat ~ Age + pregnancy_day_final,method = "nearest",caliper = c(Age = 2),data = evaluat_data0,ratio = 1,replace = F)
summary(m.out,standardize = TRUE)


##最终选择
###1:1匹配策略1：GS_week 精确匹配；Age控制在两岁
#m.out <-matchit(treat ~ Age + GS_week2,method = "nearest",exact= "GS_week2",caliper = c(Age = 0),data = evaluat_data0,ratio = 1,replace = F)
#m.out <-matchit(treat ~ Age + pregnancy_day_final,exact= "pregnancy_day_final", method = "nearest",caliper = c(Age = 2),data = evaluat_data0,ratio = 1,replace = F)
#m.out <-matchit(treat ~ Age + pregnancy_day_final,exact= "pregnancy_day_final", method = "nearest",caliper = c(Age = 1),data = evaluat_data0,ratio = 1,replace = F)
#m.out <-matchit(treat ~ Age + pregnancy_day_final,exact= "pregnancy_day_final", method = "nearest",caliper = c(Age = 0),data = evaluat_data0,ratio = 1,replace = F)
#m.out <-matchit(treat ~ Age + pregnancy_day_final, method = "nearest",data = evaluat_data0,ratio = 1,replace = F)
#m.out <-matchit(treat ~ Age + GS_week2, method = "nearest",data = evaluat_data0,ratio = 1,replace = F)
set.seed(19921010)
m.out <-matchit(treat ~ Age + GS_week2, method = "nearest",data = evaluat_data0,ratio = 1,replace = F,exact = ~ Age + GS_week2)
#m.out <-matchit(treat ~ Age + pregnancy_day_final, method = "nearest",data = evaluat_data0,ratio = 1,replace = F,exact = ~ Age + pregnancy_day_final)
m.out
summary(m.out,standardize = TRUE)
#          Control Treated
#All           394     231
#Matched       168     168
#Unmatched     226      63
mdata <- match.data(m.out)
mdata<-mdata[order(mdata$subclass,mdata$OPR_Year,decreasing = F),]

mdata2<-as.data.frame(reshape2::dcast(mdata, subclass ~ OPR_Year,value.var="Age"))
t.test(mdata2$'2022',mdata2$'2023',paired = TRUE)
#p-value = NA
mdata2<-as.data.frame(reshape2::dcast(mdata, subclass ~ OPR_Year,value.var="GS_week2"))
t.test(mdata2$'2022',mdata2$'2023',paired = TRUE)
#p-value = NA
mdata2<-as.data.frame(reshape2::dcast(mdata, subclass ~ OPR_Year,value.var="pregnancy_day_final"))
t.test(mdata2$'2022',mdata2$'2023',paired = TRUE)
##p-value = 0.02307

##可通过以下方法获得算法估计的PS值：
eps <- m.out$distance
length(eps)
## [1] 625
head(eps)
##    1         2         3         4         5         6 
#0.3867541 0.3651132 0.2789373 0.2743218 0.3758717 0.3544877 

##通过m.out$match.matrix获取配好的对子：
head(m.out$match.matrix)
#第一列是干预组的序号，第二列是和干预组配好对的，对照组的序号。
m.out$discarded##查看某个样本是否被丢弃：

##匹配后数据的平衡性检验
##检查匹配后的数据，主要是看协变量在不同组间是否已经均衡了（是不是没有差异了）。
##关于这个倾向性评分匹配后数据的平衡性检验，文献中比较推荐使用SMD和VR(variance ratio)，
##SMD<0.25说明均衡了，VR>2.0或者VR<0.5说明很不均衡（越接近1越均衡）！
##但其实也可以用假设检验，比如t检验、卡方检验等，也是没有统一的标准！

##通过summary()查看匹配前后，不同组间协变量的各种统计量。通常建议选择standardize = TRUE查看标准后的各协变量的平衡性指标：
summary(m.out,standardize = TRUE)
##结果主要是3个部分：
#Summary of Balance for All Data：原始数据中干预组和对照组的平均PS值和平均协变量，SMD,VR，每个协变量和PS的CDF（cumulative distribution functions）的均值和最大值
#Summary of Balance for Matched Data：匹配后数据的指标
#Sample Sizes：样本数量

##通过观察比较匹配前后的数据指标可知，Age均衡（-0.0578 <0.1），但是GS_week2也均衡(0.0000 <0.1)！
##这个默认的函数在计算SMD的时候会把分类变量按照连续性变量进行计算，所以计算结果是有一些问题的。
##在一开始计算匹配前数据的SMD时我们用的是tableone这个包，匹配后数据的SMD理论上也是可以用这个包的：
# 首先提取匹配后的数据
mdata <- match.data(m.out)
library(tableone)
table5 <- CreateTableOne(vars = c('Age', 'pregnancy_day_final'),data = mdata,strata = 'treat',smd=TRUE)
table5 <- print(table5, smd=TRUE, showAllLevels = TRUE,noSpaces = TRUE, printToggle = FALSE)
table5
##这个包显示SMD都小于0.25，因此均被均衡了
##tableone这个包计算的SMD也是有一些问题的，推荐大家使用cobalt包进行平衡性指标的计算，专门处理这类匹配问题的
library(cobalt)
# m.threshold表示SMD的阈值，小于这个阈值的协变量是平衡的
bal.tab(m.out, m.threshold = 0.1, un = TRUE)
#计算VR，而且根据VR来看二者均衡。
bal.tab(m.out, v.threshold = 2)

##统计检验衡量均衡性
###除了SMD和VR之外，传统的统计检验也可以用于检查匹配后的数据有没有均衡！
#首先取出匹配好的数据：
mdata <- match.data(m.out)
head(mdata);dim(mdata)
#其中distance是估计的PS，weights是权重，因为我们用的是1:1无放回匹配，所以全都是1。
##下面用t检验看看匹配后干预组和对照组有没有差异：
t.test(Age ~ OPR_Year, data = mdata)## p-value = 1
t.test(GS_week2  ~ OPR_Year, data = mdata)##p-value =1
t.test(pregnancy_day_final  ~ OPR_Year, data = mdata)##p-value = 0.5389

##结果可视化
#默认提供3种图形，但是美观性太差，就不放图了，大家感兴趣的可以自己试试看。
plot(m.out) # 默认QQ图
plot(m.out, type = 'jitter') # 散点图
plot(m.out, type = 'hist') # 直方图
##默认的不好看，还是用cobalt包进行结果的可视化。
Mutiple_plot<-cowplot::plot_grid(bal.plot(m.out, var.name = 'Age', which = 'both', grid=TRUE),
                   bal.plot(m.out, var.name = 'GS_week2', which = 'both', grid=TRUE),
                   bal.plot(m.out, var.name = 'Age', which = 'both', grid=TRUE, type="ecdf"),
                   # 还有很多参数可调整
                   love.plot(bal.tab(m.out, m.threshold=0.1),stat = "mean.diffs",
                             grid=TRUE,stars="raw",abs = F)
)
Mutiple_plot
##上面两幅图展示的是协变量在匹配前（unadjusted sample）和匹配后（adjusted sample）的数据中的分布情况，连续型变量默认是画密度图，分类变量默认是画柱状图。
##左下图是累计密度图。右下的love plot图可视化匹配前后协变量的SMD，两条竖线是0.1阈值线，匹配后x.Age在两条竖线之间，说明平衡，x.Gender不在两条竖线之间，说明还是没平衡。
ggsave(Mutiple_plot,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Matched_sample_infect_or_not_Mutiple_plot.pdf",width = 12, height =12)

#####
mdata <- match.data(m.out)
match_data<- analysis_select_2022_23[which(analysis_select_2022_23$sample %in% unique(mdata$sample)),]
head(match_data);dim(match_data)
mdata$OPR_Year<-factor(mdata$OPR_Year,levels = c(2022,2023))
match_data$OPR_Year<-factor(match_data$OPR_Year,levels = c(2022,2023))

mdata<-mdata[order(mdata$OPR_Year,mdata$subclass,decreasing = F),]
stat_data0<-compare_means(Age~OPR_Year, data=mdata, paired = TRUE)
stat_data0[which(stat_data0$p.signif != "ns"),]

global_Age_plot <- ggpaired(mdata, x = "OPR_Year", y = "Age",color = "OPR_Year", palette = ppCor,line.color = "gray", line.size = 0.4,
                             add = "jitter")+xlab("") +ylab("Ages")+ 
  labs(title =paste0("Age comparison between matched sample in two years"))+ylim(c(0,50))+
  stat_summary(fun=mean, geom="point", shape=20, size=5, color="black", fill="black") +
  stat_compare_means(method = "wilcox.test",paired = TRUE,aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 45)
global_Age_plot
global_Age_plot2<-ggwithinstats( data = mdata,x = OPR_Year, y = Age,type = "nonparametric", centrality.plotting = FALSE )

stat_data0<-compare_means(pregnancy_day_final~OPR_Year, data=mdata, paired = TRUE)
stat_data0[which(stat_data0$p.signif != "ns"),]
# 1 pregnancy_day_final 2022   2023   0.0302  0.03 0.03     *        Wilcoxon


#test <- wilcox.test(mdata2$pregnancy_day_final ~ mdata2$OPR_Year, paired = TRUE)
#test
range(mdata$pregnancy_day_final)##37 75
global_pregnancy_day_plot <- ggpaired(mdata, x = "OPR_Year", y = "pregnancy_day_final",color = "OPR_Year", palette = ppCor,line.color = "gray", line.size = 0.4,
                            add = "jitter")+xlab("") +ylab("Ages")+ 
  labs(title =paste0("pregnancy_day comparison between matched sample in two years"))+ylim(c(0,90))+
  stat_summary(fun=mean, geom="point", shape=20, size=5, color="black", fill="black") +
  stat_compare_means(method = "wilcox.test",paired = TRUE,aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 80)
global_pregnancy_day_plot
global_pregnancy_day_plot2 <-ggwithinstats( data = mdata,x = OPR_Year, y = pregnancy_day_final,type = "nonparametric", centrality.plotting = FALSE )

stat_data0<-compare_means(GS_week2~OPR_Year, data=mdata, paired = TRUE)
stat_data0[which(stat_data0$p.signif != "ns"),]

#test <- wilcox.test(mdata2$pregnancy_day_final ~ mdata2$OPR_Year, paired = TRUE)
#test
range(mdata$GS_week2)#5 10
global_GS_week_day_plot <- ggpaired(mdata, x = "OPR_Year", y = "GS_week2",color = "OPR_Year", palette = ppCor,line.color = "gray", line.size = 0.4,
                                      add = "jitter")+xlab("") +ylab("Ages")+ 
  labs(title =paste0("GS_week comparison between matched sample in two years"))+ylim(c(0,12))+
  stat_summary(fun=mean, geom="point", shape=20, size=5, color="black", fill="black") +
  stat_compare_means(method = "wilcox.test",paired = TRUE,aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 8)
global_GS_week_day_plot
global_GS_week_day_plot2 <-ggwithinstats(data = mdata,x = OPR_Year, y = GS_week2,type = "nonparametric", centrality.plotting = FALSE )


ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Age_matched_sample_compare_stat_box_two_year_compare.pdf",global_Age_plot,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Age_matched_sample_compare_stat_box_two_year_compare2.pdf",global_Age_plot2,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/pregnancy_day_matched_sample_compare_stat_box_two_year_compare.pdf",global_pregnancy_day_plot,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/pregnancy_day_matched_sample_compare_stat_box_two_year_compare2.pdf",global_pregnancy_day_plot2,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/GS_week2_matched_sample_compare_stat_box_two_year_compare.pdf",global_GS_week_day_plot,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/GS_week2_matched_sample_compare_stat_box_two_year_compare2.pdf",global_GS_week_day_plot2,width=6, height=6)


##因子在感染与非感染组的组间比较
head(match_data)
head(mdata)
write.table(as.data.frame(mdata), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/matched_sample_in_two_year_information.txt",sep = "\t",row.names=F) 
write.table(as.data.frame(match_data), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/Index_information_matched_sample_in_two_year.txt",sep = "\t",row.names=F) 

match_data2<-merge(mdata[,c("sample","subclass")],match_data,by="sample")
head(match_data2)
match_data3<-match_data2[,c("sample","Year_month","Age","subclass","OPR_Year","E_name","value")]
head(match_data3)

index_order<-c("WBC","RBC","HGB","HCT","MCV","MCH","MCHC","PLT","RDW_CV",
               "PCT","LYMPH_pct","LYMPH_av","Neut_pct","Neut_av","EO_pct","EO_av","BASO_pct","BASO_av","MONO_pct","MONO_av",
               "MPV","PDW","P-LCR","ALT","AST","TBIL","γ-GT","CK","CK-MB","BUN_Urea","Cr","UA","T-CHO","TG",
               "HDL-C","LDL-C","GLU","PT","A","INR","Fib","APTT","APTT_R","TT","R")
length(index_order)
match_data3$E_name <- factor(match_data3$E_name,levels=index_order)
match_data3$OPR_Year<-factor(match_data3$OPR_Year,levels = c(2022,2023))
match_data3<-match_data3[order(match_data3$E_name,match_data3$OPR_Year,match_data3$subclass,decreasing = F),]
match_data3<-na.omit(match_data3)
head(match_data3)
##配对比较仅保留配对数据
match_data3$type<-paste0(match_data3$E_name,":",match_data3$subclass)
number_table<-as.data.frame(table(match_data3$type))
remain_names<-as.character(number_table[which(number_table$Freq ==2),]$Var1)
match_data4<-match_data3[which(match_data3$type %in% remain_names),]
dim(match_data3);dim(match_data4)#14466 13876  

##添加分组
head(analysis_collect2)
sample_infor<-distinct(analysis_collect2[,c("sample","Infect_state","LMP_infect")])
head(sample_infor);dim(sample_infor)#249  3

class_infor<-distinct(match_data4[which(match_data4$OPR_Year==2023),c("sample","subclass")])
head(class_infor);dim(class_infor)
class_infor2<-merge(class_infor,sample_infor,by="sample")
class_infor2<-class_infor2[,-1]
head(class_infor2)

match_data5<-merge(match_data4,class_infor2,by="subclass",all.x  = T)
head(match_data5)
match_data5$group_day<-ifelse(match_data5$LMP_infect <= 7,"D7",
                                     ifelse(match_data5$LMP_infect <= 21,"D21",
                                            ifelse(match_data5$LMP_infect <= 35,"D35",
                                                   ifelse(match_data5$LMP_infect <= 49,"D49",
                                                          ifelse(match_data5$LMP_infect <= 63,"D63",
                                                                 ifelse(match_data5$LMP_infect <= 77,"D77",
                                                                        ifelse(match_data5$LMP_infect <= 91,"D91","Dmore91")))))))
table(match_data5$group_day)
match_data5$group_day<-factor(match_data5$group_day,levels = c("D7","D21","D35","D49","D63","D77","D91","Dmore91"))

match_data5<-match_data5[order(match_data5$E_name,match_data5$OPR_Year,match_data5$subclass,decreasing = F),]
head(match_data5)
table(match_data5$E_name)

dim(match_data5);dim(match_data4)# 13876    11 13876     8
#match_data4<-match_data4[order(match_data4$E_name,match_data4$OPR_Year,match_data4$subclass,decreasing = F),]
#match_data4<-match_data4[,-c(ncol(match_data4))]
#head(match_data4)
#table(match_data4$E_name)

plot_list<-list();pvalue<-data.frame()
for ( index_name in index_order){
  # index_name<-"γ-GT"#test line
  print(index_name)
  match_data6<-match_data5[which(match_data5$E_name ==index_name),]
 # match_data6<-match_data5[which(match_data4$OPR_Year ==2023),]
  #160
 # match_data6[which(match_data6$subclass ==160),]
  
  match_data6$subclass<-as.numeric(match_data6$subclass)
  #match_data5<-match_data5[,-c(3)]
  match_data6<-match_data6[order(match_data6$OPR_Year,match_data6$subclass,decreasing = F),]
  
  stat_data<-compare_means(value ~ OPR_Year, data=match_data6,paired = TRUE)
  stat_data1<-as.data.frame(stat_data)
  stat_data1$E_name<-index_name
  pvalue<-rbind(pvalue,stat_data1)
  
  compare_plot<- ggpaired(match_data6, x = "OPR_Year", y = "value",color = "OPR_Year", palette = ppCor,line.color = "gray", line.size = 0.4,add = "jitter")+
    xlab("Year") +ylab(paste0(index_name," value"))+ labs(title =index_name)+
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="black", fill="black") +
    stat_compare_means(method = "wilcox.test",paired = TRUE,aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.5)
  plot_list<-c(plot_list,list(compare_plot))
}
names(plot_list)<-index_order
pvalue
write.table(as.data.frame(pvalue), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/Index_matched_sample_compare_stat_in_two_year_compare.txt",sep = "\t",row.names=F) 

rownames(pvalue[which(pvalue$p.signif != "ns"),])
## "8"  "10" "22" "36" "37" "38" "39" "40" "44" "45"

hplot_merge<-ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],
          plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],
          plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],
          plot_list[[16]],plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]],
          plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]],
          plot_list[[26]],plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],
          plot_list[[31]],plot_list[[32]],plot_list[[33]],plot_list[[34]],plot_list[[35]],
          plot_list[[36]],plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]],
          plot_list[[41]],plot_list[[42]],plot_list[[43]],plot_list[[44]],plot_list[[45]],ncol = 7, nrow = 7)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Index_matched_sample_compare_stat_box_two_year_compare.pdf",hplot_merge,width=30, height=30)

##仅展示significant elements
hplot_merge2<-ggarrange(plot_list[[8]],plot_list[[10]],plot_list[[22]],plot_list[[36]],
                       plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]],
                       plot_list[[44]],plot_list[[45]],ncol = 4, nrow = 3)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Significant_Index_matched_sample_compare_stat_box_two_year_compare.pdf",hplot_merge2,width=16, height=13)
####################################
plot_list<-list();pvalue<-data.frame()
for ( index_name in index_order){
  # index_name<-"TT"#test line
  print(index_name)
  match_data6<-match_data5[which(match_data5$E_name ==index_name),]
  # match_data6<-match_data5[which(match_data4$OPR_Year ==2023),]
  #160
  # match_data6[which(match_data6$subclass ==160),]
  
  match_data6$subclass<-as.numeric(match_data6$subclass)
  #match_data5<-match_data5[,-c(3)]
  match_data6<-match_data6[order(match_data6$OPR_Year,match_data6$subclass,decreasing = F),]
  match_data6$group_day<-factor(match_data6$group_day,levels = c("D7","D21","D35","D49","D63","D77","D91","Dmore91"))
  
  stat_data<-compare_means(value ~ OPR_Year,data=match_data6,group.by = "group_day",paired = TRUE)
  stat_data1<-as.data.frame(stat_data)
  stat_data1$E_name<-index_name
  pvalue<-rbind(pvalue,stat_data1)
  
  compare_plot<- ggpaired(match_data6, x = "OPR_Year", y = "value",color = "OPR_Year",palette = ppCor,line.color = "gray", line.size = 0.4,add = "jitter")+
    xlab("Year") +ylab(paste0(index_name," value"))+ labs(title =index_name)+
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="black", fill="black") +
    facet_wrap(~group_day,scales = "fixed",ncol = 9)+
    stat_compare_means(method = "wilcox.test",paired = TRUE,aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.5)
  
  plot_list<-c(plot_list,list(compare_plot))
}
names(plot_list)<-index_order
pvalue
write.table(as.data.frame(pvalue), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/Index_matched_sample_compare_stat_in_two_year_compare_group_by_group_day.txt",sep = "\t",row.names=F) 

rownames(pvalue[which(pvalue$p.signif != "ns"),])

hplot_merge<-ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],
                       plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],
                       plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],
                       plot_list[[16]],plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]],
                       plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]],
                       plot_list[[26]],plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],
                       plot_list[[31]],plot_list[[32]],plot_list[[33]],plot_list[[34]],plot_list[[35]],
                       plot_list[[36]],plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]],
                       plot_list[[41]],plot_list[[42]],plot_list[[43]],plot_list[[44]],plot_list[[45]],ncol = 5, nrow = 9)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Index_matched_sample_compare_stat_box_two_year_compare_group_by_group_day.pdf",hplot_merge,width=50, height=50,limitsize = F)

##仅展示significant elements
unique(pvalue[which(pvalue$p.signif != "ns"),]$E_name)
which(index_order %in% unique(pvalue[which(pvalue$p.signif != "ns"),]$E_name))
#6  8  9 10 15 16 18 21 22 24 28 31 32 33 36 37 38 39 40 44 45

hplot_merge2<-ggarrange(plot_list[[8]],plot_list[[10]],plot_list[[22]],plot_list[[36]],
                        plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]],
                        plot_list[[44]],plot_list[[45]],plot_list[[6]],plot_list[[9]],
                        plot_list[[15]],plot_list[[16]],plot_list[[18]],plot_list[[21]],
                        plot_list[[24]],plot_list[[28]],plot_list[[31]],plot_list[[32]],
                        plot_list[[33]],plot_list[[36]],ncol = 2, nrow =11)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Significant_Index_matched_sample_compare_stat_box_two_year_compare_group_by_group_day.pdf",hplot_merge2,width=30, height=50,limitsize = F)
####################################

##########################################################################
####################################################以下不跑#####20230611
###用nonrandom进行匹配
##ref::https://www.dxy.cn/bbs/newweb/pc/post/38561237?from=recommend
library(nonrandom)
###匹配前评估
mydata.ps <- pscore(data= evaluat_data0,  formula = treat~Age+pregnancy_day_final)
plot.pscore(mydata.ps,with.legend=T,legend.cex=1,
            main = "PS distribution",
            par.1=list(col="red",lwd=2),             
            par.0=list(col="blue",lwd=2,lty=2),xlim=c(-0.5,1))
##结果：#显示处理组和非处理倾向评分密度分布,par.1处理组作图参数,#par.0非处理组作图参数。
##输入：
mydata.ps$data                    #显示倾向评分数据
summary(mydata.ps)            #显示logistic回归模型
mydata.match <- ps.match(object= mydata.ps, who.treated =1,ratio = 1,
                         caliper= "logit",x=0.2,givenTmatchingC = T,
                         matched.by = "pscore", setseed = 19921010)
#利用ps.match函数进行倾向评分匹配，who.treated =1,表示1代表处理组，ratio确定匹配比例，x代表得分容差，givenTmatchingC = T表示用未处理记录去匹配处理组记录。
##输入：summary(mydata.match)
pair<- mydata.match$data.matched   #生成匹配后的数据
#输入：
head(pair);dim(pair)
stable1 <- CreateTableOne(vars=c("Age","pregnancy_day_final"),strata="treat", data=pair)
print(stable1,showAllLevels = TRUE)

#对匹配后的数据进行统计分析
#显示各协变量差异无统计学意义，匹配效果较理想。

#输入：
dist.plot(object=mydata.match,sel=c("Age"), compare=T,lable.match=c("original data","matched sample"))
dist.plot(object=mydata.match, sel=c("Age"),compare=T,plot.type=2, lable.match=c("original data","matched sample"))
#输出分类变量匹配前后的直条图，可更换其它分类变量以显示结果。

diffm <- ps.balance(object=mydata.match,sel=c("Age","pregnancy_day_final"),cat.levels = 2,method="classical",alpha=5)
diffm  
#对匹配后的数据进行均衡性假设检验， cat.levels = 2,表示数值大于2的变量的定量变量，method="classical",
#默认对定量和分类变量进行ttest和chisq.test，alpha=5设定显著性水平为0.05。

diffm1 <- ps.balance(object=mydata.match,sel=c("Age","pregnancy_day_final"),cat.levels = 2,method="stand.diff",alpha=20)
diffm1  
#  method="stand.diff",表示采用标准差异法进行均衡性检验， alpha=20设定差异的临界值为0.2。
plot.stdf(x = diffm1,sel = c("Age","pregnancy_day_final"),
          main = "standardized differences of covariates",
          pch.p=c(20,10),col.p=c("red","blue"),
          legend.cex=1,legend.xy=c(35,6.5))

#图示化匹配前后协变量标准化差异值。

##############################
##########手动提取
###1:1匹配策略2：GS_week 精确匹配；Age差±2岁
####https://www.zhihu.com/question/37665493
evaluat_data2<-data.frame(a1=evaluat_data$Age+2,a2=evaluat_data$Age+1,a3=evaluat_data$Age,a4=evaluat_data$Age-1,a5=evaluat_data$Age-2)
evaluat_data3<-cbind(evaluat_data,evaluat_data2)
dim(evaluat_data3)##619  10
evaluat_data3<-na.omit(evaluat_data3)
dim(evaluat_data3)##615  10

evaluat_data_2022<-evaluat_data3[which(evaluat_data3$OPR_Year==2022),]
evaluat_data_2023<-evaluat_data3[which(evaluat_data3$OPR_Year==2023),]
library(data.table)
head(evaluat_data_2022)
seek_table<-as.data.table(evaluat_data_2022,keep.rownames=TRUE)
head(seek_table)
set_table<-as.data.table(evaluat_data_2023,keep.rownames=TRUE)
head(set_table)

set_table2<-reshape(data=set_table,varying = c("a1","a2","a3","a4","a5"),v.names = "merge_value",timevar = "Age_type",direction = "long")
setkey(set_table2,GS_week2,merge_value)

seek_table[,merge_value:=Age];
setkey(seek_table,GS_week2,merge_value)

Big_table <-merge(seek_table,set_table2[,list(names.b=sample,GS_week2,merge_value,age.b=Age,Age_type)],by=c("GS_week2","merge_value"),allow.cartesian = T)
Big_table <-Big_table[order(Big_table$names.b,Big_table$Age_type),]
head(Big_table)
#要想办法解决Age和pregnancy_day_final两个变量在两组间的差异，达到基线可比的目的。
# --------PSM 1:N----------
#加载数据
head(analysis_select_2022_23)
#index_order<-c("WBC","RBC","HGB","HCT","MCV","MCH","MCHC","PLT","RDW_CV",
#"PCT","LYMPH_pct","LYMPH_av","Neut_pct","Neut_av","EO_pct","EO_av","BASO_pct","BASO_av","MONO_pct","MONO_av",
#"MPV","PDW","P-LCR","ALT","AST","TBIL","γ-GT","CK","CK-MB","BUN_Urea","Cr","UA","T-CHO","TG",
#"HDL-C","LDL-C","GLU","PT","A","INR","Fib","APTT","APTT_R","TT","R")
index<-"WBC"
data<-distinct(analysis_select_2022_23[which(analysis_select_2022_23$E_name == index),c("sample","Age","OPR_Year","day_detected","GS_week")])
head(data);dim(data)#657   5
data<-na.omit(data)
table(data$OPR_Year)
#2022 2023 
# 378  229 
names(data)
#多变量因子化
str(data)
##data[,4:8]<-lapply(data[,4:8],as.factor)
#

matchlist <- matchit(OPR_Year ~ Age+GS_week, data=data,
                     method   = "nearest",
                     distance = "glm",      #
                     caliper  = 0.05,       # 卡钳值
                     ratio    = 1,          # 1:N 匹配
                     replace  = F)          # 不替换
summary(matchlist)

# 提取匹配后的病例对照
matchdata<- match.data(matchlist,
                       group = "all",
                       distance = "distance",
                       weights = "weights",
                       subclass = "subclass",
                       data = NULL,
                       include.s.weights = TRUE,
                       drop.unmatched = TRUE)
table(duplicated(matchdata$id))
table(data$OPR_Year)      # 匹配前
#2022 2023 
# 378  229 
table(matchdata$OPR_Year) # 匹配后
#2022 2023 
#218  218 

#排序配对子
matchdata$subclass <- sort(matchdata$subclass)
#查看配对子情况，有的病例只有一个对照
a=table(matchdata$subclass);table(a)

match2data <- data.frame()
for (i in seq_along(table(matchdata$subclass))) {
  if(table(matchdata$subclass)[i]==2){
    match2data[i,1]=i}
}
match2data <- tidyr::drop_na(match2data,V1)
names(match2data)[1] <- "subclass"
match2data$subclass <- as.factor(match2data$subclass)
match2data2 <- dplyr::left_join(match2data,matchdata)


#-------制作表格-----------
library(tableone)

vars <- c("sex","age","smoke","drink","hypertension",
          "diabetes" )
catVars <- c("sex","drink","hypertension",
             "diabetes","smoke")
nonVars <- c("age")
#建立表格
tabMatched <- CreateTableOne(vars = vars, 
                             strata = "Group", 
                             factorVars = catVars,
                             data = matchdata, 
                             addOverall = TRUE,#添加Overall列的分析结果 
                             test = TRUE)
#最终表格
tab5 <- print(tabMatched, 
              showAllLevels = TRUE, #显示所有水平,不折叠
              cramVars = catVars, 
              nonnormal = nonVars, 
              #exact = "M", 
              catDigits=2,        #分类变量保留小数位数
              contDigits=2,       #连续变量保留小数位数
              quote = FALSE,      #不显示引号 
              noSpaces = TRUE,    #是否删除为对齐而添加的空间。
              #如果您希望自己在其他软件中对齐数字，请使用此选项。
              printToggle = FALSE #输出matrix
) 
tab5
write.csv(tab5, "Table_One2.csv", quote=TRUE, row.names=TRUE) 

#------PSM的条件logistic回归----------
library(reportReg);library(survival)
model <- clogit(Group~ hypertension+age+ strata(subclass), 
                data=matchdata)
reportReg(model)
summary(model)
