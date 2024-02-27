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


month_select1<-c("20221","20222","20223","20224")
analysis_used <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/clinical_index_data_for_analysis_used_final.txt",header = T,sep = "\t")
dim(analysis_used)# 78529    14
index_order<-c("WBC","RBC","HGB","HCT","MCV","MCH","MCHC","PLT","RDW_CV",
               "PCT","LYMPH_pct","LYMPH_av","Neut_pct","Neut_av","EO_pct","EO_av","BASO_pct","BASO_av","MONO_pct","MONO_av",
               "MPV","PDW","P-LCR","ALT","AST","TBIL","γ-GT","CK","CK-MB","BUN_Urea","Cr","UA","T-CHO","TG",
               "HDL-C","LDL-C","GLU","PT","A","INR","Fib","APTT","APTT_R","TT","R")
length(index_order)#45
analysis_used$E_name <- factor(analysis_used$E_name,levels=index_order)
analysis_used$value<-as.numeric(analysis_used$value)

#提取2022年度前四个月信息
analysis_select_2022<-analysis_used[which(analysis_used$Year_month %in%month_select1),]
head(analysis_select_2022)


##读取被记录的康复信息完整的样本内容
collect_sample_infor <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/collect_sample_infor.txt",header = T,sep = "\t")
collect_sample_infor$ID_full<-str_pad(collect_sample_infor$ID,width =8 ,side = c("left"),pad = "0")
head(collect_sample_infor);dim(collect_sample_infor)#257   6
length(unique(collect_sample_infor$ID_full))#257

analysis_collect<-analysis_used[which(analysis_used$ID_full %in% unique(collect_sample_infor$ID_full)),]
table(analysis_collect$Year_month)
analysis_collect<-analysis_collect[which(analysis_collect$Year_month %in% c("20232","20233","20234")),]
dim(analysis_collect)# 11078    14

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

##提取所有的2023年明确感染和非感染的信息数据
head(analysis_collect2)
analysis_collect3<-analysis_collect2#[which(analysis_collect2$Infect_state =="yes"),]
analysis_select_2023<-analysis_used[which(analysis_used$sample %in% unique(analysis_collect3$sample)),]
head(analysis_select_2023)
table(analysis_select_2023$Year_month)

##合并两个年度信息数据
analysis_select_2022_23<-rbind(analysis_select_2022,analysis_select_2023)
head(analysis_select_2022_23)
analysis_select_2022_23$GS_week<-floor(analysis_select_2022_23$day_detected /7)
analysis_select_2022_23$GS_week2<-floor(analysis_select_2022_23$pregnancy_day_final/7)
##
analysis_select_2022_23$OPR_Year <- factor(analysis_select_2022_23$OPR_Year,levels=c("2022","2023"))
table(analysis_select_2022_23$OPR_Year)
head(analysis_select_2022_23)

write.table(as.data.frame(analysis_select_2022_23), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/All_collected_sample_Infect_or_not_clinical_index_4month_all.txt",sep = "\t",row.names=F) 

####数据查看
analysis_select_2022_23 <- read.table(file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/All_collected_sample_Infect_or_not_clinical_index_4month_all.txt",header = T,sep = "\t")

###进行两组间样本匹配分析：重点参考：https://www.cnblogs.com/ayueme/articles/16844697.html
##首先可以看一下原始数据的基线资料表，用的是tableone这个包，它能计算SMD（后面会介绍这个SMD的作用）
library(tableone)
evaluat_data<-distinct(analysis_select_2022_23[,c("sample","Age","pregnancy_day_final","GS_week2","OPR_Year")])## 619   2

head(evaluat_data);dim(evaluat_data)#645    5
table2 <- CreateTableOne(vars = c('Age', 'pregnancy_day_final'),data = evaluat_data, strata = 'OPR_Year', smd=TRUE)
# factorVars = c('x.Gender', 'CVD'),
table2 <- print(table2,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table2

evaluat_data0<-na.omit(evaluat_data)
evaluat_data0$treat<-ifelse(evaluat_data0$OPR_Year==2022,0,1)
table(evaluat_data0$treat)
#0   1 
#394 248

##最终选择
###1:1匹配策略1：Age和GS_week 精确匹配
m.out <-matchit(treat ~ Age + GS_week2, method = "nearest",data = evaluat_data0,ratio = 1,replace = F,exact = ~ Age + GS_week2)
#m.out <-matchit(treat ~ Age + pregnancy_day_final, method = "nearest",data = evaluat_data0,ratio = 1,replace = F,exact = ~ Age + pregnancy_day_final)
m.out
summary(m.out,standardize = TRUE)
#          Control   Treated
#All           394      248
#Matched        176     176
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
##p-value =  0.01316

##可通过以下方法获得算法估计的PS值：
eps <- m.out$distance
length(eps)
## [1] 642
head(eps)
##    1         2         3         4         5         6 
#0.4126789 0.3759879 0.2884757 0.2806527 0.3941843 0.3581350 

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
#                                 level 0              1              p       test SMD     
#n                               ""    "176"          "176"          ""      ""   ""      
#Age (mean (SD))                 ""    "31.43 (5.85)" "31.43 (5.85)" "1.000" ""   "<0.001"
#pregnancy_day_final (mean (SD)) ""    "51.52 (6.94)" "51.03 (7.03)" "0.512" ""   "0.070" 

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
t.test(pregnancy_day_final  ~ OPR_Year, data = mdata)##p-value = 0.5122

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
##上面两幅图展示的是协变量在匹配前（unadjusted sample）和匹配后（adjusted sample）的数据中的分布情况，连续型变量默认是画密度图，分类变量默认是画柱状图。
##左下图是累计密度图。右下的love plot图可视化匹配前后协变量的SMD，两条竖线是0.1阈值线，匹配后x.Age在两条竖线之间，说明平衡，x.Gender不在两条竖线之间，说明还是没平衡。
ggsave(Mutiple_plot,file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Matched_sample_for_all_collected_sample_infect_or_not_Mutiple_plot.pdf",width = 12, height =12)

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
# 1 pregnancy_day_final  2022   2023   0.0165 0.016 0.016      *        Wilcoxon


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


ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Age_matched_All_collected_sample_compare_stat_box_two_year_compare.pdf",global_Age_plot,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Age_matched_All_collected_sample_compare_stat_box_two_year_compare2.pdf",global_Age_plot2,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/pregnancy_day_matched_All_collected_sample_compare_stat_box_two_year_compare.pdf",global_pregnancy_day_plot,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/pregnancy_day_matched_All_collected_sample_compare_stat_box_two_year_compare2.pdf",global_pregnancy_day_plot2,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/GS_week2_matched_All_collected_sample_compare_stat_box_two_year_compare.pdf",global_GS_week_day_plot,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/GS_week2_matched_All_collected_sample_compare_stat_box_two_year_compare2.pdf",global_GS_week_day_plot2,width=6, height=6)


##因子在感染与非感染组的组间比较
head(match_data)
head(mdata)
write.table(as.data.frame(mdata), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/matched_All_collected_sample_in_two_year_information.txt",sep = "\t",row.names=F) 
write.table(as.data.frame(match_data), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/Index_information_matched_All_collected_sample_in_two_year.txt",sep = "\t",row.names=F) 

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
dim(match_data3);dim(match_data4)#15145   14498 

##添加分组
head(analysis_collect2)
sample_infor<-distinct(analysis_collect2[,c("sample","Infect_state","LMP_infect")])
head(sample_infor);dim(sample_infor)

class_infor<-distinct(match_data4[which(match_data4$OPR_Year==2023),c("sample","subclass")])
head(class_infor);dim(class_infor)
class_infor2<-merge(class_infor,sample_infor,by="sample")
class_infor2<-class_infor2[,-1]
head(class_infor2)

match_data5<-merge(match_data4,class_infor2,by="subclass",all.x  = T)
head(match_data5)
match_data5$group_day<-ifelse(match_data5$Infect_state =="no","CTRL",
                                    ifelse(match_data5$LMP_infect <= 7,"D7",
                                           ifelse(match_data5$LMP_infect <= 21,"D21",
                                                  ifelse(match_data5$LMP_infect <= 35,"D35",
                                                         ifelse(match_data5$LMP_infect <= 49,"D49",
                                                                ifelse(match_data5$LMP_infect <= 63,"D63",
                                                                       ifelse(match_data5$LMP_infect <= 77,"D77",
                                                                              ifelse(match_data5$LMP_infect <= 91,"D91","Dmore91"))))))))
table(match_data5$group_day)
match_data5$group_day<-factor(match_data5$group_day,levels = c("CTRL","D7","D21","D35","D49","D63","D77","D91","Dmore91"))

match_data5<-match_data5[order(match_data5$E_name,match_data5$OPR_Year,match_data5$subclass,decreasing = F),]
head(match_data5)
table(match_data5$E_name)

plot_list<-list();pvalue<-data.frame()
for ( index_name in index_order){
  # index_name<-"TT"#test line
  print(index_name)
  match_data6<-match_data5[which(match_data5$E_name ==index_name),]
 ## match_data7<-match_data6[which(match_data5$OPR_Year ==2023),]
  #160
#  match_data7[which(match_data7$subclass ==160),]
  
  match_data6$subclass<-as.numeric(match_data6$subclass)
  #match_data6<-match_data6[,-c(3)]
  match_data6<-match_data6[order(match_data6$OPR_Year,match_data6$subclass,decreasing = F),]
  
  stat_data<-compare_means(value ~ OPR_Year,data=match_data6,group.by = "Infect_state",paired = TRUE)
  stat_data1<-as.data.frame(stat_data)
  stat_data1$E_name<-index_name
  pvalue<-rbind(pvalue,stat_data1)
  
  compare_plot<- ggpaired(match_data6, x = "OPR_Year", y = "value",color = "OPR_Year", facet.by = "Infect_state",palette = ppCor,line.color = "gray", line.size = 0.4,add = "jitter")+
    xlab("Year") +ylab(paste0(index_name," value"))+ labs(title =index_name)+
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="black", fill="black") +
    stat_compare_means(method = "wilcox.test",paired = TRUE,aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.5)
  
  plot_list<-c(plot_list,list(compare_plot))
}
names(plot_list)<-index_order
pvalue
write.table(as.data.frame(pvalue), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/Index_matched_all_collected_sample_compare_stat_in_two_year_compare_group_by_infect_stage.txt",sep = "\t",row.names=F) 

rownames(pvalue[which(pvalue$p.signif != "ns"),])

hplot_merge<-ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],
                       plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],
                       plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],
                       plot_list[[16]],plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]],
                       plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]],
                       plot_list[[26]],plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],
                       plot_list[[31]],plot_list[[32]],plot_list[[33]],plot_list[[34]],plot_list[[35]],
                       plot_list[[36]],plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]],
                       plot_list[[41]],plot_list[[42]],plot_list[[43]],plot_list[[44]],plot_list[[45]],ncol = 7, nrow = 7)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Index_matched_all_collected_sample_compare_stat_box_two_year_compare_group_by_infect_stage.pdf",hplot_merge,width=30, height=30)

##仅展示significant elements
unique(pvalue[which(pvalue$p.signif != "ns"),]$E_name)
which(index_order %in% unique(pvalue[which(pvalue$p.signif != "ns"),]$E_name))
# 8  10 22 36 37 38 39 40  44 45  #多了 9 和23 和41
hplot_merge2<-ggarrange(plot_list[[8]],plot_list[[10]],plot_list[[22]],
                        plot_list[[36]],plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]],
                        plot_list[[44]],plot_list[[45]],plot_list[[9]],plot_list[[23]],plot_list[[41]],ncol = 5, nrow = 3)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Significant_Index_matched_all_collected_sample_compare_stat_box_two_year_compare_group_by_infect_stage.pdf",hplot_merge2,width=16, height=13)

####康复天数进行分组查看临床index

plot_list<-list();pvalue<-data.frame()
for ( index_name in index_order){
  # index_name<-"PCT"#test line
  print(index_name)
  match_data6<-match_data5[which(match_data5$E_name ==index_name),]
  ## match_data7<-match_data6[which(match_data5$OPR_Year ==2023),]
  #160
  #  match_data7[which(match_data7$subclass ==160),]
  
  match_data6$subclass<-as.numeric(match_data6$subclass)
  #match_data6<-match_data6[,-c(3)]
  match_data6<-match_data6[order(match_data6$OPR_Year,match_data6$subclass,decreasing = F),]
  match_data6$group_day<-factor(match_data6$group_day,levels = c("CTRL","D7","D21","D35","D49","D63","D77","D91","Dmore91"))
  
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
write.table(as.data.frame(pvalue), file="D:/PROJECT/新冠/manuscript/1.data/1.clinical_data/Index_matched_all_collected_sample_compare_stat_in_two_year_compare_group_by_group_day.txt",sep = "\t",row.names=F) 

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
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Index_matched_all_collected_sample_compare_stat_box_two_year_compare_group_by_group_day.pdf",hplot_merge,width=50, height=50,limitsize = F)

##仅展示significant elements
unique(pvalue[which(pvalue$p.signif != "ns"),]$E_name)
which(index_order %in% unique(pvalue[which(pvalue$p.signif != "ns"),]$E_name))
#  8  10 22 36 37 38 39 40 44 45# 1  6 9 14 15 16 18 24 28 31 32 34 
#hplot_merge2<-ggarrange(plot_list[[8]],plot_list[[10]],plot_list[[22]],
#                        plot_list[[36]],plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]],
#                        plot_list[[44]],plot_list[[45]],plot_list[[9]],plot_list[[23]],plot_list[[41]],ncol = 2, nrow =7)
#ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/match_group/Significant_Index_matched_all_collected_sample_compare_stat_box_two_year_compare_group_by_group_day.pdf",hplot_merge2,width=30, height=30)
