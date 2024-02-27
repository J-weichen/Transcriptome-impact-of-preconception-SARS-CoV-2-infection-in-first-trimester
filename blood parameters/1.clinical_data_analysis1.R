rm(list = ls())
library(MatchIt)

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

###第一步：处理202201至202302的手术信息
##data clean 计划生育手术室提供的手术信息导入   Note: 最终表格从LINE 323开始
clinical_data0<-as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/1.2022_计划生育手术信息.xlsx", sheet =3, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_data0)
clinical_data0_raw<-clinical_data0[,-c(5,7:10)]
head(clinical_data0_raw)

######导入操作时间与超声结果
##手术时间
clinical_data1<-as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/2.xinxike_data_3.xlsx", sheet =6, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_data1)
clinical_data1_raw<-clinical_data1[,-c(2:5)]
head(clinical_data1_raw)
##超声数据
clinical_data2<-as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/2.xinxike_data_3.xlsx", sheet =8, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_data2);dim(clinical_data2)
clinical_data2_raw<-clinical_data2[,c(1,9:11)]
head(clinical_data2_raw)
dim(clinical_data1_raw);dim(clinical_data2_raw) #1660   11  1554    4
##获得每个样本手术时间与超声数据
clinical_data_xinxike_raw<-merge(clinical_data1_raw,clinical_data2_raw,by="序号",all = T)
dim(clinical_data_xinxike_raw)#1730   14
head(clinical_data_xinxike_raw)

#########与计划生育手术信息单子合并
dim(clinical_data0_raw);dim(clinical_data_xinxike_raw)#1774   16 1730   14

clinical_data_all<-merge(clinical_data0_raw,clinical_data_xinxike_raw,by="序号",all = T)

head(clinical_data_all);dim(clinical_data_all)#1774   29
clinical_data_all$Gap_time_qual<-ifelse(clinical_data_all$JHSY_pregnancy_day == clinical_data_all$pregnancy_day,"True","False")
clinical_data_all[which(clinical_data_all$Gap_time_qual=="False"),c("JHSY_LMP_day","LMP_days","JHSY_pregnancy_day","pregnancy_day")]

clinical_data_all$GS_size_qual<-ifelse(clinical_data_all$JHSY_gestational_sac_size == clinical_data_all$gestational_sac_size,"True","False")
clinical_data_all[which(clinical_data_all$GS_size_qual=="False"),c("JHSY_gestational_sac_size","gestational_sac_size")]

###根据信息科查询以及手术室本身记录的妊娠天数，胎囊大小确定待人工审查的信息条目，确定LMP日期以及最终天数和胎囊大小
table(clinical_data_all$Gap_time_qual)
#False  True 
#  139  1517 
table(clinical_data_all$GS_size_qual)
#False  True 
#  268  1106
write.table(as.data.frame(clinical_data_all), file="D:/PROJECT/新冠/2022_data/clinical_data_all.txt",sep = "\t",row.names=F) 
######此表格在添加另外的2023年3月到4月的信息后完毕后改名为clinical_data_all_final.xlsx

#operation_date <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/xinxike_data_3.xlsx", sheet =3, col_names = T, col_types = NULL, na = "", skip = 0))
#head(operation_date)
#LMP_date <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/xinxike_data_3.xlsx", sheet =4, col_names = T, col_types = NULL, na = "", skip = 0))
#head(LMP_date)

#remain_sample<-raw_data[which(raw_data$序号 %in%  setdiff(operation_date$序号,LMP_date$序号)),]
#write.table(as.data.frame(remain_sample[,c("序号","出生日期","手术日期","门诊就诊日期","手术间隔时间","就诊科室","现病史、既往史")]), file="D:/PROJECT/新冠/2022_data/remain_sample.txt",sep = "\t",row.names=F) 

#LMP_date <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/xinxike_data_3.xlsx", sheet =2, col_names = T, col_types = NULL, na = "", skip = 0))

###########整合两部分检验检测信息###################

##导入添加了就诊科室的检验指标表格
clinical_data <-as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/3.无痛人流术前检验数据-新增就诊科室.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_data);dim(clinical_data)#88334    11
##导入先前已整理完毕的检验指标表格
clinical_data3 <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/2.xinxike_data_3.xlsx", sheet =10, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_data3);dim(clinical_data3)#85610    12

clinical_data$tag<-paste(clinical_data$序号,clinical_data$手术日期,clinical_data$检验细项,sep=":")
clinical_data3$tag<-paste(clinical_data3$序号,clinical_data3$手术日期,clinical_data3$检验细项,sep=":")
length(intersect(clinical_data$tag,clinical_data3$tag))#85610
length(setdiff(clinical_data$tag,clinical_data3$tag))#2724
rownames(clinical_data)<-clinical_data$tag
rownames(clinical_data3)<-clinical_data3$tag
##根据序号，手术日期和检查的临床指标提取先前整理的数据条目对应科室以及就诊次数
clinical_data4<-clinical_data[rownames(clinical_data3),c("对应就诊科室","对应就诊次")]
head(clinical_data4);dim(clinical_data4)#85610    2
head(clinical_data3);dim(clinical_data3)#85610    13
clinical_data05 <- cbind(clinical_data3,clinical_data4)
clinical_data5<-clinical_data05[which(!(is.na(clinical_data05$result))),]
head(clinical_data5);dim(clinical_data5)#84877    15
###确定是否在同一天对同一项目检验两次的个体==》没有
clinical_data5$tag<-paste(clinical_data5$序号,clinical_data5$E_name,clinical_data5$手术间隔时间,sep=":")
clinical_data5$tag[duplicated(clinical_data5$tag)]##character(0)
rownames(clinical_data5)<-clinical_data5$tag
head(clinical_data5);dim(clinical_data5)#84877    15

##导入202303至202304两个月的检验数据指标
clinical_data07  <-as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/4.3-4月份信息科调取数据2.xlsx", sheet =5, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_data07);dim(clinical_data07)#15016    14
clinical_data7<-clinical_data07[which(!(is.na(clinical_data07$result))),]
###确定是否在同一天对同一项目检验两次的个体==》没有
clinical_data7$tag<-paste(clinical_data7$序号,clinical_data7$E_name,clinical_data7$手术间隔时间,sep=":")
clinical_data7$tag[duplicated(clinical_data7$tag)]##character(0)
##整理格式
clinical_data8<-cbind(clinical_data7[,1:12],clinical_data7[,c(15,13,14)])
rownames(clinical_data8)<-clinical_data8$tag
head(clinical_data8);dim(clinical_data8)#15016    15

##将两次信息科调取的检验指标表格合并
clinical_data9<-rbind(clinical_data5,clinical_data8)
head(clinical_data9);dim(clinical_data9)#99893    15
range(unique(as.numeric(clinical_data9$序号)))#  1 2082
length(unique(as.numeric(clinical_data9$序号)))## 2050
#==>从信息科处直接查询到2050条序号，此时未去除仅登记手术但是未进行手术的序号（先前提交查询时候未对手术记录单子进行这部分同一个患者的去重而标记为了不同的序号）
##总的提交信息待去重文件记录为："D:/PROJECT/新冠/2022_data/病历补充核准/9.提交查询患者信息人为去重.xlsx"
write.table(as.data.frame(clinical_data9), file="D:/PROJECT/新冠/2022_data/5.clinical_data_16month_final.txt",sep = "\t",row.names=F) 

#############################计算每个检验科每个条目比例
clinical_data3 <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/5.clinical_data_16month_final.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))

##查看检测指标数目和每个检测指标的检测次数
length(unique(clinical_data3$C_name))# 51
length(unique(clinical_data3$E_name))# 51
table(clinical_data3$E_name)
###查看每个指标检测的次数
data_test<-as.data.frame(table(clinical_data3$C_name,clinical_data3$E_name))
#data_test1<-data_test[which(data_test$Freq>0 & data_test$Freq<1000),]
data_test1<-data_test[which(data_test$Freq>0),]
data_test1
###少的仅1705个，多的有2037个

#查看每个序号检测的指标数目==>少的仅8个指标，多的有61个
length(unique(clinical_data3$序号))#2050
data02<-as.data.frame(table(clinical_data3$序号))
data02<-data02[order(data02$Freq,decreasing = F),]

dim(data02[which(data02$Freq>1),])## 2050     2

##计算检测条目潜在批次
clinical_data3$term0<-paste(clinical_data3$E_name,clinical_data3$参考值范围,sep=":")
clinical_data3$term1<-paste(clinical_data3$检验细项,clinical_data3$英文缩写,sep=":")
clinical_data3$term2<-paste(clinical_data3$检验细项,clinical_data3$英文缩写,clinical_data3$参考值范围,sep=":")
head(clinical_data3)

term_list<-distinct(clinical_data3[,c("E_name","term0","term1","term2")])
dim(distinct(clinical_data3[,c("E_name","term0","term1","term2")]))#233   4
dim(distinct(clinical_data3[,c("E_name","term0")]))# 172   2
dim(distinct(clinical_data3[,c("E_name","term1")]))# 124  2
###确定不同中英文指标名字出现次数
data03<-as.data.frame(table(clinical_data3$term2))
colnames(data03) <- c("term2","Freq")
head(data03)

##合并原始数据与频率数据
term_list2<-merge(term_list,data03,by="term2")
dim(term_list2);dim(term_list)
term_list2<-term_list2[,c("E_name","term0","term1","term2","Freq")]
term_list2<-term_list2[order(term_list2$E_name,term_list2$Freq,decreasing = T),]
write.table(as.data.frame(term_list2), file="D:/PROJECT/新冠/2022_data/检验指标的频率.txt",sep = "\t") 
####对于一些单位或者检测范围特别异常的指标进行手动标记==》6.检验指标的频率.xlsx
head(term_list2)
######将每一个指标出现频率添加至检验指标总表中
clinical_data31<-merge(clinical_data3,term_list2[,c("term2","Freq")],by="term2")
head(clinical_data31)
dim(clinical_data3);dim(clinical_data31);dim(term_list2)##99893 18   19   233   5
write.table(as.data.frame(clinical_data31), file="D:/PROJECT/新冠/2022_data/clinical_data_freq_wait.txt",sep = "\t") 
length(unique(clinical_data31$序号))##2050

##查看检验科室与检验指标之间的关系
term_list30<-distinct(clinical_data31[,c("对应就诊科室","term2")])
##查看检验科室与检验序号之间的关系
term_list31<-distinct(clinical_data31[,c("对应就诊科室","序号")])
dim(term_list31)##2203    2   ===》47个序号实际上是去了不同就诊科室的人
as.data.frame(table(term_list31$对应就诊科室))

###reading revised frequency
term_list2_new <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/6.检验指标的频率.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
head(term_list2_new);dim(term_list2_new)# 233   7

######将每一个指标出现频率重新添加至检验指标总表中
head(clinical_data3);dim(clinical_data3)
index_data4<-merge(clinical_data3,term_list2_new[,c("term2","Freq","remove_tag","tag")],by="term2")
head(index_data4);dim(index_data4)
##99825    21
##去除标记为异常的指标
index_data5<-index_data4[which(index_data4$tag.y != "去除"),]
head(index_data5);dim(index_data5)
##99813    21
length(unique(index_data5$序号))##2050

##确定未查找到的人流病历信息
#clinical_data_final <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/clinical_data_all_final.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
#head(clinical_data_final)
#length(unique(clinical_data_final$序号))##2085
#clinical_data_final$ID_full<-str_pad(clinical_data_final$ID,width =8 ,side = c("left"),pad = "0")
#length(unique(clinical_data_final$ID_full))##2028

##去除重复ID患者早前的次数：那一次应该被标记为cancel
#Submit_mark<- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/病历补充核准/提交查询患者信息人为去重.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
#head(Submit_mark)
#Submit_mark2<-Submit_mark[which(Submit_mark$标注 =="R"),]
#clinical_data_final2 <-clinical_data_final[which(!(clinical_data_final$序号 %in% Submit_mark2$查询序号)),]
#head(clinical_data_final2);dim(clinical_data_final2)##2043   35

#######还差35个提交的的计划生育手术检验信息其实未找到：包括我自己收集的样本未记录的，以及确实没有找到的
#index_data6<-index_data5[which(index_data5$序号 %in% unique(clinical_data_final2$序号)),]
#dim(index_data6);length(unique(index_data6$序号))###2008

#clinical_data_remian<-clinical_data_final2[which(!(clinical_data_final2$序号 %in% unique(index_data6$序号))),]
#head(clinical_data_remian);dim(clinical_data_remian)##30 35
#write.table(as.data.frame(clinical_data_remian), file="D:/PROJECT/新冠/2022_data/clinical_data_docimology_remain2.txt",sep = "\t") 

##确定未查找到的收样样本的病历信息
#collect_sample_final <- as.data.frame(read_excel("D:/PROJECT/新冠/病例记录20230503-终版.xlsx", sheet =2, col_names = T, col_types = NULL, na = "", skip = 0))
#head(collect_sample_final);dim(collect_sample_final)# 320  23
#collect_sample_final$ID_full<-str_pad(collect_sample_final$病历号,width =8 ,side = c("left"),pad = "0")
#length(intersect(unique(collect_sample_final$ID_full),unique(clinical_data_final$ID_full)))## 312
###提取未被记录的样本
#collect_sample_remian<-collect_sample_final[which(!(collect_sample_final$ID_full %in% unique(clinical_data_final$ID_full))),]
#dim(collect_sample_remian)#4 24
#intersect(unique(collect_sample_remian$姓名),unique(clinical_data_final$姓名))##"丁心穰" "程瑞红" "常亮亮"
#setdiff(unique(collect_sample_remian$姓名),unique(clinical_data_final$姓名))## "李雨桐"   ##########"王娟"     "黄小倩"   "纵迪迪"   "孙茹燕33" "吴萌"     
#write.table(as.data.frame(collect_sample_remian), file="D:/PROJECT/新冠/2022_data/collect_sample_docimology_remain2.txt",sep = "\t") 
#intersect(unique(collect_sample_remian$ID_full),unique(clinical_data_remian$ID_full))##0

####################################################
##病历表格整理
###1) 统一三个表格中ID_full对应人名：重点于以上三人##"丁心穰" "程瑞红" "常亮亮"
###2) 手动修整病历记录终版信息且修整病历记录中没有在提交查询版本中的收样人名称归入提交版本表格中:## "王娟"     "黄小倩"   "纵迪迪"   "孙茹燕33" "吴萌"     "李雨桐" 

##导入尽可能手动补充提交查询手术单子中未在信息科查询到的样本信息
clinical_data_final <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/病历补充核准/8.clinical_data_all_final.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_data_final)
length(unique(clinical_data_final$序号)) ##2085==》相较于信息科查询到的数据，此处补充了35条手术序号信息
clinical_data_final$ID_full<-str_pad(clinical_data_final$ID,width =8 ,side = c("left"),pad = "0")
length(unique(clinical_data_final$ID_full))##2028==》表明我提交查询的手术患者实际上有2028个

##去除重复ID患者早前的次数
Submit_mark<- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/病历补充核准/9.提交查询患者信息人为去重.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
head(Submit_mark)
Submit_mark2<-Submit_mark[which(Submit_mark$标注 =="R"),]
clinical_data_final2 <-clinical_data_final[which(!(clinical_data_final$序号 %in% Submit_mark2$查询序号)),]
head(clinical_data_final2);dim(clinical_data_final2)##2043   35
length(unique(clinical_data_final2$序号))#2043
length(unique(clinical_data_final2$ID_full))#2028

##导入收集病历的信息
collect_sample_final <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/病历补充核准/7.病例记录20230503-终版.xlsx", sheet =2, col_names = T, col_types = NULL, na = "", skip = 0))
head(collect_sample_final);dim(collect_sample_final)# 320  23
collect_sample_final$ID_full<-str_pad(collect_sample_final$病历号,width =8 ,side = c("left"),pad = "0")

#个人收集的样本中未在目前检测指标表格中被记录的样本
collect_sample_remian<-collect_sample_final[which(!(collect_sample_final$ID_full %in% unique(clinical_data_final$ID_full))),]
dim(collect_sample_remian)#1 ## "李雨桐" 
##保留个人收集样本中可查询病历号且被具体记录个体信息
collect_sample_final2<-collect_sample_final[which(collect_sample_final$ID_full %in% unique(clinical_data_final$ID_full)),]
dim(collect_sample_final2)#319  24
##去除重复记录
collect_sample_final2<-collect_sample_final2[which(!(collect_sample_final2$样本编号 %in% c("Y24","Y35","B29(X9)"))),]
dim(collect_sample_final2)#316  24
##write.table(as.data.frame(collect_sample_final2), file="D:/PROJECT/新冠/2022_data/病历补充核准/collect_sample_final2.txt",sep = "\t") 
length(unique(collect_sample_final2$ID_full))#316

#######
##生化表格整理

###将收样病历中没有在提交查询版本中的收样人和提交查询表中查找到的可补充生化指标合并
###整理生化数据格式，并且与两组补充的指标进行合并
head(index_data5)
index_data6<-index_data5[which(index_data5$序号 %in% unique(clinical_data_final2$序号)),]
dim(index_data6);length(unique(index_data6$序号))###97781    21 2008
head(index_data6)
index_data_final0<-index_data6[,c("序号","手术间隔时间","E_name","result")]
head(index_data_final0);dim(index_data_final0)#97781     4

clinical_data_add <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/病历补充核准/10.生化指标最后补充.xlsx", sheet =4, col_names = T, col_types = NULL, na = "", skip = 0))
clinical_data_add<-clinical_data_add[,-c(2:4)]
head(clinical_data_add);dim(clinical_data_add)##11 53
index_data_add <-reshape2::melt(clinical_data_add,id=colnames(clinical_data_add[,c(1,2)]),variable.name="evaluation_index",value.name="value")
head(index_data_add)

##将手动查询的检测指标与信息科查询到的指标信息进行合并
colnames(index_data_add)<-colnames(index_data_final0)<-c("sample","Day","E_name","value")
index_data_final1<-rbind(index_data_final0,index_data_add)
head(index_data_final1);dim(index_data_final1)##98342     4
range(index_data_final1$Day)##0 29 ##手术间隔时间最长达29天
table(index_data_final1$Day)

############整合病人病历信息与检测指标信息
length(unique(clinical_data_final2$序号))##2043
length(unique(index_data_final1$sample))##2019
length(intersect(unique(clinical_data_final2$序号),unique(index_data_final1$sample)))##2019
length(unique(clinical_data_final2[which(clinical_data_final2$序号 %in% unique(index_data_final1$sample)),]$ID))#2004
head(clinical_data_final2)

######查看去除5月份到12月份这个患者的信息后人数
clinical_data_final2_rm<-clinical_data_final2[which(clinical_data_final2$OPR_month %in% c(1:4)),]
length(intersect(unique(clinical_data_final2_rm$序号),unique(index_data_final1$sample)))##946
length(unique(clinical_data_final2_rm[which(clinical_data_final2_rm$序号 %in% unique(index_data_final1$sample)),]$ID))#943

#########查看疾病异常序号
clinical_data_abnormal<-clinical_data_final2[which(clinical_data_final2$排除 != "NA"),]
head(clinical_data_abnormal)
clinical_data_abnormal[,c("姓名","ID_full","排除")]
unique(clinical_data_abnormal$排除)
##查看这些异常个体多少查到检测指标
dim(clinical_data_abnormal)#202  35
index_data_abnormal<-index_data_final1[which(index_data_final1$sample %in% unique(clinical_data_abnormal$序号)),]
dim(index_data_abnormal)
head(index_data_abnormal)
length(unique(index_data_abnormal$sample))##198
clinical_data_abnormal2<-clinical_data_abnormal[which(clinical_data_abnormal$序号 %in% unique(index_data_abnormal$sample)),]
clinical_data_abnormal2[,c("姓名","ID_full","排除")]
length(unique(clinical_data_abnormal2$ID_full))## 198
######查看去除5月份到12月份这个患者的信息后人数
length(unique(clinical_data_abnormal2[which(clinical_data_abnormal2$OPR_month %in% c(1:4)),]$ID_full))## 97

###########仅保留手术信息单中正常个体
clinical_data_final3<-clinical_data_final2[which(clinical_data_final2$排除 == "NA"),]
dim(clinical_data_final3)# 1841   35
head(clinical_data_final3)
length(unique(clinical_data_final3$序号))##1841
length(unique(clinical_data_final3$ID))##1827
write.table(as.data.frame(clinical_data_final3), file="D:/PROJECT/新冠/2022_data/病历补充核准/normal_characteristics_for_2022_2023_year_final.txt",sep = "\t") 

clinical_data_final4<-clinical_data_final3[,c("序号","年龄","ID_full","pregnancy_day_final","Year_month","OPR_Year","OPR_month","OPR_day","LMP_year","LMP_month","LMP_day")]
colnames(clinical_data_final4)<-c("sample","Age","ID_full","pregnancy_day_final","Year_month","OPR_Year","OPR_month","OPR_day","LMP_year","LMP_month","LMP_day")
head(clinical_data_final4);dim(clinical_data_final4)

##提取查询到检测信息的样本条目
index_data_final2<-index_data_final1[which(index_data_final1$sample %in% unique(clinical_data_final4$sample)),]
head(index_data_final2);dim(index_data_final2)#88664     4
length(unique(index_data_final2$sample))##1821
dim(distinct(index_data_final2[,1:2]))#2095    2

##根据条目合并检测信息和病人信息
evaluate_data<-merge(index_data_final2,clinical_data_final4[,c("sample","Age","ID_full","pregnancy_day_final","Year_month","OPR_Year","OPR_month","OPR_day")])
dim(evaluate_data)# 88664    11
head(evaluate_data)
evaluate_data$pregnancy_day_final<-as.numeric(evaluate_data$pregnancy_day_final)
evaluate_data$value <-as.numeric(evaluate_data$value )
str(evaluate_data)

evaluate_data$day_detected<- evaluate_data$pregnancy_day_final-evaluate_data$Day
range(na.omit(evaluate_data$pregnancy_day_final))##10 108
range(na.omit(evaluate_data$day_detected))##9 98
head(evaluate_data);dim(evaluate_data)# 88664    12

##去除同一ID患者同一指标(指标名字+查询号)早前检测的一次
evaluate_data$index<-paste0(evaluate_data$E_name,sep=":",evaluate_data$sample)
num_index<-as.data.frame(table(evaluate_data$index))
num_index[which(num_index$Freq > 1),]
head(evaluate_data)
evaluate_data<-evaluate_data[order(evaluate_data$index,evaluate_data$Day,decreasing = F),]
evaluate_data[which(evaluate_data$index =="γ-GT:2060"),]
head(evaluate_data)
evaluate_data2 <-evaluate_data[!duplicated(evaluate_data$index),]
evaluate_data2[which(evaluate_data2$index =="γ-GT:2060"),]
head(evaluate_data2);dim(evaluate_data2) ##88624    13

length(unique(evaluate_data2$sample))#1821
length(unique(evaluate_data2$ID_full))#1807

write.table(as.data.frame(evaluate_data2), file="D:/PROJECT/新冠/2022_data/病历补充核准/clinical_index_data_for_analysis.txt",sep = "\t") 

###########各个月份数据指标分析#############提取有生化，血检指标的数据
#clinical_data_final <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/clinical_data_all_final.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
analysis_data <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/clinical_index_data_for_analysis.txt",header = T,sep = "\t")
analysis_data$Year_month2<-ifelse(analysis_data$OPR_day>15,paste0(analysis_data$Year_month,"_post"),paste0(analysis_data$Year_month,"_before"))
head(analysis_data)
##排除非关注指标
analysis_used<-analysis_data[which(!(analysis_data$E_name %in% c("HBsAg","Anti-HCV","Anti-HBc","HBeAg","HBsAb","HBeAb"))),]
analysis_used$age_group<-ifelse(analysis_used$Age>=35,"AMA","YMA")
write.table(as.data.frame(analysis_used), file="D:/PROJECT/新冠/2022_data/病历补充核准/clinical_index_data_for_analysis_used_final.txt",sep = "\t") 

################################################
##此处去除了一个1800序号的样本，其仅有病毒水平数据
 length(unique(analysis_used$sample))#1820
 length(unique(analysis_used$ID_full))#1806
head(collect_sample_final2)
length(intersect(unique(analysis_used$ID_full),unique(collect_sample_final2$ID_full)))#288
collect_sample_infor<-collect_sample_final2[which(collect_sample_final2$ID_full %in% intersect(unique(analysis_used$ID_full),unique(collect_sample_final2$ID_full))),]
dim(collect_sample_infor)#288  24

###去除5到12月份之后的数目
length(unique(analysis_used[which(analysis_used$OPR_month %in% 1:4),]$sample))#848
length(unique(analysis_used[which(analysis_used$OPR_month %in% 1:4),]$ID_full))#846

###绘制具体人数分布
sample_num<-distinct(analysis_used[,c("sample","Age","Year_month","Year_month2")])
dim(sample_num)#1820   4
head(sample_num)

##计算每月人数
month_order<-c("20221","20222","20223","20224","20225","20226","20227","20228","20229","202210","202211","202212","20231","20232","20233","20234")
month_order2<-c("20221_before","20221_post","20222_before","20222_post","20223_before","20223_post","20224_before","20224_post","20225_before","20225_post","20226_before","20226_post","20227_before","20227_post","20228_before","20228_post","20229_before","20229_post","202210_before","202210_post","202211_before","202211_post","202212_before","202212_post","20231_before","20231_post","20232_before","20232_post","20233_before","20233_post","20234_before","20234_post")
sample_num$age_group<-ifelse(sample_num$Age>=35,"AMA","YMA")
sample_num$Year_month <- factor(sample_num$Year_month,levels=month_order)
sample_num$Year_month2 <- factor(sample_num$Year_month2,levels=month_order2)
table(sample_num$Year_month)
#20221  20222  20223  20224  20225  20226  20227  20228  20229 202210 202211 202212  20231  20232  20233  20234 
#  74     80    135    108    117    139    105    120    129     97    135    130     76    111    132    132 

table(sample_num$Year_month2)
#20221_before    20221_post  20222_before    20222_post  20223_before    20223_post  20224_before    20224_post  20225_before    20225_post  20226_before    20226_post 
#       33            41            33            47            65            70            52            56            49            68            82            57 
#20227_before    20227_post  20228_before    20228_post  20229_before    20229_post 202210_before   202210_post 202211_before   202211_post 202212_before   202212_post 
#       54            51            54            66            60            69            38            59            72            64            70            60 
#20231_before    20231_post  20232_before    20232_post  20233_before    20233_post  20234_before    20234_post 
#       46            30            65            46            64            69            69            63 

number_plot0<-ggplot(data=sample_num, mapping=aes(x=Year_month))+
  geom_bar(aes(fill = Year_month),stat="count",width=0.5)+coord_polar(start = 0) +
  scale_fill_manual(values=ppCor_all2)+
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_stack(0.9))+
  theme_minimal()+theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.title = element_blank(),legend.position = "none",
    axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 0),legend.title = element_text(size = 9))

number_plot00<-ggplot(data=sample_num, mapping=aes(x=Year_month,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='stack')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_stack(0.5))+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
y_max <- max(aggregate(sample~Year_month+age_group,data=sample_num,length)$sample)
number_plot002<-ggplot(data=sample_num, mapping=aes(x=Year_month,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='dodge')+
  scale_fill_manual(values=c('#999999','#E69F00'))+  ylim(0,y_max+5)+
  geom_text(stat='count',aes(label=..count..), color="black", size=3.5,position=position_dodge(0.5),vjust=-0.5)+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/number_all_sample_0.pdf",number_plot0,width=6, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/number_all_sample_1.pdf",number_plot00,width=8, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/number_all_sample_2.pdf",number_plot002,width=12, height=6)

number_plot01<-ggplot(data=sample_num, mapping=aes(x=Year_month2))+
  geom_bar(aes(fill = Year_month),stat="count",width=0.5)+coord_polar(start = 0) +
  scale_fill_manual(values=ppCor_all2)+
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_stack(0.9))+
  theme_minimal()+theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.title = element_blank(),legend.position = "none",
                        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 0),legend.title = element_text(size = 9))

number_plot001<-ggplot(data=sample_num, mapping=aes(x=Year_month2,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='stack')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_stack(0.5))+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
y_max <- max(aggregate(sample~Year_month2+age_group,data=sample_num,length)$sample)
number_plot003<-ggplot(data=sample_num, mapping=aes(x=Year_month2,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='dodge')+
  scale_fill_manual(values=c('#999999','#E69F00'))+  ylim(0,y_max+5)+
  geom_text(stat='count',aes(label=..count..), color="black", size=3.5,position=position_dodge(0.5),vjust=-0.5)+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/number_all_sample_3.pdf",number_plot01,width=8, height=8)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/number_all_sample_4.pdf",number_plot001,width=12, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/number_all_sample_5.pdf",number_plot003,width=16, height=6)

##整月##plot pecentage sperated by YMA and AMA
cell_prop0<-as.data.frame(table(sample_num$Year_month,sample_num$age_group))
cell_prop1 <- reshape2::dcast(cell_prop0,Var1 ~ Var2, value.var = "Freq")
cell_prop1$sum<-cell_prop1$AMA+cell_prop1$YMA
cell_prop1$AMA<-cell_prop1$AMA/cell_prop1$sum*100
cell_prop1$YMA<-cell_prop1$YMA/cell_prop1$sum*100

cell_prop2 <- reshape2::melt(cell_prop1[,-ncol(cell_prop1)],variable.name="Group",value.name = "value")
colnames(cell_prop2)<-c( "Year_month","Group","Ratio")
cell_prop2$Year_month<-factor(cell_prop2$Year_month,levels=month_order)

plot_cell_prop1<-ggplot(data=cell_prop2, mapping=aes(x= Year_month,y=Ratio,fill=Group))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+scale_fill_manual(name="age_group",values=c('#999999','#E69F00'))+
  labs(x="Year_month",y="Ratio(%)",title="The patient number for each mounth")+
  geom_text(aes(label=round(Ratio,2)), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
plot_cell_prop1

###半月##plot pecentage sperated by YMA and AMA
cell_prop0<-as.data.frame(table(sample_num$Year_month2,sample_num$age_group))
cell_prop1 <- reshape2::dcast(cell_prop0,Var1 ~ Var2, value.var = "Freq")
cell_prop1$sum<-cell_prop1$AMA+cell_prop1$YMA
cell_prop1$AMA<-cell_prop1$AMA/cell_prop1$sum*100
cell_prop1$YMA<-cell_prop1$YMA/cell_prop1$sum*100

cell_prop2 <- reshape2::melt(cell_prop1[,-ncol(cell_prop1)],variable.name="Group",value.name = "value")
colnames(cell_prop2)<-c( "Year_month","Group","Ratio")
cell_prop2$Year_month<-factor(cell_prop2$Year_month,levels=month_order2)

plot_cell_prop2<-ggplot(data=cell_prop2, mapping=aes(x= Year_month,y=Ratio,fill=Group))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+scale_fill_manual(name="age_group",values=c('#999999','#E69F00'))+
  labs(x="Year_month",y="Ratio(%)",title="The patient number for each mounth")+
  geom_text(aes(label=round(Ratio,2)), position = position_stack(reverse =F,vjust=0.5),size=4)+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
plot_cell_prop2
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/number_used_sample_6.pdf",plot_cell_prop1,width=8, height=7)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/number_used_sample_7.pdf",plot_cell_prop2,width=15, height=7)

########对目标指标进行比较
head(analysis_used);dim(analysis_used)##78529    14
length(unique(analysis_used$sample))#1820
head(analysis_used)
analysis_used$age_group<-ifelse(analysis_used$Age>=35,"AMA","YMA")

length(unique(analysis_used$sample))#1820
head(analysis_used)

##计算每月人数
month_order<-c("20221","20222","20223","20224","20225","20226","20227","20228","20229","202210","202211","202212","20231","20232","20233","20234")
month_order2<-c("20221_before","20221_post","20222_before","20222_post","20223_before","20223_post","20224_before","20224_post","20225_before","20225_post","20226_before","20226_post","20227_before","20227_post","20228_before","20228_post","20229_before","20229_post","202210_before","202210_post","202211_before","202211_post","202212_before","202212_post","20231_before","20231_post","20232_before","20232_post","20233_before","20233_post","20234_before","20234_post")

analysis_used$Year_month <- factor(analysis_used$Year_month,levels=month_order)
analysis_used$Year_month2 <- factor(analysis_used$Year_month2,levels=month_order2)

index_order<-c("WBC","RBC","HGB","HCT","MCV","MCH","MCHC","PLT","RDW_CV",
"PCT","LYMPH_pct","LYMPH_av","Neut_pct","Neut_av","EO_pct","EO_av","BASO_pct","BASO_av","MONO_pct","MONO_av",
"MPV","PDW","P-LCR","ALT","AST","TBIL","γ-GT","CK","CK-MB","BUN_Urea","Cr","UA","T-CHO","TG",
"HDL-C","LDL-C","GLU","PT","A","INR","Fib","APTT","APTT_R","TT","R")
length(index_order)#45
analysis_used$E_name <- factor(analysis_used$E_name,levels=index_order)

##指标的月份比较
head(analysis_used)
analysis_used$value<-as.numeric(analysis_used$value)
range(na.omit(analysis_used$value))# 0 1292
stat_data<-compare_means(value~Year_month, data=analysis_used, group.by = "E_name")
stat_data[which(stat_data$p.signif != "ns"),]
month_stat_boxplot<-ggboxplot(analysis_used, x="Year_month", y="value", color="Year_month") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "20222",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ E_name, scales = "free",ncol =7)
month_stat_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_target_index_all_used_sample_compare_stat_box1.pdf",month_stat_boxplot,width=40, height=28)

month_stat_boxplot2<-ggboxplot(analysis_used, x="Year_month2", y="value", color="Year_month2") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "20222_before",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 5,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ E_name, scales = "free",ncol =7)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_target_index_all_used_sample_compare_stat_box2.pdf",month_stat_boxplot2,width=50, height=28,limitsize = F)

head(analysis_used)
analysis_used$age_group <-factor(analysis_used$age_group,levels=c("AMA","YMA"))
month_stat_boxplot3<-ggboxplot(analysis_used,x="Year_month", y="value",fill="age_group",color="Year_month") +  
  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
 scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "20222",hide.ns = TRUE)+
  stat_summary(aes(group= age_group), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.8)) +
  #stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =7)
month_stat_boxplot3
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_target_index_all_used_sample_compare_stat_box1_age_split.pdf",month_stat_boxplot3,width=46, height=28)

month_stat_boxplot4<-ggboxplot(analysis_used, x="Year_month2", y="value", fill="age_group",color="Year_month2") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "20222_before",hide.ns = TRUE)+
  stat_summary(aes(group= age_group),fun=mean,geom="point", shape=23,size=3,fill="white",position=position_dodge(0.8)) +
  theme(axis.text.x = element_text(size = 5,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ E_name, scales = "free",ncol =7)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_target_index_all_used_sample_compare_stat_box2_age_split.pdf",month_stat_boxplot4,width=56, height=28,limitsize = F)

##########进行年龄组比较
month_stat_boxplot05<-ggplot(analysis_used, aes(x=age_group,y=value,fill=age_group))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= age_group), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =10)
month_stat_boxplot05

month_stat_boxplot5<-ggplot(analysis_used,aes(x=Year_month,y=value,fill=age_group))+
   geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE)+
  stat_summary(aes(group= age_group), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =10)
month_stat_boxplot5

ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_target_index_all_used_sample_compare_stat_box_age_group1.pdf",month_stat_boxplot05,width=30, height=17)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_target_index_all_used_sample_compare_stat_box_age_group2.pdf",month_stat_boxplot5,width=56, height=28,limitsize = F)


month_stat_boxplot6<-ggplot(analysis_used,aes(x=Year_month2,y=value,fill=age_group))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE)+
  stat_summary(aes(group= age_group), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =10)
month_stat_boxplot6

ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_target_index_all_used_sample_compare_stat_box_age_group3.pdf",month_stat_boxplot6,width=56, height=28,limitsize = F)


##绘制1到4月份中2022年与2023年各个因子组间差异
##指标的月份比较
head(analysis_used)
analysis_used$value<-as.numeric(analysis_used$value)
range(na.omit(analysis_used$value))# 0 1292

month_select<-c("20221","20231","20222","20232","20223","20233","20224","20234")
analysis_select<-analysis_used[which(analysis_used$Year_month %in%month_select),]
head(analysis_select)
analysis_select$OPR_month <- factor(analysis_select$OPR_month,levels=1:4)
analysis_select$OPR_Year <- factor(analysis_select$OPR_Year,levels=c(2022,2023))

head(analysis_select)
analysis_select$E_name_mounth<-paste0(analysis_select$E_name,":",analysis_select$OPR_month)
stat_data<-compare_means(value~OPR_Year,data=analysis_select,group.by = "E_name_mounth")
stat_data[which(stat_data$p.signif != "ns"),]
# 35 x 9
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Clinical_target_index_all_used_sample_compare_four_months_year_compare.txt",row.names=T, col.names=T) 

month_stat_boxplot08<-ggplot(analysis_select, aes(x=OPR_month,y=value,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =7)
month_stat_boxplot08
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Clinical_target_index_all_used_sample_compare_stat_box_four_months_year_compare.pdf",month_stat_boxplot08,width=30, height=30)

stat_data_sig<-stat_data[which(stat_data$p.signif != "ns"),]
stat_data_sig$E_name <- as.character(unlist(lapply(strsplit(as.character(stat_data_sig$E_name_mounth),":"), function(x) x[1])))
analysis_sig<-analysis_select[which(analysis_select$E_name %in% unique(stat_data_sig$E_name)),]
length(unique(stat_data_sig$E_name))##23

month_stat_boxplot18<-ggplot(analysis_sig, aes(x=OPR_month,y=value,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =5)
month_stat_boxplot18
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Clinical_target_index_significant_used_sample_compare_stat_box_four_months_year_compare.pdf",month_stat_boxplot18,width=20, height=20)


##age comparison
stat_data<-compare_means(Age~OPR_Year,data=analysis_select,group.by = "E_name_mounth")
stat_data[which(stat_data$p.signif != "ns"),]
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Age_all_used_sample_compare_stat_box_four_months_year_compare.txt",row.names=T, col.names=T) 

month_stat_boxplot09<-ggplot(analysis_select, aes(x=OPR_month,y=Age,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =7)
month_stat_boxplot09
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Age_all_used_sample_compare_stat_box_four_months_year_compare.pdf",month_stat_boxplot09,width=30, height=30)

##半月比较
analysis_select$OPR_month2<-ifelse(analysis_select$OPR_day>15,paste0(analysis_select$OPR_month,"_post"),paste0(analysis_select$OPR_month,"_before"))
group_set<-c("1_before","1_post","2_before","2_post","3_before","3_post","4_before","4_post")
analysis_select$OPR_month2 <- factor(analysis_select$OPR_month2,levels=group_set)
analysis_select$E_name_mounth<-paste0(analysis_select$E_name,":",analysis_select$OPR_month2)
stat_data<-compare_means(value~OPR_Year,data=analysis_select,group.by = "E_name_mounth")
stat_data[which(stat_data$p.signif != "ns"),]
write.table(as.data.frame(stat_data), file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Clinical_target_index_all_used_sample_compare_four_months_half_year_compare.txt",row.names=T, col.names=T) 

month_stat_boxplot10<-ggplot(analysis_select, aes(x=OPR_month2,y=value,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =7)
month_stat_boxplot10
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Clinical_target_index_all_used_sample_compare_stat_box_four_months_before_post_year_compare.pdf",month_stat_boxplot10,width=30, height=30)
##for significant

stat_data_sig<-stat_data[which(stat_data$p.signif != "ns"),]
stat_data_sig$E_name <- as.character(unlist(lapply(strsplit(as.character(stat_data_sig$E_name_mounth),":"), function(x) x[1])))
analysis_sig<-analysis_select[which(analysis_select$E_name %in% unique(stat_data_sig$E_name)),]
length(unique(analysis_sig$E_name))##28
month_stat_boxplot20<-ggplot(analysis_sig, aes(x=OPR_month2,y=value,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =6)
month_stat_boxplot20
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Clinical_target_index_significant_used_sample_compare_stat_box_four_months_before_post_year_compare.pdf",month_stat_boxplot10,width=30, height=25)


##compare age
month_stat_boxplot11<-ggplot(analysis_select, aes(x=OPR_month2,y=Age,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =7)
month_stat_boxplot11
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Age_all_used_sample_compare_stat_box_four_months_before_post_year_compare.pdf",month_stat_boxplot11,width=30, height=30)

###不区分各个因子进行两组年龄比较
head(analysis_select)
sample_num<-distinct(analysis_select[,c("sample","Age","OPR_month","OPR_month2","OPR_Year")])
month_stat_boxplot0<-ggplot(sample_num, aes(x=OPR_month,y=Age,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))
month_stat_boxplot0

month_stat_boxplot1<-ggplot(sample_num, aes(x=OPR_month2,y=Age,fill=OPR_Year))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= OPR_Year), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))
month_stat_boxplot1
month_stat_boxplot11<-grid.arrange(month_stat_boxplot0,month_stat_boxplot1,ncol=1)

ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/four_mounth/Age_all_used_sample_compare_stat_box_four_and_half_months_compare.pdf",month_stat_boxplot11,width=8, height=8)


##################绘制目标因子
##凝血相关指标
#c(凝血酶原时间PT,凝血酶原活动度,国际标准比率,纤维蛋白原,活化部分凝血活酶,APTT比率,凝血酶时间,TT比率)
NX_index<-c("PT","A","INR","Fib","APTT","APTT_R","TT","R")
analysis_NX<-analysis_used[which(analysis_used$E_name %in% NX_index),]

analysis_NX$Year_month <- factor(analysis_NX$Year_month,levels=month_order)
analysis_NX$Year_month2 <- factor(analysis_NX$Year_month2,levels=month_order2)

##进行月份间比较
month_stat_boxplot<-ggboxplot(analysis_NX, x="Year_month", y="value", color="Year_month") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "20222",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ E_name, scales = "free",ncol =4)
month_stat_boxplot
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_NX_index_all_used_sample_compare_stat_box1.pdf",month_stat_boxplot,width=25, height=10,limitsize = F)

month_stat_boxplot2<-ggboxplot(analysis_NX, x="Year_month2", y="value", color="Year_month2") +  scale_color_manual(values=ppCor_all)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "20222_before",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 5,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ E_name, scales = "free",ncol =4)
month_stat_boxplot2
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_NX_index_all_used_sample_compare_stat_box2.pdf",month_stat_boxplot2,width=35, height=10,limitsize = F)

##########进行年龄组比较
month_stat_boxplot07<-ggplot(analysis_NX, aes(x=age_group,y=value,fill=age_group))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= age_group), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =4)
month_stat_boxplot07

month_stat_boxplot7<-ggplot(analysis_NX,aes(x=Year_month,y=value,fill=age_group))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE)+
  stat_summary(aes(group= age_group), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =4)
month_stat_boxplot7

ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_NX_index_all_used_sample_compare_stat_box_age_group1.pdf",month_stat_boxplot07,width=10, height=6)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_NX_index_all_used_sample_compare_stat_box_age_group2.pdf",month_stat_boxplot7,width=50, height=10,limitsize = F)

month_stat_boxplot8<-ggplot(analysis_NX,aes(x=Year_month2,y=value,fill=age_group))+
  geom_boxplot(position=position_dodge(),width=0.5)+#geom_jitter(width = 0.2,color="grey",alpha=0.7)+ 
  scale_color_manual(values=ppCor_all)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE)+
  stat_summary(aes(group= age_group), fun = "mean", geom = "point",shape=23,size=1,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),
        legend.title = element_text(size = 9))+ facet_wrap(~ E_name, scales = "free",ncol =4)
month_stat_boxplot8

ggsave("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/Clinical_NX_index_all_used_sample_compare_stat_box_age_group3.pdf",month_stat_boxplot8,width=56, height=10,limitsize = F)

###########对指标随着新冠感染康复日期变化的动态绘制
head(collect_sample_final2);dim(collect_sample_final2)##316  24
collect_sample_final20<-collect_sample_final2[which(collect_sample_final2$异常与否 == "OK"),]
dim(collect_sample_final20)#287  24
#排除不确定感染还是未感染的
collect_sample_final3<-collect_sample_final20[which(collect_sample_final20$初定新冠阳性日期 != "na"),]
dim(collect_sample_final3)#251  24
head(collect_sample_final3)
collect_sample_info<-collect_sample_final3[,c("ID_full","末次月经","初定新冠阳性日期","末次月经相较于手术日_T","手术日相对于感染的时间_T","末次月经相较于感染日_T" )]
head(collect_sample_info)
colnames(collect_sample_info)<-c("ID_full","LMP","Infect_state","LMP_operate","operate_infect","LMP_infect")
collect_sample_info$Infect_state <-ifelse(collect_sample_info$Infect_state == "no","no","yes")
head(collect_sample_info);dim(collect_sample_info)

###排除确认是否在总的患者信息记录表格中存在以及疾病状态是否对应
head(analysis_used);head(clinical_data_final4)
length(intersect(unique(collect_sample_info$ID_full),unique(clinical_data_final4$ID_full)))##251==>疾病状态对应OK
length(intersect(unique(collect_sample_info$ID_full),unique(analysis_used$ID_full)))##249 ==>两个患者的临床检测指标信息没查到

setdiff(unique(collect_sample_info$ID_full),unique(analysis_used$ID_full))
##"19894361"  ：A89 未感染患者 
##"22751341"  : C93 仅查到尿蛋白信息，故总的检测指标信息中未纳入

#将这两个患者也排除
collect_sample_info_final<-collect_sample_info[which(collect_sample_info$ID_full %in% intersect(unique(collect_sample_info$ID_full),unique(analysis_used$ID_full))),]
dim(collect_sample_info_final)#249   6
write.table(as.data.frame(collect_sample_info_final), file="D:/PROJECT/新冠/2022_data/病历补充核准/collect_sample_infor.txt",sep = "\t") 
