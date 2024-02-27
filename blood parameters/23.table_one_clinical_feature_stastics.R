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

#####################
analysis_used <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/clinical_index_data_for_analysis_used_final.txt",header = T,sep = "\t")
head(analysis_used)
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
clinical_data<-distinct(analysis_select[,c("sample","ID_full","Age","Year_month")])
clinical_data[which(clinical_data$ID_full %in% clinical_data$ID_full[duplicated(clinical_data$ID_full)]),]
analysis_select[which(analysis_select$sample=="1694"),]$ID_full<-"22361716_B"
clinical_data<-distinct(analysis_select[,c("sample","ID_full","Year_month","OPR_Year","Age","pregnancy_day_final","age_group")])


head(clinical_data)
dim(clinical_data)#848   7
length(unique(clinical_data$ID_full))##847
table(clinical_data$OPR_Year)
#2022 2023 
# 397  451 
clinical_data$GS_week2<-floor(clinical_data$pregnancy_day_final/7)
clinical_data$GS_week <-clinical_data$pregnancy_day_final/7

##首先可以看一下原始数据的基线资料表，用的是tableone这个包，它能计算SMD（后面会介绍这个SMD的作用）
library(tableone)
head(clinical_data);dim(clinical_data)#848   8
table1 <- CreateTableOne(vars = c('Age', 'GS_week','age_group'),data = clinical_data, smd=TRUE)
table1 <- print(table1,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table1

table2 <- CreateTableOne(vars = c('Age', 'GS_week','age_group'),data = clinical_data, strata = 'OPR_Year', smd=TRUE)#ttest、wilcox.test、chisq.test
table2 <- print(table2,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table2
###############进行分布评估
library(rstatix)
library(car)
##GS_week
tem_data <-clinical_data[,c("OPR_Year","GS_week")]
head(tem_data)
shapiro_data<-tem_data %>% group_by(OPR_Year) %>%  shapiro_test(GS_week)
head(shapiro_data)
leveneTest_result <- leveneTest(GS_week~OPR_Year, data = tem_data)
print(leveneTest_result)

compare_means(GS_week~OPR_Year, data=tem_data,method = "wilcox.test")
compare_means(GS_week~OPR_Year, data=tem_data,method = "t.test")
#.y.      group1 group2         p    p.adj p.format p.signif method  
# GS_week 2022   2023   0.00000224 0.0000022 2.2e-06  ****     Wilcoxon
# GS_week 2022   2023   0.0000106 0.000011 1.1e-05  ****     T-test

##Age
tem_data <-clinical_data[,c("OPR_Year","Age")]
head(tem_data)
shapiro_data<-tem_data %>% group_by(OPR_Year) %>%  shapiro_test(Age)
head(shapiro_data)
leveneTest_result <- leveneTest(Age~OPR_Year, data = tem_data)
print(leveneTest_result)

compare_means(Age~OPR_Year, data=tem_data,method = "wilcox.test")
compare_means(Age~OPR_Year, data=tem_data,method = "t.test")
#.y.      group1 group2         p    p.adj p.format p.signif method  
# Age       2022   2023   0.325  0.32 0.32     ns       Wilcoxon
# Age       2022   2023   0.240  0.24 0.24     ns       T-test

##Age_state
tem_data <-clinical_data[,c("OPR_Year","age_group")]
# 创建列联表
table <- table(tem_data$OPR_Year, tem_data$age_group)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)


########################
clinical_data$ID_full[duplicated(clinical_data$ID_full)]
clinical_data[which(clinical_data$ID_full=="22361716"),]
clinical_data_final <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/normal_characteristics_for_2022_2023_year_final.txt",header = T,sep = "\t")
head(clinical_data_final)
clinical_data_final[which(clinical_data_final$ID_full=="22361716"),]
clinical_data_final2<-clinical_data_final[which(clinical_data_final$ID_full %in% unique(clinical_data$ID_full)),]
dim(clinical_data_final2)#859  35
clinical_data_final3<-clinical_data_final2[which(clinical_data_final2$'序号' %in% unique(clinical_data$sample)),]
dim(clinical_data_final3)#848  35
clinical_data_final3[which(clinical_data_final3$ID_full=="22361716"),]
table(clinical_data_final3$'产次')
write.table(as.data.frame(clinical_data_final3), file="D:/PROJECT/新冠/2022_data/病历补充核准/normal_characteristics_for_848_cycle_final.txt",row.names=T, col.names=T) 

############孕产情况评估############
GP_data_final <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/normal_characteristics_for_848_cycle_final3.txt",header = T,sep = "\t")
head(GP_data_final)

GP_data_final$multigravida<-ifelse(GP_data_final$G_state !="G0",ifelse(is.na(GP_data_final$G_state),NA,"Yes"),"No")
GP_data_final$multiparous<-ifelse(GP_data_final$P_state !="P0",ifelse(is.na(GP_data_final$P_state),NA,"Yes"),"No")
GP_data_final2<-GP_data_final[which(!(is.na(GP_data_final$G_P_state))),]

table(GP_data_final2$multigravida)
table(GP_data_final2$multiparous)

table4 <- CreateTableOne(vars = c('multigravida', 'multiparous','Abortion_state'),data = GP_data_final2,smd=TRUE)#ttest、wilcox.test、chisq.test
table4 <- print(table4,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table4

table3 <- CreateTableOne(vars = c('multigravida', 'multiparous','Abortion_state'),data = GP_data_final2, strata = 'OPR_Year', smd=TRUE)#ttest、wilcox.test、chisq.test
table3 <- print(table3,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table3

##multigravida
tem_data <-GP_data_final2[,c("OPR_Year","multigravida")]
# 创建列联表
table <- table(tem_data$OPR_Year, tem_data$multigravida)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##Fisher’s exact test
result_fisher <- fisher.test(table)
print(result_fisher)

##multiparous
tem_data <-GP_data_final2[,c("OPR_Year","multiparous")]
# 创建列联表
table <- table(tem_data$OPR_Year, tem_data$multiparous)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##chisq exact test
result_chi2 <- chisq.test(table)
print(result_chi2)

##Abortion_state
tem_data <-GP_data_final2[,c("OPR_Year","Abortion_state")]
# 创建列联表
table <- table(tem_data$OPR_Year, tem_data$Abortion_state)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##Fisher’s exact test
result_fisher <- fisher.test(table)
print(result_fisher)


####################################
####for comparison the infection and no infection############
##读取被记录的康复信息完整的样本内容
collect_sample_info_final <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/collect_sample_infor.txt",header = T,sep = "\t")
collect_sample_info_final$ID_full<-str_pad(collect_sample_info_final$ID,width =8 ,side = c("left"),pad = "0")
head(collect_sample_info_final);dim(collect_sample_info_final)#249   6
length(unique(collect_sample_info_final$ID_full))#257

dim(collect_sample_info_final)#249   6
table(collect_sample_info_final$Infect_state)
# no yes 
# 17 232 
collect_sample_info_final2<-collect_sample_info_final[which(!(is.na(collect_sample_info_final$LMP_infect))),]
table5 <- CreateTableOne(vars = c('LMP_infect'),data = collect_sample_info_final2,smd=TRUE)#ttest、wilcox.test、chisq.test
table5 <- print(table5,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table5
#LMP_infect (mean (SD)) ""    "49.07 (25.66)";n="231"     

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
head(analysis_collect)
length(unique(analysis_collect$sample))#249
clinical_data<-distinct(analysis_collect[,c("sample","ID_full","Year_month","OPR_Year","Age","pregnancy_day_final","age_group")])

head(clinical_data)
dim(clinical_data)#249   7
length(unique(clinical_data$ID_full))##847
table(clinical_data$OPR_Year)
#2023 
#249
clinical_data$GS_week2<-floor(clinical_data$pregnancy_day_final/7)
clinical_data$GS_week <-clinical_data$pregnancy_day_final/7

clinical_data2<-merge(collect_sample_info_final[,c("ID_full","Infect_state")],clinical_data,by="ID_full")
dim(clinical_data2)
head(clinical_data2)

##首先可以看一下原始数据的基线资料表，用的是tableone这个包，它能计算SMD（后面会介绍这个SMD的作用）
library(tableone)
head(clinical_data);dim(clinical_data)#848   8
table1 <- CreateTableOne(vars = c('Age','GS_week','age_group'),data = clinical_data2, smd=TRUE)
table1 <- print(table1,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table1

table2 <- CreateTableOne(vars = c('Age', 'GS_week','age_group'),data = clinical_data2, strata = 'Infect_state', smd=TRUE)#ttest、wilcox.test、chisq.test
table2 <- print(table2,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table2
###############进行分布评估
library(rstatix)
library(car)
##GS_week
tem_data <-clinical_data2[,c("Infect_state","GS_week")]
head(tem_data)
shapiro_data<-tem_data %>% group_by(Infect_state) %>%  shapiro_test(GS_week)
head(shapiro_data)
leveneTest_result <- leveneTest(GS_week~Infect_state, data = tem_data)
print(leveneTest_result)

compare_means(GS_week~Infect_state, data=tem_data,method = "wilcox.test")
compare_means(GS_week~Infect_state, data=tem_data,method = "t.test")
#.y.      group1 group2         p    p.adj p.format p.signif method  
# GS_week yes    no     0.825  0.83 0.83     ns       Wilcoxon
#GS_week yes    no     0.654  0.65 0.65     ns       T-test

##Age
tem_data <-clinical_data2[,c("Infect_state","Age")]
head(tem_data)
shapiro_data<-tem_data %>% group_by(Infect_state) %>%  shapiro_test(Age)
head(shapiro_data)
leveneTest_result <- leveneTest(Age~Infect_state, data = tem_data)
print(leveneTest_result)

compare_means(Age~Infect_state, data=tem_data,method = "wilcox.test")
compare_means(Age~Infect_state, data=tem_data,method = "t.test")
#.y.      group1 group2         p    p.adj p.format p.signif method  
# Age        yes    no     0.00742 0.0074 0.0074   **       Wilcoxon
# Age        yes    no      0.0124 0.012 0.012    *        T-test

##Age_state
tem_data <-na.omit(clinical_data2[,c("Infect_state","age_group")])
# 创建列联表
table <- table(tem_data$Infect_state, tem_data$age_group)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)

result_chi2 <- chisq.test(table)
result_chi2

############孕产情况评估############
GP_data_final <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/normal_characteristics_for_848_cycle_final3.txt",header = T,sep = "\t")
head(GP_data_final)

GP_data_final$multigravida<-ifelse(GP_data_final$G_state !="G0",ifelse(is.na(GP_data_final$G_state),NA,"Yes"),"No")
GP_data_final$multiparous<-ifelse(GP_data_final$P_state !="P0",ifelse(is.na(GP_data_final$P_state),NA,"Yes"),"No")
GP_data_final2<-GP_data_final[which(!(is.na(GP_data_final$G_P_state))),]
head(GP_data_final2);dim(GP_data_final2)

table(GP_data_final2$multigravida)
table(GP_data_final2$multiparous)

head(clinical_data2)
GP_data_final3<-merge(clinical_data2[,c("sample","Infect_state")],GP_data_final2,by="sample")
dim(GP_data_final3)
head(GP_data_final3)

table4 <- CreateTableOne(vars = c('multigravida', 'multiparous','Abortion_state'),data = GP_data_final3,smd=TRUE)#ttest、wilcox.test、chisq.test
table4 <- print(table4,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table4

table3 <- CreateTableOne(vars = c('multigravida', 'multiparous','Abortion_state'),data = GP_data_final3, strata = 'Infect_state', smd=TRUE)#ttest、wilcox.test、chisq.test
table3 <- print(table3,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table3

##multigravida
tem_data <-GP_data_final3[,c("Infect_state","multigravida")]
# 创建列联表
table <- table(tem_data$Infect_state, tem_data$multigravida)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##Fisher’s exact test
result_fisher <- fisher.test(table)
print(result_fisher)

##multiparous
tem_data <-GP_data_final3[,c("Infect_state","multiparous")]
# 创建列联表
table <- table(tem_data$Infect_state, tem_data$multiparous)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##chisq exact test
result_fisher <- fisher.test(table)
print(result_fisher)

##Abortion_state
tem_data <-GP_data_final3[,c("Infect_state","Abortion_state")]
# 创建列联表
table <- table(tem_data$Infect_state, tem_data$Abortion_state)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##Fisher’s exact test
result_chi2 <- chisq.test(table)
print(result_chi2)

##########################################################################
##########for transcriptome sample###############
###读取最终用于计算的转录组信息
#RNA_meta <- read.table(file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/Final_fielt_Decidua_Vill_analysis_metadata.txt",sep="\t",header = T)
RNA_meta <- read.table(file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/Final_fielt_Decidua_Vill_analysis_metadata_genderadd.txt",sep="\t",header = T)
head(RNA_meta);dim(RNA_meta)#245  13
range(na.omit(RNA_meta$LMP_infect))# -2 110
RNA_meta$age_group<-ifelse(RNA_meta$Age>=35,"AMA","YMA")

RNA_sample_info_final <- read.table(file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/Final_fielt_Decidua_Vill_analysis_metadata_clinic.txt",header = T,sep = "\t")
head(RNA_sample_info_final);dim(RNA_sample_info_final)##132   9

RNA_meta_information<-distinct(RNA_meta[,c("Age","LMP_infect","LMP_operate","infect_state","gender","age_group")])
RNA_meta_information$GS_week2<-floor(RNA_meta_information$LMP_operate/7)
RNA_meta_information$GS_week <-RNA_meta_information$LMP_operate/7
RNA_meta_information1<-RNA_meta_information[which(!(is.na(RNA_meta_information$LMP_infect))),]
dim(RNA_meta_information)# 132   8
dim(RNA_meta_information1)# 120  8

head(RNA_meta_information)

table7 <- CreateTableOne(vars = c('LMP_infect'),data = RNA_meta_information1,smd=TRUE)#ttest、wilcox.test、chisq.test
table7 <- print(table7,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table7
#LMP_infect (mean (SD)) ""   "49.50 (25.79)";n="120"     

##首先可以看一下原始数据的基线资料表，用的是tableone这个包，它能计算SMD（后面会介绍这个SMD的作用）
head(RNA_meta_information);dim(RNA_meta_information)#848   8
table1 <- CreateTableOne(vars = c('Age','GS_week','age_group',"gender"),data = RNA_meta_information, smd=TRUE)
table1 <- print(table1,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table1

table2 <- CreateTableOne(vars = c('Age','GS_week','age_group',"gender"),data = RNA_meta_information, strata = 'infect_state', smd=TRUE)#ttest、wilcox.test、chisq.test
table2 <- print(table2,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table2

###############进行分布评估
library(rstatix)
library(car)
##GS_week
tem_data <-RNA_meta_information[,c("infect_state","GS_week")]
head(tem_data)
shapiro_data<-tem_data %>% group_by(infect_state) %>%  shapiro_test(GS_week)
head(shapiro_data)
leveneTest_result <- leveneTest(GS_week~infect_state, data = tem_data)
print(leveneTest_result)

compare_means(GS_week~infect_state, data=tem_data,method = "wilcox.test")
compare_means(GS_week~infect_state, data=tem_data,method = "t.test")
#.y.      group1 group2         p    p.adj p.format p.signif method  
# GS_week Infect No     0.956  0.96 0.96     ns       Wilcoxon
# GS_week Infect No     0.510  0.51 0.51     ns       T-test

##Age
tem_data <-RNA_meta_information[,c("infect_state","Age")]
head(tem_data)
shapiro_data<-tem_data %>% group_by(infect_state) %>%  shapiro_test(Age)
head(shapiro_data)
leveneTest_result <- leveneTest(Age~infect_state, data = tem_data)
print(leveneTest_result)

compare_means(Age~infect_state, data=tem_data,method = "wilcox.test")
compare_means(Age~infect_state, data=tem_data,method = "t.test")
#.y.      group1 group2         p    p.adj p.format p.signif method  
# Age   Infect No     0.156  0.16 0.16     ns       Wilcoxon
# Age   Infect No     0.196   0.2 0.2      ns       T-test

##Age_state
tem_data <-na.omit(RNA_meta_information[,c("infect_state","age_group")])
# 创建列联表
table <- table(tem_data$infect_state, tem_data$age_group)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##Fisher’s exact test
result_fisher <- fisher.test(table)
print(result_fisher)

##Age_state
tem_data <-na.omit(RNA_meta_information[,c("infect_state","gender")])
# 创建列联表
table <- table(tem_data$infect_state, tem_data$gender)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)

result_chi2 <- chisq.test(table)
print(result_chi2)


############孕产情况评估############
RNA_sample_info_Vill<-RNA_sample_info_final[,c("Vill_code","ID_full")]
RNA_sample_info_Decidua<-RNA_sample_info_final[,c("Decidua_code","ID_full")]
colnames(RNA_sample_info_Decidua)<-colnames(RNA_sample_info_Vill)<-c("sample_code","ID_full")
RNA_sample_info_final2<-rbind(RNA_sample_info_Vill,RNA_sample_info_Decidua)
head(RNA_sample_info_final2)
RNA_sample_info_final3<-merge(RNA_sample_info_final2,RNA_meta,by="sample_code")
head(RNA_sample_info_final3)
dim(RNA_sample_info_final3)#245  15

GP_data_final <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/normal_characteristics_for_848_cycle_final3.txt",header = T,sep = "\t")
head(GP_data_final)
GP_data_final$ID_full<-str_pad(GP_data_final$ID_full,width =8 ,side = c("left"),pad = "0")

GP_data_final$multigravida<-ifelse(GP_data_final$G_state !="G0",ifelse(is.na(GP_data_final$G_state),NA,"Yes"),"No")
GP_data_final$multiparous<-ifelse(GP_data_final$P_state !="P0",ifelse(is.na(GP_data_final$P_state),NA,"Yes"),"No")

GP_data_final2<-GP_data_final[which(!(is.na(GP_data_final$G_P_state))),]
table(GP_data_final2$multigravida)
table(GP_data_final2$multiparous)
head(GP_data_final2);dim(GP_data_final2)#842  20
dim(GP_data_final)#848  20

head(RNA_sample_info_final3)
RNA_sample_info_final3$ID_full<-str_pad(RNA_sample_info_final3$ID_full,width =8 ,side = c("left"),pad = "0")
RNA_sample_info_final4<-distinct(RNA_sample_info_final3[,c("ID_full","infect_state")])
dim(RNA_sample_info_final4)#132   2

GP_data_final4<-merge(RNA_sample_info_final4,GP_data_final2,by="ID_full")
dim(GP_data_final4)#130  21
head(GP_data_final4)

table(GP_data_final4$infect_state)
#Infect     No 
#   118     12 

setdiff(RNA_sample_info_final4$ID_full,GP_data_final2$ID_full) ## "22516163" "22506670"
setdiff(RNA_sample_info_final4$ID_full,GP_data_final$ID_full) 

table4 <- CreateTableOne(vars = c('multigravida', 'multiparous','Abortion_state'),data = GP_data_final4,smd=TRUE)#ttest、wilcox.test、chisq.test
table4 <- print(table4,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table4

table3 <- CreateTableOne(vars = c('multigravida', 'multiparous','Abortion_state'),data = GP_data_final4, strata = 'infect_state', smd=TRUE)#ttest、wilcox.test、chisq.test
table3 <- print(table3,smd=TRUE,showAllLevels = TRUE,noSpaces = TRUE,printToggle = FALSE)
table3

##multigravida
tem_data <-GP_data_final4[,c("infect_state","multigravida")]
# 创建列联表
table <- table(tem_data$infect_state, tem_data$multigravida)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##Fisher’s exact test
result_fisher <- fisher.test(table)
print(result_fisher)

##multiparous
tem_data <-GP_data_final4[,c("infect_state","multiparous")]
# 创建列联表
table <- table(tem_data$infect_state, tem_data$multiparous)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##chisq exact test
result_chi2 <- chisq.test(table)
print(result_chi2)

##Abortion_state
tem_data <-GP_data_final4[,c("infect_state","Abortion_state")]
# 创建列联表
table <- table(tem_data$infect_state, tem_data$Abortion_state)
print(table)
# 执行卡方检验
result_chi2 <- chisq.test(table)
# 评估每个单元格的预期频数是否至少为5
print(result_chi2$expected)
# 打印观察频数与预期频数之间的差异:小于2
print(result_chi2$residuals)
# 打印拟合优度检验结果：>0.05
print(result_chi2)
##Fisher’s exact test
result_fisher <- fisher.test(table)
print(result_fisher)
