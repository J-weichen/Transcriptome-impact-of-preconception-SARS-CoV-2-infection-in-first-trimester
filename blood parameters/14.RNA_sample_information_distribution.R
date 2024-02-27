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


####纳入转录测序候选样本信息
sample_data <- as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/测序样本选择-20230516.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
sample_data<-sample_data[,c(2,4,5,6,13)]
colnames(sample_data)<-c("Exe_code","ID","sample","Age","day_collect")
sample_data$ID_full<-str_pad(sample_data$ID,width =8 ,side = c("left"),pad = "0")
head(sample_data);dim(sample_data)# 320   4
length(unique(sample_data$ID_full))#317

sample_data2 <- as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/测序样本选择-20230516.xlsx", sheet =6, col_names = T, col_types = NULL, na = "", skip = 0))
sample_data2<-sample_data2[,c(1,2,4,5)]
colnames(sample_data2)<-c("indi_code","Exe_code","Vill_code","Decidua_code")
head(sample_data2);dim(sample_data2)# 136  4
length(unique(sample_data2$Exe_code))#134

length(intersect(unique(sample_data$Exe_code),unique(sample_data2$Exe_code))) #134
sample_data3<-merge(sample_data2,sample_data,by="Exe_code")
head(sample_data3);dim(sample_data3)#136 7
length(unique(sample_data3$Exe_code))#134

###读取最终用于计算的转录组信息
RNA_meta <- read.table(file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/Final_fielt_Decidua_Vill_analysis_metadata.txt",sep="\t",header = T)
head(RNA_meta);dim(RNA_meta)#245   8
range(na.omit(RNA_meta$LMP_infect))# -2 110
##提取样本信息获取最终用于计算的样本
R_code<-sample_data3[which(sample_data3$Vill_code %in% unique(as.character(RNA_meta$sample_code))),]$indi_code
D_code<-sample_data3[which(sample_data3$Decidua_code %in% unique(as.character(RNA_meta$sample_code))),]$indi_code
length(R_code);length(D_code)#120 125
code_final<-unique(c(R_code,D_code))
length(code_final)#132

sample_data4<-sample_data3[which(sample_data3$indi_code %in% code_final),]
dim(sample_data4)#132   9
head(sample_data4)
write.table(as.data.frame(sample_data4), file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/Final_fielt_Decidua_Vill_analysis_metadata_clinic.txt",sep = "\t",row.names=F, col.names=T) 

##family remove
setdiff(sample_data3$indi_code,code_final)# "E63"   "E79-2" "E106"  "E43-2"
setdiff(code_final,sample_data3$indi_code)#

length(unique(as.character(na.omit(c(sample_data4$Vill_code,sample_data4$Decidua_code)))))#264
library_rm<-setdiff(unique(as.character(na.omit(c(sample_data4$Vill_code,sample_data4$Decidua_code)))), unique(as.character(RNA_meta$sample_code)))#266
sample_data4[which(sample_data4$Vill_code %in% library_rm),]$Vill_code<- NA
sample_data4[which(sample_data4$Decidua_code %in% library_rm),]$Decidua_code<- NA

##########################################
##读取临床表型数据
analysis_used <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/clinical_index_data_for_analysis_used_final.txt",header = T,sep = "\t")
head(analysis_used)
analysis_collect<-analysis_used[which(analysis_used$ID_full %in% unique(sample_data4$ID_full)),]
dim(analysis_collect)# 5911   15
analysis_collect<-analysis_collect[which(analysis_collect$Year_month %in% c("20232","20233","20234")),]
dim(analysis_collect)# 5866   15
table(analysis_collect$Year_month)
#20232 20233 20234 
#1147  2149  2570 
length(unique(analysis_collect$ID_full))## 134
analysis_index_wide <- reshape2::dcast(analysis_collect,ID_full ~ E_name, value.var = "value")
head(analysis_index_wide);dim(analysis_index_wide)

##转录检测数据与临床检测数据合并
analysis_index_wide2 <-merge(sample_data4[,c("ID_full","Vill_code","Decidua_code")],analysis_index_wide,by="ID_full")
head(analysis_index_wide2);dim(analysis_index_wide2)## 132  48
which(is.na(analysis_index_wide2$Decidua_code))#  51  65  70  79  99 105 127
which(is.na(analysis_index_wide2$Vill_code))   #   7   9  11  19  20  34  72  93 100 124 125 130
names(analysis_index_wide2)[4] <- "Gamma-GT"

write.table(as.data.frame(analysis_index_wide2), file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/analysis_index_wide.txt",sep = "\t",row.names=F, col.names=T) 

##提取检测样本临床数据
##读取被记录的康复信息完整的样本内容
collect_sample_infor <- read.table(file="D:/PROJECT/新冠/2022_data/病历补充核准/collect_sample_infor.txt",header = T,sep = "\t")
collect_sample_infor$ID_full<-str_pad(collect_sample_infor$ID,width =8 ,side = c("left"),pad = "0")
head(collect_sample_infor);dim(collect_sample_infor)#249   6
length(unique(collect_sample_infor$ID_full))#249
table(collect_sample_infor$Infect_state)
# no yes 
# 17 232

collect_RNA_infor_final <-distinct(merge(sample_data4[,c("ID_full","Age")],collect_sample_infor,by="ID_full"))
collect_RNA_infor_final$age_group<-ifelse(collect_RNA_infor_final$Age>=35,"AMA","YMA")
head(collect_RNA_infor_final);dim(collect_RNA_infor_final)#132   8
table(collect_RNA_infor_final$Infect_state)
#no yes 
#12 120 
rownames(collect_RNA_infor_final)<-collect_RNA_infor_final$ID_full
collect_RNA_infor_final<-collect_RNA_infor_final[order(collect_RNA_infor_final$Infect_state,collect_RNA_infor_final$Age,decreasing = F),]
head(collect_RNA_infor_final)

####绘图####################
heat_data<-analysis_index_wide2#[-which(analysis_index_wide2$Vill_code %in% c("R43.2","R79.2")),]
head(heat_data);dim(heat_data)# 132  48
rownames(heat_data)<-heat_data$ID_full
heat_data<-heat_data[,-1]
sample_names<-rownames(heat_data)

heat_data[is.na(heat_data)] <- 0 # =="na"
heat_data[heat_data !=0] <- 1
heat_data<-as.data.frame(lapply(heat_data,as.numeric))
rownames(heat_data)<-sample_names
heat_data[1:4,1:5];dim(heat_data)
#03246814 15134013
index_order<-c("WBC","RBC","HGB","HCT","MCV","MCH","MCHC","PLT","RDW_CV","PCT","LYMPH_pct","LYMPH_av","Neut_pct","Neut_av","EO_pct","EO_av","BASO_pct","BASO_av","MONO_pct","MONO_av",
               "MPV","PDW","P.LCR",
               "ALT","AST","TBIL","Gamma.GT","CK","CK.MB","BUN_Urea","Cr","UA","T.CHO","TG","HDL.C","LDL.C","GLU",
               "PT","A","INR","Fib","APTT","APTT_R","TT","R")

heat_data<-heat_data[collect_RNA_infor_final$ID_full,c("Vill_code","Decidua_code",index_order)]

heat_data2<-t(heat_data)
dim(heat_data2) ## 47 132
range(heat_data2)# 0 1
heat_data2[1:4,1:5]

head(collect_RNA_infor_final)
annotation_col<-collect_RNA_infor_final[,c("age_group","Infect_state")]
annotation_row <-data.frame(Class=factor(c(rep("RNA",2),rep("Blood_routine_index",23),rep("Serum_biochemical_parameter",14),rep("Coagulation_profile",8))))
rownames(annotation_row) = c("Vill_code","Decidua_code",index_order)

anno_colors = list(
  Class =c(RNA=ppCor[1],Blood_routine_index=ppCor[2],Serum_biochemical_parameter=ppCor[3],Coagulation_profile=ppCor[4]),
  Infect_state=c(yes=ppCor[1],no=ppCor[2]),
  age_group=c(AMA=ppCor[6],YMA=ppCor[3]))

#labels_col = c("")
p1<-pheatmap(heat_data2,cluster_rows=F, cluster_cols =F,show_colnames = F,
             annotation_col = annotation_col,
             annotation_row=annotation_row,
             annotation_colors = anno_colors, 
             gaps_row = c(2,25,39),
             main = "Sample distribution",
             legend_breaks = c(0.2,0.8),legend_labels = c("none","detected"),
             color = colorRampPalette(colors = c("purple","gold"))(2))

pdf("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/no_rm_long/RNA_index_detected_in_bulk_RNA_used_sample_heatmap.pdf",width=8, height=8)
print(p1)
dev.off()

##绘制每个index的样本数目
index_sample_num<-as.data.frame(rowSums(heat_data2))
head(index_sample_num)
index_sample_num$index<-rownames(index_sample_num)
colnames(index_sample_num) <-c("sample_number","Index")
index_sample_num$Index<-factor(index_sample_num$Index,levels = c("Vill_code","Decidua_code",index_order))
index_sample_num<-index_sample_num[order(index_sample_num$Index,decreasing = F),]
index_sample_num$Class<-c(rep("RNA",2),rep("Blood_routine_index",23),rep("Serum_biochemical_parameter",14),rep("Coagulation_profile",8))
index_sample_num$Class<-factor(index_sample_num$Class,levels = c("RNA","Blood_routine_index","Serum_biochemical_parameter","Coagulation_profile"))

numplot0<-ggplot(data=index_sample_num, mapping=aes(x= Index,y=sample_number,fill=Class))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+
  geom_text(aes(label=sample_number),size=3,position=position_stack(vjust=1.01))+
  scale_fill_manual(name="Index_class",values=ppCor)+
  theme_classic()+labs(x="Index class",y="Detected sample number",title="The sample number for each index")+
  scale_y_continuous(breaks = seq(0,160,by=20))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
numplot0
##ONLY RNA
numplot01<-ggplot(data=index_sample_num[which(index_sample_num$Class =="RNA"),], mapping=aes(x= Index,y=sample_number))+
  geom_bar(stat="identity",width=0.8,position= 'stack',fill=c("blue","darkred"))+
  geom_text(aes(label=sample_number),size=5,position=position_stack(vjust=1.01))+
  theme_classic()+labs(x="Index class",y="Detected sample number",title="The sample number for two Tissue")+
  scale_y_continuous(breaks = seq(0,160,by=20))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
numplot01

ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/no_rm_long/number_in_used_RNA_library_index_all.pdf",numplot0,width=12, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/no_rm_long/number_in_used_RNA_library_index_RNA.pdf",numplot01,width=6, height=6)

###sample number
collect_RNA_infor_final$Age<-as.numeric(collect_RNA_infor_final$Age)
collect_RNA_infor_final$LMP_infect2<-collect_RNA_infor_final$LMP_infect
collect_RNA_infor_final[is.na(collect_RNA_infor_final$LMP_infect),]$LMP_infect2<- 160
collect_RNA_infor_final$GS_week <-floor(collect_RNA_infor_final$LMP_operate/7)
write.table(as.data.frame(collect_RNA_infor_final), file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/collect_RNA_infor_final.txt",sep = "\t",row.names=F, col.names=T) 

range(collect_RNA_infor_final$Age)##18 46
range(collect_RNA_infor_final$LMP_operate)##38 83
range(na.omit(collect_RNA_infor_final$LMP_infect))##-31 156
range(collect_RNA_infor_final$LMP_infect2)##-31 160
range(collect_RNA_infor_final$GS_week)#5 11

numplot1<-ggplot(collect_RNA_infor_final,aes(x=Age))+geom_histogram(binwidth=1,color='black',fill='#03B0AB',cex=1)+
  stat_bin(binwidth=1, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(18,46,1))+
  labs(x = "Age(years)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
      axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
      axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))

numplot2<-ggplot(collect_RNA_infor_final,aes(x=LMP_infect2))+geom_histogram(binwidth=5,color='black',fill='pink',cex=1)+
  stat_bin(binwidth=5, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(-30,160,5))+
  labs(x = "LMP_infect(days)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
numplot3<-ggplot(collect_RNA_infor_final,aes(x=LMP_infect2))+geom_histogram(binwidth=1,color='black',fill='pink',cex=1)+
  stat_bin(binwidth=1, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(-30,160,1))+
  labs(x = "LMP_infect(days)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))

numplot4<-ggplot(collect_RNA_infor_final,aes(x=LMP_operate))+geom_histogram(binwidth=1,color='black',fill='lightblue',cex=1)+
  stat_bin(binwidth=1, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(38,83,1))+
  labs(x = "LMP_operate(days)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
numplot5<-ggplot(collect_RNA_infor_final,aes(x=GS_week))+geom_histogram(binwidth=1,color='black',fill='lightblue',cex=1)+
  stat_bin(binwidth=1, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(5,11,1))+
  labs(x = "GS_week(weeks)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
numplot2
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/no_rm_long/number_in_used_RNA_library_Age.pdf",numplot1,width=8, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/no_rm_long/number_in_used_RNA_library_infect_day.pdf",numplot2,width=8, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/no_rm_long/number_in_used_RNA_library_infect_day2.pdf",numplot3,width=16, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/no_rm_long/number_in_used_RNA_library_pregnancy_day.pdf",numplot4,width=8, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/no_rm_long/number_in_used_RNA_library_GS_week.pdf",numplot5,width=6, height=5)



####################rm two sample with recover day in -31 and 156#####本次未跑20230907
#######################################################################
collect_RNA_infor_final_rm<-collect_RNA_infor_final[which(!(collect_RNA_infor_final$LMP_infect %in% c(-31,156))),]
dim(collect_RNA_infor_final_rm)#130  10
collect_RNA_infor_final[which(collect_RNA_infor_final$LMP_infect %in% c(-31,156)),]$ID_full#"17394215" "22682021"

heat_data<-analysis_index_wide2[-which(analysis_index_wide2$Vill_code %in% c("R43.2","R79.2")),]
heat_data_rm<-heat_data[which(!(heat_data$ID_full %in% c("17394215","22682021"))),]

head(heat_data_rm);dim(heat_data_rm)#130  48
rownames(heat_data_rm)<-heat_data_rm$ID_full
heat_data_rm<-heat_data_rm[,-1]
sample_names<-rownames(heat_data_rm)

###############
heat_data_rm[is.na(heat_data_rm)] <- 0 # =="na"
heat_data_rm[heat_data_rm !=0] <- 1
heat_data_rm<-as.data.frame(lapply(heat_data_rm,as.numeric))
rownames(heat_data_rm)<-sample_names
heat_data_rm[1:4,1:5];dim(heat_data_rm)
#03246814 15134013
index_order<-c("WBC","RBC","HGB","HCT","MCV","MCH","MCHC","PLT","RDW_CV","PCT","LYMPH_pct","LYMPH_av","Neut_pct","Neut_av","EO_pct","EO_av","BASO_pct","BASO_av","MONO_pct","MONO_av",
               "MPV","PDW","P.LCR",
               "ALT","AST","TBIL","γ.GT","CK","CK.MB","BUN_Urea","Cr","UA","T.CHO","TG","HDL.C","LDL.C","GLU",
               "PT","A","INR","Fib","APTT","APTT_R","TT","R")

heat_data_rm<-heat_data_rm[collect_RNA_infor_final_rm$ID_full,c("Vill_code","Decidua_code",index_order)]

heat_data_rm2<-t(heat_data_rm)
dim(heat_data_rm2) ## 47 130
range(heat_data_rm2)# 0 1
heat_data_rm2[1:4,1:5]

head(collect_RNA_infor_final_rm)
annotation_col<-collect_RNA_infor_final_rm[,c("age_group","Infect_state")]
annotation_row <-data.frame(Class=factor(c(rep("RNA",2),rep("Blood_routine_index",23),rep("Serum_biochemical_parameter",14),rep("Coagulation_profile",8))))
rownames(annotation_row) = c("Vill_code","Decidua_code",index_order)

anno_colors = list(
  Class =c(RNA=ppCor[1],Blood_routine_index=ppCor[2],Serum_biochemical_parameter=ppCor[3],Coagulation_profile=ppCor[4]),
  Infect_state=c(yes=ppCor[1],no=ppCor[2]),
  age_group=c(AMA=ppCor[6],YMA=ppCor[3]))

#labels_col = c("")
p1<-pheatmap(heat_data_rm2,cluster_rows=F, cluster_cols =F,show_colnames = F,
             annotation_col = annotation_col,
             annotation_row=annotation_row,
             annotation_colors = anno_colors, 
             gaps_row = c(2,25,39),
             main = "Sample distribution",
             legend_breaks = c(0.2,0.8),legend_labels = c("none","detected"),
             color = colorRampPalette(colors = c("purple","gold"))(2))

pdf("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/rm_long/RNA_index_detected_in_bulk_RNA_used_sample_heatmap.pdf",width=8, height=8)
print(p1)
dev.off()


##绘制每个index的样本数目
index_sample_num<-as.data.frame(rowSums(heat_data_rm2))
head(index_sample_num)
index_sample_num$index<-rownames(index_sample_num)
colnames(index_sample_num) <-c("sample_number","Index")
index_sample_num$Index<-factor(index_sample_num$Index,levels = c("Vill_code","Decidua_code",index_order))
index_sample_num<-index_sample_num[order(index_sample_num$Index,decreasing = F),]
index_sample_num$Class<-c(rep("RNA",2),rep("Blood_routine_index",23),rep("Serum_biochemical_parameter",14),rep("Coagulation_profile",8))
index_sample_num$Class<-factor(index_sample_num$Class,levels = c("RNA","Blood_routine_index","Serum_biochemical_parameter","Coagulation_profile"))

numplot0<-ggplot(data=index_sample_num, mapping=aes(x= Index,y=sample_number,fill=Class))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+
  geom_text(aes(label=sample_number),size=3,position=position_stack(vjust=1.01))+
  scale_fill_manual(name="Index_class",values=ppCor)+
  theme_classic()+labs(x="Index class",y="Detected sample number",title="The sample number for each index")+
  scale_y_continuous(breaks = seq(0,160,by=20))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
numplot0
##ONLY RNA
numplot01<-ggplot(data=index_sample_num[which(index_sample_num$Class =="RNA"),], mapping=aes(x= Index,y=sample_number))+
  geom_bar(stat="identity",width=0.8,position= 'stack',fill=c("blue","darkred"))+
  geom_text(aes(label=sample_number),size=5,position=position_stack(vjust=1.01))+
  theme_classic()+labs(x="Index class",y="Detected sample number",title="The sample number for two Tissue")+
  scale_y_continuous(breaks = seq(0,160,by=20))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
numplot01

ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/rm_long/number_in_used_RNA_library_index_all.pdf",numplot0,width=12, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/rm_long/number_in_used_RNA_library_index_RNA.pdf",numplot01,width=6, height=6)

###sample number
collect_RNA_infor_final_rm$Age<-as.numeric(collect_RNA_infor_final_rm$Age)
collect_RNA_infor_final_rm$LMP_infect2<-collect_RNA_infor_final_rm$LMP_infect
collect_RNA_infor_final_rm[is.na(collect_RNA_infor_final_rm$LMP_infect),]$LMP_infect2<- 120
collect_RNA_infor_final_rm$GS_week <-floor(collect_RNA_infor_final_rm$LMP_operate/7)

range(collect_RNA_infor_final_rm$Age)##18 46
range(collect_RNA_infor_final_rm$LMP_operate)##38 83
range(na.omit(collect_RNA_infor_final_rm$LMP_infect))## -2 110
range(collect_RNA_infor_final_rm$LMP_infect2)## -2 120
range(collect_RNA_infor_final_rm$GS_week)#5 11

numplot1<-ggplot(collect_RNA_infor_final_rm,aes(x=Age))+geom_histogram(binwidth=1,color='black',fill='#03B0AB',cex=1)+
  stat_bin(binwidth=1, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(18,46,1))+
  labs(x = "Age(years)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))

numplot2<-ggplot(collect_RNA_infor_final_rm,aes(x=LMP_infect2))+geom_histogram(binwidth=5,color='black',fill='pink',cex=1)+
  stat_bin(binwidth=5, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(-5,120,5))+
  labs(x = "LMP_infect(days)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
numplot3<-ggplot(collect_RNA_infor_final_rm,aes(x=LMP_infect2))+geom_histogram(binwidth=1,color='black',fill='pink',cex=1)+
  stat_bin(binwidth=1, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(-2,120,1))+
  labs(x = "LMP_infect(days)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))

numplot4<-ggplot(collect_RNA_infor_final_rm,aes(x=LMP_operate))+geom_histogram(binwidth=1,color='black',fill='lightblue',cex=1)+
  stat_bin(binwidth=1, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(38,83,1))+
  labs(x = "LMP_operate(days)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
numplot5<-ggplot(collect_RNA_infor_final_rm,aes(x=GS_week))+geom_histogram(binwidth=1,color='black',fill='lightblue',cex=1)+
  stat_bin(binwidth=1, geom='text', color='black', size=4,aes(label=..count..), position=position_stack(vjust=1.05))+
  theme_classic(base_size = 20)+scale_x_continuous(breaks = seq(5,11,1))+
  labs(x = "GS_week(weeks)", y = "sample_num(n)", title ="All used  sample number in RNAs")+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
numplot2
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/rm_long/number_in_used_RNA_library_Age.pdf",numplot1,width=8, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/rm_long/number_in_used_RNA_library_infect_day.pdf",numplot2,width=8, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/rm_long/number_in_used_RNA_library_infect_day2.pdf",numplot3,width=16, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/rm_long/number_in_used_RNA_library_pregnancy_day.pdf",numplot4,width=8, height=5)
ggsave("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/fielt_data/rm_long/number_in_used_RNA_library_GS_week.pdf",numplot5,width=6, height=5)
