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
#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

###########计算信息中心人流信息数据指标分析
clinical_data_final <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/信息中心人流信息_所有.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_data_final)
clinical_data_final$Year_month <-clinical_data_final$年月
clinical_data_final$Year_month2<-ifelse(clinical_data_final$日期...8>15,paste0(clinical_data_final$Year_month,"_post"),paste0(clinical_data_final$Year_month,"_before"))

head(clinical_data_final)
##计算每月人数
table(clinical_data_final$Year_month)
##20221 202210 202211 202212  20222  20223  20224  20225  20226  20227  20228  20229  20231  20232  20233  20234 
##   91    107    158    168     86    150    119    128    157    124    130    140     92    124    159    148 
table(clinical_data_final$Year_month2)
##20221_before    20221_post 202210_before   202210_post 202211_before   202211_post 202212_before   202212_post  20222_before 
##     40            51            43            64            82            76            90            78            35 
##20222_post  20223_before    20223_post  20224_before    20224_post  20225_before    20225_post  20226_before    20226_post 
##     51            73            77            58            61            54            74            88            69 
##20227_before    20227_post  20228_before    20228_post  20229_before    20229_post  20231_before    20231_post  20232_before 
##     67            57            57            73            64            76            56            36            70 
##20232_post  20233_before    20233_post  20234_before    20234_post 
##     54            75            84            73            75 
month_order<-c("20221","20222","20223","20224","20225","20226","20227","20228","20229","202210","202211","202212","20231","20232","20233","20234")
month_order2<-c("20221_before","20221_post","20222_before","20222_post","20223_before","20223_post","20224_before","20224_post","20225_before","20225_post","20226_before","20226_post","20227_before","20227_post","20228_before","20228_post","20229_before","20229_post","202210_before","202210_post","202211_before","202211_post","202212_before","202212_post","20231_before","20231_post","20232_before","20232_post","20233_before","20233_post","20234_before","20234_post")
clinical_data_final$age_group<-ifelse(clinical_data_final$年龄>=35,"AMA","YMA")
patient_ID0<-distinct(clinical_data_final[,c("ID","age_group","Year_month","Year_month2")])
patient_ID0$Year_month <- factor(patient_ID0$Year_month,levels=month_order)
patient_ID0$Year_month2 <- factor(patient_ID0$Year_month2,levels=month_order2)
head(patient_ID0)
number_plot00<-ggplot(data=patient_ID0, mapping=aes(x=Year_month,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='stack')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_stack(0.5))+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
y_max <- max(aggregate(ID~Year_month+age_group,data=patient_ID0,length)$ID)
number_plot002<-ggplot(data=patient_ID0, mapping=aes(x=Year_month,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='dodge')+
  scale_fill_manual(values=c('#999999','#E69F00'))+  ylim(0,y_max+5)+
  geom_text(stat='count',aes(label=..count..), color="black", size=3.5,position=position_dodge(0.5),vjust=-0.5)+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
ggsave("D:/PROJECT/新冠/2022_data/Final_number_all_sample_1.pdf",number_plot00,width=8, height=6)
ggsave("D:/PROJECT/新冠/2022_data/Final_number_all_sample_2.pdf",number_plot002,width=12, height=6)

number_plot001<-ggplot(data=patient_ID0, mapping=aes(x=Year_month2,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='stack')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_stack(0.5))+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
y_max <- max(aggregate(ID~Year_month2+age_group,data=patient_ID0,length)$ID)
number_plot003<-ggplot(data=patient_ID0, mapping=aes(x=Year_month2,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='dodge')+
  scale_fill_manual(values=c('#999999','#E69F00'))+  ylim(0,y_max+5)+
  geom_text(stat='count',aes(label=..count..), color="black", size=3.5,position=position_dodge(0.5),vjust=-0.5)+
  theme_minimal()+theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0,angle = 90),legend.title = element_text(size = 9))
ggsave("D:/PROJECT/新冠/2022_data/Final_number_all_sample_3.pdf",number_plot001,width=12, height=6)
ggsave("D:/PROJECT/新冠/2022_data/Final_number_all_sample_4.pdf",number_plot003,width=16, height=6)

####计算收样的统计
sample_data <- as.data.frame(read_excel("D:/PROJECT/新冠/2022_data/病历补充核准/7.病例记录20230503-终版.xlsx", sheet =2, col_names = T, col_types = NULL, na = "", skip = 0))
sample_data<-sample_data[,-(1:2)]
head(sample_data)
dim(sample_data)#320  17
sample_data1<-sample_data[which(sample_data$异常与否 =="OK"),]
dim(sample_data1)#301  17
head(sample_data1)

###
###sample number
sample_number<-sample_data1[,c(1:7)]
colnames(sample_number)<-c("sample","Age","Vill_80","Decidual_80","Vill_later","Decidual_later","Plasma")
head(sample_number)

sample_number2<- melt(sample_number,id.vars=c("sample","Age"),variable.name="class",value.name = "number")
head(sample_number2)
sample_number2<-sample_number2[which(sample_number2$number != "na"),]
sample_number2$age_group<-ifelse(sample_number2$Age>=35,"AMA","YMA")
head(sample_number2);dim(sample_number2)##1157    4

number_plot001<-ggplot(data=sample_number2, mapping=aes(x=class,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='stack')+
  scale_fill_manual(values=c('#999999','#E69F00',"blue"))+
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_stack(0.5))+
  theme_classic() + labs(x = "Age(year)", y = "sample_num(n)", title ="All normal sample number in each class")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 6,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
number_plot001
ggsave(file="D:/PROJECT/重点研发/2022年年度总结/All_normal_sample_number_in_each_class.pdf",number_plot001,width = 6, height =5)

##re图展示样本分布
heat_data<-sample_number
rownames(heat_data)<-heat_data$sample
heat_data<-heat_data[,-c(1:2)]
heat_data[heat_data =="na"] <- 0
heat_data[heat_data !=0] <- 1

heat_data<-as.data.frame(lapply(heat_data,as.numeric))
rownames(heat_data)<-sample_number$sample
heat_data[1:4,1:5]

heat_data2<-t(heat_data)
dim(heat_data2) ##  5 291
range(heat_data2)# 0 1
heat_data2[1:4,1:5]
col_data<-sample_number[,1:2]
col_data$age_group<-ifelse(col_data$Age>=35,"AMA","YMA")
rownames(col_data) = col_data$sample

annotation_col<-data.frame(age_group = factor(col_data$age_group))
rownames(annotation_col) = colnames(heat_data2)
annotation_row <-data.frame(Class=factor(c("Vill_80","Decidual_80","Vill_later","Decidual_later","Plasma")))
rownames(annotation_row) = annotation_row$Class

anno_colors = list(
  Class =c(Vill_80=ppCor[1],Decidual_80=ppCor[2],Vill_later=ppCor[3],Decidual_later=ppCor[4],Plasma=ppCor[5]),
  age_group=c(AMA=ppCor[6],YMA=ppCor[3]))
labels_col = c("")
library(pheatmap)
heat_plot<-pheatmap(heat_data2,cluster_rows=F, cluster_cols =F, na_col = "black",
             annotation_col = annotation_col, annotation_row=annotation_row,
             annotation_colors = anno_colors, 
             labels_col = labels_col,
             #gaps_row = c(18),
             main = "Sample distribution",
             legend_breaks = c(0.2,0.8),legend_labels = c("none","Yes"),
             color = colorRampPalette(colors = c("purple","gold"))(2))

pdf(file="D:/PROJECT/重点研发/2022年年度总结/heatmap_distribution_All_normal_sample_number_in_each_class.pdf",width = 6, height =4)
print(heat_plot)
dev.off()
#######################age distribution#############
##all normal sample 
All_sample<-as.data.frame(table(sample_data1$年龄))
All_sample<-All_sample[which(All_sample$Var1 !="NA"),]
age_line<-ggplot(All_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "Orange") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Age(year)", y = "Tissue_num(n)", title ="All normal sample number in each age")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(All_sample$Freq)+5)
age_line
ggsave(file="D:/PROJECT/重点研发/2022年年度总结/Frequency_distribution_all_normal_sample_number.pdf",age_line,width = 6, height =3)

##only Tissue
sample_data1$type <-"Tissue"
sample_data1[which(sample_data1$`绒毛-80存储`=="na" & sample_data1$`蜕膜-80存储`=="na" & sample_data1$`绒毛-RNAlater`=="na" & sample_data1$`蜕膜-RNAlater`=="na" ),]$type <-"Plasma_only"
table(sample_data1$type)
#Plasma_only      Tissue 
#       26         265 
sample_data2<-sample_data1[which(sample_data1$type =="Tissue"),]
head(sample_data2)
Tissue_sample<-as.data.frame(table(sample_data2$年龄))
Tissue_sample<-Tissue_sample[which(Tissue_sample$Var1 !="NA"),]
Tissue_age_line<-ggplot(Tissue_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "blue") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Age(year)", y = "Tissue_num(n)", title ="Tissue sample number in each age")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(Tissue_sample$Freq)+5)
ggsave(file="D:/PROJECT/重点研发/2022年年度总结/Frequency_distribution_tissue_normal_sample_number.pdf",Tissue_age_line,width = 6, height =3)

head(sample_data2)

#######################pregnancy days distribution#############
##all normal sample 
All_sample<-as.data.frame(table(sample_data1$末次月经相较于手术日_T))
All_sample$day<-as.numeric(as.character(All_sample$Var1))
All_sample<-All_sample[which(All_sample$day >0),]
ggplot(All_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "Red") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="All normal sample number in each gestation days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(All_sample$Freq)+5)

##only Tissue
sample_data1$type <-"Tissue"
sample_data1[which(sample_data1$`绒毛-80存储`=="na" & sample_data1$`蜕膜-80存储`=="na" & sample_data1$`绒毛-RNAlater`=="na" & sample_data1$`蜕膜-RNAlater`=="na" ),]$type <-"Plasma_only"
table(sample_data1$type)
#Plasma_only      Tissue 
#       32         269 
sample_data2<-sample_data1[which(sample_data1$type =="Tissue"),]
head(sample_data2)
Tissue_sample<-as.data.frame(table(sample_data2$末次月经相较于手术日_T))
Tissue_sample$day<-as.numeric(as.character(Tissue_sample$Var1))
Tissue_sample<-Tissue_sample[which(Tissue_sample$day >0),]
ggplot(Tissue_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "green") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "Tissue_num(n)", title ="Tissue sample number in each gestation days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(Tissue_sample$Freq)+5)

#######################infection days distribution#############
sample_data1_raw<-sample_data1[which(sample_data1$手术日相对于感染的时间_T >0 & sample_data1$手术日相对于感染的时间_T <200),]
sample_data1_raw2<-sample_data1_raw[which(sample_data1_raw$手术日相对于感染的时间_T != "NA"),];dim(sample_data1_raw2)

ggplot(data=sample_data1_raw2,aes(x=手术日相对于感染的时间_T)) + 
  geom_histogram(binwidth = 15,colour="black", fill="Yellow")+
  stat_bin(binwidth=15,geom="text",aes(label=..count..),vjust = -1) +ylim(0,55)+
  scale_x_continuous(breaks = seq(0,max(sample_data1_raw2$手术日相对于感染的时间_T),15),labels = comma)+
theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="All normal sample number in each infection days")
  

##all normal sample 
All_sample<-as.data.frame(table(sample_data1$手术日相对于感染的时间_T))
All_sample$day<-as.numeric(as.character(All_sample$Var1))
All_sample<-All_sample[which(All_sample$day >0 & All_sample$day <200),]
ggplot(All_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "Yellow") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="All normal sample number in each infection days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(All_sample$Freq)+5)+
  theme(axis.text.x = element_text(size =8,colour = 'black',vjust=0.5,hjust=1,angle = 90))

##only Tissue
sample_data1$type <-"Tissue"
sample_data1[which(sample_data1$`绒毛-80存储`=="na" & sample_data1$`蜕膜-80存储`=="na" & sample_data1$`绒毛-RNAlater`=="na" & sample_data1$`蜕膜-RNAlater`=="na" ),]$type <-"Plasma_only"
table(sample_data1$type)
#Plasma_only      Tissue 
#       32         269 
sample_data2<-sample_data1[which(sample_data1$type =="Tissue"),]
head(sample_data2)

sample_data2_raw<-sample_data2[which(sample_data2$手术日相对于感染的时间_T >0 & sample_data2$手术日相对于感染的时间_T <200),]
sample_data2_raw2<-sample_data2_raw[which(sample_data2_raw$手术日相对于感染的时间_T != "NA"),];dim(sample_data2_raw2)

ggplot(data=sample_data2_raw2,aes(x=手术日相对于感染的时间_T)) + 
  geom_histogram(binwidth = 15,colour="black", fill="purple")+
  stat_bin(binwidth=15,geom="text",aes(label=..count..),vjust = -1) +ylim(0,55)+
  scale_x_continuous(breaks = seq(0,max(sample_data1_raw2$手术日相对于感染的时间_T),15),labels = comma)+
  theme_classic() + labs(x = "Time(days)", y = "Tissue_num(n)", title ="Tissue sample number in each infection days")

Tissue_sample<-as.data.frame(table(sample_data2$手术日相对于感染的时间_T))
Tissue_sample$day<-as.numeric(as.character(Tissue_sample$Var1))
Tissue_sample<-Tissue_sample[which(Tissue_sample$day >0 & Tissue_sample$day <200),]
ggplot(Tissue_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "purple") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "Tissue_num(n)", title ="Tissue sample number in each infection days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(Tissue_sample$Freq)+5)+
  theme(axis.text.x = element_text(size =8,colour = 'black',vjust=0.5,hjust=1,angle = 90))

#######################infect_lmp days distribution#############
sample_data1_raw2<-sample_data1[which(sample_data1$末次月经相较于感染日_T != "NA"),];dim(sample_data1_raw2)

ggplot(data=sample_data1_raw2,aes(x=末次月经相较于感染日_T)) + 
  geom_histogram(binwidth = 15,colour="black", fill="green")+
  stat_bin(binwidth=15,geom="text",aes(label=..count..),vjust = -1) +ylim(0,65)+xlim(-30,150)+
  scale_x_continuous(breaks = seq(-30,max(sample_data1_raw2$末次月经相较于感染日_T),15),labels = comma)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="All normal sample number in LMP to infection days")


##all normal sample 
All_sample<-as.data.frame(table(sample_data1$末次月经相较于感染日_T))
All_sample$day<-as.numeric(as.character(All_sample$Var1))

ggplot(All_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "green") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="All normal sample number in  LMP to infection days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(All_sample$Freq)+5)+
  theme(axis.text.x = element_text(size =8,colour = 'black',vjust=0.5,hjust=1,angle = 90))

##only Tissue
sample_data1$type <-"Tissue"
sample_data1[which(sample_data1$`绒毛-80存储`=="na" & sample_data1$`蜕膜-80存储`=="na" & sample_data1$`绒毛-RNAlater`=="na" & sample_data1$`蜕膜-RNAlater`=="na" ),]$type <-"Plasma_only"
table(sample_data1$type)
#Plasma_only      Tissue 
#       32         269 
sample_data2<-sample_data1[which(sample_data1$type =="Tissue"),]
head(sample_data2)
sample_data2_raw2<-sample_data2[which(sample_data2$末次月经相较于感染日_T != "NA"),];dim(sample_data2)

ggplot(data=sample_data2_raw2,aes(x=末次月经相较于感染日_T)) + 
  geom_histogram(binwidth = 15,colour="black", fill="red")+
  stat_bin(binwidth=15,geom="text",aes(label=..count..),vjust = -1) +ylim(0,60)+
  scale_x_continuous(breaks = seq(-30,max(sample_data1_raw2$末次月经相较于感染日_T),15),labels = comma)+
  theme_classic() + labs(x = "Time(days)", y = "Tissue_num(n)", title ="Tissue sample number in LMP to infection days")

Tissue_sample<-as.data.frame(table(sample_data2$末次月经相较于感染日_T))
Tissue_sample$day<-as.numeric(as.character(Tissue_sample$Var1))
ggplot(Tissue_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "red") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "Tissue_num(n)", title ="Tissue sample number in LMP to infection days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(Tissue_sample$Freq)+5)+
  theme(axis.text.x = element_text(size =8,colour = 'black',vjust=0.5,hjust=1,angle = 90))

#####################################################################
#####################################################################
##选取同时RNAlater 存有绒毛和蜕膜
head(sample_number)
sample_number1<-sample_number[which(sample_number$Vill_later == "Y" & sample_number$Decidual_later =="Y"),]
dim(sample_number1)##220   7
sample_number2<- melt(sample_number1,id.vars=c("sample","Age"),variable.name="class",value.name = "number")
head(sample_number2)
sample_number2<-sample_number2[which(sample_number2$number != "na"),]
sample_number2$age_group<-ifelse(sample_number2$Age>=35,"AMA","YMA")
head(sample_number2);dim(sample_number2)##1157    4

number_plot001<-ggplot(data=sample_number2, mapping=aes(x=class,fill=age_group))+
  geom_bar(stat="count",width=0.5,position='stack')+
  scale_fill_manual(values=c('#999999','#E69F00',"blue"))+
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_stack(0.5))+
  theme_classic() + labs(x = "Age(year)", y = "sample_num(n)", title ="Vill and Decidual both in RNA later sample number in each class")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 6,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
number_plot001

#######################age distribution#############
head(sample_data1)
sample_data1_both<-sample_data1[which(sample_data1$'绒毛-RNAlater' == "Y" & sample_data1$'蜕膜-RNAlater'  =="Y"),]
head(sample_data1_both)
##all normal sample 
All_sample<-as.data.frame(table(sample_data1_both$年龄))
All_sample<-All_sample[which(All_sample$Var1 !="NA"),]
ggplot(All_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "Orange") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Age(year)", y = "Tissue_num(n)", title ="Vill and Decidual both in RNA later sample number in each age")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(All_sample$Freq)+5)

#######################pregnancy days distribution#############
##all normal sample 
All_sample<-as.data.frame(table(sample_data1_both$末次月经相较于手术日_T))
All_sample$day<-as.numeric(as.character(All_sample$Var1))
All_sample<-All_sample[which(All_sample$day >0),]
ggplot(All_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "Red") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="Vill and Decidual both in RNA later sample number in each gestation days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(All_sample$Freq)+5)

#######################infection days distribution#############
range(na.omit(sample_data1_both$手术日相对于感染的时间_T))##7 45044
sample_data1_raw<-sample_data1_both[which(sample_data1_both$手术日相对于感染的时间_T >0 & sample_data1_both$手术日相对于感染的时间_T <200),]
sample_data1_raw2<-sample_data1_raw[which(sample_data1_raw$手术日相对于感染的时间_T != "NA"),];dim(sample_data1_raw2)
range(sample_data1_raw2$手术日相对于感染的时间_T)##7 154

ggplot(data=sample_data1_raw2,aes(x=手术日相对于感染的时间_T)) + 
  geom_histogram(binwidth = 15,colour="black", fill="Yellow")+
  stat_bin(binwidth=15,geom="text",aes(label=..count..),vjust = -1) +ylim(0,55)+
  scale_x_continuous(breaks = seq(0,max(sample_data1_raw2$手术日相对于感染的时间_T),15),labels = comma)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="Vill and Decidual both in RNA later sample number in each infection days")


All_sample<-as.data.frame(table(sample_data1_both$手术日相对于感染的时间_T))
All_sample$day<-as.numeric(as.character(All_sample$Var1))
All_sample<-All_sample[which(All_sample$day >0 & All_sample$day <200),]
ggplot(All_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "Yellow") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="Vill and Decidual both in RNA later sample number in each infection days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(All_sample$Freq)+2)+
  theme(axis.text.x = element_text(size =8,colour = 'black',vjust=0.5,hjust=1,angle = 90))

#######################infect_lmp days distribution#############
sample_data1_raw2<-sample_data1_both[which(sample_data1_both$末次月经相较于感染日_T != "NA"),]
dim(sample_data1_raw2)##191  21

ggplot(data=sample_data1_raw2,aes(x=末次月经相较于感染日_T)) + 
  geom_histogram(binwidth = 15,colour="black", fill="green")+
  stat_bin(binwidth=15,geom="text",aes(label=..count..),vjust = -1) +ylim(0,65)+xlim(-30,155)+
  scale_x_continuous(breaks = seq(-30,max(sample_data1_raw2$末次月经相较于感染日_T),15),labels = comma)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="Vill and Decidual both in RNA later sample number in LMP to infection days")


All_sample<-as.data.frame(table(sample_data1_raw2$末次月经相较于感染日_T))
All_sample$day<-as.numeric(as.character(All_sample$Var1))

ggplot(All_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "green") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="Vill and Decidual both in RNA later sample number in  LMP to infection days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(All_sample$Freq)+5)+
  theme(axis.text.x = element_text(size =8,colour = 'black',vjust=0.5,hjust=1,angle = 90))

#######################infect_lmp days distribution#############
sample_data1_raw2<-sample_data1_both[which(sample_data1_both$末次月经相较于感染日_T != "NA"),];dim(sample_data1_raw2)
sample_data1_raw2<-sample_data1_raw2[which(sample_data1_raw2$末次月经相较于感染日_T %in% c(-31,-2,0)),];dim(sample_data1_raw2)
range(sample_data1_raw2$末次月经相较于感染日_T)
ggplot(data=sample_data1_raw2,aes(x=末次月经相较于感染日_T)) + 
  geom_histogram(binwidth = 15,colour="black", fill="green")+
  stat_bin(binwidth=15,geom="text",aes(label=..count..),vjust = -1) +ylim(0,65)+xlim(0,150)+
  scale_x_continuous(breaks = seq(-30,max(sample_data1_raw2$末次月经相较于感染日_T),15),labels = comma)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="Vill and Decidual both in RNA later sample number in LMP to infection days")


All_sample<-as.data.frame(table(sample_data1_raw2$末次月经相较于感染日_T))
All_sample$day<-as.numeric(as.character(All_sample$Var1))

ggplot(All_sample, aes(x=factor(Var1), y=Freq)) + geom_line(group = 1,linewidth=1,color= "green") + geom_point(size=2, shape=15)+
  theme_classic() + labs(x = "Time(days)", y = "sample_num(n)", title ="Vill and Decidual both in RNA later sample number in  LMP to infection days")+
  geom_text(aes(y = Freq, label= Freq), hjust= 0.5, vjust = -2, size= 4, color= "black") +ylim(0,max(All_sample$Freq)+5)+
  theme(axis.text.x = element_text(size =8,colour = 'black',vjust=0.5,hjust=1,angle = 90))
