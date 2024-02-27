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


###########evaluation between FACS and SC data####################
##for Decidua sample
FACS_sc_Decidua_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/absolute_value/final_relationship_between_FACS_SC_three_sample.txt",header = T,sep = "\t")

head(FACS_sc_Decidua_cell_ratio)
FACS_immune<-FACS_sc_Decidua_cell_ratio[,c(1:5)]
colnames(FACS_immune)<-c("Bulk_code","Ts","NKs","Macs","DCs")
head(FACS_immune)

SC_immune<-FACS_sc_Decidua_cell_ratio[,c(1,6:9)]
colnames(SC_immune)<-c("Bulk_code","NKs","Ts","Macs","DCs")
head(SC_immune)

FACS_data <-melt(FACS_immune,id=c("Bulk_code"),variable.name="Cell",value.name="ratio_FACS")
head(FACS_data)
SC_data <-melt(SC_immune,id=c("Bulk_code"),variable.name="Cell",value.name="ratio_SC")
head(SC_data)

FACS_data$sample_cell<-paste0(FACS_data$Bulk_code,"-",FACS_data$Cell)
SC_data$sample_cell<-paste0(SC_data$Bulk_code,"-",SC_data$Cell)
merge_data<-merge(FACS_data[,c(4,3)],SC_data[,c(4,3)],by="sample_cell")

merge_data$sample<-as.character(unlist(lapply(strsplit(merge_data$sample_cell,"-"), function(x) x[1])))
merge_data$cell<-as.character(unlist(lapply(strsplit(merge_data$sample_cell,"-"), function(x) x[2])))
head(merge_data)

##输出相关性统计结果
correlation_results <- merge_data %>%  group_by(cell) %>% summarize(correlation = cor(ratio_FACS, ratio_SC,method = "spearman"), p_value = cor.test(ratio_FACS, ratio_SC,method = "spearman")$p.value)

write.table(as.data.frame(merge_data),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/absolute_value/file1_Immune_Cell_ratio_decidual_facs_CsX_data_",index,".txt"),sep = "\t",row.names=F) 
write.table(as.data.frame(correlation_results),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/absolute_value/file2_Immune_Cell_ratio_decidual_facs_CsX_correlationship_data_",index,".txt"),sep = "\t",row.names=F) 
correlation_results[which(correlation_results$p_value<0.05),]

index_plot<-ggplot(merge_data, aes(x = ratio_FACS, y = ratio_SC)) +
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(method = "spearman")+
  xlab("Facs_Cell ratio(%)")+ylab("scRNA_Cell ratio(%)") +labs(title = "Decidua_correlationship")+
  facet_wrap(~ cell, scales = "free",ncol =2)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/absolute_value/file3_Immune_Cell_ratio_decidual_facs_CsX_correlationship_",index,".pdf"),width = 8, height =8)

###########evaluation between Facs and psudo bulk RNA####################
index<-"manu_sc_1k_abs"
##for Decidua sample
FACS_Decidua_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/absolute_value/Final_relationship_between_FACS_sc_bulk_decode.txt",header = T,sep = "\t")
CIBERSORTx_pesudo_Decidua <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/absolute_value/CIBERSORTx_cell_ratio_manu_significnat_1K_adjust.txt",header = T,sep = "\t")

head(FACS_Decidua_cell_ratio)
FACS_immune<-FACS_Decidua_cell_ratio[,c(9,3:8)]
colnames(FACS_immune)<-c("Bulk_code","Bs","Ts","NKs", "Myloids","Macs","DCs")
head(FACS_immune)

final_CsX_value_data<-CIBERSORTx_pesudo_Decidua[which(CIBERSORTx_pesudo_Decidua$Mixture %in% FACS_Decidua_cell_ratio$Bulk_code),]
head(final_CsX_value_data)
table_immune<-final_CsX_value_data[,c(1,8:12)]
table_immune2<-as.data.frame(t(apply(table_immune[,2:ncol(table_immune)],1,function(x) round(prop.table(x)*100,2))))
table_immune2$Bulk_code<-as.character(table_immune$Mixture)
head(table_immune2)

FACS_data <-melt(FACS_immune,id=c("Bulk_code"),variable.name="Cell",value.name="ratio_FACS")
head(FACS_data)
Cx_data <-melt(table_immune2,id=c("Bulk_code"),variable.name="Cell",value.name="ratio_CX")
head(Cx_data)
FACS_data$sample_cell<-paste0(FACS_data$Bulk_code,"-",FACS_data$Cell)
Cx_data$sample_cell<-paste0(Cx_data$Bulk_code,"-",Cx_data$Cell)
merge_data<-merge(FACS_data[,c(4,3)],Cx_data[,c(4,3)],by="sample_cell")

merge_data$sample<-as.character(unlist(lapply(strsplit(merge_data$sample_cell,"-"), function(x) x[1])))
merge_data$cell<-as.character(unlist(lapply(strsplit(merge_data$sample_cell,"-"), function(x) x[2])))
head(merge_data)

##输出相关性统计结果
correlation_results <- merge_data %>%  group_by(cell) %>% summarize(correlation = cor(ratio_FACS, ratio_CX,method = "spearman"), p_value = cor.test(ratio_FACS, ratio_CX,method = "spearman")$p.value)

write.table(as.data.frame(merge_data),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/absolute_value/file1_Immune_Cell_ratio_decidual_facs_CsX_data_",index,".txt"),sep = "\t",row.names=F) 
write.table(as.data.frame(correlation_results),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/absolute_value/file2_Immune_Cell_ratio_decidual_facs_CsX_correlationship_data_",index,".txt"),sep = "\t",row.names=F) 
correlation_results[which(correlation_results$p_value<0.05),]

index_plot<-ggplot(merge_data, aes(x = ratio_FACS, y = ratio_CX)) +
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(method = "spearman")+
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "Decidua_correlationship")+
  facet_wrap(~ cell, scales = "free",ncol =2)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/absolute_value/file3_Immune_Cell_ratio_decidual_facs_CsX_correlationship_",index,".pdf"),width = 8, height =8)

