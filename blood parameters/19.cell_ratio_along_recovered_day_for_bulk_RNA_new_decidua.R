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

###########evaluation using psudo bulk RNA####################
True_sc_Decidua_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/True_sc_Decidua_cell_ratio.txt",header = T,sep = "\t")
#
flag_tag<-"manu_sig_no_1K"
CIBERSORTx_pesudo_Decidua <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_no_sc_adjust/CIBERSORTx_pesudo_Mixture_Decidua.txt",header = T,sep = "\t")

head(True_sc_Decidua_cell_ratio)
head(CIBERSORTx_pesudo_Decidua)

True_data <-melt(True_sc_Decidua_cell_ratio,id=c("sample_code"),variable.name="Cell",value.name="ratio_true")
head(True_data)
CB_data <-melt(CIBERSORTx_pesudo_Decidua,id=c("Mixture"),variable.name="Cell",value.name="ratio_CB")
head(CB_data)
CB_data$ratio_CB<-CB_data$ratio_CB*100
colnames(True_data)<-c("Mixture","Cell","ratio_true")
com_data<-merge(True_data,CB_data,by=c("Mixture","Cell"))
head(com_data)
table(com_data$Cell)
##绘制相关性图
##输出相关性统计结果
correlation_results <- com_data %>%  group_by(Cell) %>% summarize(correlation = cor(ratio_true, ratio_CB,method = "spearman"), p_value = cor.test(ratio_true, ratio_CB,method = "spearman")$p.value)

write.table(as.data.frame(com_data),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_no_sc_adjust/file1_True_SC_Cell_ratio_CsX_Ratio_value_",flag_tag,".txt"),sep = "\t",row.names=F) 
write.table(as.data.frame(correlation_results),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_no_sc_adjust/file2_True_SC_Cell_ratio_CsX_correlationship_data_",flag_tag,".txt"),sep = "\t",row.names=F) 
correlation_results[which(correlation_results$p_value<0.05),]

index_plot<-ggplot(com_data, aes(x = ratio_true, y =  ratio_CB)) +
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(method = "spearman")+
  xlab("sc_RNA_Cell ratio(%)")+ylab("CsX_Cell ratio(%)") +labs(title = "Decidua_correlationship")+
  facet_wrap(~ Cell, scales = "free",ncol =4)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_no_sc_adjust/file3_True_SC_Cell_ratio_CsX_correlationship_",flag_tag,".pdf"),width = 12, height =10)

##线性相关性分析1）计算相关性 reference:https://zhuanlan.zhihu.com/p/543768987
#ref:https://www.coder.work/article/6652294#google_vignette
#ref:https://www.saoniuhuo.com/question/detail-2441094.html
#https://rdrr.io/cran/ggpmisc/man/stat_poly_eq.html
cor.test(com_data$ratio_true,com_data$ratio_CB,data=com_data)

###根据结果可以知道相关性系数为0.9374051，且呈正相关，P值为2.2e-16 <0.05
###回归系数的95% 置信区间为[ 0.9072704 0.9579627]
df_cor<-lm(ratio_CB~ratio_true,data=com_data)
summary(df_cor)
###回归方程为ratio_CB=0.9934360*ratio_CB+0.0003455 
###拟合优度R^2=0.8787,即拟合度很好 
##模型 p 值，而不是斜率 p 值(通常) p-value: < 2.2e-16

range(com_data$ratio_true)# 0.0000 74.7855
range(com_data$ratio_CB)#  0.00000 83.21935
formula <- y ~ x
index_plot<-ggplot(com_data, aes(x = ratio_true, y = ratio_CB)) + #color = group,
  xlab("True SC Decidua Cell ratio")+ylab("CIBERSORTx Decidua Cell ratio") +labs(title = "Decidua correlationship(spearman)")+
  scale_x_continuous(breaks = seq(0,100,10)) + 
  scale_color_manual(values=ppCor)+
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  #geom_text(aes(x=55,y=124,label="R^2=0.88\ny=0.75x+58.96"),color="red",family = "serif",fontface = "plain",size = 5)+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+##method one
  # stat_cor(aes(label = after_stat(rr.label)), color = "red", geom = "label")+##method two
  #stat_fit_glance(method = 'lm', method.args = list(formula = formula),geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.35, size = 3)+
  geom_point(alpha = 0.5,aes(color = Mixture,size =2)) +  
  ggpubr::stat_cor(method = "spearman",size =4,label.x = 10,label.y = 80)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,..p.value.label..,sep ="~~~")),size =4,colour = "black", formula = formula, label.x = 0.1,label.y = 0.7, parse = TRUE)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_no_sc_adjust/file4_True_SC_Cell_ratio_CsX_correlationship_no_split_",flag_tag,".pdf"),width = 12, height =10)

#
flag_tag<-"manu_sig_1K"
CIBERSORTx_pesudo_Decidua <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_sc_adjust/CIBERSORTx_pesudo_Mixture_Decidua_manu_1K.txt",header = T,sep = "\t")

head(True_sc_Decidua_cell_ratio)
head(CIBERSORTx_pesudo_Decidua)

True_data <-melt(True_sc_Decidua_cell_ratio,id=c("sample_code"),variable.name="Cell",value.name="ratio_true")
head(True_data)
CB_data <-melt(CIBERSORTx_pesudo_Decidua,id=c("Mixture"),variable.name="Cell",value.name="ratio_CB")
head(CB_data)
CB_data$ratio_CB<-CB_data$ratio_CB*100
colnames(True_data)<-c("Mixture","Cell","ratio_true")
com_data<-merge(True_data,CB_data,by=c("Mixture","Cell"))
head(com_data)
table(com_data$Cell)
##绘制相关性图
##输出相关性统计结果
correlation_results <- com_data %>%  group_by(Cell) %>% summarize(correlation = cor(ratio_true, ratio_CB,method = "spearman"), p_value = cor.test(ratio_true, ratio_CB,method = "spearman")$p.value)

write.table(as.data.frame(com_data),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_sc_adjust/file1_True_SC_Cell_ratio_CsX_Ratio_value_",flag_tag,".txt"),sep = "\t",row.names=F) 
write.table(as.data.frame(correlation_results),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_sc_adjust/file2_True_SC_Cell_ratio_CsX_correlationship_data_",flag_tag,".txt"),sep = "\t",row.names=F) 
correlation_results[which(correlation_results$p_value<0.05),]

index_plot<-ggplot(com_data, aes(x = ratio_true, y =  ratio_CB)) +
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(method = "spearman")+
  xlab("sc_RNA_Cell ratio(%)")+ylab("CsX_Cell ratio(%)") +labs(title = "Decidua_correlationship")+
  facet_wrap(~ Cell, scales = "free",ncol =4)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_sc_adjust/file3_True_SC_Cell_ratio_CsX_correlationship_",flag_tag,".pdf"),width = 12, height =10)

cor.test(com_data$ratio_true,com_data$ratio_CB,data=com_data)

###根据结果可以知道相关性系数为0.9374051，且呈正相关，P值为2.2e-16 <0.05
###回归系数的95% 置信区间为[ 0.9072704 0.9579627]
df_cor<-lm(ratio_CB~ratio_true,data=com_data)
summary(df_cor)
###回归方程为ratio_CB=0.9934360*ratio_CB+0.0003455 
###拟合优度R^2=0.8787,即拟合度很好 
##模型 p 值，而不是斜率 p 值(通常) p-value: < 2.2e-16

range(com_data$ratio_true)# 0.0000 74.7855
range(com_data$ratio_CB)#  0.00000 83.21935
formula <- y ~ x
index_plot<-ggplot(com_data, aes(x = ratio_true, y = ratio_CB)) + #color = group,
  xlab("True SC Decidua Cell ratio")+ylab("CIBERSORTx Decidua Cell ratio") +labs(title = "Decidua correlationship(spearman)")+
  scale_x_continuous(breaks = seq(0,100,10)) + 
  scale_color_manual(values=ppCor)+
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  #geom_text(aes(x=55,y=124,label="R^2=0.88\ny=0.75x+58.96"),color="red",family = "serif",fontface = "plain",size = 5)+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+##method one
  # stat_cor(aes(label = after_stat(rr.label)), color = "red", geom = "label")+##method two
  #stat_fit_glance(method = 'lm', method.args = list(formula = formula),geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.35, size = 3)+
  geom_point(alpha = 0.5,aes(color = Mixture,size =2)) +  
  ggpubr::stat_cor(method = "spearman",size =4,label.x = 10,label.y = 80)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,..p.value.label..,sep ="~~~")),size =4,colour = "black", formula = formula, label.x = 0.1,label.y = 0.7, parse = TRUE)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_sc_adjust/file4_True_SC_Cell_ratio_CsX_correlationship_no_split_",flag_tag,".pdf"),width = 12, height =10)

#
flag_tag<-"manu_NEW_no_adj"
CIBERSORTx_pesudo_Decidua <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_pesudo_Mixture_Decidua_manu_NEW.txt",header = T,sep = "\t")

head(True_sc_Decidua_cell_ratio)
head(CIBERSORTx_pesudo_Decidua)

True_data <-melt(True_sc_Decidua_cell_ratio,id=c("sample_code"),variable.name="Cell",value.name="ratio_true")
head(True_data)
CB_data <-melt(CIBERSORTx_pesudo_Decidua,id=c("Mixture"),variable.name="Cell",value.name="ratio_CB")
head(CB_data)
CB_data$ratio_CB<-CB_data$ratio_CB*100
colnames(True_data)<-c("Mixture","Cell","ratio_true")
com_data<-merge(True_data,CB_data,by=c("Mixture","Cell"))
head(com_data)
table(com_data$Cell)
##绘制相关性图
##输出相关性统计结果
correlation_results <- com_data %>%  group_by(Cell) %>% summarize(correlation = cor(ratio_true, ratio_CB,method = "spearman"), p_value = cor.test(ratio_true, ratio_CB,method = "spearman")$p.value)
correlation_results[which(correlation_results$p_value<0.05),]

write.table(as.data.frame(com_data),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_no_adjust/file1_True_SC_Cell_ratio_CsX_Ratio_value_",flag_tag,".txt"),sep = "\t",row.names=F) 
write.table(as.data.frame(correlation_results),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_no_adjust/file2_True_SC_Cell_ratio_CsX_correlationship_data_",flag_tag,".txt"),sep = "\t",row.names=F) 

index_plot<-ggplot(com_data, aes(x = ratio_true, y =  ratio_CB)) +
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(method = "spearman")+
  xlab("sc_RNA_Cell ratio(%)")+ylab("CsX_Cell ratio(%)") +labs(title = "Decidua_correlationship")+
  facet_wrap(~ Cell, scales = "free",ncol =4)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_no_adjust/file3_True_SC_Cell_ratio_CsX_correlationship_",flag_tag,".pdf"),width = 12, height =10)

cor.test(com_data$ratio_true,com_data$ratio_CB,data=com_data)
###根据结果可以知道相关性系数为0.9374051，且呈正相关，P值为2.2e-16 <0.05
###回归系数的95% 置信区间为[ 0.9072704 0.9579627]
df_cor<-lm(ratio_CB~ratio_true,data=com_data)
summary(df_cor)
###回归方程为ratio_CB=0.9934360*ratio_CB+0.0003455 
###拟合优度R^2=0.8787,即拟合度很好 
##模型 p 值，而不是斜率 p 值(通常) p-value: < 2.2e-16

range(com_data$ratio_true)# 0.0000 74.7855
range(com_data$ratio_CB)#  0.00000 83.21935
formula <- y ~ x
index_plot<-ggplot(com_data, aes(x = ratio_true, y = ratio_CB)) + #color = group,
  xlab("True SC Decidua Cell ratio")+ylab("CIBERSORTx Decidua Cell ratio") +labs(title = "Decidua correlationship(spearman)")+
  scale_x_continuous(breaks = seq(0,100,10)) + 
  scale_color_manual(values=ppCor)+
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  #geom_text(aes(x=55,y=124,label="R^2=0.88\ny=0.75x+58.96"),color="red",family = "serif",fontface = "plain",size = 5)+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+##method one
  # stat_cor(aes(label = after_stat(rr.label)), color = "red", geom = "label")+##method two
  #stat_fit_glance(method = 'lm', method.args = list(formula = formula),geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.35, size = 3)+
  geom_point(alpha = 0.5,aes(color = Mixture,size =2)) +  
  ggpubr::stat_cor(method = "spearman",size =4,label.x = 10,label.y = 80)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,..p.value.label..,sep ="~~~")),size =4,colour = "black", formula = formula, label.x = 0.1,label.y = 0.7, parse = TRUE)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_no_adjust/file4_True_SC_Cell_ratio_CsX_correlationship_no_split_",flag_tag,".pdf"),width = 12, height =10)

#
flag_tag<-"manu_NEW_sig_1K"
CIBERSORTx_pesudo_Decidua <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_adjust/CIBERSORTx_pesudo_Mixture_Decidua_manu_NEW_1K.txt",header = T,sep = "\t")

head(True_sc_Decidua_cell_ratio)
head(CIBERSORTx_pesudo_Decidua)

True_data <-melt(True_sc_Decidua_cell_ratio,id=c("sample_code"),variable.name="Cell",value.name="ratio_true")
head(True_data)
CB_data <-melt(CIBERSORTx_pesudo_Decidua,id=c("Mixture"),variable.name="Cell",value.name="ratio_CB")
head(CB_data)
CB_data$ratio_CB<-CB_data$ratio_CB*100
colnames(True_data)<-c("Mixture","Cell","ratio_true")
com_data<-merge(True_data,CB_data,by=c("Mixture","Cell"))
head(com_data)
table(com_data$Cell)
##绘制相关性图
##输出相关性统计结果
correlation_results <- com_data %>%  group_by(Cell) %>% summarize(correlation = cor(ratio_true, ratio_CB,method = "spearman"), p_value = cor.test(ratio_true, ratio_CB,method = "spearman")$p.value)

write.table(as.data.frame(com_data),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_adjust/file1_True_SC_Cell_ratio_CsX_Ratio_value_",flag_tag,".txt"),sep = "\t",row.names=F) 
write.table(as.data.frame(correlation_results),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_adjust/file2_True_SC_Cell_ratio_CsX_correlationship_data_",flag_tag,".txt"),sep = "\t",row.names=F) 
correlation_results[which(correlation_results$p_value<0.05),]

index_plot<-ggplot(com_data, aes(x = ratio_true, y =  ratio_CB)) +
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(method = "spearman")+
  xlab("sc_RNA_Cell ratio(%)")+ylab("CsX_Cell ratio(%)") +labs(title = "Decidua_correlationship")+
  facet_wrap(~ Cell, scales = "free",ncol =4)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_adjust/file3_True_SC_Cell_ratio_CsX_correlationship_",flag_tag,".pdf"),width = 12, height =10)

##线性相关性分析1）计算相关性 reference:https://zhuanlan.zhihu.com/p/543768987
#ref:https://www.coder.work/article/6652294#google_vignette
#ref:https://www.saoniuhuo.com/question/detail-2441094.html
#https://rdrr.io/cran/ggpmisc/man/stat_poly_eq.html
cor.test(com_data$ratio_true,com_data$ratio_CB,data=com_data)

###根据结果可以知道相关性系数为0.9374051，且呈正相关，P值为2.2e-16 <0.05
###回归系数的95% 置信区间为[ 0.9072704 0.9579627]
df_cor<-lm(ratio_CB~ratio_true,data=com_data)
summary(df_cor)
###回归方程为ratio_CB=0.9934360*ratio_CB+0.0003455 
###拟合优度R^2=0.8787,即拟合度很好 
##模型 p 值，而不是斜率 p 值(通常) p-value: < 2.2e-16

range(com_data$ratio_true)# 0.0000 74.7855
range(com_data$ratio_CB)#  0.00000 83.21935
formula <- y ~ x
index_plot<-ggplot(com_data, aes(x = ratio_true, y = ratio_CB)) + #color = group,
  xlab("True SC Decidua Cell ratio")+ylab("CIBERSORTx Decidua Cell ratio") +labs(title = "Decidua correlationship(spearman)")+
  scale_x_continuous(breaks = seq(0,100,10)) + 
  scale_color_manual(values=ppCor)+
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  #geom_text(aes(x=55,y=124,label="R^2=0.88\ny=0.75x+58.96"),color="red",family = "serif",fontface = "plain",size = 5)+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+##method one
  # stat_cor(aes(label = after_stat(rr.label)), color = "red", geom = "label")+##method two
  #stat_fit_glance(method = 'lm', method.args = list(formula = formula),geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.35, size = 3)+
  geom_point(alpha = 0.5,aes(color = Mixture,size =2)) +  
  ggpubr::stat_cor(method = "spearman",size =4,label.x = 10,label.y = 80)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,..p.value.label..,sep ="~~~")),size =4,colour = "black", formula = formula, label.x = 0.1,label.y = 0.7, parse = TRUE)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_manu_NEW_sc_adjust/file4_True_SC_Cell_ratio_CsX_correlationship_no_split_",flag_tag,".pdf"),width = 12, height =10)

#
flag_tag<-"CsX_sig_1K"
CIBERSORTx_pesudo_Decidua <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_adjust/CIBERSORTx_pesudo_Mixture_Decidua_CsX_1K.txt",header = T,sep = "\t")

head(True_sc_Decidua_cell_ratio)
head(CIBERSORTx_pesudo_Decidua)

True_data <-melt(True_sc_Decidua_cell_ratio,id=c("sample_code"),variable.name="Cell",value.name="ratio_true")
head(True_data)
CB_data <-melt(CIBERSORTx_pesudo_Decidua,id=c("Mixture"),variable.name="Cell",value.name="ratio_CB")
head(CB_data)
CB_data$ratio_CB<-CB_data$ratio_CB*100
colnames(True_data)<-c("Mixture","Cell","ratio_true")
com_data<-merge(True_data,CB_data,by=c("Mixture","Cell"))
head(com_data)
table(com_data$Cell)
##绘制相关性图
##输出相关性统计结果
correlation_results <- com_data %>%  group_by(Cell) %>% summarize(correlation = cor(ratio_true, ratio_CB,method = "spearman"), p_value = cor.test(ratio_true, ratio_CB,method = "spearman")$p.value)
correlation_results[which(correlation_results$p_value<0.05),]

write.table(as.data.frame(com_data),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_adjust/file1_True_SC_Cell_ratio_CsX_Ratio_value_",flag_tag,".txt"),sep = "\t",row.names=F) 
write.table(as.data.frame(correlation_results),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_adjust/file2_True_SC_Cell_ratio_CsX_correlationship_data_",flag_tag,".txt"),sep = "\t",row.names=F) 

index_plot<-ggplot(com_data, aes(x = ratio_true, y =  ratio_CB)) +
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(method = "spearman")+
  xlab("sc_RNA_Cell ratio(%)")+ylab("CsX_Cell ratio(%)") +labs(title = "Decidua_correlationship")+
  facet_wrap(~ Cell, scales = "free",ncol =4)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_adjust/file3_True_SC_Cell_ratio_CsX_correlationship_",flag_tag,".pdf"),width = 12, height =10)

cor.test(com_data$ratio_true,com_data$ratio_CB,data=com_data)
###根据结果可以知道相关性系数为0.9374051，且呈正相关，P值为2.2e-16 <0.05
###回归系数的95% 置信区间为[ 0.9072704 0.9579627]
df_cor<-lm(ratio_CB~ratio_true,data=com_data)
summary(df_cor)
###回归方程为ratio_CB=0.9934360*ratio_CB+0.0003455 
###拟合优度R^2=0.8787,即拟合度很好 
##模型 p 值，而不是斜率 p 值(通常) p-value: < 2.2e-16

range(com_data$ratio_true)# 0.0000 74.7855
range(com_data$ratio_CB)#  0.00000 83.21935
formula <- y ~ x
index_plot<-ggplot(com_data, aes(x = ratio_true, y = ratio_CB)) + #color = group,
  xlab("True SC Decidua Cell ratio")+ylab("CIBERSORTx Decidua Cell ratio") +labs(title = "Decidua correlationship(spearman)")+
  scale_x_continuous(breaks = seq(0,100,10)) + 
  scale_color_manual(values=ppCor)+
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  #geom_text(aes(x=55,y=124,label="R^2=0.88\ny=0.75x+58.96"),color="red",family = "serif",fontface = "plain",size = 5)+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+##method one
  # stat_cor(aes(label = after_stat(rr.label)), color = "red", geom = "label")+##method two
  #stat_fit_glance(method = 'lm', method.args = list(formula = formula),geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.35, size = 3)+
  geom_point(alpha = 0.5,aes(color = Mixture,size =2)) +  
  ggpubr::stat_cor(method = "spearman",size =4,label.x = 10,label.y = 80)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,..p.value.label..,sep ="~~~")),size =4,colour = "black", formula = formula, label.x = 0.1,label.y = 0.7, parse = TRUE)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_adjust/file4_True_SC_Cell_ratio_CsX_correlationship_no_split_",flag_tag,".pdf"),width = 12, height =10)

#
flag_tag<-"CsX_no_adj"
CIBERSORTx_pesudo_Decidua <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_no_adjust/CIBERSORTx_pesudo_Mixture_Decidua_CsX_no_adj.txt",header = T,sep = "\t")

head(True_sc_Decidua_cell_ratio)
head(CIBERSORTx_pesudo_Decidua)

True_data <-melt(True_sc_Decidua_cell_ratio,id=c("sample_code"),variable.name="Cell",value.name="ratio_true")
head(True_data)
CB_data <-melt(CIBERSORTx_pesudo_Decidua,id=c("Mixture"),variable.name="Cell",value.name="ratio_CB")
head(CB_data)
CB_data$ratio_CB<-CB_data$ratio_CB*100
colnames(True_data)<-c("Mixture","Cell","ratio_true")
com_data<-merge(True_data,CB_data,by=c("Mixture","Cell"))
head(com_data)
table(com_data$Cell)
##绘制相关性图
##输出相关性统计结果
correlation_results <- com_data %>%  group_by(Cell) %>% summarize(correlation = cor(ratio_true, ratio_CB,method = "spearman"), p_value = cor.test(ratio_true, ratio_CB,method = "spearman")$p.value)
correlation_results[which(correlation_results$p_value<0.05),]

write.table(as.data.frame(com_data),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_no_adjust/file1_True_SC_Cell_ratio_CsX_Ratio_value_",flag_tag,".txt"),sep = "\t",row.names=F) 
write.table(as.data.frame(correlation_results),file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_no_adjust/file2_True_SC_Cell_ratio_CsX_correlationship_data_",flag_tag,".txt"),sep = "\t",row.names=F) 

index_plot<-ggplot(com_data, aes(x = ratio_true, y =  ratio_CB)) +
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  geom_point(alpha = 1,size =1.5,color = "black") +  ggpubr::stat_cor(method = "spearman")+
  xlab("sc_RNA_Cell ratio(%)")+ylab("CsX_Cell ratio(%)") +labs(title = "Decidua_correlationship")+
  facet_wrap(~ Cell, scales = "free",ncol =4)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_no_adjust/file3_True_SC_Cell_ratio_CsX_correlationship_",flag_tag,".pdf"),width = 12, height =10)

cor.test(com_data$ratio_true,com_data$ratio_CB,data=com_data)
###根据结果可以知道相关性系数为0.9374051，且呈正相关，P值为2.2e-16 <0.05
###回归系数的95% 置信区间为[ 0.9072704 0.9579627]
df_cor<-lm(ratio_CB~ratio_true,data=com_data)
summary(df_cor)
###回归方程为ratio_CB=0.9934360*ratio_CB+0.0003455 
###拟合优度R^2=0.8787,即拟合度很好 
##模型 p 值，而不是斜率 p 值(通常) p-value: < 2.2e-16

range(com_data$ratio_true)# 0.0000 74.7855
range(com_data$ratio_CB)#  0.00000 83.21935
formula <- y ~ x
index_plot<-ggplot(com_data, aes(x = ratio_true, y = ratio_CB)) + #color = group,
  xlab("True SC Decidua Cell ratio")+ylab("CIBERSORTx Decidua Cell ratio") +labs(title = "Decidua correlationship(spearman)")+
  scale_x_continuous(breaks = seq(0,100,10)) + 
  scale_color_manual(values=ppCor)+
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  #geom_text(aes(x=55,y=124,label="R^2=0.88\ny=0.75x+58.96"),color="red",family = "serif",fontface = "plain",size = 5)+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+##method one
  # stat_cor(aes(label = after_stat(rr.label)), color = "red", geom = "label")+##method two
  #stat_fit_glance(method = 'lm', method.args = list(formula = formula),geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.35, size = 3)+
  geom_point(alpha = 0.5,aes(color = Mixture,size =2)) +  
  ggpubr::stat_cor(method = "spearman",size =4,label.x = 10,label.y = 80)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,..p.value.label..,sep ="~~~")),size =4,colour = "black", formula = formula, label.x = 0.1,label.y = 0.7, parse = TRUE)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/seurat_cluster/pesudo_Decidua_cell_ratio_CsX_sc_no_adjust/file4_True_SC_Cell_ratio_CsX_correlationship_no_split_",flag_tag,".pdf"),width = 12, height =10)

##########################################################
##for Decidua sample
Decidua_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_manu.txt",header = T,sep = "\t")
Decidua_cell_ratio2<-Decidua_cell_ratio[,-((ncol(Decidua_cell_ratio)-2):ncol(Decidua_cell_ratio))]
rowSums(Decidua_cell_ratio2[,-1])
analysis_data <-melt(Decidua_cell_ratio2,id=c("Mixture"),variable.name="Cell",value.name="ratio")
head(analysis_data)
colnames(analysis_data)<-c("sample","Cell","ratio")
#analysis_data$sample <- factor(analysis_data$sample,levels=paste0("T"))

Ratio_0<-ggplot(analysis_data,aes(sample,ratio,fill=Cell))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors)+guides(fill=guide_legend(title=NULL))
Ratio_0
ggsave(Ratio_0,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_barplot.pdf",width = 16, height =8)

##合并康复日期信息
colData_used <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/Final_fielt_Decidua_Vill_analysis_metadata_genderadd.txt",sep="\t",header =T)
rownames(colData_used)<-as.character(colData_used$sample_code)
dim(colData_used)#245  13
head(colData_used)

colData_Decidua<-colData_used[which(colData_used$sample_code %in% as.character(Decidua_cell_ratio2$Mixture)),]
head(colData_Decidua);dim(colData_Decidua)#125  13
###############################
colData_Decidua$group_day2<-ifelse(colData_Decidua$infect_state =="No","CTRL",
                                   ifelse(colData_Decidua$LMP_infect <= 7,"D7",
                                          ifelse(colData_Decidua$LMP_infect <= 21,"D21",
                                                 ifelse(colData_Decidua$LMP_infect <= 35,"D35",
                                                        ifelse(colData_Decidua$LMP_infect <= 49,"D49",
                                                               ifelse(colData_Decidua$LMP_infect <= 63,"D63",
                                                                      ifelse(colData_Decidua$LMP_infect <= 77,"D77","Dmore77")))))))
table(colData_Decidua$group_day2)
class_order<-c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77")
colData_Decidua$group_day2<-factor(colData_Decidua$group_day2,levels = class_order)

Decidua_cell_ratio2$sample_code <-Decidua_cell_ratio2$Mixture
analysis_Decidua<-merge(colData_Decidua,Decidua_cell_ratio2,by="sample_code")
analysis_Decidua[1:5,1:16]
analysis_Decidua$infect_state <- factor(analysis_Decidua$infect_state,levels=c("Infect","No"))
analysis_Decidua$LMP_infect<-as.numeric(analysis_Decidua$LMP_infect)

table(analysis_Decidua$Year_month)
head(analysis_Decidua)
dim(distinct(analysis_Decidua[,c("sample","infect_state")]))# 125  2

##感染与未感染的组间比较
#因子在感染与非感染组的组间比较
analysis_collect0<-analysis_Decidua[,c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2",colnames(Decidua_cell_ratio2[,-c(1,ncol(Decidua_cell_ratio2))]))]
head(analysis_collect0)
data_anlysis0 <- melt(analysis_collect0,variable.name="Cell",value.name = "Cell_ratio",id.vars = c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2"))
data_anlysis0$age_group<-ifelse(data_anlysis0$Age>=35,"AMA","YMA")
head(data_anlysis0)
#data_anlysis0<-data_anlysis0[which(data_anlysis0$Cell_ratio>0),]

stat_data<-compare_means(Cell_ratio~infect_state, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]

stat_boxplot<-ggplot(data_anlysis0, aes(x=infect_state,y=Cell_ratio,fill=infect_state))+
  geom_boxplot(position=position_dodge(),width=0.5)+geom_jitter(width = 0.2,size=1,color="grey",alpha=0.8)+ 
  scale_color_manual(values=ppCor)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_compare_means(label="p.format", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= infect_state), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ Cell, scales = "free",ncol =4)
stat_boxplot
#ggsave(stat_boxplot,file=paste0("/mnt/data/chenwei/covid19/bulk_RNA_result/1.rm_long/Decidua_result/Decidua_compare_stat_box_infect_noremove_",index_name,".pdf"),width = 16, height =8)
ggsave(stat_boxplot,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_compare_stat_box_infect_state.pdf",width = 16, height =12)

################# 比较不同区间分组
head(data_anlysis0)
stat_data<-compare_means(Cell_ratio~group_day2, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]
# A tibble:  12 x 9
#  Cell  .y.        group1 group2        p p.adj p.format p.signif method  
# 1 LECs  Cell_ratio D21    D63     0.0334      1 0.0334   *        Wilcoxon
# 2 Macs  Cell_ratio D77    Dmore77 0.0271      1 0.0271   *        Wilcoxon
# 3 DCs   Cell_ratio CTRL   D7      0.0365      1 0.0365   *        Wilcoxon
# 4 DCs   Cell_ratio D7     D21     0.0103      1 0.0103   *        Wilcoxon
# 5 DCs   Cell_ratio D7     D49     0.0183      1 0.0183   *        Wilcoxon
# 6 DCs   Cell_ratio D7     D63     0.0139      1 0.0139   *        Wilcoxon
# 7 DCs   Cell_ratio D7     Dmore77 0.00608     1 0.0061   **       Wilcoxon
# 8 NKs   Cell_ratio CTRL   D7      0.00711     1 0.0071   **       Wilcoxon
# 9 NKs   Cell_ratio D7     D21     0.0365      1 0.0365   *        Wilcoxon
#10 NKs   Cell_ratio D7     D35     0.00941     1 0.0094   **       Wilcoxon
#11 NKs   Cell_ratio D7     D49     0.00942     1 0.0094   **       Wilcoxon
#12 NKs   Cell_ratio D35    Dmore77 0.0322      1 0.0322   *        Wilcoxon

head(data_anlysis0)
Daygroup_stat_boxplot1 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =4)
Daygroup_stat_boxplot1
Daygroup_stat_boxplot10<-Daygroup_stat_boxplot1+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

Daygroup_stat_boxplot2 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =4)
Daygroup_stat_boxplot2
Daygroup_stat_boxplot20<-Daygroup_stat_boxplot2+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

ggsave(Daygroup_stat_boxplot1,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot2,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_D7_compare.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare2.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot20,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_D7_compare2.pdf",width = 25, height =22)

##绘制康复曲线
data_anlysis1<-data_anlysis0[which(data_anlysis0$infect_state == "Infect"),]
range(na.omit(data_anlysis1$LMP_infect))##-2 110
##输出相关性统计结果
correlation_results <- data_anlysis1 %>%  group_by(Cell) %>% summarize(correlation = cor(LMP_infect, Cell_ratio,method = "spearman"), p_value = cor.test(LMP_infect, Cell_ratio,method = "spearman")$p.value)
write.table(as.data.frame(correlation_results), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/Spearman_pvalue_each_cell_decidua_data.txt",sep = "\t",row.names=F) 

index_plot<-ggplot(data_anlysis1, aes(x = LMP_infect, y = Cell_ratio)) + #color = group,
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "Decidua_cell_ratio")+
  geom_rug(position = "jitter", size = 0.1, color = "black")+
  scale_x_continuous(breaks = seq(-5,120,10)) + facet_wrap(~ Cell, scales = "free",ncol =4)+
  scale_color_manual(values=my_morandi_colors)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot10<-index_plot+  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)+ geom_point(alpha = 0.5,aes(color = age_group,size =2))
index_plot30<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 0.5,aes(color = age_group,size =2)) +  ggpubr::stat_cor(label.x = 3,method = "spearman")
index_plot31<-index_plot+ geom_smooth(aes(color = age_group,fill=age_group), method = lm,se = TRUE, linewidth=2, alpha = 0.2, fullrange = TRUE)+ geom_point(alpha = 0.5,aes(color = age_group),size = 2) +   ggpubr::stat_cor(aes(color = age_group), label.x = 0,method = "spearman")

index_plot32<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+geom_smooth(method = lm, se = TRUE, fill = "grey",size=2, alpha = 0.5, fullrange = TRUE) +  ggpubr::stat_cor(label.x = 0,method = "spearman")
index_plot33<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+ geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)

#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
ggsave(index_plot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_sm_trend_plot.pdf",width = 16, height =12)
ggsave(index_plot30,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_lm_single_plot.pdf",width = 16, height =12)
ggsave(index_plot31,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_lm_age_split_plot.pdf",width = 16, height =12)
ggsave(index_plot32,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_lm_no_age_split_plot.pdf",width = 16, height =12)
ggsave(index_plot33,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_sm_no_age_split_plot.pdf",width = 16, height =12)

#########################
##########################################################
##for Decidua sample:manu_new_no_adjust
Decidua_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_manu_NEW.txt",header = T,sep = "\t")
Decidua_cell_ratio2<-Decidua_cell_ratio[,-((ncol(Decidua_cell_ratio)-2):ncol(Decidua_cell_ratio))]
rowSums(Decidua_cell_ratio2[,-1])
analysis_data <-melt(Decidua_cell_ratio2,id=c("Mixture"),variable.name="Cell",value.name="ratio")
head(analysis_data)
colnames(analysis_data)<-c("sample","Cell","ratio")
#analysis_data$sample <- factor(analysis_data$sample,levels=paste0("T"))

Ratio_0<-ggplot(analysis_data,aes(sample,ratio,fill=Cell))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors)+guides(fill=guide_legend(title=NULL))
Ratio_0
ggsave(Ratio_0,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_barplot.pdf",width = 16, height =8)

##合并康复日期信息
colData_used <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/Final_fielt_Decidua_Vill_analysis_metadata_genderadd.txt",sep="\t",header =T)
rownames(colData_used)<-as.character(colData_used$sample_code)
dim(colData_used)#245  13
head(colData_used)

colData_Decidua<-colData_used[which(colData_used$sample_code %in% as.character(Decidua_cell_ratio2$Mixture)),]
head(colData_Decidua);dim(colData_Decidua)#125  13
###############################
colData_Decidua$group_day2<-ifelse(colData_Decidua$infect_state =="No","CTRL",
                                   ifelse(colData_Decidua$LMP_infect <= 7,"D7",
                                          ifelse(colData_Decidua$LMP_infect <= 21,"D21",
                                                 ifelse(colData_Decidua$LMP_infect <= 35,"D35",
                                                        ifelse(colData_Decidua$LMP_infect <= 49,"D49",
                                                               ifelse(colData_Decidua$LMP_infect <= 63,"D63",
                                                                      ifelse(colData_Decidua$LMP_infect <= 77,"D77","Dmore77")))))))
table(colData_Decidua$group_day2)
class_order<-c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77")
colData_Decidua$group_day2<-factor(colData_Decidua$group_day2,levels = class_order)

Decidua_cell_ratio2$sample_code <-Decidua_cell_ratio2$Mixture
analysis_Decidua<-merge(colData_Decidua,Decidua_cell_ratio2,by="sample_code")
analysis_Decidua[1:5,1:16]
analysis_Decidua$infect_state <- factor(analysis_Decidua$infect_state,levels=c("Infect","No"))
analysis_Decidua$LMP_infect<-as.numeric(analysis_Decidua$LMP_infect)

table(analysis_Decidua$Year_month)
head(analysis_Decidua)
dim(distinct(analysis_Decidua[,c("sample","infect_state")]))# 125  2

##感染与未感染的组间比较
#因子在感染与非感染组的组间比较
analysis_collect0<-analysis_Decidua[,c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2",colnames(Decidua_cell_ratio2[,-c(1,ncol(Decidua_cell_ratio2))]))]
head(analysis_collect0)
data_anlysis0 <- melt(analysis_collect0,variable.name="Cell",value.name = "Cell_ratio",id.vars = c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2"))
data_anlysis0$age_group<-ifelse(data_anlysis0$Age>=35,"AMA","YMA")
head(data_anlysis0)
#data_anlysis0<-data_anlysis0[which(data_anlysis0$Cell_ratio>0),]

stat_data<-compare_means(Cell_ratio~infect_state, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]

stat_boxplot<-ggplot(data_anlysis0, aes(x=infect_state,y=Cell_ratio,fill=infect_state))+
  geom_boxplot(position=position_dodge(),width=0.5)+geom_jitter(width = 0.2,size=1,color="grey",alpha=0.8)+ 
  scale_color_manual(values=ppCor)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_compare_means(label="p.format", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= infect_state), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ Cell, scales = "free",ncol =4)
stat_boxplot
#ggsave(stat_boxplot,file=paste0("/mnt/data/chenwei/covid19/bulk_RNA_result/1.rm_long/Decidua_result/Decidua_compare_stat_box_infect_noremove_",index_name,".pdf"),width = 16, height =8)
ggsave(stat_boxplot,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_compare_stat_box_infect_state.pdf",width = 16, height =12)

################# 比较不同区间分组
head(data_anlysis0)
stat_data<-compare_means(Cell_ratio~group_day2, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]
# A tibble:  12 x 9
#   Cell  .y.        group1 group2       p p.adj p.format p.signif method  
#1 Trbs  Cell_ratio CTRL   D21    0.0473      1 0.0473   *        Wilcoxon
#2 Trbs  Cell_ratio D21    D49    0.0494      1 0.0494   *        Wilcoxon
#3 Trbs  Cell_ratio D21    D77    0.0339      1 0.0339   *        Wilcoxon
#4 FBs_1 Cell_ratio CTRL   D21    0.0473      1 0.0473   *        Wilcoxon
#5 NKs   Cell_ratio CTRL   D7     0.0103      1 0.0103   *        Wilcoxon
#6 NKs   Cell_ratio D7     D21    0.0103      1 0.0103   *        Wilcoxon
#7 NKs   Cell_ratio D7     D35    0.00776     1 0.0078   **       Wilcoxon
#8 NKs   Cell_ratio D7     D49    0.0183      1 0.0183   *        Wilcoxon
#9 NKs   Cell_ratio D7     D63    0.0208      1 0.0208   *        Wilcoxon
#10 Masts Cell_ratio D7     D35    0.00584     1 0.0058   **       Wilcoxon
#11 Masts Cell_ratio D7     D63    0.0124      1 0.0124   *        Wilcoxon
#12 Masts Cell_ratio D7     D77    0.0481      1 0.0481   *        Wilcoxon

head(data_anlysis0)
Daygroup_stat_boxplot1 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =4)
Daygroup_stat_boxplot1
Daygroup_stat_boxplot10<-Daygroup_stat_boxplot1+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

Daygroup_stat_boxplot2 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =4)
Daygroup_stat_boxplot2
Daygroup_stat_boxplot20<-Daygroup_stat_boxplot2+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

ggsave(Daygroup_stat_boxplot1,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot2,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_D7_compare.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare2.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot20,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_D7_compare2.pdf",width = 25, height =22)

##绘制康复曲线
data_anlysis1<-data_anlysis0[which(data_anlysis0$infect_state == "Infect"),]
range(na.omit(data_anlysis1$LMP_infect))##-2 110
##输出相关性统计结果
correlation_results <- data_anlysis1 %>%  group_by(Cell) %>% summarize(correlation = cor(LMP_infect, Cell_ratio,method = "spearman"), p_value = cor.test(LMP_infect, Cell_ratio,method = "spearman")$p.value)
write.table(as.data.frame(correlation_results), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/Spearman_pvalue_each_cell_decidua_data.txt",sep = "\t",row.names=F) 

index_plot<-ggplot(data_anlysis1, aes(x = LMP_infect, y = Cell_ratio)) + #color = group,
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "Decidua_cell_ratio")+
  geom_rug(position = "jitter", size = 0.1, color = "black")+
  scale_x_continuous(breaks = seq(-5,120,10)) + facet_wrap(~ Cell, scales = "free",ncol =4)+
  scale_color_manual(values=my_morandi_colors)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot10<-index_plot+  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)+ geom_point(alpha = 0.5,aes(color = age_group,size =2))
index_plot30<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 0.5,aes(color = age_group,size =2)) +  ggpubr::stat_cor(label.x = 3,method = "spearman")
index_plot31<-index_plot+ geom_smooth(aes(color = age_group,fill=age_group), method = lm,se = TRUE, linewidth=2, alpha = 0.2, fullrange = TRUE)+ geom_point(alpha = 0.5,aes(color = age_group),size = 2) +   ggpubr::stat_cor(aes(color = age_group), label.x = 0,method = "spearman")

index_plot32<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+geom_smooth(method = lm, se = TRUE, fill = "grey",size=2, alpha = 0.5, fullrange = TRUE) +  ggpubr::stat_cor(label.x = 0,method = "spearman")
index_plot33<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+ geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)

#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
ggsave(index_plot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_sm_trend_plot.pdf",width = 16, height =12)
ggsave(index_plot30,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_lm_single_plot.pdf",width = 16, height =12)
ggsave(index_plot31,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_lm_age_split_plot.pdf",width = 16, height =12)
ggsave(index_plot32,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_lm_no_age_split_plot.pdf",width = 16, height =12)
ggsave(index_plot33,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_NEW_sc_no_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_sm_no_age_split_plot.pdf",width = 16, height =12)

##########################################################
##for Decidua sample:manu_1K_adjust
Decidua_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_manu_1K_adjust.txt",header = T,sep = "\t")
Decidua_cell_ratio2<-Decidua_cell_ratio[,-((ncol(Decidua_cell_ratio)-2):ncol(Decidua_cell_ratio))]
rowSums(Decidua_cell_ratio2[,-1])
analysis_data <-melt(Decidua_cell_ratio2,id=c("Mixture"),variable.name="Cell",value.name="ratio")
head(analysis_data)
colnames(analysis_data)<-c("sample","Cell","ratio")
#analysis_data$sample <- factor(analysis_data$sample,levels=paste0("T"))

Ratio_0<-ggplot(analysis_data,aes(sample,ratio,fill=Cell))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+theme_bw()+
  theme(panel.grid=element_blank(),axis.ticks.length=unit(0.5,'cm'),legend.position="right",axis.text.x = element_text(angle=90,hjust=1, vjust=0.5,size=10), axis.text.y = element_text(size=10))+
  scale_fill_manual(values=my_morandi_colors)+guides(fill=guide_legend(title=NULL))
Ratio_0
ggsave(Ratio_0,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_barplot.pdf",width = 16, height =8)

##合并康复日期信息
colData_used <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/Final_fielt_Decidua_Vill_analysis_metadata_genderadd.txt",sep="\t",header =T)
rownames(colData_used)<-as.character(colData_used$sample_code)
dim(colData_used)#245  13
head(colData_used)

colData_Decidua<-colData_used[which(colData_used$sample_code %in% as.character(Decidua_cell_ratio2$Mixture)),]
head(colData_Decidua);dim(colData_Decidua)#125  13
###############################
colData_Decidua$group_day2<-ifelse(colData_Decidua$infect_state =="No","CTRL",
                                   ifelse(colData_Decidua$LMP_infect <= 7,"D7",
                                          ifelse(colData_Decidua$LMP_infect <= 21,"D21",
                                                 ifelse(colData_Decidua$LMP_infect <= 35,"D35",
                                                        ifelse(colData_Decidua$LMP_infect <= 49,"D49",
                                                               ifelse(colData_Decidua$LMP_infect <= 63,"D63",
                                                                      ifelse(colData_Decidua$LMP_infect <= 77,"D77","Dmore77")))))))
table(colData_Decidua$group_day2)
class_order<-c("CTRL","D7","D21","D35","D49","D63","D77","Dmore77")
colData_Decidua$group_day2<-factor(colData_Decidua$group_day2,levels = class_order)

Decidua_cell_ratio2$sample_code <-Decidua_cell_ratio2$Mixture
analysis_Decidua<-merge(colData_Decidua,Decidua_cell_ratio2,by="sample_code")
analysis_Decidua[1:5,1:16]
analysis_Decidua$infect_state <- factor(analysis_Decidua$infect_state,levels=c("Infect","No"))
analysis_Decidua$LMP_infect<-as.numeric(analysis_Decidua$LMP_infect)

table(analysis_Decidua$Year_month)
head(analysis_Decidua)
dim(distinct(analysis_Decidua[,c("sample","infect_state")]))# 125  2

##感染与未感染的组间比较
#因子在感染与非感染组的组间比较
analysis_collect0<-analysis_Decidua[,c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2",colnames(Decidua_cell_ratio2[,-c(1,ncol(Decidua_cell_ratio2))]))]
head(analysis_collect0)
data_anlysis0 <- melt(analysis_collect0,variable.name="Cell",value.name = "Cell_ratio",id.vars = c("sample_code","Age","LMP_infect","LMP_operate","infect_state","gender","group_day2"))
data_anlysis0$age_group<-ifelse(data_anlysis0$Age>=35,"AMA","YMA")
head(data_anlysis0)
#data_anlysis0<-data_anlysis0[which(data_anlysis0$Cell_ratio>0),]

stat_data<-compare_means(Cell_ratio~infect_state, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]

stat_boxplot<-ggplot(data_anlysis0, aes(x=infect_state,y=Cell_ratio,fill=infect_state))+
  geom_boxplot(position=position_dodge(),width=0.5)+geom_jitter(width = 0.2,size=1,color="grey",alpha=0.8)+ 
  scale_color_manual(values=ppCor)+ scale_fill_manual(values=c("blue","red"))+
  stat_compare_means(label="p.signif", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_compare_means(label="p.format", method="wilcox.test",hide.ns = TRUE,label.x=1.5)+
  stat_summary(aes(group= infect_state), fun = "mean", geom = "point",shape=23,size=3,fill="white",position=position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 9))+ facet_wrap(~ Cell, scales = "free",ncol =4)
stat_boxplot
#ggsave(stat_boxplot,file=paste0("/mnt/data/chenwei/covid19/bulk_RNA_result/1.rm_long/Decidua_result/Decidua_compare_stat_box_infect_noremove_",index_name,".pdf"),width = 16, height =8)
ggsave(stat_boxplot,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_compare_stat_box_infect_state.pdf",width = 16, height =12)

################# 比较不同区间分组
head(data_anlysis0)
stat_data<-compare_means(Cell_ratio~group_day2, data=data_anlysis0,group.by = "Cell")
stat_data[which(stat_data$p.signif != "ns"),]
# A tibble:  10 x 9
#   Cell  .y.        group1 group2       p p.adj p.format p.signif method  
#1 FBs_1 Cell_ratio CTRL   D21     0.0104      1 0.0104   *        Wilcoxon
#2 FBs_1 Cell_ratio D21    D49     0.0275      1 0.0275   *        Wilcoxon
#3 VECs  Cell_ratio D21    D35     0.0272      1 0.0272   *        Wilcoxon
#4 LECs  Cell_ratio D7     D49     0.0330      1 0.0330   *        Wilcoxon
#5 LECs  Cell_ratio D7     Dmore77 0.0328      1 0.0328   *        Wilcoxon
#6 DCs   Cell_ratio D7     Dmore77 0.0284      1 0.0284   *        Wilcoxon
#7 DCs   Cell_ratio D21    Dmore77 0.00957     1 0.0096   **       Wilcoxon
#8 DCs   Cell_ratio D63    Dmore77 0.0357      1 0.0357   *        Wilcoxon
#9 NKs   Cell_ratio D7     D49     0.0472      1 0.0472   *        Wilcoxon
#10 Masts Cell_ratio D7     D35     0.0369      1 0.0369   *        Wilcoxon

head(data_anlysis0)
Daygroup_stat_boxplot1 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "CTRL",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =4)
Daygroup_stat_boxplot1
Daygroup_stat_boxplot10<-Daygroup_stat_boxplot1+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

Daygroup_stat_boxplot2 <-ggboxplot(data_anlysis0, x="group_day2", y="Cell_ratio", color="group_day2",add = "jitter") +  
  scale_color_manual(values=ppCor)+ #, palette = "jco"
  # stat_compare_means(method="anova", label.y=40) + 
  stat_compare_means(label="p.signif", method="wilcox.test",ref.group = "D7",hide.ns = TRUE)+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black") +
  theme(axis.text.x = element_text(size = 15,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))+
  facet_wrap(~ Cell, scales = "free",ncol =4)
Daygroup_stat_boxplot2
Daygroup_stat_boxplot20<-Daygroup_stat_boxplot2+ylim(0,max(data_anlysis0$Cell_ratio)+0.1)

ggsave(Daygroup_stat_boxplot1,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot2,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_D7_compare.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_CTRL_compare2.pdf",width = 25, height =22)
ggsave(Daygroup_stat_boxplot20,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_box_static_compare_group_day_D7_compare2.pdf",width = 25, height =22)

##绘制康复曲线
data_anlysis1<-data_anlysis0[which(data_anlysis0$infect_state == "Infect"),]
range(na.omit(data_anlysis1$LMP_infect))##-2 110
##输出相关性统计结果
correlation_results <- data_anlysis1 %>%  group_by(Cell) %>% summarize(correlation = cor(LMP_infect, Cell_ratio,method = "spearman"), p_value = cor.test(LMP_infect, Cell_ratio,method = "spearman")$p.value)
write.table(as.data.frame(correlation_results), file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/Spearman_pvalue_each_cell_decidua_data.txt",sep = "\t",row.names=F) 

index_plot<-ggplot(data_anlysis1, aes(x = LMP_infect, y = Cell_ratio)) + #color = group,
  xlab("infect day to last menstrual period")+ylab("Cell ratio(%)") +labs(title = "Decidua_cell_ratio")+
  geom_rug(position = "jitter", size = 0.1, color = "black")+
  scale_x_continuous(breaks = seq(-5,120,10)) + facet_wrap(~ Cell, scales = "free",ncol =4)+
  scale_color_manual(values=my_morandi_colors)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot10<-index_plot+  geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)+ geom_point(alpha = 0.5,aes(color = age_group,size =2))
index_plot30<-index_plot+  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  geom_point(alpha = 0.5,aes(color = age_group,size =2)) +  ggpubr::stat_cor(label.x = 3,method = "spearman")
index_plot31<-index_plot+ geom_smooth(aes(color = age_group,fill=age_group), method = lm,se = TRUE, linewidth=2, alpha = 0.2, fullrange = TRUE)+ geom_point(alpha = 0.5,aes(color = age_group),size = 2) +   ggpubr::stat_cor(aes(color = age_group), label.x = 0,method = "spearman")

index_plot32<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+geom_smooth(method = lm, se = TRUE, fill = "grey",size=2, alpha = 0.5, fullrange = TRUE) +  ggpubr::stat_cor(label.x = 0,method = "spearman")
index_plot33<-index_plot+  geom_point(alpha = 0.5,colour ="pink",size =2,position=position_dodge(0.1))+ geom_smooth(method = "auto", se = TRUE, fill = "grey",linewidth=2, alpha = 0.5)

#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
ggsave(index_plot10,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_sm_trend_plot.pdf",width = 16, height =12)
ggsave(index_plot30,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_lm_single_plot.pdf",width = 16, height =12)
ggsave(index_plot31,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_lm_age_split_plot.pdf",width = 16, height =12)
ggsave(index_plot32,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_lm_no_age_split_plot.pdf",width = 16, height =12)
ggsave(index_plot33,file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result1203/Bulk_result/Decidua_cell_ratio_manu_sc_1K_adjust/CIBERSORTx_Decidua_bulk_cell_ratio_sm_no_age_split_plot.pdf",width = 16, height =12)
