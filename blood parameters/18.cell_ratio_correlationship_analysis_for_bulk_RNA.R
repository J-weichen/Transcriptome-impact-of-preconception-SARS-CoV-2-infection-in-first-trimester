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
#set colors
pal1<-pal_nejm("default",alpha = 1)(8)
pal2<-pal_jama("default",alpha = 1)(7)
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)
pal5 <- pal_npg("nrc", alpha=0.5)(10)
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor_all)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]

##for Decidua sample
True_sc_Decidua_cell_ratio <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/True_sc_Decidua_cell_ratio.txt",header = T,sep = "\t")
CIBERSORTx_pesudo_Decidua <- read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_pesudo_Mixture_Decidua.txt",header = T,sep = "\t")

head(True_sc_Decidua_cell_ratio)
head(CIBERSORTx_pesudo_Decidua)
CB_data <-melt(CIBERSORTx_pesudo_Decidua,id=c("Mixture"),variable.name="Cell",value.name="ratio_CB")
head(CB_data)
True_data<-True_sc_Decidua_cell_ratio[,c(2,1,3)]
colnames(True_data)<-c("Mixture","Cell","ratio_true")
com_data<-merge(True_data,CB_data,by=c("Mixture","Cell"))

##绘制相关性图
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

range(com_data$ratio_true)# 0.0000000 0.7373995
range(com_data$ratio_CB)# 0.0000000 0.7621307
formula <- y ~ x
index_plot<-ggplot(com_data, aes(x = ratio_true, y = ratio_CB)) + #color = group,
  xlab("True SC Decidua Cell ratio")+ylab("CIBERSORTx Decidua Cell ratio") +labs(title = "Decidua correlationship(pearson)")+
  scale_x_continuous(breaks = seq(0,1,0.1)) + 
  scale_color_manual(values=ppCor_all)+
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  #geom_text(aes(x=55,y=124,label="R^2=0.88\ny=0.75x+58.96"),color="red",family = "serif",fontface = "plain",size = 5)+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+##method one
  # stat_cor(aes(label = after_stat(rr.label)), color = "red", geom = "label")+##method two
  #stat_fit_glance(method = 'lm', method.args = list(formula = formula),geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.35, size = 3)+
  geom_point(alpha = 0.5,aes(color = Mixture,size =2)) +  
  ggpubr::stat_cor(method = "pearson",size =4,label.x = 0.1,label.y = 0.8)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,..p.value.label..,sep ="~~~")),size =4,colour = "black", formula = formula, label.x = 0.1,label.y = 0.85, parse = TRUE)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file="D:/高龄-卵巢项目/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/pesudo_Decidua_cell_ratio_correlationship.pdf",width = 8, height =8)

##for Vill sample
True_sc_Vill_cell_ratio <- read.table(file="D:/高龄-卵巢项目/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/True_sc_Vill_cell_ratio.txt",header = T,sep = "\t")
CIBERSORTx_pesudo_Vill <- read.table(file="D:/高龄-卵巢项目/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/CIBERSORTx_pesudo_Mixture_Vill.txt",header = T,sep = "\t")

head(True_sc_Vill_cell_ratio)
head(CIBERSORTx_pesudo_Vill)
CB_data <-melt(CIBERSORTx_pesudo_Vill,id=c("Mixture"),variable.name="Cell",value.name="ratio_CB")
head(CB_data)
True_data<-True_sc_Vill_cell_ratio[,c(2,1,3)]
colnames(True_data)<-c("Mixture","Cell","ratio_true")
com_data<-merge(True_data,CB_data,by=c("Mixture","Cell"))
table(com_data$Mixture)
table(com_data$Cell)

##绘制相关性图
##线性相关性分析1）计算相关性 reference:https://zhuanlan.zhihu.com/p/543768987
#ref:https://www.coder.work/article/6652294#google_vignette
#ref:https://www.saoniuhuo.com/question/detail-2441094.html
#https://rdrr.io/cran/ggpmisc/man/stat_poly_eq.html
cor.test(com_data$ratio_true,com_data$ratio_CB,data=com_data)
###根据结果可以知道相关性系数为0.7687359，且呈正相关，P值为7.84e-15<0.05
###回归系数的95% 置信区间为[0.6514320 0.8501456]
df_cor<-lm(ratio_CB~ratio_true,data=com_data)
summary(df_cor)
###回归方程为ratio_CB=0.927231*ratio_CB+0.005198  
###拟合优度R^2= 7.84e-15,即拟合度一般 
##模型 p 值，而不是斜率 p 值(通常) p-value: < 7.84e-15

range(com_data$ratio_true)# 0.0000000  0.704165
range(com_data$ratio_CB)# 0.0000000 0.5912516
formula <- y ~ x
index_plot<-ggplot(com_data, aes(x = ratio_true, y = ratio_CB)) + #color = group,
  xlab("True SC Vill Cell ratio")+ylab("CIBERSORTx Vill Cell ratio") +labs(title = "Vill correlationship(pearson)")+
  scale_x_continuous(breaks = seq(0,1,0.1)) + 
  scale_color_manual(values=my_morandi_colors)+
  geom_smooth(method = lm, se = TRUE, fill = "grey",size=1, alpha = 0.5, fullrange = TRUE) +  
  #geom_text(aes(x=55,y=124,label="R^2=0.88\ny=0.75x+58.96"),color="red",family = "serif",fontface = "plain",size = 5)+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+##method one
  # stat_cor(aes(label = after_stat(rr.label)), color = "red", geom = "label")+##method two
  #stat_fit_glance(method = 'lm', method.args = list(formula = formula),geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.35, size = 3)+
  geom_point(alpha = 0.5,aes(color = Mixture,size =2)) +  
  ggpubr::stat_cor(method = "pearson",size =4,label.x = 0.1,label.y = 0.8)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label..,..p.value.label..,sep ="~~~")),size =4,colour = "black", formula = formula, label.x = 0.1,label.y = 0.85, parse = TRUE)+
  theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 16), 
                    plot.caption = element_text(size = 12), axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 15)) 
index_plot
ggsave(index_plot,file="D:/高龄-卵巢项目/新冠/manuscript/2.Figure/转录组结果/bulk_RNA_result/4.cellratio_result/pesudo_Vill_cell_ratio_correlationship.pdf",width = 8, height =8)

#########plot the formal data
