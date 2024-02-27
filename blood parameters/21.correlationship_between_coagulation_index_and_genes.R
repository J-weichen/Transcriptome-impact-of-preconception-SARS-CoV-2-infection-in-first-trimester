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
library(VennDiagram)

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
analysis_index_wide <- read.table(file= "D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/analysis_index_wide.txt",sep="\t",header = T)
analysis_index_wide$ID_full<-str_pad(analysis_index_wide$ID_full,width =8 ,side = c("left"),pad = "0")
head(analysis_index_wide);dim(analysis_index_wide)# 132  48
###选择凝血相关指标
Coagulation_factor<-c( "PT","APTT","APTT_R","A","INR","Fib","TT","R")
analysis_index_wide<-analysis_index_wide[,c("ID_full", "Vill_code", "Decidua_code",Coagulation_factor)]
head(analysis_index_wide)

###读取最终用于计算的转录组信息
RNA_meta <- read.table(file="D:/PROJECT/新冠/manuscript/3.Table/bulk_RNAseq/Final_fielt_Decidua_Vill_analysis_metadata.txt",sep="\t",header = T)
head(RNA_meta);dim(RNA_meta)#245   8

code_data <-melt(analysis_index_wide[,1:3],id=c("ID_full"),variable.name="type",value.name="sample_code")
code_data2 <-na.omit(code_data[,c(1,3)])
RNA_meta2<-merge(code_data2,RNA_meta,by="sample_code")
RNA_meta3<-distinct(RNA_meta2[,c(2,4:6)])
head(RNA_meta3);dim(RNA_meta3)#132   4

analysis_index_wide2<-merge(analysis_index_wide,RNA_meta3,by="ID_full")
head(analysis_index_wide2);dim(analysis_index_wide2)#132   51

#####输入表达矩阵
count_table_all<-read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/转录组结果/1.rm_long/Decidual_Vill_final_245_sample_gene_expression_matrix_log2_normalized_counts.txt",sep="\t", row.names=1,header =T)
count_table_all[1:4,1:4]

#####输入凝血因子编码基因
##读取目标基因
clinical_index<-as.data.frame(read_excel("D:/PROJECT/新冠/manuscript/3.Table/Coagulation/1.Blood_coagulation_gene.xlsx", sheet =1, col_names = T, col_types = NULL, na = "", skip = 0))
head(clinical_index)
table(clinical_index$Source)
# GO_term_GO_0007596 HALLMARK_COAGULATION 
#     42                  138 
index_maker0<-unique(clinical_index$Gene_symbols)
index_maker1<-index_maker0[which(index_maker0 %in% colnames(count_table_all))]
length(index_maker0);length(index_maker1)## 159 ##158

##读取分泌蛋白列表
secretome_index1 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_SPOCTOPUS.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))
secretome_index2 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_SignalP.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))
secretome_index3 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_Secreted.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))
secretome_index4 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_Phobius.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))
secretome_index5 <- as.data.frame(read.table("D:/PROJECT/新冠/manuscript/3.Table/secretome/protein_class_Predicted.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE))

index_maker0<-unique(c(secretome_index1$Gene,secretome_index2$Gene,secretome_index3$Gene,secretome_index4$Gene,secretome_index5$Gene))
index_maker2<-index_maker0[which(index_maker0 %in% colnames(count_table_all))]
length(index_maker0);length(index_maker2)## 4741 ##4143

index_maker3<-Reduce(intersect, list(index_maker1,index_maker2))
length(index_maker3)#107

######for all coagulation_gene detected############
#######提取Vill 和decidua 各自的表达矩阵###########
count_table_all[1:4,1:4]

Decidua_count <-count_table_all[as.character(na.omit(analysis_index_wide2$Decidua_code)),index_maker1]
Decidua_count<-Decidua_count[,which(colSums(Decidua_count)>0)]
dim(Decidua_count)#125 158

Vill_count <-count_table_all[as.character(na.omit(analysis_index_wide2$Vill_code)),index_maker1]
Vill_count<-Vill_count[,which(colSums(Vill_count)>0)]
dim(Vill_count)#120 158
Vill_count[1:4,1:4]

###############合并凝血相关基因与外周凝血指标数据#############
Vill_index<-analysis_index_wide2[which(!is.na(analysis_index_wide2$Vill_code)),-c(1,3)]
decidua_index<-analysis_index_wide2[which(!is.na(analysis_index_wide2$Decidua_code)),-c(1,2)]
rownames(Vill_index)<-Vill_index$Vill_code;Vill_index<-Vill_index[,-1]
rownames(decidua_index)<-decidua_index$Decidua_code;decidua_index<-decidua_index[,-1]
head(decidua_index)
Vill_marix<-merge(Vill_index,Vill_count,by=0)
Decidua_marix<-merge(decidua_index,Decidua_count,by=0)
dim(Decidua_marix)

###############相关性分析####################
matrix_list<-list(Vill_marix,Decidua_marix)
matrix_name<-c("Vill","Decidua")
index_order<-colnames(Vill_index)

for (list_num in 1:2){
  #list_num<-1
  marix_data<-matrix_list[[list_num]]
  name_matrx<-matrix_name[list_num]
  print(name_matrx)
  
  stat_data<-data.frame()
  for (gene_select in index_maker1){
    #gene_select ="A2M"
   # print(gene_select)
    index_name<-c();cor_r<-c();cor_pvalue <- c()
    for ( index_select in index_order){
      # index_select<-"PT"
      #print(index_select)
      data_anlysis<-na.omit(marix_data[,c(gene_select,index_select)])
      colnames(data_anlysis)<-c("gene","index")
      
     # if (sd(data_anlysis$index)==0) {next}
     # if (sd(data_anlysis$gene)==0) {next}
      
      c_r   <- cor(as.numeric(data_anlysis$gene),as.numeric(data_anlysis$index),method="spearman")
      p_vlue<- cor.test(as.numeric(data_anlysis$gene),as.numeric(data_anlysis$index),method="spearman")[[3]]
      index_name<-c(index_name,index_select)
      cor_r<-c(cor_r,c_r)
      cor_pvalue<-c(cor_pvalue,p_vlue)
    }
    
    length(cor_r);length(cor_pvalue);length(index_name)
    cor_data_df<-data.frame(index=index_name,spearman_cor=cor_r,pvalue=cor_pvalue)
    cor_data_df$gene<-gene_select
    stat_data<-rbind(stat_data,cor_data_df)
  }
  
  
  #names(cor_data_df) <- c("DMR_region","correlation","pvalue")
  head(stat_data);dim(stat_data)##1738    4
  which(stat_data$pvalue < 0.05)
  
  stat_data$group<-"no_sig"
  #stat_data[which(abs(stat_data$spearman_cor) >= 0.600 & stat_data$pvalue < 0.05),]$group <- "candidate_region"
  stat_data[which(stat_data$pvalue < 0.05),]$group <- "candidate_region"
  table(stat_data$group)
  #candidate_region    no_sig 
  #    205             1533 
  
  head(stat_data)
  stat_data$Adjusted_P_Value <- p.adjust(stat_data$pvalue, method = "BH")
  #dim(stat_data[which(stat_data$Adjusted_P_Value<0.05),])# 0 6
  write.table(stat_data, file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_correlated_pvalue_for_coagulation_index_gene_all.txt"),row.names=T, col.names=T,sep="\t") 
}

###################绘制相关性热图#######################
###for villi
#matrix_name<-c(,"Decidua")
name_matrx<-"Vill"
stat_data<-read.table(file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_correlated_pvalue_for_coagulation_index_gene_all.txt"),sep="\t", row.names=1,header =T)
index_order<-c("PT","APTT","APTT_R","A","INR","Fib","TT","R","Age","LMP_infect","LMP_operate")
head(stat_data)
corr2 <- dcast(stat_data[,c("gene","index","spearman_cor")], gene ~ index , value.var = "spearman_cor")
p.mat2 <- dcast(stat_data[,c("gene","index","pvalue")],  gene ~ index , value.var = "pvalue")
rownames(corr2)<-corr2$gene;corr2<-corr2[,-1]
rownames(p.mat2)<-p.mat2$gene;p.mat2<-p.mat2[,-1]
corr2<-corr2[,index_order]
p.mat2<-p.mat2[,index_order]

corheatmap_oneside<-ggcorrplot(t(corr2),outline.color = "grey",ggtheme = ggplot2::theme_bw,
                               colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman",lab = F)+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))

corheatmap_up<-corheatmap_oneside+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                        axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                        axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up

corr_oneside<-t(corr2);pmat_oneside<-t(p.mat2)
corheatmap_oneside2<-ggcorrplot(corr_oneside,outline.color = "grey",ggtheme = ggplot2::theme_bw,p.mat = pmat_oneside,insig = "blank",sig.level = 0.05,
                                colors = c("#6D9EC1", "white", "#E46726"),lab = F,title = "Corralation:spearman")+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_up2<-corheatmap_oneside2+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                          axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                          axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up2

ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_all.pdf"),corheatmap_up,width=10, height=100,limitsize = FALSE)
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_P005.pdf"),corheatmap_up2,width=10, height=100,limitsize = FALSE)

##plot after removing null object and index with significant value
rw_column_fielt <-pmat_oneside
rw_column_fielt[rw_column_fielt > 0.05] <-NA
rw_column_fielt2<-rw_column_fielt[rowSums(!is.na(rw_column_fielt)) !=0,colSums(!is.na(rw_column_fielt)) !=0]
head(rw_column_fielt2)

corr_oneside2<-corr_oneside[rownames(rw_column_fielt2),colnames(rw_column_fielt2)]
pmat_oneside2<-pmat_oneside[rownames(rw_column_fielt2),colnames(rw_column_fielt2)]
corheatmap_filt_lable <-ggcorrplot(round(corr_oneside2,2),outline.color = "grey", p.mat = pmat_oneside2,insig = "blank",
                                   ggtheme = ggplot2::theme_bw,lab = T,lab_size =2, #method = "circle",
                                   colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+                      
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_filt_nolable <-ggcorrplot(corr_oneside2,outline.color = "grey",ggtheme = ggplot2::theme_bw,lab = F, p.mat = pmat_oneside2,insig = "blank",
                                     colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))

corheatmap_up3<-corheatmap_filt_lable+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                            axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                            axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up4<-corheatmap_filt_nolable+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                              axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                              axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up3
corheatmap_up4
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_lable.pdf"),corheatmap_up3,width=10, height=45,limitsize = FALSE)
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_NOlable.pdf"),corheatmap_up4,width=10, height=45,limitsize = FALSE)

###for only eight factors####################
index_order<-c("PT","APTT","APTT_R","A","INR","Fib","TT","R")
head(stat_data)
stat_data2<-stat_data[which(stat_data$index %in% index_order),]
corr2 <- dcast(stat_data2[,c("gene","index","spearman_cor")], gene ~ index , value.var = "spearman_cor")
p.mat2 <- dcast(stat_data2[,c("gene","index","pvalue")],  gene ~ index , value.var = "pvalue")
rownames(corr2)<-corr2$gene;corr2<-corr2[,-1]
rownames(p.mat2)<-p.mat2$gene;p.mat2<-p.mat2[,-1]
corr2<-corr2[,index_order]
p.mat2<-p.mat2[,index_order]

corheatmap_oneside<-ggcorrplot(t(corr2),outline.color = "grey",ggtheme = ggplot2::theme_bw,
                               colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman",lab = F)+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))

corheatmap_up<-corheatmap_oneside+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                        axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                        axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up

corr_oneside<-t(corr2);pmat_oneside<-t(p.mat2)
corheatmap_oneside2<-ggcorrplot(corr_oneside,outline.color = "grey",ggtheme = ggplot2::theme_bw,p.mat = pmat_oneside,insig = "blank",sig.level = 0.05,
                                colors = c("#6D9EC1", "white", "#E46726"),lab = F,title = "Corralation:spearman")+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_up2<-corheatmap_oneside2+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                          axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                          axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up2

ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_all_eight_index.pdf"),corheatmap_up,width=10, height=100,limitsize = FALSE)
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_P005_eight_index.pdf"),corheatmap_up2,width=10, height=100,limitsize = FALSE)

##plot after removing null object and index with significant value
rw_column_fielt <-pmat_oneside
rw_column_fielt[rw_column_fielt > 0.05] <-NA
rw_column_fielt2<-rw_column_fielt[rowSums(!is.na(rw_column_fielt)) !=0,colSums(!is.na(rw_column_fielt)) !=0]
head(rw_column_fielt2)

corr_oneside2<-corr_oneside[rownames(rw_column_fielt2),colnames(rw_column_fielt2)]
pmat_oneside2<-pmat_oneside[rownames(rw_column_fielt2),colnames(rw_column_fielt2)]
corheatmap_filt_lable <-ggcorrplot(round(corr_oneside2,2),outline.color = "grey", p.mat = pmat_oneside2,insig = "blank",
                                   ggtheme = ggplot2::theme_bw,lab = T,lab_size =2, #method = "circle",
                                   colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+                      
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_filt_nolable <-ggcorrplot(corr_oneside2,outline.color = "grey",ggtheme = ggplot2::theme_bw,lab = F, p.mat = pmat_oneside2,insig = "blank",
                                     colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))

corheatmap_up3<-corheatmap_filt_lable+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                            axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                            axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up4<-corheatmap_filt_nolable+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                              axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                              axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up3
corheatmap_up4
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_lable_eight_index.pdf"),corheatmap_up3,width=10, height=12,limitsize = FALSE)
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_NOlable_eight_index.pdf"),corheatmap_up4,width=10, height=12,limitsize = FALSE)


###order by cell number
pmat_oneside3<-as.data.frame(pmat_oneside2)
head(pmat_oneside3)
pmat_oneside3$index<-rownames(pmat_oneside3)
pmat_oneside4 <-reshape2::melt(pmat_oneside3,id="index",variable.name="gene",value.name="value")
head(pmat_oneside4)
pmat_oneside4<-pmat_oneside4[which(pmat_oneside4$value < 0.05),]
pmat_oneside4$index_gene<-paste0(pmat_oneside4$index,"-",pmat_oneside4$gene)
head(pmat_oneside4)


corr_oneside2_num<-as.data.frame(corr_oneside2)
#corr_oneside2_num[corr_oneside2_num > 0] <- NA
corr_oneside3<-as.data.frame(corr_oneside2) #as.data.frame(corr_oneside2_num[rowSums(!is.na(corr_oneside2_num)) !=0,colSums(!is.na(corr_oneside2_num)) !=0])
corr_oneside3[1:4,1:4]
corr_oneside3$index<-as.character(rownames(corr_oneside3))
corr_oneside4 <-reshape2::melt(corr_oneside3,id="index",variable.name="gene",value.name="value")
head(corr_oneside4)

corr_oneside4<-na.omit(corr_oneside4)
corr_oneside4$gene<-as.character(corr_oneside4$gene)
corr_oneside4$index_gene<-paste0(corr_oneside4$index,"-",corr_oneside4$gene)

pmat_oneside4<-pmat_oneside4[which(pmat_oneside4$index_gene %in% unique(corr_oneside4$index_gene)),]
corr_oneside4<-corr_oneside4[which(corr_oneside4$index_gene %in% unique(pmat_oneside4$index_gene)),]
dim(pmat_oneside4);dim(corr_oneside4)#88  4 88  4

corr_oneside5 <- dcast(corr_oneside4, index ~ gene , value.var = "value")
pmat_oneside5 <- dcast(pmat_oneside4, index ~ gene , value.var = "value")
rownames(corr_oneside5)<-corr_oneside5$index;corr_oneside5<-corr_oneside5[,-1]
rownames(pmat_oneside5)<-pmat_oneside5$index;pmat_oneside5<-pmat_oneside5[,-1]
dim(pmat_oneside5);dim(corr_oneside5)
head(corr_oneside5)
###order by cell number
data_cellnumber<-as.data.frame(colSums(!is.na(corr_oneside5)))
data_cellnumber$gene<-rownames(data_cellnumber)
colnames(data_cellnumber)<-c("Cell_number","gene_names")
data_cellnumber<-data_cellnumber[order(data_cellnumber$Cell_number,decreasing = T),]
data_cellnumber$gene_names<-factor(data_cellnumber$gene_names,levels=data_cellnumber$gene_names)

cell_remain<-index_order[which(index_order %in% rownames(corr_oneside5))]
corr_oneside5<-corr_oneside5[cell_remain,rev(rownames(data_cellnumber))]
head(corr_oneside5)

corheatmap_filt_lable <-ggcorrplot(round(corr_oneside5,2),outline.color = "grey", p.mat = pmat_oneside5,insig = "blank",
                                   ggtheme = ggplot2::theme_bw,lab = T,lab_size =3, #method = "circle",
                                   colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+                      
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_up3<-corheatmap_filt_lable+theme(axis.title.y = element_text(size =10, colour = "black"),
                                            axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                            axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up3
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_lable_eight_index_order.pdf"),corheatmap_up3,width=10, height=12,limitsize = FALSE)

########################
###for Decidua
name_matrx<-"Decidua"
stat_data<-read.table(file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_correlated_pvalue_for_coagulation_index_gene_all.txt"),sep="\t", row.names=1,header =T)
index_order<-c("PT","APTT","APTT_R","A","INR","Fib","TT","R","Age","LMP_infect","LMP_operate")
head(stat_data)
corr2 <- dcast(stat_data[,c("gene","index","spearman_cor")], gene ~ index , value.var = "spearman_cor")
p.mat2 <- dcast(stat_data[,c("gene","index","pvalue")],  gene ~ index , value.var = "pvalue")
rownames(corr2)<-corr2$gene;corr2<-corr2[,-1]
rownames(p.mat2)<-p.mat2$gene;p.mat2<-p.mat2[,-1]
corr2<-corr2[,index_order]
p.mat2<-p.mat2[,index_order]

corheatmap_oneside<-ggcorrplot(t(corr2),outline.color = "grey",ggtheme = ggplot2::theme_bw,
                               colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman",lab = F)+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))

corheatmap_up<-corheatmap_oneside+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                        axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                        axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up

corr_oneside<-t(corr2);pmat_oneside<-t(p.mat2)
corheatmap_oneside2<-ggcorrplot(corr_oneside,outline.color = "grey",ggtheme = ggplot2::theme_bw,p.mat = pmat_oneside,insig = "blank",sig.level = 0.05,
                                colors = c("#6D9EC1", "white", "#E46726"),lab = F,title = "Corralation:spearman")+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_up2<-corheatmap_oneside2+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                          axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                          axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up2

ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_all.pdf"),corheatmap_up,width=10, height=100,limitsize = FALSE)
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_P005.pdf"),corheatmap_up2,width=10, height=100,limitsize = FALSE)

##plot after removing null object and index with significant value
rw_column_fielt <-pmat_oneside
rw_column_fielt[rw_column_fielt > 0.05] <-NA
rw_column_fielt2<-rw_column_fielt[rowSums(!is.na(rw_column_fielt)) !=0,colSums(!is.na(rw_column_fielt)) !=0]
head(rw_column_fielt2)

corr_oneside2<-corr_oneside[rownames(rw_column_fielt2),colnames(rw_column_fielt2)]
pmat_oneside2<-pmat_oneside[rownames(rw_column_fielt2),colnames(rw_column_fielt2)]
corheatmap_filt_lable <-ggcorrplot(round(corr_oneside2,2),outline.color = "grey", p.mat = pmat_oneside2,insig = "blank",
                                   ggtheme = ggplot2::theme_bw,lab = T,lab_size =2, #method = "circle",
                                   colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+                      
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_filt_nolable <-ggcorrplot(corr_oneside2,outline.color = "grey",ggtheme = ggplot2::theme_bw,lab = F, p.mat = pmat_oneside2,insig = "blank",
                                     colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))

corheatmap_up3<-corheatmap_filt_lable+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                            axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                            axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up4<-corheatmap_filt_nolable+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                              axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                              axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up3
corheatmap_up4
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_lable.pdf"),corheatmap_up3,width=10, height=45,limitsize = FALSE)
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_NOlable.pdf"),corheatmap_up4,width=10, height=45,limitsize = FALSE)

###for only eight factors####################
index_order<-c("PT","APTT","APTT_R","A","INR","Fib","TT","R")
head(stat_data)
stat_data2<-stat_data[which(stat_data$index %in% index_order),]
corr2 <- dcast(stat_data2[,c("gene","index","spearman_cor")], gene ~ index , value.var = "spearman_cor")
p.mat2 <- dcast(stat_data2[,c("gene","index","pvalue")],  gene ~ index , value.var = "pvalue")
rownames(corr2)<-corr2$gene;corr2<-corr2[,-1]
rownames(p.mat2)<-p.mat2$gene;p.mat2<-p.mat2[,-1]
corr2<-corr2[,index_order]
p.mat2<-p.mat2[,index_order]

corheatmap_oneside<-ggcorrplot(t(corr2),outline.color = "grey",ggtheme = ggplot2::theme_bw,
                               colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman",lab = F)+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))

corheatmap_up<-corheatmap_oneside+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                        axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                        axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up

corr_oneside<-t(corr2);pmat_oneside<-t(p.mat2)
corheatmap_oneside2<-ggcorrplot(corr_oneside,outline.color = "grey",ggtheme = ggplot2::theme_bw,p.mat = pmat_oneside,insig = "blank",sig.level = 0.05,
                                colors = c("#6D9EC1", "white", "#E46726"),lab = F,title = "Corralation:spearman")+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_up2<-corheatmap_oneside2+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                          axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                          axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up2

ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_all_eight_index.pdf"),corheatmap_up,width=10, height=100,limitsize = FALSE)
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_P005_eight_index.pdf"),corheatmap_up2,width=10, height=100,limitsize = FALSE)

##plot after removing null object and index with significant value
rw_column_fielt <-pmat_oneside
rw_column_fielt[rw_column_fielt > 0.05] <-NA
rw_column_fielt2<-rw_column_fielt[rowSums(!is.na(rw_column_fielt)) !=0,colSums(!is.na(rw_column_fielt)) !=0]
head(rw_column_fielt2)

corr_oneside2<-corr_oneside[rownames(rw_column_fielt2),colnames(rw_column_fielt2)]
pmat_oneside2<-pmat_oneside[rownames(rw_column_fielt2),colnames(rw_column_fielt2)]
corheatmap_filt_lable <-ggcorrplot(round(corr_oneside2,2),outline.color = "grey", p.mat = pmat_oneside2,insig = "blank",
                                   ggtheme = ggplot2::theme_bw,lab = T,lab_size =2, #method = "circle",
                                   colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+                      
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_filt_nolable <-ggcorrplot(corr_oneside2,outline.color = "grey",ggtheme = ggplot2::theme_bw,lab = F, p.mat = pmat_oneside2,insig = "blank",
                                     colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))

corheatmap_up3<-corheatmap_filt_lable+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                            axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                            axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up4<-corheatmap_filt_nolable+theme(axis.title.y = element_text(size = 2, colour = "black"),
                                              axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                              axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up3
corheatmap_up4
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_lable_eight_index.pdf"),corheatmap_up3,width=10, height=12,limitsize = FALSE)
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_NOlable_eight_index.pdf"),corheatmap_up4,width=10, height=12,limitsize = FALSE)


###order by cell number
pmat_oneside3<-as.data.frame(pmat_oneside2)
head(pmat_oneside3)
pmat_oneside3$index<-rownames(pmat_oneside3)
pmat_oneside4 <-reshape2::melt(pmat_oneside3,id="index",variable.name="gene",value.name="value")
head(pmat_oneside4)
pmat_oneside4<-pmat_oneside4[which(pmat_oneside4$value < 0.05),]
pmat_oneside4$index_gene<-paste0(pmat_oneside4$index,"-",pmat_oneside4$gene)
head(pmat_oneside4)


corr_oneside3<-as.data.frame(corr_oneside2) #as.data.frame(corr_oneside2_num[rowSums(!is.na(corr_oneside2_num)) !=0,colSums(!is.na(corr_oneside2_num)) !=0])
corr_oneside3[1:4,1:4]
corr_oneside3$index<-as.character(rownames(corr_oneside3))
corr_oneside4 <-reshape2::melt(corr_oneside3,id="index",variable.name="gene",value.name="value")
head(corr_oneside4)

corr_oneside4<-na.omit(corr_oneside4)
corr_oneside4$gene<-as.character(corr_oneside4$gene)
corr_oneside4$index_gene<-paste0(corr_oneside4$index,"-",corr_oneside4$gene)

pmat_oneside4<-pmat_oneside4[which(pmat_oneside4$index_gene %in% unique(corr_oneside4$index_gene)),]
corr_oneside4<-corr_oneside4[which(corr_oneside4$index_gene %in% unique(pmat_oneside4$index_gene)),]
dim(pmat_oneside4);dim(corr_oneside4)#88  4 88  4

corr_oneside5 <- dcast(corr_oneside4, index ~ gene , value.var = "value")
pmat_oneside5 <- dcast(pmat_oneside4, index ~ gene , value.var = "value")
rownames(corr_oneside5)<-corr_oneside5$index;corr_oneside5<-corr_oneside5[,-1]
rownames(pmat_oneside5)<-pmat_oneside5$index;pmat_oneside5<-pmat_oneside5[,-1]
dim(pmat_oneside5);dim(corr_oneside5)
head(corr_oneside5)
###order by cell number
data_cellnumber<-as.data.frame(colSums(!is.na(corr_oneside5)))
data_cellnumber$gene<-rownames(data_cellnumber)
colnames(data_cellnumber)<-c("Cell_number","gene_names")
data_cellnumber<-data_cellnumber[order(data_cellnumber$Cell_number,decreasing = T),]
data_cellnumber$gene_names<-factor(data_cellnumber$gene_names,levels=data_cellnumber$gene_names)

cell_remain<-index_order[which(index_order %in% rownames(corr_oneside5))]
corr_oneside5<-corr_oneside5[cell_remain,rev(rownames(data_cellnumber))]
head(corr_oneside5)

corheatmap_filt_lable <-ggcorrplot(round(corr_oneside5,2),outline.color = "grey", p.mat = pmat_oneside5,insig = "blank",
                                   ggtheme = ggplot2::theme_bw,lab = T,lab_size =3, #method = "circle",
                                   colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:spearman")+                      
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0,limits=c(-1,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
corheatmap_up3<-corheatmap_filt_lable+theme(axis.title.y = element_text(size =10, colour = "black"),
                                            axis.text.x = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 90),
                                            axis.text.y = element_text(size = 7,colour = 'black',vjust=0.5,hjust=1,angle = 0))
corheatmap_up3
ggsave(paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx,"_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005_lable_eight_index_order.pdf"),corheatmap_up3,width=7, height=12,limitsize = FALSE)

###############组织间相关性讨论############
index_order<-c("PT","APTT","APTT_R","A","INR","Fib","TT","R")
name_matrx1<-"Decidua"
Decidua_stat_data<-read.table(file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx1,"_correlated_pvalue_for_coagulation_index_gene_all.txt"),sep="\t", row.names=1,header =T)
Decidua_data<-Decidua_stat_data[which(Decidua_stat_data$index %in% index_order),]

name_matrx2<-"Vill"
Vill_stat_data<-read.table(file=paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/",name_matrx2,"_correlated_pvalue_for_coagulation_index_gene_all.txt"),sep="\t", row.names=1,header =T)
Villi_data<-Vill_stat_data[which(Vill_stat_data$index %in% index_order),]

head(Decidua_data)
Decidua_data2<-Decidua_data[which(Decidua_data$pvalue < 0.05),]
Villi_data2<-Villi_data[which(Villi_data$pvalue < 0.05),]
Decidua_data2$trend<-ifelse(Decidua_data2$spearman_cor<0,"Down","Up")
Villi_data2$trend<-ifelse(Villi_data2$spearman_cor<0,"Down","Up")

Decidua_data2$term<-paste(Decidua_data2$index,Decidua_data2$gene,Decidua_data2$trend,sep = ":")
Villi_data2$term<-paste(Villi_data2$index,Villi_data2$gene,Villi_data2$trend,sep = ":")

intersect(unique(Decidua_data2$gene),unique(Villi_data2$gene))## "CASP9" "F9" "FURIN" "HTRA1" "MEP1A" "F13A1"
length(intersect(unique(Decidua_data2$term),unique(Villi_data2$term))) #0

unique(Decidua_data2$gene)
# [1] "C1R"      "C3"       "CAPN2"    "CASP9"    "CFH"      "CSRP1"    "F10"      "F9"       "FN1"      "FURIN"    "HTRA1"    "ITIH1"    "LRP1"     "MAFF"    
#[15] "MEP1A"    "MMP1"     "MMP8"     "MMP9"     "RAPGEF3"  "S100A1"   "SERPINA1" "SH2B2"    "SPARC"    "THBS1"    "USP11"    "VWF"      "ADORA2A"  "C4BPB"   
#[29] "F13A1"    "F2RL3"    "MMRN1"    "PABPC4"  

unique(Villi_data2$gene)
#[1] "ANXA1"    "C1QA"     "C2"       "C8A"      "C8G"      "C9"       "CASP9"    "CFD"      "CLU"      "CPN1"     "CTSH"     "CTSK"     "DCT"      "DPP4"    
#[15] "DUSP14"   "DUSP6"    "F11"      "F13B"     "F2"       "F2RL2"    "F9"       "FURIN"    "GP9"      "HRG"      "HTRA1"    "KLKB1"    "LAMP2"    "MBL2"    
#[29] "MEP1A"    "MMP10"    "MMP11"    "PLAT"     "PLEK"     "PROC"     "RABIF"    "RAC1"     "SERPINE1" "SIRT2"    "CD36"     "CD40LG"   "F13A1"    "GGCX"    

####relationship among DEGs and correlated coagulation genes
##for decidual genes 
tag<-"single_group_compare_vs_CTRL"
file_used<-paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/Final_",tag,"_DEGs_p001_FC2_number_merge_decidua.txt")
Decidual_gene<-read.table(file=file_used,sep="\t", row.names=1,header =T)
Decidual_gene[1:4,1:4]
Decidual_gene$day<-rownames(Decidual_gene)
Decidual_gene2 <- melt(Decidual_gene,variable.name="gene",value.name = "count",id.vars = c("day"))
Decidual_Up<-as.character(unique(Decidual_gene2[which(Decidual_gene2$count==1),]$gene))
Decidual_Down<-as.character(unique(Decidual_gene2[which(Decidual_gene2$count==-1),]$gene))

file_used<-paste0("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/Final_",tag,"_DEGs_p001_FC2_number_merge_villi.txt")
Villi_gene<-read.table(file=file_used,sep="\t", row.names=1,header =T)
Villi_gene[1:4,1:4]
Villi_gene$day<-rownames(Villi_gene)
Villi_gene2 <- melt(Villi_gene,variable.name="gene",value.name = "count",id.vars = c("day"))
Villi_Up<-as.character(unique(Villi_gene2[which(Villi_gene2$count==1),]$gene))
Villi_Down<-as.character(unique(Villi_gene2[which(Villi_gene2$count==-1),]$gene))
##relationship
length(Decidual_Up);length(Decidual_Down);length(Villi_Up);length(Villi_Down)
#700 156 1468 286

Decidua_coagulation <-unique(Decidua_data2$gene)
Villi_coagulation <-unique(Villi_data2$gene)

#Villi_Up=Villi_Up,Villi_Down=Villi_Down,

venn <-venn.diagram(list(Decidual_Up=Decidual_Up,Decidual_Down=Decidual_Down,Decidua_coagulation=Decidua_coagulation),
                    alpha=c(0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","blue","grey"), 
                    cex = 1.5,cat.col=c("red","blue","grey"), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage(); 
grid.draw(venn)

pdf("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/Venn_relationship_among_Decidual_DEGs_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005.pdf",width = 6,height = 6)
grid.draw(venn)
dev.off()

venn <-venn.diagram(list(Villi_Up=Villi_Up,Villi_Down=unique(Villi_Down),Villi_coagulation=Villi_coagulation),
                    alpha=c(0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","blue","grey"), 
                    cex = 1.5,cat.col=c("red","blue","grey"), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage(); 
grid.draw(venn)

pdf("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/Venn_relationship_among_vill_DEGs_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005.pdf",width = 6,height = 6)
grid.draw(venn)
dev.off()

unique(Reduce(intersect, list(Decidual_Up,Decidua_coagulation)))#"THBS1"
unique(Reduce(intersect, list(Villi_Up,Villi_coagulation)))     #"CFD"   "PLAT"  "CLU"   "F2"    "KLKB1"
unique(Reduce(intersect, list(Villi_Down,Villi_coagulation)))#"FURIN" "GP9" 

Decidual_gene[,c("THBS1")]#Dmore77

Villi_gene[,c("CFD","PLAT","CLU","F2","KLKB1","FURIN","GP9")]
#         CFD PLAT CLU F2 KLKB1 FURIN GP9
#D7        1    0   0  0     0    -1  -1
#D35       0    1   1  0     0     0   0
#D49       0    0   0  1     0     0   0
#D63       0    0   0  1     1     0   0
#D77       0    0   0  1     0     0   0

################relationship among glm genes and correlated coagulation genes
Vill_results<-read.table(file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/Vill_gene_Lm_model2.txt",sep="\t", row.names=1,header =T)
head(Vill_results)
Vill_results$Gene_Name<-as.character(Vill_results$Gene_Name)
Villus_RD<-Vill_results[which(Vill_results$index == "LMP_infect"),]
Villus_RD$Adjusted_P_Value <- p.adjust(Villus_RD$p_values_type2, method = "BH")
#plot_data_RD$log10_pvalue<-(-log10(plot_data_RD$Adjusted_P_Value))
Villus_RD$threshold <-ifelse(Villus_RD$Adjusted_P_Value <0.05, ifelse(Villus_RD$Coefficient> 0,'Up','Down'),"No")
Villus_RD$threshold <-factor(Villus_RD$threshold,levels = c('Up','Down',"No"))
Villus_RD_Up<-Villus_RD[which(Villus_RD$threshold == "Up"),]
Villus_RD_Down<-Villus_RD[which(Villus_RD$threshold == "Down"),]
dim(Villus_RD_Up);dim(Villus_RD_Down)#1071    7 1555    7

##############reading decidual  recovery days related glm genes 
Decidual_results<-read.table( file="D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/decidua_gene_Lm_model2.txt",sep="\t", row.names=1,header =T)
head(Decidual_results)
Decidual_results$Gene_Name<-as.character(Decidual_results$Gene_Name)
Decidual_RD<-Decidual_results[which(Decidual_results$index == "LMP_infect"),]
Decidual_RD$Adjusted_P_Value <- p.adjust(Decidual_RD$p_values_type2, method = "BH")
#plot_data_RD$log10_pvalue<-(-log10(plot_data_RD$Adjusted_P_Value))
Decidual_RD$threshold <-ifelse(Decidual_RD$Adjusted_P_Value <0.05, ifelse(Decidual_RD$Coefficient> 0,'Up','Down'),"No")
Decidual_RD$threshold <-factor(Decidual_RD$threshold,levels = c('Up','Down',"No"))
Decidual_RD_Up<-Decidual_RD[which(Decidual_RD$threshold == "Up"),]
Decidual_RD_Down<-Decidual_RD[which(Decidual_RD$threshold == "Down"),]
dim(Decidual_RD_Up);dim(Decidual_RD_Down)#211   7 320   7
Decidual_RD[which(Decidual_RD$Gene_Name %in% c("FURIN","VWF","USP11")),]
######################################################
head(Decidual_RD_Up)
Decidual_Up_glm<-unique(Decidual_RD_Up$Gene_Name)
Decidual_down_glm<-unique(Decidual_RD_Down$Gene_Name)
Vill_Up_glm<-unique(Villus_RD_Up$Gene_Name)
Vill_down_glm<-unique(Villus_RD_Down$Gene_Name)

#####################################################
Decidua_coagulation <-unique(Decidua_data2$gene)
Villi_coagulation <-unique(Villi_data2$gene)

#Villi_Up=Villi_Up,Villi_Down=Villi_Down,

venn <-venn.diagram(list(Decidual_Up=Decidual_Up_glm,Decidual_Down=Decidual_down_glm,Decidua_coagulation=Decidua_coagulation),
                    alpha=c(0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","blue","grey"), 
                    cex = 1.5,cat.col=c("red","blue","grey"), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage(); 
grid.draw(venn)

pdf("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/Venn_relationship_among_Decidual_glm_gene_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005.pdf",width = 6,height = 6)
grid.draw(venn)
dev.off()

venn <-venn.diagram(list(Villi_Up=Vill_Up_glm,Villi_Down=Vill_down_glm,Villi_coagulation=Villi_coagulation),
                    alpha=c(0.8,0.8,0.8),lwd=1,lty=1, col="white",fill=c("red","blue","grey"), 
                    cex = 1.5,cat.col=c("red","blue","grey"), cat.fontface=4,cat.cex = 1.5,    
                    main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                    rotation.degree=0,filename = NULL)
grid.newpage(); 
grid.draw(venn)

pdf("D:/PROJECT/新冠/manuscript/2.Figure/clincal_index/figure_plot/RNA_coagulation_relationship/Venn_relationship_among_vill_glm_gene_spearman_correlated_pvalue_for_coagulation_index_gene_pvalue005.pdf",width = 6,height = 6)
grid.draw(venn)
dev.off()

##############    #####
unique(Reduce(intersect,list(Vill_down_glm,Villi_coagulation)))  ### "CFD"
##############    #####
##分泌性质蛋白相关性index_maker3
unique(Reduce(intersect,list(index_maker3,Villi_coagulation)))  ### "CFD"
#[1] "ANXA1"    "C1QA"     "C2"       "C8A"      "C8G"      "C9"       "CFD"      "CLU"      "CPN1"     "CTSH"     "CTSK"     "DPP4"     "DUSP6"    "F11"     
#[15] "F13B"     "F2"       "F9"       "HRG"      "HTRA1"    "KLKB1"    "MBL2"     "MMP10"    "MMP11"    "PLAT"     "PROC"     "SERPINE1" "CD40LG"   "F13A1" 
unique(Reduce(intersect,list(index_maker3,Decidua_coagulation)))  ### "CFD"
#[1] "C1R"      "C3"       "CFH"      "F10"      "F9"       "FN1"      "HTRA1"    "ITIH1"    "LRP1"     "MMP1"     "MMP8"     "MMP9"     "SERPINA1" "SPARC"   
#[15] "THBS1"    "USP11"    "VWF"      "C4BPB"    "F13A1"    "F2RL3"    "MMRN1"    "PABPC4"  