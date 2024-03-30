### 绘制各个cluster 差异基因的散点图
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel) #用于注释文本
library(magrittr)
library(dplyr)
rm(list = ls())
setwd("D:/Projects/IFNprojects/analysis/")
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# ClusterMarker_AM=read.csv("./DFgenes/DFgenes_Alpha_vs_Mock_27ct.csv",row.names = 1)
# ClusterMarker_LM=read.csv("./DFgenes/DFgenes_Lambda_vs_Mock_27ct.csv",row.names = 1)
# 
ct=c("SUS","Bowman gland cell")
# df=rbind(ClusterMarker_LM[which(ClusterMarker_LM$Tcelltype %in% ct),],
#          ClusterMarker_AM[which(ClusterMarker_AM$Tcelltype %in% ct),])
# write.csv(df,file = "scatterDEGplot_testData.csv")
df = read.csv("scatterDEGplot_testData.csv",row.names = 1)

ClusterMarker=df[which(df$Tcelltype %in% ct),]

ClusterMarker$cluster=paste0(ClusterMarker$Tcelltype,"_",ClusterMarker$compare)
# 查看每个cluster的marker基因数量
table(ClusterMarker$cluster)

# 根据自己计算的marker基因数量确定log2FC的阈值，这里先定为0.5
ClusterMarker <- subset(ClusterMarker,p_val_adj < 0.05 & abs(avg_logFC) > 0.15)# 
ClusterMarker$thr_signi <- as.factor(ifelse(ClusterMarker$avg_logFC > 0 , 'Up', 'Down'))
dim(ClusterMarker)
table(ClusterMarker$thr_signi)
# Down   Up 
# 3044  890 

##指定类别排列顺序
orderLevels=c(paste0(rep(ct[1],2),c("_Alpha_vs_Mock","_Lambda_vs_Mock")),#两细胞类别
              paste0(rep(ct[2],2),c("_Alpha_vs_Mock","_Lambda_vs_Mock")))
ClusterMarker$cluster <- factor(ClusterMarker$cluster,
                                levels=orderLevels)
##类别变量均需转换为数字变量
ClusterMarker$cluster=plyr::mapvalues(ClusterMarker$cluster,
                                      from = orderLevels,
                                      to = seq(0,length(unique(ClusterMarker$cluster))-1))
ClusterMarker$cluster=ClusterMarker$cluster %>% as.vector(.) %>% as.numeric(.)

## 这里挑选log2FC为top5的基因进行展示
top_up_label <- ClusterMarker %>% 
  subset(., thr_signi%in%"Up" ) %>% #& avg_logFC >1.5
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_logFC) %>%
  as.data.frame()

top_down_label <- ClusterMarker %>% 
  subset(., thr_signi%in%"Down" ) %>% #& avg_logFC < -1.5
  group_by(cluster) %>%
  top_n(n = -5, wt = avg_logFC) %>%
  as.data.frame()

top_label <- rbind(top_up_label,top_down_label)
top_label$thr_signi %<>% 
  factor(., levels = c("Up","Down"))

## 开始绘图 
### 准备绘制暗灰色背景所需数据 
background_position <- ClusterMarker %>%
  dplyr::group_by(cluster) %>%
  summarise(Min = min(avg_logFC) - 0.2, Max = max(avg_logFC) + 0.2) %>%
  as.data.frame()
background_position$Max=ifelse(background_position$Max<0,0,background_position$Max)

##设置背景框的宽度
background_position$start <- background_position$cluster - 0.3
background_position$end <- background_position$cluster + 0.3

### 准备绘制中间区域cluster彩色bar所需数据
cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$cluster - 0.5
cluster_bar_position$end <- cluster_bar_position$cluster + 0.5
cluster_bar_position$cluster %<>% 
  factor(., levels = c(0:max(as.vector(.))))


## 设置填充颜色
tissueCols = c("#0a8cf7","#febe08","#ed2b2b")#new 更鲜艳配色 
cols_thr_signi <- c("Up" = tissueCols[3],
                    "Down" = tissueCols[1])
cols_cluster <- c("0" = tissueCols[2],
                  "1" = tissueCols[3],
                  "2" = tissueCols[2],
                  "3" = tissueCols[3])

## 将需要标注的基因合并至原数据框中
dat=ClusterMarker
dat$tmp1=paste0(dat$gene,"_",dat$Tcelltype,"_",dat$compare)
top_label$tmp2=paste0(top_label$gene,"_",top_label$Tcelltype,"_",top_label$compare)
intersect(dat$tmp1,top_label$tmp2)
dat$gene_label=ifelse(dat$tmp1 %in% top_label$tmp2,dat$gene,"")

library(ggrastr)
p <- ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,
                                            ymax = Max),
            fill = "grey60", alpha = 0.1) + ###添加灰色背景色
  geom_jitter(data = dat, aes(x = cluster, y = avg_logFC, 
                              colour = thr_signi,size=-log10(p_val_adj)),
              # alpha=0.8 ,
              position = position_jitter(width = 0.3,seed = 1)
  )+
  geom_point_rast(raster.dpi=600)+
  geom_text_repel(data = dat, aes(x = cluster, y = avg_logFC, label = gene_label),
                  position = position_jitter(width = 0.3, seed = 1),
                  show.legend = F, size = 4,
                  box.padding = unit(0.2, "lines")) +
  scale_size_continuous(range = c(0.1,3)) +
  scale_color_manual(values = cols_thr_signi)+
  scale_x_continuous(limits = c(-0.5, max(ClusterMarker$cluster) + 0.5),
                     breaks = seq(0, max(ClusterMarker$cluster), 1),
                     labels = NULL) + #修改坐标轴显示刻度
  geom_rect(data = cluster_bar_position, aes(xmin = start+0.1, xmax = end-0.1, ymin = -0.15,
                                             ymax = 0.15, fill = cluster), 
            color = "black", alpha = 1, show.legend = F) +
  scale_fill_manual(values = cols_cluster) +
  geom_vline(xintercept=1.5,lty=2,col="black",lwd=0.6) +
  labs(x = NULL, y = "Log2(FC)",color=NULL,size="-Log10(adjust P value)") +
  theme_bw() 
plot1 <- p + theme(panel.grid.minor = element_blank(), ##去除网格线
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(colour = 'black', size = 12),
                   axis.title = element_text(colour = 'black', size = 14),
                   # axis.text.x = element_text(colour = 'black', size = 12, vjust = 22,hjust=0.5), #调整x轴坐标,vjust的值按照最终结果稍加调整
                   panel.border = element_blank(), ## 去掉坐标轴
                   axis.ticks.x = element_blank(), ## 去掉的坐标刻度线
                   axis.ticks.y = element_line(colour = "black",size = 0.8),
                   axis.ticks.length.y = unit(0.15,"cm"),
                   legend.margin = margin(l=-0.9,unit = "cm"), #缩小图例与图间的距离
                   axis.line.y = element_line(colour = "black",size = 0.8,
                                              arrow = arrow(length = unit(0.2, 'cm'), #画箭头
                                                            ends = "both", #设置双向实心箭头
                                                            type = "closed")))+ #添加y轴坐标轴
  annotate(geom="text", x=c(0.5,2.5), y=-2.7, #调整细胞类别注释的位置
           label=c("SUS","BGC"),color="black",size=5)+
  annotate(geom="text", x=seq(0, max(ClusterMarker$cluster), 1), y=0, 
           label=c(expression("IFN-"*alpha*" vs Mock"),expression("IFN-"*lambda*" vs Mock"),
                   expression("IFN-"*alpha*" vs Mock"),expression("IFN-"*lambda*" vs Mock")),#标注每组的类别名称
           color="black",size=4)+
  guides(color=guide_legend(reverse =T,override.aes = list(size=3)))
plot1

ggsave(filename = paste0("./DFgenes/","SUS-Bowman","_combineVolcanoplot.pdf"),
       plot = plot1, width = 8, height = 5.5)