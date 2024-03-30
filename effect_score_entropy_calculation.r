library(Seurat)
library(ggplot2)
library(dplyr)
library(plyr)
library(Cairo)
library(RColorBrewer)
library(pheatmap)
library(Hmisc)
library(ggpubr)
library(ggrastr)

##---------------- compute heterogeneity (Entropy) of different group in each celltype --------------------
Sdata = FindNeighbors(Sdata, dims = 1:50,return.neighbor = T,k.param = 20)

NN = Sdata@neighbors$RNA.nn@nn.idx
head(NN)
rownames(NN) = Sdata@neighbors$RNA.nn@cell.names

colnames(NN) = paste0('nn_',seq(1,20))
NN = as.data.frame(NN)

tmp = data.frame(row.names = seq(1,27985),cellName = Sdata$orig.ident)
head(tmp)
NN_trasform = as.data.frame(apply(NN, 1:2, function(x){
  gsub(x,tmp[x,1],x = x)
}))
head(NN_trasform)

# calculate entropy based on KNN results
# 熵的公式应该是P(x)log(P(x))，p(x)就是每个细胞及其周围的邻居中来自不同组别的细胞的发生概率
myentropyFun = function(data=NULL){apply(data, 1, function(x) {
  p_all = 0
  vars_inrow = names(table(t(x)))
  for(i in 1:length(vars_inrow)){
    count_tmp = sum(x == vars_inrow[i])
    p_tmp = ( -log2(count_tmp/length(x)) * count_tmp/length(x) )
    p_all = p_all + p_tmp
  }
  # x=test[1,]
  # count1 = sum(x == 'Alpha')
  # count2 = sum(x == 'Mock')
  # count3 = sum(x == 'Lambda')
  # p_entropy = ( -log2(count1/length(x)) * count1/length(x) ) +
  #   ( -log2(count2/length(x)) * count2/length(x) ) +
  #   ( -log2(count3/length(x)) * count3/length(x) )
  return(p_all)
})}
maxEntropy = (-log2(6/20) * 6/20) +(-log2(7/20) * 7/20)+(-log2(7/20) * 7/20)

p_entropy = myentropyFun(NN_trasform)
NN_entropy = data.frame(row.names = names(p_entropy),
                        entropy = p_entropy,
                        Tcelltype = Sdata$Tcelltype[names(p_entropy)])
NN_entropy$heterogeneity = max(NN_entropy$entropy)-NN_entropy$entropy
head(NN_entropy)
#画图细胞类型按中位数排序
tgc = tapply(NN_entropy$heterogeneity, NN_entropy$Tcelltype , median)
NN_entropy$Tcelltype = factor(NN_entropy$Tcelltype,levels = names(tgc[order(tgc,decreasing = F)]))
head(NN_entropy)

# box plot of entropy
sampleoutpath = './mergedSamples_pure/mergedSamples_pure_3/mergedSamples_rmImmuEry/'
ggplot(NN_entropy, aes(x = Tcelltype , y = heterogeneity, fill= Tcelltype)) +
  stat_boxplot(geom = 'errorbar', aes(ymin=..ymax..,color=Tcelltype),width=0.25,alpha=1) + #仅添加上误差线
  stat_boxplot(geom = 'errorbar', aes(ymax=..ymin..,color=Tcelltype),width=0.25,alpha=1) + #仅添加下误差线
  geom_boxplot(linetype='dashed',aes(color=Tcelltype),width = 0.6) + #加上下竖线，改为虚线
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fill=Tcelltype),width = 0.6,color='black',alpha=1,outlier.shape = NA)+ #仅画boxplot中的box
  # stat_boxplot(geom = "errorbar",width=0.25,aes(color=Tcelltype)) +#正常箱线图
  # geom_boxplot(width = 0.6,alpha=1,aes(color=Tcelltype),outlier.shape = NA) +  #正常箱线图
  scale_fill_manual(values = TcelltypeCols) +  
  scale_color_manual(values = TcelltypeCols) +  
  theme_bw()+
  theme(panel.grid = element_blank(), 
        # axis.line = element_line(colour = 'black', size = 0.8), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.title = element_blank(),
        legend.position = 'none',
        legend.key.size = unit(1,'cm'),
        legend.text = element_text(size = 14, color = 'black'),
        axis.text = element_text(size = 13, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = NULL, y = 'Heterogeneity', title = NULL) 
ggsave(filename = paste0(sampleoutpath,'Entropy/','Heterogeneity_group_all_boxplot_20230914.pdf'),w=7,h=3.5)
ggsave(filename = paste0(sampleoutpath,'Entropy/','Heterogeneity_group_all_boxplot_20240124.pdf'),w=12,h=3.5)

# UMAP plot of entropy
NN_entropy$UMAP1 = Sdata@reductions$umap@cell.embeddings[,1]
NN_entropy$UMAP2 = Sdata@reductions$umap@cell.embeddings[,2]
NN_entropy = NN_entropy[order(NN_entropy$heterogeneity,decreasing = F),]
write.csv(NN_entropy,file = paste0(sampleoutpath,'Entropy/Heterogeneity_group_all_data_20230914.csv'))

sampleoutpath = './mergedSamples_pure/mergedSamples_pure_3/mergedSamples_rmImmuEry/'
NN_entropy = read.csv(paste0(sampleoutpath,'Entropy/Heterogeneity_group_all_data_20230914.csv'),row.names = 1)
head(NN_entropy)
library(tidydr)
library(ggrastr)
ggplot(NN_entropy,aes(x=UMAP1,y=UMAP2,color=heterogeneity))+
  geom_point_rast(size =0.1 , alpha =1 )  +  
  theme_dr()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15),axis.line = element_line(color='black'),
        legend.text = element_text(size=13,color='black'),
        legend.title = element_text(size=13,color='black'),
        legend.direction = "vertical"  #"horizontal" #
  )+
  scale_colour_gradient2(low = "#59C5E6",mid='white', high = "#E27688", na.value = NA,midpoint = 0.8) +
  # scale_colour_gradient2(low = "grey",mid='grey90', high = "red", na.value = NA,midpoint = 0.8) +
  # scale_colour_gradient(low ="grey90",high="red")+
  labs(x='UMAP 1',y='UMAP 2',title = NULL,color='Heterogeneity')
ggsave(filename = paste0(sampleoutpath,'Entropy/','Entropy_group_all_heterogeneity_UMAP.pdf'),w=5.6,h=4)

# density plot of entropy
ggplot(NN_entropy, aes(group = Tcelltype , x = heterogeneity, color= Tcelltype)) + geom_density(alpha=0.3,adjust=1.5) + 
  scale_color_manual(values = TcelltypeCols)

#