##-------------------- ISG geneset score -----------
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
##------- calculate ISG geneset score by AUCell
extractAUCFUN=function(sampleoutpath=NULL,group=NULL,Sdata=NULL,Genes=NULL,top_pct=0.05,return.res=F){
  library(AUCell)
  cells = rownames(Sdata@meta.data)[Sdata$orig.ident == group]
  subData = Sdata[["RNA"]]@data[,cells]
  metaData=Sdata@meta.data[cells,]

  cells_rankings <- AUCell_buildRankings(subData)
  geneSets <- list("ISG score"=Genes)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*top_pct)
  
  df = getAUC(cells_AUC)
  df1=data.frame("cell" = colnames(df),"score" = df[1,])
  df1=df1[cells,]
  df1$celltype=metaData[cells,]$Tcelltype
  df1$group=group
  write.csv(df1,file = paste0(sampleoutpath,group,"_AUCscore.csv"))
  if(isTRUE(return.res)){
      return(df1)
  }
}
groupList=c("Lambda","Alpha","Mock")
sampleoutpath = './mergedSamples_pure/mergedSamples_pure_3/mergedSamples_rmImmuEry/'
sampleoutpath_ISGscore=paste0(sampleoutpath,'ISGscore/')
if(!dir.exists(sampleoutpath_ISGscore)){dir.create(sampleoutpath_ISGscore)}
m=readLines("./ISG_mouse_final.txt")
for(i in 1:3){
  group=groupList[i]
  extractAUCFUN(sampleoutpath = sampleoutpath_ISGscore,Sdata = Sdata, group=group,Genes = m)
}
##------- load AUC results
sampleoutpath = './mergedSamples_pure/mergedSamples_pure_3/mergedSamples_rmImmuEry/'
sampleoutpath_ISGscore=paste0(sampleoutpath,'ISGscore/')
df=NULL
groupList=c("Lambda","Alpha","Mock")
for (i in 1:3) {
  i=1
  group=groupList[i]
  df_tmp=read.csv(paste0(sampleoutpath_ISGscore,group,"_AUCscore.csv"),row.names = 1)
  df=rbind(df,df_tmp)
}
df$group=factor(df$group,levels = c("Mock","Alpha","Lambda"))

df0=df
# olfactory
ct=c("Olfactory HBC","GBC","INP","OSN","SUS","BGC","Oligodendrocyte")
ctlabels = c("Olfactory HBC","GBC","INP","OSN","SUS","BGC","Oligodendrocyte")
# respiratory
ct=c("RBC","Ciliated cell","Fezf2+ epithelial cell",
     "Cyp4b1+ gland cell","Dmbt1+ gland cell","Bpifb5+ gland cell")
ctlabels=c("RBC","Ciliated cell","Fezf2+ epithelial cell",
     "Cyp4b1+ gland cell","Dmbt1+ gland cell","Bpifb5+ gland cell")
# ctlabels=c("RBC","Ciliated cell","Fezf2+\nepithelial cell",
#      "Cyp4b1+\ngland cell","Dmbt1+\ngland cell","Bpifb5+\ngland cell")
# common function
ct=c("Goblet cell","Club cell", "Hillock cell","Ionocyte","Tuft cell","Obp1a+ gland cell")
ctlabels = c("Goblet cell","Club cell", "Hillock cell","Ionocyte","Tuft cell","Obp1a+\ngland cell")
# other cell type
ct=c("VSN", "Astrocyte","M cell","Dental basal cell","Fibroblast","Endothelial cell")
ctlabels = c("VSN", "Astrocyte","M cell","Dental\nbasal cell","Fibroblast","Endothelial\ncell") #标题换行

df=df0[which(df0$celltype %in% ct),]
df$celltype=factor(df$celltype,levels = ct,labels = ctlabels)

### violin plot for ISG geneset score calculated by AUC
compaired <- list(c("Alpha","Mock"),c("Lambda","Mock"),c("Alpha", "Lambda"))
p=ggplot(df,aes(group,score,fill=group))+
  # stat_boxplot(geom="errorbar",width=0.35,color="black")+
  geom_violin(scale = "width",alpha=1,width=0.9,size=0.01)+
  geom_boxplot(outlier.shape = NA,width=0.15,color='black',fill='white',size=0.01)+
  # stat_summary(fun=median,geom='point',color="black",shape=18)+
  facet_wrap(.~celltype,scales = "free_y",nrow = 1)+
  scale_fill_manual(values = tissueCols,labels=c("Lambda"=expression("IFN-"*lambda),
                                                 "Alpha"=expression("IFN-"*alpha),"Mock"="Mock")) +
  # scale_color_manual(values = tissueCols) +
  labs(x=NULL,y="ISG score",title=NULL,fill="Group") +
  theme_classic() +
  theme(plot.title=element_text(size = 13,hjust = 0.5,colour = "black"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.title=element_text(size = 15,colour = "black"),
        axis.line.x.bottom = element_line(colour = "black",size=0.6),
        axis.line.y.left = element_line(colour = "black",size=0.6),
        panel.grid = element_blank(),
        legend.text = element_text(size=12,colour = "black"),
        legend.text.align = 0,
        strip.placement = "outside", #分面标签位置
        strip.text.x = element_text(angle=0,vjust=0,hjust = 0.5,size=15), #分面标签 文字倾斜; 字号
        strip.background = element_blank()#, #分面标签 不要背景
  ) +
  geom_signif(data=df,comparisons = compaired,
              # y_position = c(0.31,0.27,0.18),#横线标记的位置c(0.3,0.28,0.26,0.12,0.1)
              tip_length = 0,#连线的长度
              step_increase = 0.1,
              map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05," "=2),
              vjust = 0.5,textsize = 5,
              # map_signif_level=function(s)sprintf("p = %.2g", s),
              test = "wilcox.test")
p
ggsave(p,filename = paste0(sampleoutpath_ISGscore,"ISGgenes_score_celltype_AUC_LambdaHigh.pdf"),w=17.2,h=2.4,useDingbats = F)
ggsave(p,filename = paste0(sampleoutpath_ISGscore,"ISGgenes_score_celltype_AUC_AlphaHigh.pdf"),w=15.2,h=2.4,useDingbats = F)
ggsave(p,filename = paste0(sampleoutpath_ISGscore,"ISGgenes_score_celltype_AUC_equal.pdf"),w=12.2,h=2.6,useDingbats = F)
ggsave(p,filename = paste0(sampleoutpath_ISGscore,"ISGgenes_score_celltype_AUC_other.pdf"),w=12.2,h=2.6,useDingbats = F)

## single gene expression between 3 groups
gene= 'Tslp' #'Usp18' #'Socs1'
df = data.frame(row.names = rownames(Sdata@meta.data),
                score=Sdata[['RNA']]@counts[gene,],
                group=Sdata$orig.ident,
                celltype=Sdata$Tcelltype)
### violin plot for ISG geneset score calculated by AUC
compaired <- list(c("Alpha","Mock"),c("Lambda","Mock"),c("Alpha", "Lambda"))
p=ggplot(df,aes(group,score,fill=group))+
  stat_boxplot(geom="errorbar",width=0.35,color="black")+
  geom_boxplot()+
  # geom_violin(scale = "width",alpha=1,width=0.9,size=0.05)+
  stat_summary(fun=median,geom='point',color="black",shape=18)+
  facet_wrap(.~df$celltype,scales = "free_y",nrow = 4)+
  scale_fill_manual(values = tissueCols,labels=c("Lambda"=expression("IFN-"*lambda),
                                                 "Alpha"=expression("IFN-"*alpha),"Mock"="Mock")) +
  # scale_color_manual(values = tissueCols) +
  labs(x=NULL,y=gene,title=NULL,fill="Group") +#"ISG score"
  theme_classic() +
  theme(plot.title=element_text(size = 13,hjust = 0.5,colour = "black"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=12,colour = "black"),
        axis.title=element_text(size = 13,colour = "black"),
        axis.line.x.bottom = element_line(colour = "black",size=0.6),
        axis.line.y.left = element_line(colour = "black",size=0.6),
        panel.grid = element_blank(),
        legend.text = element_text(size=12,colour = "black"),
        legend.text.align = 0,
        strip.placement = "outside", #分面标签位置
        strip.text.x = element_text(angle=0,vjust=0,hjust = 0.5,size=13), #分面标签 文字倾斜; 字号
        strip.background = element_blank()#, #分面标签 不要背景
  ) +
  geom_signif(data=df,comparisons = compaired,
              # y_position = c(0.31,0.27,0.18),#横线标记的位置c(0.3,0.28,0.26,0.12,0.1)
              tip_length = 0.01,#连线的长度
              step_increase = 0.1,
              map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05,"ns"=2),
              vjust = 0.4,
              # map_signif_level=function(s)sprintf("p = %.2g", s),
              test = "wilcox.test")
p
# ggsave(p,filename = paste0("ISGgenes_score_celltype_AUC_LambdaHigh.pdf"),w=9.2,h=2.3,useDingbats = F)
#