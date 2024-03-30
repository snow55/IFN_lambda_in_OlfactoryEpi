library(Seurat)
library(ggplot2)
library(dplyr)
library(plyr)
library(Cairo)
library(RColorBrewer)
library(pheatmap)
library(Hmisc)
library(ggpubr)

TctCols=c("#E889BD","#67C2A3","#8EA0C9")
tissueCols = c("#0a8cf7","#febe08","#ed2b2b")#new 更鲜艳配色 

TcelltypeCols=c("#e986b7","#f6999a","#d84b71","#EC188B", "#e2a226",
                "#f9a270","#f7aa5d","#e86502",
                '#f2991a',"#c5aac3",'#cf6eb8',
                "#fdcf16",'#76c6ba',
                "#648bc9","#368b9d","#be71af",
                "#d7ce96",
                '#1a734c','#3387b5',"#9452a0","#1e77fe",'#72daf2',
                "#d4af29","#b3d78b","#9b9a9b"#,"#e57371",Erythrocyte
                # "#cab2d6" Immune cell
)

##-------------------- multi violin plot for each Tcelltype
Genes=readLines("./vlnplots_Genes.txt")

Genes = c('Krt5','Krt14','Nrcam','Vit','Ascl1','Neurod1','Neurog1','Gap43','Omp','Stoml3', #Olfactory HBC - OSN
          'Muc2','Chil6', #SUS - BGC
          'Calr4','Plp1','Slc6a11', #VSN - Astrocyte
          'Adh7','Cyp2f2','Ly6d','Dnah5','Fezf2','Cyp4b1','Dmbt1','Bpifb5',
          'Muc5b','Cyp4a12a','Tnfaip2','Gp2',
          'Krt13','Cftr','Lrmp','Klk14',
          'Ambn','Col1a1','Cdh5'
)
# Genes=readLines("./MarkerGenes.txt")
Genes=intersect(Genes,rownames(Sdata[["RNA"]]@data))
df=as.data.frame(t(as.matrix(Sdata[["RNA"]]@data[Genes,])))
df$celltype=Sdata$Tcelltype
library(reshape2)
df1=melt(df,id = "celltype")

##纵向vlnplot
p2 = ggplot(df1,aes(x=celltype,y=value,fill=celltype,color=celltype))+
  geom_violin(#color="white",#aes(color = cluster), 
    trim = TRUE,
    scale = "width",alpha=1,width=0.9,size=0.05) +
  facet_wrap(.~variable,ncol = 1, strip.position = "left",scales = "free_y") + #
  scale_color_manual(values = TcelltypeCols)+
  scale_fill_manual(values = TcelltypeCols)+
  # scale_y_continuous(expand = c(0,0),position = "left")+
  scale_y_continuous(expand = c(0, 0), position="right" #, 
    #                  labels = function(x)
    # c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")
    ) + #标记y轴数值
  theme_bw() +
  theme(
    panel.grid = element_blank(), #不要背景网格
    
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank(),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=45,vjust = 1,hjust = 1,size=13,color="black"), #y坐标刻度字号
    axis.title = element_blank(), #y标题字号
    
    legend.position = "none", #不要图例
    panel.spacing=unit(0,"cm"), #分面的间距
    
    strip.placement = "outside", #分面标签位置
    strip.text.y.left = element_text(angle = 0,size=13,color="black",hjust = 1),#调整分面变量横向靠右对齐
    # strip.text.x = element_text(angle=-90,vjust=0.5,hjust = 1,size=13,color="black"), #分面标签 文字倾斜; 字号
    strip.background = element_blank() #分面标签 不要背景
  )+
  labs(y=NULL,x=NULL)
print(p2)
ggsave(p2,filename = paste0(sampleoutpath,"./multiVlnplot_MarkerGenes_column_20230914.pdf"),w=6.8,h=12)


##横向vlnplot
df1$celltype=factor(df1$celltype,levels = rev(levels(df1$celltype)))
p2 = ggplot(df1,aes(x=value,y=celltype,fill=celltype,color=celltype))+
  geom_violin(#aes(color = celltype), 
    # trim = FALSE,
    color=NA,#"white",
    scale = "width",
    alpha=1,
    width=0.9,size=0.5) +
  facet_wrap(variable~.,nrow = 1, strip.position = "top",scale='free_x') + #
  scale_color_manual(values = rev(TcelltypeCols))+
  scale_fill_manual(values = rev(TcelltypeCols))+
  scale_x_continuous(expand = c(0,0),position = "top")+
  theme_bw() +
  theme(
    panel.grid = element_blank(), #不要背景网格
    # axis.line.x.bottom = element_line(color="black",size=0.6),
    # axis.line.y.left = element_line(color="black",size=0.6),
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13,color="black"), #y坐标刻度字号
    axis.title.x = element_blank(), #y标题字号
    
    legend.position = "none", #不要图例
    panel.spacing=unit(0,"cm"), #分面的间距
    
    strip.placement = "outside", #分面标签位置
    # strip.text.y = element_text(angle=-90,vjust=0,hjust = 0,size=15), #分面标签 文字倾斜; 字号
    strip.text.x = element_text(angle=90,vjust=0.5,hjust = 0,size=13,color="black"), #分面标签 文字倾斜; 字号
    strip.background = element_blank() #分面标签 不要背景
  )+
  labs(y=NULL,x=NULL)
print(p2)
ggsave(p2,filename = paste0(sampleoutpath,"./multiVlnplot_MarkerGenes_row.pdf"),w=12,h=4)


###--------------- violin plot for several genes * several ct  #EDFig 5,6
Genes = c("Mx1", "Isg15", "Ifit1","Ifit3", "Rsad2", "Oasl2", "Stat1","Stat2",
          "Ifi44", "Ifit3b", "Irgm1", "Gbp7")
# Genes = c("Mx1", "Isg15","Isg20", "Stat1","Stat2", "Rsad2","Oasl1", "Oasl2", "Ifit1","Ifit3","Ifi44", 
          # "Ifit3b","Ifi35","Irgm1","Gbp3", "Gbp7","Usp18","Irf7","Bst2","Trim25")
Genes = c("Ace2","Tmprss2","Furin")
Genes = intersect(Genes,rownames(Sdata[["RNA"]]@data))

df=as.data.frame(t(as.data.frame(Sdata[["RNA"]]@data[Genes,])))
df$group=Sdata$orig.ident
df$celltype=Sdata$Tcelltype
head(df)

ct=c("SUS","BGC")
# ct=c("Ciliated cell","Fezf2+ epithelial cell")
# ct=c("Goblet cell","Club cell")
tmp=data.frame(row.names = rownames(df),
               # exp=df[,gene],
               group=df$group,
               celltype=df$celltype)
tmp=cbind(tmp,df[,Genes])
tmp=tmp[which(tmp$celltype %in% ct),]
library(reshape2)
tmp=melt(tmp,id.vars = c("group","celltype"))
##vlnplot
# mycolor=TcelltypeCols[c(5:8,12,13)]
mycolor=tissueCols
compaired <- list(c("Alpha","Mock"),c("Lambda","Mock"),c("Alpha", "Lambda"))

ggplot(tmp,aes(x=group,y=value,fill=group,color=group))+
  geom_violin(scale = "width",alpha=1,width=0.9)+
  # geom_jitter(size=0.8,alpha=0.6)+
  # geom_boxplot(width=0.2,color='black',fill='white',outlier.shape = NA,linewidth=0.1)+
  facet_grid(celltype~variable,scales = "free")+
  theme_classic() +
  theme(plot.title=element_text(size = 16,hjust = 0.5,colour = "black"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.title=element_text(size = 16,colour = "black"),
        axis.line.x.bottom = element_line(colour = "black",size=0.6),
        axis.line.y.left = element_line(colour = "black",size=0.6),
        panel.grid = element_blank(),
        legend.text = element_text(size=15,colour = "black"),
        legend.text.align = 0,
        strip.placement = "outside", #分面标签位置
        strip.text.y = element_text(size=15,color="black"), #分面标签 文字倾斜; 字号
        strip.text.x = element_text(angle=0,vjust=0,hjust = 0.5,size=17,color="black"), #分面标签 文字倾斜; 字号
        strip.background.x = element_blank()#, #分面标签 不要背景
  ) +
  geom_signif(data=tmp,comparisons = compaired,
              # y_position = c(0.31,0.27,0.18),#横线标记的位置c(0.3,0.28,0.26,0.12,0.1)
              tip_length = 0,#连线的长度
              color="black",
              map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, " "=2), #将NS改为ns
              # map_signif_level = T,
              step_increase = 0.1,vjust = 0.5,textsize = 5,
              # map_signif_level=function(s)sprintf("p = %.2g", s),
              test = "wilcox.test")+
  scale_color_manual(values = mycolor,labels=c("Mock"="Mock","Alpha"=expression("IFN-"*alpha),
                                               "Lambda"=expression("IFN-"*lambda))) +
  scale_fill_manual(values = mycolor,labels=c("Mock"="Mock","Alpha"=expression("IFN-"*alpha),
                                              "Lambda"=expression("IFN-"*lambda))) +
  labs(x="",y="Expression",title=NULL,fill="Group",color="Group")

sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_3/mergedSamples_rmImmuEry/"
vlnPath=paste0(sampleoutpath,"violinplot/")
if(!dir.exists(vlnPath)){dir.create(vlnPath,recursive = T)}
ggsave(filename = paste0(vlnPath,"vlnplot_","12IsgGenes","_","2ct_SUS-BGC","_exp_20230914.pdf"),w=17,h=4.4)#
ggsave(filename = paste0(vlnPath,"vlnplot_","12IsgGenes","_","2ct_Ciliated-Fezf2","_exp_20230914.pdf"),w=17,h=4.7)#SUS-BGC

ggsave(filename = paste0(vlnPath,"vlnplot_","Ace2-Tmprss2-Furin","_","SUS-BGC","_exp_20230914.pdf"),w=7,h=4)
# ggsave(filename = paste0(vlnPath,"vlnplot_","Ace2-Tmprss2","_","Ciliated-Fezf2","_exp_20230914.pdf"),w=4.2,h=4.3)

### calculate exact p value
head(tmp)
useData = tmp[tmp$celltype=='BGC' & tmp$variable=='Ace2',]
wilcox.test(useData[useData$group==groupList[1],'value'],useData[useData$group==groupList[2],'value'],alternative='less')
wilcox.test(useData[useData$group==groupList[1],'value'],useData[useData$group==groupList[2],'value'],alternative='two.sided')
# two.sided
# SUS: Tmprss2: 
# Lambda < Mock p-value = 1.217e-13
# Alpha < Mock p-value = 1.79e-09
# Lambda < Alpha p-value = 0.009841
# BGC: Tmprss2: 
# Lambda < Alpha p-value = 0.02399
# BGC: Ace2: 
# Lambda < Mock p-value = 0.002486
# Lambda < Alpha p-value = 0.0001566

##--------------- IFN receptors expression 
Genes=c("Ifnlr1","Il10rb","Ifnar1","Ifnar2")
df=as.data.frame(t(as.data.frame(Sdata[["RNA"]]@data[Genes,])))
df$group=Sdata$orig.ident
df$celltype=Sdata$Tcelltype
head(df)

# olfactory
ct=c("Olfactory HBC","GBC","INP","OSN","SUS","BGC","Oligodendrocyte")
ctlabels = c("Olfactory\nHBC","GBC","INP","OSN","SUS","BGC","Oligodendrocyte")
# respiratory
ct=c("RBC","Ciliated cell","Fezf2+ epithelial cell",
     "Cyp4b1+ gland cell","Dmbt1+ gland cell","Bpifb5+ gland cell")
ctlabels=c("RBC","Ciliated cell","Fezf2+\nepithelial cell",
           "Cyp4b1+\ngland cell","Dmbt1+\ngland cell","Bpifb5+\ngland cell")
# # other cell type
# ct=setdiff(levels(usemetaData$Tcelltype),
#            c(c("Olfactory HBC","GBC","INP","OSN","SUS","BGC"),
#              c("Respiratory basal cell","Ciliated cell","Fezf2+ epithelial cell",
#                "Cyp4b1+ gland cell","Dmbt1+ gland cell","Bpifb5+ gland cell")))
# ct=levels(Sdata$Tcelltype)

gene=c("Ifnar1","Ifnar2")
gene=c("Ifnlr1","Il10rb")
tmp=data.frame(row.names = rownames(df),
               # exp=df[,gene],
               group=df$group,
               celltype=df$celltype)
tmp=cbind(tmp,df[,gene])
tmp=tmp[which(tmp$celltype %in% ct),]
tmp$celltype=factor(tmp$celltype,levels = ct,labels = ctlabels)
library(reshape2)
tmp=melt(tmp,id.vars = c("group","celltype"))
##vlnplot
# mycolor=TcelltypeCols[c(5:8,12,13)]
mycolor=tissueCols
compaired <- list(c("Alpha","Mock"),c("Lambda","Mock"),c("Alpha", "Lambda"))

ggplot(tmp,aes(x=group,y=value,fill=group,color=group))+
  geom_violin(scale = "width",alpha=1,width=0.9)+
  geom_jitter_rast(size=0.8,alpha=0.6)+
  facet_grid(variable~celltype,scales = "free")+
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
        strip.text.y = element_text(size=14,color="black"), #分面标签 文字倾斜; 字号
        strip.text.x = element_text(angle=0,vjust=0,hjust = 0.5,size=13,color="black"), #分面标签 文字倾斜; 字号
        strip.background.x = element_blank()#, #分面标签 不要背景
  ) +
  geom_signif(data=tmp,comparisons = compaired,
              # y_position = c(0.31,0.27,0.18),#横线标记的位置c(0.3,0.28,0.26,0.12,0.1)
              tip_length = 0,#连线的长度
              color="black",
              step_increase = 0.1,
              map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, " "=2), #将NS改为ns
              vjust = 0.5,textsize = 5,
              # map_signif_level=function(s)sprintf("p = %.2g", s),
              test = "wilcox.test")+
  scale_color_manual(values = mycolor,labels=c("Mock"="Mock","Alpha"=expression("IFN-"*alpha),
                                               "Lambda"=expression("IFN-"*lambda))) +
  scale_fill_manual(values = mycolor,labels=c("Mock"="Mock","Alpha"=expression("IFN-"*alpha),
                                              "Lambda"=expression("IFN-"*lambda))) +
  labs(x="",y="Expression",title=NULL,fill="Group",color="Group")

vlnPath=paste0(sampleoutpath,"violinplot/")
if(!dir.exists(vlnPath)){dir.create(vlnPath,recursive = T)}
ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnlr1-Il10rb","_","all","_exp_withsig.pdf"),w=17,h=3.5)#OlfactoryEpi
ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnar1-Ifnar2","_","all","_exp_withsig.pdf"),w=17,h=3.5)#OlfactoryEpi

ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnlr1-Il10rb","_","ResEpi","_exp_withsig.pdf"),w=7,h=3.5)#OlfactoryEpi
ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnar1-Ifnar2","_","ResEpi","_exp_withsig.pdf"),w=7,h=3.5)#OlfactoryEpi

ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnlr1-Il10rb","_","OlfEpi","_exp_withsig.pdf"),w=6.5,h=3.5)#OlfactoryEpi
# ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnar1-Ifnar2","_","ResEpi","_exp_withsig.pdf"),w=6.5,h=3.5)#OlfactoryEpi
# ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnlr1-Il10rb","_","OlfEpi","_exp_nosig.pdf"),w=6.5,h=3)#OlfactoryEpi
# ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnar1-Ifnar2","_","ResEpi","_exp_nosig.pdf"),w=6.5,h=3)#OlfactoryEpi


##-------------------- IFN receptors in Mock 
Genes=c("Ifnlr1","Il10rb","Ifnar1","Ifnar2")
df=as.data.frame(t(as.data.frame(Sdata[["RNA"]]@data[Genes,])))
df$group=Sdata$orig.ident
df$celltype=Sdata$Tcelltype
head(df)

# olfactory
ct=c("Olfactory HBC","GBC","INP","OSN","SUS","BGC","Oligodendrocyte")
ctlabels = c("Olfactory\nHBC","GBC","INP","OSN","SUS","BGC","Oligodendrocyte")
# respiratory
ct=c("RBC","Ciliated cell","Fezf2+ epithelial cell",
     "Cyp4b1+ gland cell","Dmbt1+ gland cell","Bpifb5+ gland cell")
ctlabels=c("RBC","Ciliated cell","Fezf2+\nepithelial cell",
           "Cyp4b1+\ngland cell","Dmbt1+\ngland cell","Bpifb5+\ngland cell")
# respiratory
ct=c("Ciliated cell","Fezf2+ epithelial cell")
ctlabels=c("Ciliated cell","Fezf2+\nepithelial cell")

gene=c("Ifnar1","Ifnlr1")
tmp=data.frame(row.names = rownames(df),
               # exp=df[,gene],
               group=df$group,
               celltype=df$celltype)
tmp=cbind(tmp,df[,gene])
tmp$group = as.character(tmp$group)
tmp=tmp[which(tmp$celltype %in% ct),]
tmp = tmp[which(tmp$group %in% 'Mock'),]
tmp$celltype=factor(tmp$celltype,levels = ct,labels = ctlabels)
library(reshape2)
tmp=melt(tmp,id.vars = c("group","celltype"))
##vlnplot
# mycolor=TcelltypeCols[c(5:8,12,13)]
mycolor=c("#febe08", "#ed2b2b")
# compaired <- list(c("Alpha","Mock"),c("Lambda","Mock"),c("Alpha", "Lambda"))
compaired <- list(gene)

ggplot(tmp,aes(x=variable,y=value,fill=variable,color=variable))+
  geom_violin(scale = "width",alpha=1,width=0.9)+
  geom_jitter_rast(size=0.8,alpha=0.6)+
  facet_grid(~celltype,scales = "free")+
  theme_classic() +
  theme(plot.title=element_text(size = 13,hjust = 0.5,colour = "black"),
        axis.text.x=element_text(size=14,colour = "black"),#,angle=45,hjust = 1,vjust = 1
        axis.text.y=element_text(size=12,colour = "black"),
        axis.title=element_text(size = 14,colour = "black"),
        axis.line.x.bottom = element_line(colour = "black",size=0.6),
        axis.line.y.left = element_line(colour = "black",size=0.6),
        panel.grid = element_blank(),
        legend.text = element_text(size=12,colour = "black"),
        legend.text.align = 0, legend.position = 'none',
        strip.placement = "outside", #分面标签位置
        strip.text.y = element_text(size=15,color="black"), #分面标签 文字倾斜; 字号
        strip.text.x = element_text(angle=0,vjust=0,hjust = 0.5,size=13,color="black"), #分面标签 文字倾斜; 字号
        strip.background.x = element_blank()#, #分面标签 不要背景
  ) +
  geom_signif(data=tmp,comparisons = compaired,
              # y_position = c(0.31,0.27,0.18),#横线标记的位置c(0.3,0.28,0.26,0.12,0.1)
              tip_length = 0,#连线的长度
              color="black",
              step_increase = 0.1,
              map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, " "=2), #将NS改为ns
              vjust = 0.5,textsize = 5,
              # map_signif_level=function(s)sprintf("p=%.2g", s),
              test = "wilcox.test")+
  scale_color_manual(values = mycolor) +
  scale_fill_manual(values = mycolor) +
  labs(x="",y="Expression",title=NULL,fill="",color="")

vlnPath=paste0(sampleoutpath,"violinplot/")
if(!dir.exists(vlnPath)){dir.create(vlnPath,recursive = T)}
ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnar1-Ifnlr1","_","ResEpi","_exp_withsig.pdf"),w=7,h=3.1)#OlfactoryEpi
ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnar1-Ifnlr1","_","OlfEpi","_exp_withsig.pdf"),w=7,h=3.1)#OlfactoryEpi

ggsave(filename = paste0(vlnPath,"vlnplot_","Ifnar1-Ifnlr1","_","Ciliated-Fezf2+Epi","_exp_withsig.pdf"),w=5,h=2.5)#OlfactoryEpi


tmp1 = tmp[which(tmp$celltype=='Ciliated cell'),]
wilcox.test(tmp1$value[tmp1$variable=='Ifnlr1'],tmp1$value[tmp1$variable=='Ifnar1'],alternative = c("less"))#, "greater""two.sided", 