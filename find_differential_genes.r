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

##---------------- find differential genes #################
### find differential genes between any two group
#匿名函数
(function(){
  group1="Alpha"
  group2="Lambda" 
  group3="Mock"
  
  df=NULL
  # celltypelist=c("Sustentacular cell","Goblet cell","Hillock cell")
  # celltypelist=c("Erythroblast")
  # celltypelist=setdiff(unique(Sdata$Tcelltype),celltypelist)
  celltypelist=unique(Sdata$Tcelltype)
  for (ct in celltypelist) {
    # for (ct in levels(Sdata$Tcelltype)) {
    # ct=levels(Sdata$Tcelltype)[1]
    subobj=subset(Sdata,cells=rownames(Sdata@meta.data)[Sdata$Tcelltype==ct & Sdata$orig.ident %in% c(group1,group2)])
    Idents(subobj)="orig.ident"
    df0=FindMarkers(subobj,ident.1 = group2,ident.2 = group1,logfc.threshold = 0.001,min.pct = 0.0005)
    df0$Tcelltype=ct
    df0$gene=rownames(df0)
    df0=df0[df0$p_val<0.05,]
    df=rbind(df,df0)
  }
  df$compare=paste0(group2,"_vs_",group1)
  # write.csv(df,file=paste0(sampleoutpath,"DFgenes_",group1,"_vs_",group2,"_3ct.csv"))
  write.csv(df,file=paste0(sampleoutpath,"DFgenes_",group2,"_vs_",group1,"_25ct.csv"))
  
  
  df=NULL
  for (ct in celltypelist) {
    # ct=levels(Sdata$Tcelltype)[1]
    subobj=subset(Sdata,cells=rownames(Sdata@meta.data)[Sdata$Tcelltype==ct & Sdata$orig.ident %in% c(group1,group3)])
    Idents(subobj)="orig.ident"
    df0=FindMarkers(subobj,ident.1 = group1,ident.2 = group3,logfc.threshold = 0.001,min.pct = 0.0005)
    df0$Tcelltype=ct
    df0$gene=rownames(df0)
    df0=df0[df0$p_val<0.05,]
    df=rbind(df0,df)
  }
  df$compare=paste0(group1,"_vs_",group3)
  write.csv(df,file=paste0(sampleoutpath,"DFgenes_",group1,"_vs_",group3,"_25ct.csv"))
  
  df=NULL
  for (ct in celltypelist) {
    # ct=levels(Sdata$Tcelltype)[1]
    subobj=subset(Sdata,cells=rownames(Sdata@meta.data)[Sdata$Tcelltype==ct & Sdata$orig.ident %in% c(group2,group3)])
    Idents(subobj)="orig.ident"
    df0=FindMarkers(subobj,ident.1 = group2,ident.2 = group3,logfc.threshold = 0.001,min.pct = 0.0005)
    df0$Tcelltype=ct
    df0$gene=rownames(df0)
    df0=df0[df0$p_val<0.05,]
    df=rbind(df0,df)
  }
  df$compare=paste0(group2,"_vs_",group3)
  write.csv(df,file=paste0(sampleoutpath,"DFgenes_",group2,"_vs_",group3,"_25ct.csv"))
})()
#

##---------------- 将差异基因写成metascape中可直接读入的多list格式 ---------------
comparisonList = c("Lambda_vs_Mock","Alpha_vs_Mock","Lambda_vs_Alpha")
group=comparisonList[3]
DFgenes=read.csv(paste0(sampleoutpath,"DFgenes/DFgenes_",group,"_25ct.csv"),row.names = 1)
DFgenes$compare = group

##--------- single group high genelist
MaxLength=max(as.data.frame(table(DFgenes[which(DFgenes$avg_log2FC>0 & DFgenes$p_val_adj<0.05),"Tcelltype"]))[2])#最多的基因数
DFgenelist <- c()
for(i in unique(DFgenes$Tcelltype)){
  # i="Bpifp5+ gland cell"
  genes = DFgenes$gene[which(DFgenes$Tcelltype==i & DFgenes$avg_log2FC>0 & DFgenes$p_val_adj<0.05)]
  length(genes)=MaxLength  ##合并不等长的向量：先让其长度均等于最大向量长度再合并
  DFgenelist = cbind(DFgenelist,genes)
}
colnames(DFgenelist) <- paste0('Lambda_high_',unique(DFgenes$Tcelltype))
write.table(DFgenelist,file=paste0(sampleoutpath,'DFgenes/',group,"_DFgenelist_logFC0_padj0.05_Lambda_high_25ct.csv"),sep = ",",col.names = T,row.names = F)

##--------- combined both genelist
MaxLength=max(as.data.frame(table(DFgenes[which(DFgenes$avg_log2FC<0 & DFgenes$p_val_adj<0.05),"Tcelltype"]))[2])#最多的基因数
DFgenelist <- c()
for(i in unique(DFgenes$Tcelltype)){
  genes = DFgenes$gene[which(DFgenes$Tcelltype==i & DFgenes$avg_log2FC>0 & DFgenes$p_val_adj<0.05)]
  length(genes)=MaxLength  
  DFgenelist = cbind(DFgenelist,genes)
}
colnames(DFgenelist) <- paste0('Lambda_high_',unique(DFgenes$Tcelltype))
for(i in unique(DFgenes$Tcelltype)){
  genes = DFgenes$gene[which(DFgenes$Tcelltype==i & DFgenes$avg_log2FC<0 & DFgenes$p_val_adj<0.05)]
  length(genes)=MaxLength  
  DFgenelist = cbind(DFgenelist,genes)
}
colnames(DFgenelist)[26:50] <- paste0('Alpha_high_',unique(DFgenes$Tcelltype))
dim(DFgenelist)

DFgenelist_clean = DFgenelist
rmCol=''
for(i in colnames(DFgenelist)){
  if(table(is.na(DFgenelist[,i]))[1]<5){
    rmCol = c(rmCol,i)
  }
}    
DFgenelist_clean = DFgenelist[,-which(colnames(DFgenelist) %in% rmCol)]
write.table(DFgenelist,file=paste0(sampleoutpath,'DFgenes/',group,"_DFgenelist_logFC0_padj0.05_Lambda_vs_Alpha_clean_25ct.csv"),sep = ",",col.names = T,row.names = F)


### 只画一个对比的差异基因数量图
comparisonList = c("Lambda_vs_Mock","Alpha_vs_Mock","Lambda_vs_Alpha")
i=3
group=comparisonList[i]
DFgenes=read.csv(paste0("./DFgenes_",group,"_27ct.csv"),row.names = 1)
df0=as.data.frame(table(DFgenes$Tcelltype[DFgenes$p_val_adj<0.05 & DFgenes$avg_logFC<0]))
colnames(df0)=c("Tcelltype","number")
df0=df0[order(df0$number,decreasing = T),]
df0$Tcelltype=factor(df0$Tcelltype,levels = unique(df0$Tcelltype))
ggplot(df0,aes(x=Tcelltype,y=number,fill=group))+geom_bar(stat = "identity") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 12,color = "black"),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.title.y = element_text(size = 14,color = "black"),
        legend.text=element_text(size = 12,color = "black"),
        legend.title=element_text(size = 12,color = "black"),
        legend.position = c(0.8,0.8))+
  scale_fill_manual(values = c("Lambda_vs_Mock"=tissueCols[3],
                               "Alpha_vs_Mock"=tissueCols[2],
                               "Lambda_vs_Alpha"=tissueCols[1]),
                    labels = c("Lambda_vs_Mock"="Lambda vs Mock",
                               "Alpha_vs_Mock"="Alpha vs Mock",
                               "Lambda_vs_Alpha"="Lambda vs Alpha"))+
  labs(x=NULL,y="Differential gene number",fill="Comparison")
ggsave(paste0(comparisonList[i],"_DFgenes_number_padj0.05_logFCless0.pdf"),w=6.5,h=5.2)

#