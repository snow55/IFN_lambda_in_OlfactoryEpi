#################### 过滤完重新分析后的数据的 densityplot, Dotplot, Geneplot and Vlnplot ###############
####### Dotplot 
Genes <- read.table("./mGenesDot_Genes.txt",header = F)##wholeMarkers
Genes <- read.table("./celltypeMarker_selected_2.txt",header = F)##wholeMarkers
# Genes <- read.table("./wholeMarkers_20190727.txt",header = F)##wholeMarkers
Genes <- Genes$V1 ##wholeMarkers
# Genes <- capitalize(tolower(Genes))
Genes <- intersect(Genes,rownames(Sdata[["RNA"]]@counts))
dotplot_outpath <- paste0(sampleoutpath,sample,"_Dotplot/")#_recluster
if(!dir.exists(dotplot_outpath)){
  dir.create(dotplot_outpath)
}
pdf(paste0(dotplot_outpath, sample,"_celltypeGenes_2_reorder_labeled",".pdf"),
    w=15+ceiling(length(Genes)/4),h=2+ceiling(length(Genes)/8),useDingbats = T)
# pdf(paste0(dotplot_outpath, sample,"_wholeGenes",".pdf"),
#     w=15+ceiling(length(Genes)/4),h=2+ceiling(length(Genes)/6),useDingbats = T)
marker_dot_plot <- DotPlot(object = Sdata, features = Genes,group.by = "celltype",
                           cols = c("lightgrey","red"),col.min = 0, dot.scale = 8)
marker_dot_plot <- marker_dot_plot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                           axis.text=element_text(size=20))
print(marker_dot_plot)
dev.off()

######### Geneplot and Vlnplot
pdf(paste0(sampleoutpath,"Geneplot_filtered_",sample,".pdf"),width = 12,height = 6)
plot1 <- FeatureScatter(Sdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Sdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(CombinePlots(plots = list(plot1, plot2)))
dev.off()
# pdf(paste0(sampleoutpath,"VlnPlot_batch_Smallpoint_",sample,".pdf"),width = 8,height = 6)
# vp <- VlnPlot(Sdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#               group.by = "orig.ident", ncol = 3,pt.size=0)
# print(vp)
# dev.off()

pdf(paste0(sampleoutpath,"VlnPlot_cluster_",sample,".pdf"),width = 8,height = 10)
vp <- VlnPlot(Sdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "seurat_clusters", ncol = 1,pt.size=0.1)
print(vp)
dev.off()

######### boxplot of nFeatures
# m <- aggregate(Sdata@meta.data$nFeature_RNA,list(Sdata@meta.data$seurat_clusters),summary)
pdf(paste0(sampleoutpath, sample,"_res0.8_nFeature_boxplot",".pdf") , w=12,h=6) #_reorder
p <- ggplot(NULL,aes(y=Sdata$nFeature_RNA,x=Sdata@active.ident)) + geom_boxplot() +  
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) +
  labs(x="cluster",y="nFeature_RNA")+
  geom_hline(aes(yintercept=c(200)),color="red", linetype="dashed", size=0.5) +  
  geom_hline(aes(yintercept=c(300)),color="red", linetype="dashed", size=0.5) +  
  geom_hline(aes(yintercept=c(500)),color="red", linetype="dashed", size=0.5)
print(p)
dev.off()

######## densityplot
sampleoutpath1 <- paste0(sampleoutpath,"densityplot_cluster/")#_filtered_50
if(!dir.exists(sampleoutpath1)){dir.create(sampleoutpath1)}
library(grid)
for(f in 1:3){
  # Sdata[["umipergene"]] <- Sdata[["nCount_RNA"]]/Sdata[["nFeature_RNA"]]
  variable <- Sdata@meta.data[,c("nFeature_RNA","nCount_RNA","percent.mt")]
  pdf(paste0(sampleoutpath1,"multplot_",colnames(variable)[f],".pdf"),width =10+max(as.numeric(levels(Idents(Sdata))))/4,height = 6+max(as.numeric(levels(Idents(Sdata))))/4)
  grid.newpage()
  n=max(as.numeric(levels(Idents(Sdata))))%/%5+1
  pushViewport(viewport(layout = grid.layout(n,5)))
  vplayout=function(x,y){
    viewport(layout.pos.row = x,layout.pos.col = y)
  }
  
  for(g in 0:max(as.numeric(levels(Idents(Sdata))))){
    r1=g%/%5+1
    c1=g%%5+1
    # p <- ggplot(diamonds[1:1000,],aes(x=price,y=carat)) + geom_point()+labs(title = g)
    useData=Sdata@meta.data[Sdata@meta.data$seurat_clusters==g,]
    p <- ggplot(useData, aes(x=variable[,f][Sdata@meta.data$seurat_clusters==g])) +
      geom_histogram(aes(y=..density..),colour="black", fill="white") +
      geom_vline(aes(xintercept=c(100)),   # Ignore NA values for mean
                 color="red", linetype="dashed", size=0.5) +
      geom_vline(aes(xintercept=c(200)),   # Ignore NA values for mean
                 color="green", linetype="dashed", size=0.5) +
      geom_vline(aes(xintercept=c(300)),   # Ignore NA values for mean
                 color="blue", linetype="dashed", size=0.5) +
      geom_vline(aes(xintercept=c(50)),   # Ignore NA values for mean
                 color="orange", linetype="dashed", size=0.5) +
      geom_vline(aes(xintercept=quantile(variable[,f][Sdata@meta.data$seurat_clusters==g],probs = c(0.99))),   # Ignore NA values for mean
                 color="red", linetype="dashed", size=0.5) +
      geom_density(alpha=.2, fill="#FF6666") + # 重叠部分采用透明设置
      labs(x=colnames(variable)[f],title = paste0("Cluster",g))
    print(p,vp=vplayout(r1,c1))
  } 
  dev.off()
}
###### remove 上限P99
for(f in 1:3){
  variable <- Sdata@meta.data[,c("nFeature_RNA","nCount_RNA","percent.mt")]
  pdf(paste0(sampleoutpath1,"multplot_removeP99_",colnames(variable)[f],".pdf"),width =6+max(as.numeric(levels(Idents(Sdata))))/4,height = 6+max(as.numeric(levels(Idents(Sdata))))/4)
  grid.newpage()
  n=max(as.numeric(levels(Idents(Sdata))))%/%5+1
  pushViewport(viewport(layout = grid.layout(n,5)))
  vplayout=function(x,y){
    viewport(layout.pos.row = x,layout.pos.col = y)
  }
  
  for(g in 0:max(as.numeric(levels(Idents(Sdata))))){
    r1=g%/%5+1
    c1=g%%5+1
    # p <- ggplot(diamonds[1:1000,],aes(x=price,y=carat)) + geom_point()+labs(title = g)
    useData=Sdata@meta.data[Sdata@meta.data$seurat_clusters==g,]
    p99=quantile(variable[,f][Sdata@meta.data$seurat_clusters==g],probs = c(0.99))
    vx = variable[,f][Sdata@meta.data$seurat_clusters==g]
    vx = subset(vx,vx < p99)
    p <- ggplot(NULL, aes(x=vx)) +
      geom_histogram(aes(y=..density..),colour="black", fill="white") +
      geom_vline(aes(xintercept=c(100)),   # Ignore NA values for mean
                 color="red", linetype="dashed", size=0.5) +
      geom_vline(aes(xintercept=c(200)),   # Ignore NA values for mean
                 color="green", linetype="dashed", size=0.5) +
      geom_vline(aes(xintercept=c(300)),   # Ignore NA values for mean
                 color="blue", linetype="dashed", size=0.5) +
      geom_vline(aes(xintercept=c(50)),   # Ignore NA values for mean
                 color="orange", linetype="dashed", size=0.5) +
      # geom_vline(aes(xintercept=quantile(variable[,f][Sdata@meta.data$seurat_clusters==g],probs = c(0.99))),   # Ignore NA values for mean
      #            color="red", linetype="dashed", size=0.5) +
      geom_density(alpha=.2, fill="#FF6666") + # 重叠部分采用透明设置
      labs(x=colnames(variable)[f],title = paste0("Cluster",g))
    print(p,vp=vplayout(r1,c1))
  } 
  dev.off()
}

for(i in 1:3){
  for(j in 1:3){
    if(i==1 & j==3){
      sample = paste0(sampleset1[i],sampleset2[j],"_100genes_2")
    }else{
      sample = paste0(sampleset1[i],sampleset2[j],"_100genes")
    }
    sampleoutpath <- paste0(anaPath,sampleset1[i],"/",paste0(sampleset1[i],sampleset2[j]),"/test_20genes/")
    load(paste0(sampleoutpath,sample,".Rda"))
    ## 写出expression matrix，采用scrublet find doublet cells
    sample = paste0(sampleset1[i],sampleset2[j])
    write.csv(as.data.frame(t(Sdata[["RNA"]]@counts)),file = paste0(sampleoutpath,sample,"_matrix.csv"),row.names = T)
  }
}

#################### DotPlot ###########
Idents(Sdata)="celltype"
df1=df[df$cellType=="Muc2+ Cyp2g1+ Olfactory Epithelial Cells",]
top10 =df1 %>% top_n(n=25, wt=avg_logFC)
DotPlot(Sdata,features = top10$gene,cols = c("grey","red"))+RotatedAxis()
DotPlot(Sdata,features = df1$gene[1:15],cols = c("grey","red"))+RotatedAxis()

DotPlot(Sdata,features = c("Epcam","Ptprc","Muc5b","Bpifa1","Bpifb2","Rgs5", "Prss33", "Crabp1",
                           "Scgb3a2","Scgb1a1","Cxcl17","Cxcl1","Sftpc","Scgb1c1","Bpifb4",
                           "Sox9","Icam1","Krt5","Trp63","Ascl1","Notch1","Hes1","Hes6","Kit","Omp","Cftr",
                           "Gabrp", "Steap4","Msln",
                           "Ace2","Tmprss2","Furin","Stat1",
                           "Mx1","Isg15"),cols = c("lightgrey","red"))+RotatedAxis()

Genes=DFgenes$gene[DFgenes$celltype=="Ionocyte 1" & DFgenes$pct.2<0.02]
DotPlot(Sdata,features =Genes, cols = c("grey","red"),group.by="celltype")+RotatedAxis()

DotPlot(Sdata,features = c("Atat1","Lipf","Tff2","Muc5ac","Ace2","Tmprss2","Stat1","Furin",
                           "Gp2","Marcks","Spib","Tnfsf11","Cyp4a12a","Adgrg6","Trpm5"),
        cols = c("white","red"),group.by="celltype")+RotatedAxis()

FeaturePlot(Sdata,features = c("Elane","Ighd","Vpreb1","Csf1r","Camp","Cxcr2"))

Genes=c("Epcam","Gfi1b","Gp2","Gjb2", "Ubd", "Ctsh", "Anxa5", "Ccl20", "Tnfaip2", "Msln", 
        "Il4i1", "Adgrd1", "Ctsd", "Aif1", "H2-M2", "Serpinb6a", 
        "Ccl6", "Ccl9", "Marcksl1", "Pglyrp1", "Rac2", "Spib",
        "Sox8", "Mier3", "Mycn", "Tcf15", "Etohi1", "Tcf15", 
        # "Atat1","Lipf","Tff2","Muc5ac","Ace2","Tmprss2","Stat1","Furin",
        # "Gp2","Marcks","Spib","Tnfsf11","Cyp4a12a","Adgrg6","Trpm5"
        "Scgb1a1","Krt13","Cftr",
        "Scgb1c1", "Bpifb4", "Furin","Cyp2a5","Rgs5","Svopl","Crlf1","Prss33","Hes6","Dclk1","Il25",
        "Lipf","Dcpp3","Dcpp1","Muc5ac","Tff2",
        "Scgb1a1", "Scgb3a1", "Il13ra1", "Reg3g", "Lgr6", "Bpifb1","Ascl1","Ngn2","Cyp2f2","Calca"
)
Genes=intersect(Genes,rownames(Sdata[["RNA"]]@data))
DotPlot(Epi,features = Genes,
        cols = c("white","red"))+RotatedAxis()

Genes <- readLines("./mGenesDot_Genes.txt")##cell type Markers
Genes <- readLines("./featureplot_Genes.txt")
top10 =DFgenes %>% group_by(newcluster) %>% top_n(n=5, wt=avg_logFC)
Genes <- intersect(top10$gene,rownames(Sdata[["RNA"]]@counts))

Genes <- readLines("./wholeMarkers_20190727.txt")##wholeMarkers
Genes <- capitalize(tolower(Genes))

celltypeMarkers=readLines("./celltypeMarker_selected.txt")
Genes <- intersect(celltypeMarkers,rownames(Sdata[["RNA"]]@counts))

Genes = c("Ace2","Tmprss2","Furin","Stat1")

Genes <- intersect(Genes,rownames(Sdata[["RNA"]]@counts))
dotplot_outpath <- paste0(sampleoutpath,sample,"_Dotplot/")
if(!dir.exists(dotplot_outpath)){
  dir.create(dotplot_outpath)
}
# pdf(paste0(dotplot_outpath, sample,"_celltypeGenes",".pdf"),
#     w=15+ceiling(length(Genes)/4),h=2+ceiling(length(Genes)/8),useDingbats = T)
# pdf(paste0(dotplot_outpath, sample,"_celltypeGenes2",".pdf"),
#     w=6+ceiling(length(Genes)/4),h=4+ceiling(length(Genes)/4),useDingbats = T)
# pdf(paste0(dotplot_outpath, sample,"_wholeGenes",".pdf"),
#     w=15+ceiling(length(Genes)/4),h=2+ceiling(length(Genes)/6),useDingbats = T)

pdf(paste0(dotplot_outpath, sample,"_res1.2_celltypeMarkers_selected",".pdf"),
    w=18+ceiling(length(Genes)/4),h=6+ceiling(length(Genes)/6),useDingbats = T)

# pdf(paste0(dotplot_outpath, sample,"_top5DFGenes_celltype_reorder",".pdf"),
#     w=18+ceiling(length(Genes)/4),h=2+ceiling(length(Genes)/8),useDingbats = T)
# pdf(paste0(dotplot_outpath, sample,"_covGenes_celltype_reorder",".pdf"),
#     w=9,h=10,useDingbats = T)
marker_dot_plot <- DotPlot(object = Sdata, features = Genes,#group.by = "celltype",
                           cols = c("grey","red"),col.min = 0, dot.scale = 8) 
marker_dot_plot <- marker_dot_plot + theme_bw() +
  theme(panel.grid = element_line(linetype = "dashed",color = "lightgrey"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text=element_text(size=16)) 
print(marker_dot_plot)
dev.off()

###### for each cluster
dotplot_outpath <- paste0(sampleoutpath,sample,"_Dotplot_cluster/")
if(!dir.exists(dotplot_outpath)){dir.create(dotplot_outpath)}
for(i in 0:(length(unique(Sdata$seurat_clusters))-1)){
  # i=0
  df = DFgenes[DFgenes$cluster==i,]
  Genes = df$gene[1:50]
  pdf(paste0(dotplot_outpath, "Cluster",i,"_top50DFGenes",".pdf"),
      w=6+ceiling(length(Genes)/4),h=2+ceiling(length(Genes)/4),useDingbats = T)
  marker_dot_plot <- DotPlot(object = Sdata, features = Genes,
                             cols = c("lightgrey","red"),col.min = 0, dot.scale = 8)
  marker_dot_plot <- marker_dot_plot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                             axis.text=element_text(size=20))
  print(marker_dot_plot)
  dev.off()
}

#
Genes = c("Gng2","5430417L22Rik","Tsnax","Stbd1","Trpc2","Gm7534","Vil1","Clptm1", "Abhd2",
          "Calr4","Atp2b2","Pla2g16","Tmprss6","Ankrd63") #Vomeronasal Receptor Expressing Neurons 
Genes = c("Krt5","Aqp3","Ltf","Azgp1","Pifo","Muc5ac","Lypd2","Jak1","Cdhr3","Ace2","Tmprss2",
          "Npy","Gria2","Arx","Fabp7", "Apoe", "Atp1a2")
DotPlot(Sdata,features = Genes,cols = c("grey","red"))+RotatedAxis()
#