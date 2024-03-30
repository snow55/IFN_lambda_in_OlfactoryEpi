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
mycol_41 = c('#fdc086','#f2991a',  '#cc051a', "#e986b7", '#d94c1a',
             '#e6cc4c', '#b2d199','#1a734c', '#72daf2', '#3387b5',
             '#b7a39f', "#d8e755", "#c5aac3","#ea7425", '#dbb972',
             "#3abca5","#b49674","#e2a226",'#76c6ba',"#e58f70",
             "#648bc9", "#368b9d","#e76da2",'#7fc97f',"#be71af",
             "#d7ce96","#be2583","#42a16c","#9b9a9b",  "#fdcf16",
             "#bf6067", "#d84b71", "#f6999a", "#e57371", "#9452a0",
             "#2fa147",'#beaed4',"#1e77fe","#a6cfe4",  '#cf6eb8',
             "#EC188B","#1a79b5", "#fdbf6d", "#b3d78b")
             #,, "#d4af29","#cab2d6","#ff00ff"
show_col(mycol_41)
mycol_49=c("#e986b7","#e76da2","#ff00ff","#bf6067", "#f6999a","#d84b71","#EC188B", "#e2a226", #basal cell 8
          "#f9a270","#f7aa5d","#d66551","#e86502","#e58f70","#ea7425", "#ed2b2b",'#d94c1a', '#cc051a', #Neuron 9
          '#f2991a',"#c5aac3",'#cf6eb8' , #VSN/Olig/Astrocyte 3
          "#fdcf16","#d8e755",'#b2d199', #SUS 3
          '#76c6ba',"#42a16c", #Bowman gland cell 2
          "#648bc9", "#368b9d","#be71af",
          "#d7ce96","#b49674", #Ionocyte 2
          '#1a734c', '#3387b5',"#9452a0",'#beaed4',"#1e77fe",'#72daf2',"#1a79b5", #special Epithelial cell 7
          "#d4af29",'#dbb972',#Dental epithelial cell 2
          "#e57371","#be2583", #Erythroblast 2
          "#b3d78b",'#7fc97f',"#3abca5","#2fa147", #Fibroblast 4
          "#9b9a9b", #Endothelial 1
          "#cab2d6", "#a6cfe4",'#b7a39f'#Immune cell 3
)
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
library(scales)
show_col(mycol_49)

set.seed(2022)
setwd("/home/xfwang/data/projects/Cov19/analysis_new/")

##------------------- analysis from cellranger filtered matrix -------------

filterMatrixClusterFun = function(i=NULL,sampleset=NULL,nfeatures=3000,resolution=0.8){
  sample = sampleset[i]
  sampleoutpath=paste0(sampleset[i],"/filteredMatrix/")
  if(!dir.exists(sampleoutpath)){dir.create(sampleoutpath,recursive = T)}
  
  inputdir=paste0("../cellranger_results/",sampleset[i],"_CD45neg_2_out/outs/filtered_feature_bc_matrix/")
  ftmat = Read10X(data.dir = inputdir)
  colnames(ftmat) <- paste0(sample,"_",colnames(ftmat))
  
  Sdata<- CreateSeuratObject(ftmat,project = sample, min.cells = 0, min.features = 0)
  dim(Sdata[["RNA"]]@counts)
  Sdata[["percent.mt"]] <- PercentageFeatureSet(Sdata, pattern = "^mt-")
  head(Sdata@meta.data, 5)
  
  Sdata <- NormalizeData(Sdata, normalization.method = "LogNormalize", scale.factor = 10000)
  Sdata <- FindVariableFeatures(Sdata, selection.method = "vst", nfeatures = nfeatures)
  top10 <- head(VariableFeatures(Sdata), 10)
  pdf(paste0(sampleoutpath,"VarFeatureplot_filtered_",sample,".pdf"),width = 8,height = 4)
  plot1 <- VariableFeaturePlot(Sdata)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(CombinePlots(plots = list(plot1, plot2)))
  dev.off()
  all.genes <- rownames(Sdata)
  Sdata <- ScaleData(Sdata, features = all.genes)
  Sdata <- RunPCA(Sdata, features = VariableFeatures(object = Sdata))
  # print(Sdata[["pca"]], dims = 1:5, nfeatures = 5)
  
  PCAplotFun = function(sample){
    pdf(paste0(sampleoutpath,"PCAplot_filtered_",sample,"_0.pdf"),width = 5.5,height = 5)
    p1 <- VizDimLoadings(Sdata, dims = 1:2, reduction = "pca")
    p2 <- DimPlot(Sdata, reduction = "pca", dims = c(1, 2))
    p3 <- DimPlot(Sdata, reduction = "pca", dims = c(2, 3))
    p4 <- DimHeatmap(Sdata, dims = 1, cells = 500, balanced = TRUE)
    # DimHeatmap(Sdata, dims = 1:15, cells = 500, balanced = TRUE)
    p5 <- ElbowPlot(Sdata,ndims = 50)
    print(list(p1, p2,p3,p4,p5))
    dev.off()
  }
  PCAplotFun(sample)
  
  Sdata <- FindNeighbors(Sdata, dims = 1:50)
  Sdata <- FindClusters(Sdata, resolution = resolution)
  
  Sdata <- RunTSNE(Sdata, dims = 1:50,check_duplicates = FALSE)
  pdf(paste0(sampleoutpath,"tSNEplot_filtered_",sample,".pdf"),width = 5,height = 5)
  p <- DimPlot(Sdata, reduction = "tsne",label = TRUE,repel = TRUE)+ NoLegend()
  print(p)
  dev.off()
  save(Sdata,file = paste0(sampleoutpath,sample,".Rda"))
  Sdata <- RunUMAP(Sdata, reduction = "pca", dims = 1:50)
  pdf(paste0(sampleoutpath,"UMAPplot_filtered_batch_",sample,".pdf"),width = 5,height = 5)
  p <- DimPlot(Sdata, reduction = "umap",label = TRUE,repel = TRUE) + NoLegend()
  print(p)
  dev.off()
  
  save(Sdata,file = paste0(sampleoutpath,sample,".Rda"))
}
for(i in 1:3){
  filterMatrixClusterFun(i=i,sampleset=c("Alpha","Lambda","Mock"),resolution=3)
}

i=2
sampleset = c("Alpha","Lambda","Mock")
sample=sampleset[i]
sampleoutpath = paste0(sampleset[i],"/filteredMatrix/")
load(paste0(sampleset[i],"/filteredMatrix/",sampleset[i],".Rda"))

###################### The 1th Step：过滤掉gene<100,percent.mt>50的细胞，剩余细胞数
## Alpha: 11136 cells #remove cluster 3,9,40 with low gene number and high mt.pct
## Lambda: 10776 cells #remove cluster 2,14 with low gene number and high mt.pct
## Mock: 10404 cells #remove cluster 4,8,23,31,5 with low gene number and high mt.pct
ftcells = rownames(Sdata@meta.data)[(!Sdata$seurat_clusters %in% c(2,14)) & 
                                      Sdata$nFeature_RNA>100 & Sdata$percent.mt<50]

UMAPPlot(Sdata,label=T,repel=T,cells.highlight=ftcells,sizes.highlight =0.1)+NoLegend()
ggsave(paste0(sampleoutpath,sample,"_filterGreycells_100gene50mt.pdf"),w=5,h=5)

UMAPPlot(Sdata,label=T,repel=T,cells.highlight=rownames(Sdata@meta.data)[Sdata$seurat_clusters==54],sizes.highlight =0.1)+NoLegend()

## 写出expression matrix，采用scrublet find doublet cells
write.csv(t(as.data.frame(Sdata[["RNA"]]@counts[,ftcells])),file = paste0(sampleoutpath,sample,"_filteredMatrix.csv"),row.names = T)


###################### The 2th Step：用Scrublet过滤doublet cells
# doublet cell:
## Alpha: 203 cells
## Lambda: 0 cells
## Mock: 0 cells
# outliers: Mock: ncounts>90000
outliers=rownames(Sdata@meta.data)[Sdata$nCount_RNA>90000]

#Alpha:有doublets,无outliers
i=1
sampleset = c("Alpha","Lambda","Mock")
sample=sampleset[i]
doubCells = readLines(paste0("/data/jinwf/wangxf/jupyter_notebook/Python/Scrublet/Cov19/",sample,"/",sample,"_double_cell_name"))
keepCells = setdiff(ftcells,doubCells)
writeLines(keepCells,paste0(sampleoutpath,sample,"_filtered_pureCells.txt"))

#Lambda:无doublets and outliers
keepCells=ftcells
writeLines(ftcells,paste0(sampleoutpath,sample,"_filtered_pureCells.txt"))

#Mock:无doublets,有outliers
keepCells=setdiff(ftcells,outliers)
writeLines(keepCells,paste0(sampleoutpath,sample,"_filtered_pureCells.txt"))

subData = subset(Sdata,cells = keepCells)#ftcells for lambda
subData_Alpha = subData #10933
subData_Lambda = subData #10776
subData_Mock = subData #10401

###################### The 3th Step：合并三个matrix
mergeData = merge(x=subData_Alpha,y=c(subData_Lambda,subData_Mock))
write.csv(mergeData[["RNA"]]@counts,file=paste0(sampleoutpath,sample,"_countsMatrix.csv"))

mergeWorkFlowFun = function(Sdata=NULL,sampleoutpath=NULL,sample=NULL,
                            nfeatures = 3000,resolution = 0.5,dims = 1:50,
                            min.cells = 3, min.features = 0){
  if(!dir.exists(sampleoutpath)){dir.create(sampleoutpath,recursive = T)}

  Sdata <- CreateSeuratObject(Sdata[["RNA"]]@counts,project = sample, min.cells = min.cells, min.features = min.features)
  Sdata[["percent.mt"]] <- PercentageFeatureSet(Sdata, pattern = "^mt-")
  Sdata <- NormalizeData(Sdata, normalization.method = "LogNormalize", scale.factor = 10000)
  Sdata <- FindVariableFeatures(Sdata, selection.method = "vst", nfeatures = nfeatures)
  top10 <- head(VariableFeatures(Sdata), 10)
  pdf(paste0(sampleoutpath,"VarFeatureplot_filtered_",sample,".pdf"),width = 12,height = 6)
  plot1 <- VariableFeaturePlot(Sdata)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(CombinePlots(plots = list(plot1, plot2)))
  dev.off()
  all.genes <- rownames(Sdata)
  Sdata <- ScaleData(Sdata, features = all.genes)
  Sdata <- RunPCA(Sdata, features = VariableFeatures(object = Sdata))
  # print(Sdata[["pca"]], dims = 1:5, nfeatures = 5)
  
  PCAplotFun = function(sample){
    pdf(paste0(sampleoutpath,"PCAplot_filtered_",sample,".pdf"),width = 5.5,height = 5)
    p1 <- VizDimLoadings(Sdata, dims = 1:2, reduction = "pca")
    p2 <- DimPlot(Sdata, reduction = "pca", dims = c(1, 2))
    p3 <- DimPlot(Sdata, reduction = "pca", dims = c(2, 3))
    p4 <- DimHeatmap(Sdata, dims = 1, cells = 500, balanced = TRUE)
    # DimHeatmap(Sdata, dims = 1:15, cells = 500, balanced = TRUE)
    p5 <- ElbowPlot(Sdata,ndims = 50)
    print(list(p1, p2,p3,p4,p5))
    dev.off()
  }
  PCAplotFun(sample)
  
  Sdata <- FindNeighbors(Sdata, dims = 1:50)
  Sdata <- FindClusters(Sdata, resolution = resolution)
  # head(Idents(Sdata), 5)
  Sdata <- RunUMAP(Sdata, dims = 1:50)
  pdf(paste0(sampleoutpath,"UMAPplot_filtered_batch_",sample,".pdf"),width = 6,height = 5)
  p1 <- DimPlot(Sdata, reduction = "umap",group.by = "orig.ident")
  p2 <- DimPlot(Sdata, reduction = "umap",label = TRUE,repel = TRUE) + NoLegend()
  print(list(p1,p2))
  dev.off()
  
  Sdata <- RunTSNE(Sdata, dims = 1:50)
  pdf(paste0(sampleoutpath,"tSNEplot_filtered_batch_",sample,".pdf"),width = 5.5,height = 5)
  p1 <- DimPlot(Sdata, reduction = "tsne",group.by = "orig.ident")
  p2 <- DimPlot(Sdata, reduction = "tsne",label = TRUE,repel = TRUE) + NoLegend()
  print(list(p1,p2))
  dev.off()
  
  pdf(paste0(sampleoutpath,"PCAplot_filtered_batch_",sample,".pdf"),width = 5.5,height = 5)
  p1 <- DimPlot(Sdata, reduction = "pca", dims = c(1, 2),group.by = "orig.ident")
  p2 <- DimPlot(Sdata, reduction = "pca", dims = c(1, 2))
  print(list(p1,p2))
  dev.off()
  
  pdf(paste0(sampleoutpath,"histgram_",sample,".pdf"))
  p1 <- hist(Sdata@meta.data$percent.mt,breaks = 50,xlab = "percent.mito",main = "Histgram of percent.mito",col = "cornsilk")
  p2 <- hist(log10(Sdata@meta.data$nCount_RNA),breaks = 50,xlab = "nCount_RNA",main = "Histgram of log10(nCount_RNA)",col = "cornsilk")
  p3 <- hist(Sdata@meta.data$nFeature_RNA,breaks = 50,xlab = "nFeature_RNA",main = "Histgram of nFeature_RNA",col = "cornsilk")
  print(list(p1,p2,p3))
  dev.off()
  
  pdf(paste0(sampleoutpath,"Geneplot_filtered_",sample,".pdf"),width = 12,height = 6)
  plot1 <- FeatureScatter(Sdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Sdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(CombinePlots(plots = list(plot1, plot2)))
  dev.off()
  pdf(paste0(sampleoutpath,"Geneplot_filtered_batch_",sample,".pdf"),width = 12,height = 6)
  plot1 <- FeatureScatter(Sdata, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
  plot2 <- FeatureScatter(Sdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")
  print(CombinePlots(plots = list(plot1, plot2)))
  dev.off()
  
  save(Sdata,file=paste0(sampleoutpath,sample,".Rda"))
}
mergeWorkFlowFun(Sdata=mergeData,sample="mergedSamples",
                 sampleoutpath = "./mergedSamples/",
                 nfeatures = 3000,resolution = 1.2)

sample="mergedSamples"
sampleoutpath = "./mergedSamples/"
load(paste0(sampleoutpath,sample,".Rda"))

Sdata = FindClusters(Sdata,resolution = 3.5)
UMAPPlot(Sdata,label=T,group.by="RNA_snn_res.3.5",cols=rep(mycol_41,2))+NoLegend()

##------------------- The 4th Step：从整体层面-人工识别和去除doublet cell -------------
##### remove doublet cells since co-expression marker genes of different cell types
#remove cluster 38,56,59,64,69,71
cellRemove=rownames(Sdata@meta.data)[Sdata$seurat_clusters %in% c(38,56,59,64,69,71)]#resolution=3.5
UMAPPlot(Sdata,label=T,repel=T,cells.highlight=cellRemove,sizes.highlight =0.1)+NoLegend()

# DFgenes_59=FindMarkers(Sdata,ident.1 = "59",only.pos = T)
# DFgenes_59$gene=rownames(DFgenes_59)
# DFgenes_38=FindMarkers(Sdata,ident.1 = "38",only.pos = T)
# DFgenes_38$gene=rownames(DFgenes_38)
# save(DFgenes_59,file = paste0(sampleoutpath,"DFgenes_cluster59_resolution3.5.csv"))
# save(DFgenes_38,file = paste0(sampleoutpath,"DFgenes_cluster38_resolution3.5.csv"))

cluster62=rownames(Sdata@meta.data)[Sdata$seurat_clusters==62]#resolution=3.5
cluster47=rownames(Sdata@meta.data)[Sdata$seurat_clusters==47]#resolution=3.5
writeLines(cluster62,paste0(sampleoutpath,sample,"_CellsInCluster62_res.3.5_beforeCellremove.txt"))

top10=DFgenes_41 %>% top_n(20,avg_logFC)
DotPlot(Sdata,features = rownames(top10)[1:20],cols = c("lightgrey","red"))+RotatedAxis()

############### 去除后重新聚类 
subData=subset(Sdata,cells = cellRemove,invert=T)#22902 genes x 31361 cells
mergeWorkFlowFun(Sdata=subData,sample="mergedSamples",sampleoutpath = "./mergedSamples_pure/",
                 nfeatures = 3000,resolution = 1.2)

sample="mergedSamples"
sampleoutpath = "./mergedSamples_pure/"
load(paste0(sampleoutpath,sample,".Rda"))

Sdata = FindClusters(Sdata,resolution = 1)#final use:resolution=1
# Sdata = RunUMAP(Sdata,dims=1:50,n.neighbors = 25,min.dist = 0.3)#final use: n.neighbors = 25,min.dist = 0.3 - 默认参数
UMAPPlot(Sdata,label=T,cols=rep(mycol_41,2))+NoLegend()
ggsave(filename = paste0(sampleoutpath,"UMAPplot_cluster.pdf"),w=5,h=4.8)
# save(Sdata,file = paste0(sampleoutpath,sample,".Rda"))

clustree(Sdata@meta.data,prefix = "RNA_snn_res.")

############## redefine cluster 
# add a new cluster to separate cluster30 which distributed in two places
cluster30=rownames(Sdata@meta.data)[Sdata$seurat_clusters==30]
cellUse=setdiff(cluster30,cluster47)
UMAPPlot(Sdata,label=T,repel=T,cells.highlight=cellUse,sizes.highlight =0.1)+NoLegend()
FeaturePlot(Sdata,features = c("Krt5","Ascl1","Notch1"))

Sdata$redefineCluster=as.character(Sdata$seurat_clusters)
Sdata$redefineCluster[cellUse]="40"
Sdata$redefineCluster=factor(Sdata$redefineCluster,levels = seq(0,40))
UMAPPlot(Sdata,label=T,group.by="redefineCluster",cols=mycol_41)+NoLegend()
table(Sdata$orig.ident,Sdata$redefineCluster)
Idents(Sdata)="redefineCluster"

############# reorder and define cell types for merged data
celltypeFile = read.csv("./reorder_celltype.csv",header = T)
Sdata$newcluster=plyr::mapvalues(Sdata$redefineCluster,from=as.character(celltypeFile$oldcluster),to=as.character(celltypeFile$newcluster))
Sdata$newcluster=factor(Sdata$newcluster,levels = as.character(celltypeFile$newcluster))
Sdata$celltype=plyr::mapvalues(Sdata$newcluster,from = as.character(celltypeFile$newcluster),to = as.character(celltypeFile$celltype))
Sdata$Tcelltype=plyr::mapvalues(Sdata$newcluster,from = as.character(celltypeFile$newcluster),to = as.character(celltypeFile$Tcelltype))
Idents(Sdata)="newcluster"
Sdata$celltypeColor=plyr::mapvalues(Sdata$newcluster,from = as.character(celltypeFile$newcluster),to = mycol_41)
Sdata$orig.ident=factor(Sdata$orig.ident,levels=c("Mock","Alpha","Lambda"))
Sdata$tissueColor=plyr::mapvalues(Sdata$orig.ident,from = levels(Sdata$orig.ident),to = tissueCols)

celltypeFile$color=mycol_41
write.csv(celltypeFile,file = paste0("./","reorder_celltype_1.csv"),row.names = F)
save(Sdata,file = paste0(sampleoutpath,sample,".Rda"))

FeaturePlot(Sdata,features = c("Epcam","Ptprc"))

UMAPPlot(Sdata,label=F,repel=T,cols=levels(Sdata$celltypeColor),group.by="Tcelltype")+NoLegend()
TSNEPlot(Sdata,label=T,repel=T,cols=mycol_41)+NoLegend()

##------------------- The 5th Step：分别取出不同模块的细胞进行精细化-人工识别和去除doublet cell -------------
############### subset Epithelial cells ##########
##--- 纯化Epithelial cell ---
# 1.取出所有Epithial cell：cluster8:20，重聚类
subData=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$Tcelltype %in% 
                                                             c(levels(Sdata$Tcelltype)[8:20]))])
mergeWorkFlowFun(Sdata=subData,sample="Epithelial",
                 sampleoutpath = "./mergedSamples_pure/Epithelial/",
                 nfeatures = 3000,resolution = 1.2)

sample="Epithelial"
sampleoutpath = "./mergedSamples_pure/Epithelial/"
load(paste0(sampleoutpath,sample,".Rda"))

EpiAll=Sdata
# 2.第一次过滤：去除res=2时，cluster 7,9,27,39,47 - 表达神经元相关的基因
subOSN=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$RNA_snn_res.2 %in% c(7,9,27,39,47))])
subEpi=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$RNA_snn_res.2 %in% c(7,9,27,39,47))],invert=T)

mergeWorkFlowFun(Sdata=subEpi,sample="Epithelial_pure",
                 sampleoutpath = "./mergedSamples_pure/Epithelial/Epithelial_pure/",
                 nfeatures = 3000,resolution = 1.2)

sample="Epithelial_pure"
sampleoutpath = "./mergedSamples_pure/Epithelial/Epithelial_pure/"
load(paste0(sampleoutpath,sample,".Rda"))

EpiAll_2=Sdata
# 2.第二次过滤：去除res=3时，cluster 7,9,27,39,47 - 表达其他细胞的基因
subEpi=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$RNA_snn_res.3 %in% c(52,53,41))],invert=T)
mergeWorkFlowFun(Sdata=subEpi,sample="Epithelial_pure",
                 sampleoutpath = "./mergedSamples_pure/Epithelial/Epithelial_pure_2/",
                 nfeatures = 3000,resolution = 1.2)

sample="Epithelial_pure"
sampleoutpath = "./mergedSamples_pure/Epithelial/Epithelial_pure_2/"
load(paste0(sampleoutpath,sample,".Rda"))

Sdata = FindClusters(Sdata,resolution = 1.5)
UMAPPlot(Sdata,label=T,group.by="RNA_snn_res.1.5",cols=rep(mycol_41,2))+NoLegend()
UMAPPlot(Sdata,cells.highlight=rownames(Sdata@meta.data)[which(Sdata@active.ident=="8")])

EpiAll_3=Sdata

Sdata$celltype_0=subData$celltype[rownames(Sdata@meta.data)]

genes=intersect(c("Omp","Hes6","Neurog1","Neurod1","Sox2","Tp63", #neuron
                  "Cdhr3","Foxj1","Ccdc153","Ccdc113","Gap2","Pifo",#Ciliated cell
                  "Scgb1a1","Krt15","Cyp2f2","Lypd2","Cbr2",#Club cell
                  "Ascl3","Foxi1","Cftr","Coch", #Ionocyte
                  "Ascl2","Rgs13","Pou2f3","Trpm5","Lrmp", #Tuft cell
                  "Krt4","Trp63","Krt5","Krt8","Bcam","Dcn","Krt14","Col17a1","Igfbp3", # Basal cell
                  "Krt13","Ecm1", # Hillock cell
                  "Gp2","Spib","Sox2","Tff1","Tff2","Tff3","Muc5ac","Lipf","Dcpp1","Dcpp2","Dcpp3", #Goblet cell
                  "Muc5b","Bpifa1","Bpifb2","Scgb1c1","Muc2","Sec14l3","Cyp2g1",
                  "Ltf","Azgp1","Crlf1","Rgs5","Cspg4",
                  "Bpifb4","Bpifb6","Bpifb5","Bpifb9a",
                  "Il25","Tslp",
                  "Gnb3","Gng13", "Fxyd6", "Ovol3","Gnat3",
                  "Alox5ap","Mgst3","Ptgs1","Dclk1","Msln","Gpx2","Mgst2","Pla2g4a","Sdc4","Il13ra1",
                  "Gfi1b","Spib","Sox9","Ascl1","Neurog1","Tuj1","Gap43","Sox8","Mier3"
),rownames(Sdata[["RNA"]]@data))


DotPlot(Sdata,features = genes,cols = c("lightgrey","red"))+RotatedAxis()#,group.by = "celltype"
# Club cell: Nfia, Scgb1a1
# Neuroendocrine cell: Ascl1
# Tuft cell: Ascl2, Rgs13, Pou2f3
# Ionocyte: Ascl3
# Goblet cell: Muc5ac
# Ciliated cell: Cdhr3,  Foxj1
# Basal cell: Krt5, Trp63
# Hillock cell: Krt13, Ecm1
Idents(Sdata)="celltype"
table(Sdata@active.ident,Sdata$celltype)

## markers from paper
DFgenes_NatureNeu <- read.csv("DFgenes_NatureNeu.csv",header=T)
colnames(DFgenes_NatureNeu)=c("gene",
                              "p_val",
                              "avg_logFC",
                              "pct.1",
                              "pct.2",
                              "p_val_adj",
                              "cluster")
DFgenes_NatureNeu=DFgenes_NatureNeu[-1,]
table(DFgenes_NatureNeu$cluster)
ctlist=c("Olfactory Horizontal Basal Cells",
         "Respiratory Horizontal Basal Cells",
         # "Fibroblasts/Stromal Cells",
         # "Pericytes",
         "Sustentacular Cells",
         # "Vascular Smooth Muscle Cells",
         "Respiratory Secretory cells",
         "Bowman's Gland",
         "Respiratory Columnar Cells",
         "Respiratory Ciliated Cells",
         "Immature Neurons",
         "Olfactory Ensheathing Glia",
         "Mature Neurons",
         "Respiratory Epithelial Cells",
         "Olfactory Microvillar Cells",
         "Globose Basal Cells")
DFgenes_NatureNeu=DFgenes_NatureNeu[which(DFgenes_NatureNeu$cluster %in% ctlist),]

top5 <- DFgenes_NatureNeu[DFgenes_NatureNeu$cluster==ctlist[5],] %>% group_by(cluster) %>% top_n(n=30, wt=avg_logFC)
# top5 <- DFgenes_NatureNeu %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
Genes <- intersect(human2mouse(top5$gene)$mouseGene,rownames(Sdata[["RNA"]]@counts))
Genes <- intersect(capitalize(tolower(top5$gene)),rownames(Sdata[["RNA"]]@counts))

DotPlot(Sdata,features = Genes,cols = c("lightgrey","red"))+RotatedAxis()#,group.by = "celltype"

library(homologene)
human2mouse(top5$gene)$mouseGene

Genes = c("Epcam","Krt17","Krt16","Krt6a","Krt6b","Krt15","Krt19","Krt7", #basal cell
          "Krt18","Krt8",
          # "Tgm2","Loricrin","Flg","Ivl",
          "Krt10","Krt1", 
          "Krt5","Krt14","Trp63","Icam1","Sox2", #HBC
          "Kit", "Ascl1","Cxcl14","Neurod1","Neurog1","Hes6","Ezh2","Hmga1","Mki67","Stmn1","Top2a",#GBC
          "Adh7", #respiratory basal cell/HBC
          "Ascl3","Cftr", #Ionocyte
          "Trpm5","Sox9", #Tuft cell
          "Cyp2g1","Muc2", "Sec14l3","Cyp1a2","Notch2"#, #SUS
)

Genes = intersect(Genes, rownames(Sdata[["RNA"]]@data))

DotPlot(Sdata,features = Genes,
        cols = c("lightgrey","red"))+RotatedAxis()#,group.by = "celltype"

############### subset Neurons ##########
##--- 纯化Neuron cell ---
# 1.取出所有Neuron cell：cluster1:7，与Epithelial cell中高表达Neuron相关基因的细胞合并后，重聚类
subData=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$Tcelltype %in% 
                                                             c(levels(Sdata$Tcelltype)[1:7]))])
subData=merge(subData,subOSN)#subOSN from Epithelial

mergeWorkFlowFun(Sdata=subData,sample="Neuron",
                 sampleoutpath = "./mergedSamples_pure/Neuron/",
                 nfeatures = 3000,resolution = 1.2)

sample="Neuron"
sampleoutpath = "./mergedSamples_pure/Neuron/"
load(paste0(sampleoutpath,sample,".Rda"))

Sdata = FindClusters(Sdata,resolution = 0.3)
UMAPPlot(Sdata,label=T,group.by="RNA_snn_res.0.3",cols=rep(mycol_41,2))+NoLegend()
UMAPPlot(Sdata,cells.highlight=rownames(Sdata@meta.data)[which(Sdata@active.ident=="8")])


NeuronAll=Sdata
# 2.第一次过滤：去除res=1.2时，cluster 10, 16,19,22, 25,26,27 - 同时表达神经元与其他细胞相关的基因
subData=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$RNA_snn_res.1.2 %in% c(10, 16,19,22, 25,26,27))],invert=T)

mergeWorkFlowFun(Sdata=subData,sample="Neuron_pure",
                 sampleoutpath = "./mergedSamples_pure/Neuron/Neuron_pure/",
                 nfeatures = 3000,resolution = 0.8)
length(intersect(rownames(subOSN@meta.data),rownames(subData@meta.data)))

sample="Neuron_pure"
sampleoutpath = "./mergedSamples_pure/Neuron/Neuron_pure/"
load(paste0(sampleoutpath,sample,".Rda"))

ctlist=c("Olfactory Horizontal Basal Cells",
         "Respiratory Horizontal Basal Cells",
         # "Fibroblasts/Stromal Cells",
         # "Pericytes",
         "Sustentacular Cells",
         # "Vascular Smooth Muscle Cells",
         "Respiratory Secretory cells",
         "Bowman's Gland",
         "Respiratory Columnar Cells",
         "Respiratory Ciliated Cells",
         "Immature Neurons",
         "Olfactory Ensheathing Glia",
         "Mature Neurons",
         "Respiratory Epithelial Cells",
         "Olfactory Microvillar Cells",
         "Globose Basal Cells")
DFgenes_NatureNeu=DFgenes_NatureNeu[which(DFgenes_NatureNeu$cluster %in% ctlist),]

top5 <- DFgenes_NatureNeu[DFgenes_NatureNeu$cluster==ctlist[9],] %>% group_by(cluster) %>% top_n(n=30, wt=avg_logFC)
# top5 <- DFgenes_NatureNeu %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
Genes <- intersect(human2mouse(top5$gene)$mouseGene,rownames(Sdata[["RNA"]]@counts))
Genes <- intersect(capitalize(tolower(top5$gene)),rownames(Sdata[["RNA"]]@counts))

DotPlot(Sdata,features = Genes,cols = c("lightgrey","red"))+RotatedAxis()#,group.by = "celltype"

Genes = c("Epcam","Krt17","Krt16","Krt6a","Krt6b","Krt15","Krt19","Krt7", #basal cell
          "Krt18","Krt8",
          # "Tgm2","Loricrin","Flg","Ivl",
          "Krt10","Krt1", 
          "Krt5","Krt14","Trp63","Icam1","Sox2", #HBC
          "Kit", "Ascl1","Cxcl14","Neurod1","Neurog1","Hes6","Ezh2","Hmga1","Mki67","Stmn1","Top2a",#GBC
          "Adh7", #respiratory basal cell
          "Ascl3","Cftr", #Ionocyte
          "Trpm5","Sox9", #Tuft cell
          "Cyp2g1","Muc2", "Sec14l3","Cyp1a2","Notch2", #SUS
          "Omp","Gng13", #mature OSN
          "Gng8","Gap43","Tuj1", #immature OSN
          "Sox2", "Pax6", "Hes2", "Sox9",
          "Rbm24", "Neurod1", "Elavl4", "Scg2","Lhx2", "Tex15", "Crabp1", #intermediate neural precursor (INP)
          "Calr4","Bmp6", #vomeronasal sensory neurons, VSN
          "Plp1", "Npy", "Pmp22", "Mpz", "Ttyh1", #Olfactory Ensheathing Glia
          "Slc6a11", "Gria2", "Plpp3", "Mt3", "Bcan", #Gria2 Neuron
          "Sox11","Insm1","Cdyl2","Thsd7b"
)
Genes = intersect(Genes, rownames(Sdata[["RNA"]]@data))

DotPlot(Sdata,features = Genes,
        cols = c("lightgrey","red"))+RotatedAxis()

Neuron_pure=Sdata

#

############### subset ImmuneFibroblast cells ############
##--- 纯化ImmuneFibroblast cell ---
# 1.取出所有ImmuneFibroblast cell：cluster21:25，重聚类
subData=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$Tcelltype %in% 
                                                             c(levels(Sdata$Tcelltype)[21:25]))])

mergeWorkFlowFun(Sdata=subData,sample="ImmuneFibroblast",
                 sampleoutpath = "./mergedSamples_pure/ImmuneFibroblast/",
                 nfeatures = 3000,resolution = 0.8)

sample="ImmuneFibroblast"
sampleoutpath = "./mergedSamples_pure/ImmuneFibroblast/"
load(paste0(sampleoutpath,sample,".Rda"))

Sdata = FindClusters(Sdata,resolution = 0.5)
UMAPPlot(Sdata,label=T,group.by="RNA_snn_res.0.5",cols=rep(mycol_41,2))+NoLegend()
UMAPPlot(Sdata,cells.highlight=rownames(Sdata@meta.data)[which(Sdata@active.ident=="8")])


ImmuneFibroblastAll=Sdata
# 2.第一次过滤：去除res=0.8时，cluster 5,16,17 - 同时表达神经元与其他细胞相关的基因
subData=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$RNA_snn_res.0.8 %in% c(5,16,17))],invert=T)

mergeWorkFlowFun(Sdata=subData,sample="ImmuneFibroblast_pure",
                 sampleoutpath = "./mergedSamples_pure/ImmuneFibroblast/ImmuneFibroblast_pure/",
                 nfeatures = 3000,resolution = 0.8)

sample="ImmuneFibroblast_pure"
sampleoutpath = "./mergedSamples_pure/ImmuneFibroblast/ImmuneFibroblast_pure/"
load(paste0(sampleoutpath,sample,".Rda"))


top5 <- DFgenes %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
# top5 <- DFgenes_NatureNeu %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
Genes <- intersect(human2mouse(top5$gene)$mouseGene,rownames(Sdata[["RNA"]]@counts))
Genes <- intersect(capitalize(tolower(top5$gene)),rownames(Sdata[["RNA"]]@counts))


Genes = c("Pf4","Gp9","Ppbp", #megakaryocyte
          "Cdh5", "Kdr","Sema3g","Efnb2","Dll4", #Vascular endothelial cell 
          "Lyve1","Prox1","Tbx1", #Lymphatic endothelial cell
          "Rgs5","Cspg4", #Pericyte
          "Omp","Nrxn1", #Neuron
          "Gypa","Gypb","Tfrc","Alas2","Hbb−bt","Hba−a1","Bpgm", #erythroblast
          "Cpa3","Ms4a2","Fcer1a", #Mast cell
          "Elane","Mpo","Lyz2", # GMP
          "Cd68","Cd14","Csf1r", "Fcgr1", #Monocyte
          "Csf3r","Cxcr2", "Ngp", "S100a8", "S100a9", #Neutrophil
          "Cd79a", "Vpreb3", "Pax5", #B cell  
          "Pdgfra","Dcn","Col1a1", #Fibroblast
          "Cthrc1", "Acta2", "Postn", "Adam12",
          "Col15a1", "Col4a1", "Hspg2",
          "Pi16", "Dpp4", "Ly6h",
          "Gp2","Tnfaip2","Ly6d"
)
DotPlot(Sdata,features = Genes,cols = c("lightgrey","red"))+RotatedAxis()

############### merge pure cell together ######################
setwd("/home/wxf/data/projects/Cov19/analysis_new/")

ctNameList=c("Epithelial","Neuron","ImmuneFibroblast")
# i=1
#Epithelial cell过滤了2次，其他的过滤1次
if(i=1){
  sample=paste0(ctNameList[i],"_pure")
  sampleoutpath = paste0("./mergedSamples_pure/",ctNameList[i],"/",ctNameList[i],"_pure_2/")
  load(paste0(sampleoutpath,sample,".Rda"))
}else{
  sample=paste0(ctNameList[i],"_pure")
  sampleoutpath = paste0("./mergedSamples_pure/",ctNameList[i],"/",ctNameList[i],"_pure/")
  load(paste0(sampleoutpath,sample,".Rda"))
}

subData=merge(x=EpiAll_3,y=c(Neuron_pure,Sdata))

mergeWorkFlowFun(Sdata=subData,sample="mergedSamples_pure",
                 sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_2/",
                 nfeatures = 3000,resolution = 1.5)

sample="mergedSamples_pure"
sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_2/"
load(paste0(sampleoutpath,sample,".Rda"))

Sdata = FindClusters(Sdata,resolution = 1.5)
UMAPPlot(Sdata,label=T,group.by="RNA_snn_res.1.2",cols=rep(mycol_41,3))+NoLegend()
UMAPPlot(Sdata,cells.highlight=rownames(EpiAll_3@meta.data)[which(EpiAll_3$RNA_snn_res.1.5=="37")])


Genes = c(
  # "Pf4","Gp9","Ppbp", #megakaryocyte
  # "Cdh5", "Kdr","Sema3g","Efnb2","Dll4", #Vascular endothelial cell
  # "Lyve1","Prox1","Tbx1", #Lymphatic endothelial cell
  # "Rgs5","Cspg4", #Pericyte
  # "Omp","Nrxn1", #Neuron
  # "Gypa","Gypb","Tfrc","Alas2","Hbb−bt","Hba−a1","Bpgm", #erythroblast
  # "Cpa3","Ms4a2","Fcer1a", #Mast cell
  # "Elane","Mpo","Lyz2", # GMP
  # "Cd68","Cd14","Csf1r", "Fcgr1", #Monocyte
  # "Csf3r","Cxcr2", "Ngp", "S100a8", "S100a9", #Neutrophil
  # "Cd79a", "Vpreb1","Vpreb3", "Pax5","Ighm","Ighd","Sell", #B cell
  # "Pdgfra","Dcn","Col1a1", #Fibroblast
  # "Cthrc1", "Acta2", "Postn", "Adam12",
  # "Col15a1", "Col4a1", "Hspg2",
  # "Pi16", "Dpp4", "Ly6h",
  # "Epcam"
  # ,
  "Epcam","Krt17","Krt16","Krt6a","Krt6b","Krt15","Krt19","Krt7", #basal cell
  "Krt18","Krt8",
  # "Tgm2","Loricrin","Flg","Ivl",
  "Krt10","Krt1","Krt2","Krt13","Krt4",
  "Krt5","Krt14","Trp63","Icam1", "Aqp4", "Ly6d","Sox2","Pax6","Hes1","Sox9", #HBC
  "Adh7", #respiratory basal cell/HBC
  "Sprr1a", #Activated/Cycling HBC
  "Kit", "Ascl1","Cxcl14","Neurod1","Neurog1","Hes6","Ezh2","Hmga1","Mki67","Stmn1","Top2a",#GBC
  "Ascl3","Cftr","Fgf10","Fgfr2b", #Ionocyte
  "Trpm5","Avil","Sh2d7","Lrmp",#Tuft cell
  "Dnah5","Tmem212","Cdhr3", #Ciliated cell
  "Tff2","Muc5ac","Muc5b","Bpifa1", #Goblet cell
  # "Lipf","Dcpp1","Dcpp2","Dcpp3","Notch3",#Goblet cell2
  "Scgb1a1","Scgb3a2","Cyp2f2","Lypd2","Cbr2","Nfia", #Club cell
  "Foxj1", "Pifo", "Ift57", "Lrrc23",
  "Nrcam","Vit", #olfactory HBC
  "Gp2","Tnfaip2","Ly6d", #M cell
  "Gsto1","Csta1","Asprv1","Psca","Prss22", #Hillock cell
  "Sox9","Rgs5","Msln","Vmo1","Chil6","Prss33","Scgb1c1","Ldhd","Gpx6","Bpifb6", #Bowman gland cell
  "Cyp2g1","Muc2", "Sec14l3","Cyp1a2","Notch2","Ermn", "Cyp2a13","Cyp2j2", #SUS
  "Ephx1", #SUS ventral
  "Sult1c1","Aqp5", #SUS dorsal
  "Rbm24", "Neurod1", "Elavl4", "Scg2","Lhx2", "Tex15", "Crabp1", #intermediate neural precursor (INP)
  "Gng8","Gap43","Tuj1", #immature OSN
  "Omp","Gng13","Cd36", "Cnga2", "Adcy3" #mature OSN
  # "Calr4","Bmp6", #vomeronasal sensory neurons, VSN
  # "Plp1", "Npy", "Pmp22", "Mpz","Mog","Mbp","Klk6","Apod",#"Rip","Olig1","Olig2","Itgam","Iba1", #Npy Neuron - Oligodendrocyte
  # "Slc6a11","Ttyh1","Gria2","Mt3", "Bcan","Slc1a3", "Slc4a4","Aqp4","Gjal","Slc1a2","Gfap", #Gria2 Neuron - Astrocyte
  # "Fezf2","Meg3","Car9",#Fezf2+ Epithelial cell
  # "Scgb2b27","Car6","Pax3","Dnase1","Klk14","Obp1a","Obp1b","Obp2a","Obp2b", #Obp1a+Obp1b+ Epithelial cell
  # "Cxcl17","Ccl9","Ltf","Bglap3","Aox3","Cyp4b1","Tac1","Dmbt1","Slc39a8", #Cxcl17+Ccl9+ Epithelial cell
  # "Bpifb9a","Bpifb9b","Bpifb5","Siglcc8", #Bpifb9a+Bpifb5+ Epithelial cell
  # "Fxyd2","Klk1", "Krt7", "Krt19","Tfcp2l1","Scnn1b", #Krt19, Krt7, and Tcfcp2l1 are ductal maturation markers
  # "Pip", # markers for salivary gland maturation
  # "Aqp5","Bhlha15"#,  #are markers for terminal differentiation of SMG epithelial cells (acinar and duct, respectively)
  # "Gstt1", 
  # "Gfra3","Kit",
  # "Etv4", "Fgfr2","Mmp2","Esp1","Esp22"
)
Genes <- intersect(Genes,rownames(Sdata[["RNA"]]@counts))

DotPlot(EpiAll_3,group.by = "RNA_snn_res.1.5",features = Genes,cols = c("lightgrey","red"))+RotatedAxis()
DotPlot(allCell_pure,features = Genes,cols = c("lightgrey","red"))+RotatedAxis()#group.by = "RNA_snn_res.3.5",
DotPlot(Sdata,features = Genes,cols = c("lightgrey","red"))+RotatedAxis()#group.by = "RNA_snn_res.3.5",

UMAPPlot(Sdata,group.by="RNA_snn_res.1.5",label=T,cols=rep(mycol_41,3))+NoLegend()

UMAPPlot(Sdata,cells.highlight=rownames(EpiAll_3@meta.data)[which(EpiAll_3$RNA_snn_res.1.5=="35")])
allCell_pure=Sdata

############### subset special epithelial cell ################
sample="mergedSamples_pure"
sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_2/"
load(paste0(sampleoutpath,sample,".Rda"))

# 1.取出所有special epithelial cell：cluster16,23,26,53,35,45,59,44，重聚类
subData=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata@active.ident %in% c(16,23,26,53,35,45,59,44))])

mergeWorkFlowFun(Sdata=subData,sample="specialEpithelial",
                 sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_2/specialEpithelial/",
                 nfeatures = 3000,resolution = 0.8)

sample="specialEpithelial"
sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_2/specialEpithelial/"
load(paste0(sampleoutpath,sample,".Rda"))

Sdata = FindClusters(Sdata,resolution = 0.5)
UMAPPlot(Sdata,label=T,group.by="RNA_snn_res.0.8",cols=rep(mycol_41,2))+NoLegend()
UMAPPlot(Sdata,cells.highlight=rownames(Sdata@meta.data)[which(Sdata@active.ident=="8")])

DFgenes_Neuron=read.csv("./mergedSamples_pure/mergedSamples_pure_2/specialEpithelial/41593_2021_851_MOESM2_ESM.csv",header = T)
# > unique(DFgenes_Neuron$cluster)
# [1] 01 Dentate Gyrus Granule Cells  02 Dentate Gyrus Granule Cells  03 Dorsal CA1 Pyramids          04 Ventral CA1 Pyramids        
# [5] 05 CA1 Neurons                  06 CA2/CA3 Pyramids             07 CA2/CA3 Neurons              08 Mossy Cells                 
# [9] 09 Subiculum/Entorhinal Neurons 10 PV/SST Interneurons          11 VIP Interneurons             12 RELN Interneurons           
# [13] 13 Cajal Retzius Cells          14 Interneurons                 15 Neurons, Unresolved          16 Neurons, Unresolved         
# [17] 17 Astrocytes                   18 Astrocytes                   19 OPCs                         20 Maturing OPCs               
# [21] 21 Mature Oligodendrocytes      22 Oligodendrocytes             23 Oligodendrocytes             24 Microglia/Macrophage        
# [25] 25 Endothelial                  26 Fibroblast-like              27 Choroid Plexus   
i=1

top5 <- DFgenes_Neuron[DFgenes_Neuron$cluster==unique(DFgenes_Neuron$cluster)[i],] %>% group_by(cluster) %>% top_n(n=50, wt=avg_logFC)#p_val_adj
# top5 <- DFgenes_NatureNeu %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
# Genes <- intersect(human2mouse(top5$gene)$mouseGene,rownames(Sdata[["RNA"]]@counts))
Genes <- intersect(top5$gene,rownames(Sdata[["RNA"]]@counts))


DotPlot(Sdata,features = Genes,cols = c("lightgrey","red"))+RotatedAxis()#,group.by = "celltype"
DotPlot(allCell_pure,features = Genes,cols = c("lightgrey","red"))+RotatedAxis()#group.by = "RNA_snn_res.3.5",

############### subset basal cell ##############
sample="mergedSamples_pure"
sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_2/"
load(paste0(sampleoutpath,sample,".Rda"))

#reorder RNA_snn_res3.5 of allCell_pure
current.id=c(9, 17,60,37,49,15,28,55,47,33,52,
             7,62,42,27,25,54,  10,18,58,  
             63,  43,29,14,61,40,6,8,3,1,2,12,4,5,0,24,50, 
             11,21,13,48, 22, 32,56, 38,66, 44, 16,23, 45,59, 26,53,35, 
             36,46,19, 20,30,34,65,57, 39,41,64,31,51)
current.id=c(11,15,26,39, 13,34, 20,28, 8,47,25,14, 6,44, 33,9,36,7,1,3,5,2,0, 46, 40,22,  4,12, 21, 19, 31, 35, 10,37,45, 23,42,30, 32,38, 16, 17,27,29,48,43, 18,24,41
)
Sdata$seurat_clusters=factor(Sdata$RNA_snn_res.1.5,levels = current.id)

#
# subset basal cell
subData=subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$RNA_snn_res.3.5 %in% c(9, 17,60,37,49,15,28,55,47,33,52,63))])
mergeWorkFlowFun(Sdata=subData,sample="basalCell",
                 sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_2/basalCell/",
                 nfeatures = 3000,resolution = 0.5)

sample="basalCell"
sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_2/basalCell/"
load(paste0(sampleoutpath,sample,".Rda"))

Sdata = FindClusters(Sdata,resolution = 0.8)
UMAPPlot(Sdata,label=T,group.by="RNA_snn_res.0.8",cols=rep(mycol_41,2))+NoLegend()
UMAPPlot(Sdata,cells.highlight=rownames(allCell_pure@meta.data)[which(allCell_pure$RNA_snn_res.3.5=="63")])

# Mcell=rownames(EpiAll_3@meta.data)[EpiAll_3$RNA_snn_res.3==22]
# c63=rownames(allCell_pure@meta.data)[which(allCell_pure$RNA_snn_res.3.5=="63")]

### special basal cell: cluster 10,5,8,11 -- dental epithelial cell
# cluster13 in res0.8-basal cell need to be removed since it has high gene number, both expressed markers of basal cell and mature OSN, and dispersed in many cluster of allCell_pure
cluster13=rownames(Sdata@meta.data)[Sdata$RNA_snn_res.0.8==13]
writeLines(cluster13,paste0(sampleoutpath,"doublet_basalCell_needremove.txt"))
Genes=c("Krt14","Cdh1", #epithelial cell
        "Gli1", "Igfbp5","Sfrp5","Acta2","Shh","Vwa2","Enpp2","Cdh6","Cpne5", "Col22a1","Vwde",#cluster10 Dental stratum intermedium progenitor
        "Gpx3","Cpe","Tfrc","Npl","Pmch",  #cluster 5  Dental stratum intermedium cell
        "Enam","Klk4","Calb1","Gpr155","Gad1","Bmp4","Bmp7",#cluster 8 Ameloblast
        "Slc1a1","Ptchd4","Slc39a2","Slc24a4","Gm17660","Ptpn22","Stmn2" #cluster11  Ameloblast
)
Genes <- intersect(Genes,rownames(Sdata[["RNA"]]@counts))

DotPlot(Sdata,features = Genes,cols = c("lightgrey","red"))+RotatedAxis()#group.by = "RNA_snn_res.3.5",

# dental epithelial stem cells: Sox2, Lrig1, Bmi1, Gli1, Igfbp5, and Lgr5
# ameloblast differentiation, including 
# preameloblasts (Shh+ cluster 11), 
# secretory (Enam+ cluster 5), 
# maturation (Klk4+ cluster 10) and 
# postmaturation (Gm17660+ cluster 6) 

###

##------------------ The 6th Step：merge all pure cell again - final -------------
setwd("/home/wxf/data/projects/Cov19/analysis_new/")

sample="mergedSamples_pure"
sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_2/"
load(paste0(sampleoutpath,sample,".Rda"))

cellRM=readLines(paste0(sampleoutpath,"doublet_basalCell_needremove.txt"))

subData=subset(Sdata,cells=cellRM,invert=T)

mergeWorkFlowFun(Sdata=subData,sample="mergedSamples_pure",
                 sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_3/",
                 nfeatures = 3000,resolution = 1.5)

sample="mergedSamples_pure"
sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_3/"
load(paste0(sampleoutpath,sample,".Rda"))

Sdata = FindClusters(Sdata,resolution = 1.5)
UMAPPlot(Sdata,label=T,group.by="RNA_snn_res.1.5",cols=rep(mycol_41,3))+NoLegend()
UMAPPlot(Sdata,cells.highlight=rownames(allCell_pure@meta.data)[which(allCell_pure$RNA_snn_res.3.5=="43")])

FeaturePlot(Sdata,features = c("Mki67"))

Genes = c(
  ## blood cell
  # "Pf4","Gp9","Ppbp", #megakaryocyte
  # "Cdh5", "Kdr","Sema3g","Efnb2","Dll4", #Vascular endothelial cell
  # "Lyve1","Prox1","Tbx1", #Lymphatic endothelial cell
  # "Rgs5","Cspg4", #Pericyte
  # "Omp","Nrxn1", #Neuron
  # "Gypa","Gypb","Tfrc","Alas2","Hbb−bt","Hba−a1","Bpgm", #erythroblast
  # "Cpa3","Ms4a2","Fcer1a", #Mast cell
  # "Elane","Mpo","Lyz2", # GMP
  # "Cd68","Cd14","Csf1r", "Fcgr1", #Monocyte
  # "Csf3r","Cxcr2", "Ngp", "S100a8", "S100a9", #Neutrophil
  # "Cd79a", "Vpreb1","Vpreb3", "Pax5","Ighm","Ighd","Sell", #B cell
  # "Pdgfra","Dcn","Col1a1", #Fibroblast
  # "Cthrc1", "Acta2", "Postn", "Adam12",
  # "Col15a1", "Col4a1", "Hspg2",
  # "Pi16", "Dpp4", "Ly6h",
  # "Epcam"
  # ,
  "Epcam","Krt17","Krt16","Krt6a","Krt6b","Krt15","Krt19","Krt7", #basal cell
  "Krt18","Krt8",
  # "Tgm2","Loricrin","Flg","Ivl",
  "Krt10","Krt1","Krt2","Krt13","Krt4",
  "Krt5","Krt14","Trp63","Icam1", "Aqp4", "Ly6d","Sox2","Pax6","Hes1","Sox9", #HBC
  "Adh7", #respiratory basal cell/HBC
  "Sprr1a", #Activated/Cycling HBC
  "Kit", "Ascl1","Cxcl14","Neurod1","Neurog1","Hes6","Ezh2","Hmga1","Mki67","Stmn1","Top2a",#GBC
  # "Ascl3","Cftr","Fgf10","Fgfr2b", #Ionocyte
  # "Trpm5","Avil","Sh2d7","Lrmp",#Tuft cell
  # "Dnah5","Tmem212","Cdhr3", #Ciliated cell
  # "Tff2","Muc5ac","Muc5b","Bpifa1", #Goblet cell
  # # "Lipf","Dcpp1","Dcpp2","Dcpp3","Notch3",#Goblet cell2
  # "Scgb1a1","Scgb3a2","Cyp2f2","Lypd2","Cbr2","Nfia", #Club cell
  # "Foxj1", "Pifo", "Ift57", "Lrrc23",
  # "Nrcam","Vit", #olfactory HBC
  # "Gp2","Tnfaip2","Ly6d", #M cell
  # "Gsto1","Csta1","Asprv1","Psca","Prss22", #Hillock cell
  # "Sox9","Rgs5","Msln","Vmo1","Chil6","Prss33","Scgb1c1","Ldhd","Gpx6","Bpifb6", #Bowman gland cell
  # "Cyp2g1","Muc2", "Sec14l3","Cyp1a2","Notch2","Ermn", "Cyp2a13","Cyp2j2", #SUS
  # "Ephx1", #SUS ventral
  # "Sult1c1","Aqp5", #SUS dorsal
  "Rbm24", "Neurod1", "Elavl4", "Scg2","Lhx2", "Tex15", "Crabp1", #intermediate neural precursor (INP)
  "Gng8","Gap43","Tuj1", #immature OSN
  "Omp","Gng13","Cd36", "Cnga2", "Adcy3", #mature OSN
  "Calr4","Bmp6", #vomeronasal sensory neurons, VSN
  "Plp1", "Npy", "Pmp22", "Mpz","Mog","Mbp","Klk6","Apod",#"Rip","Olig1","Olig2","Itgam","Iba1", #Npy Neuron - Oligodendrocyte
  "Slc6a11","Ttyh1","Gria2","Mt3", "Bcan","Slc1a3", "Slc4a4","Aqp4","Gjal","Slc1a2","Gfap" #Gria2 Neuron - Astrocyte
  # "Fezf2","Meg3","Car9",#Fezf2+ Epithelial cell
  # "Scgb2b27","Car6","Pax3","Dnase1","Klk14","Obp1a","Obp1b","Obp2a","Obp2b", #Obp1a+Obp1b+ Epithelial cell
  
  # "Cxcl17","Ccl9","Ltf","Bglap3","Aox3","Cyp4b1","Tac1","Dmbt1","Slc39a8", #Cxcl17+Ccl9+ Epithelial cell
  # "Bpifb9a","Bpifb9b","Bpifb5","Siglcc8", #Bpifb9a+Bpifb5+ Epithelial cell
  # "Fxyd2","Klk1", "Krt7", "Krt19","Tfcp2l1","Scnn1b", #Krt19, Krt7, and Tcfcp2l1 are ductal maturation markers
  # "Pip", # markers for salivary gland maturation
  # "Aqp5","Bhlha15"#,  #are markers for terminal differentiation of SMG epithelial cells (acinar and duct, respectively)
  # "Gstt1", 
  # "Gfra3","Kit",
  # "Etv4", "Fgfr2","Mmp2","Esp1","Esp22"
)
Genes <- intersect(Genes,rownames(Sdata[["RNA"]]@counts))

DotPlot(Sdata,features = Genes,group.by = "celltype",cols = c("lightgrey","red"))+RotatedAxis()#group.by = "RNA_snn_res.3.5",

UMAPPlot(Sdata,group.by="RNA_snn_res.1.5",label=T,cols=rep(mycol_41,3))+NoLegend()

##------------------ The 7th Step：redefine cluster ------------
#cluster 68 - allCell res4
cyclingRHBC=rownames(Sdata@meta.data)[which(Sdata$RNA_snn_res.4=="68")]
#cluster 6 - BasalCell res0.8
Mcell=rownames(BasalCell@meta.data)[which(BasalCell$RNA_snn_res.0.8=="6")]
UMAPPlot(Sdata,cells.highlight=rownames(Sdata@meta.data)[which(Sdata$RNA_snn_res.4=="68")])
UMAPPlot(Sdata,cells.highlight=rownames(BasalCell@meta.data)[which(BasalCell$RNA_snn_res.0.8=="6")],sizes.highlight = 0.1)

## redefine two define cluster: add cyclingRHBC and Mcell
Sdata$redefineCluster=as.character(Sdata$RNA_snn_res.1.5)
Sdata$redefineCluster[cyclingRHBC]="47" #73 cells
Sdata$redefineCluster[Mcell]="48" #242 cells
UMAPPlot(Sdata,group.by="redefineCluster",label=T,cols=rep(mycol_41,3))+NoLegend()

table(Sdata$orig.ident,Sdata$redefineCluster)
Idents(Sdata)="redefineCluster"

##------------------ The 8th Step：reorder and define cell types for merged data -------------
celltypeFile = read.csv("./mergedSamples_pure/mergedSamples_pure_3/reorder_celltype_new.csv")
Sdata$newcluster=plyr::mapvalues(Sdata$redefineCluster,from=as.character(celltypeFile$oldcluster),to=as.character(celltypeFile$newcluster))
Sdata$newcluster=factor(Sdata$newcluster,levels = as.character(celltypeFile$newcluster))
Sdata$celltype=plyr::mapvalues(Sdata$newcluster,from = as.character(celltypeFile$newcluster),to = as.character(celltypeFile$celltype))
Sdata$Tcelltype=plyr::mapvalues(Sdata$newcluster,from = as.character(celltypeFile$newcluster),to = as.character(celltypeFile$Tcelltype))
Idents(Sdata)="newcluster"
Sdata$celltypeColor=plyr::mapvalues(Sdata$newcluster,from = as.character(celltypeFile$newcluster),to = mycol_49)
Sdata$TcelltypeColor=plyr::mapvalues(Sdata$Tcelltype,from = as.character(unique(Sdata$Tcelltype)),to = TcelltypeCols)
Sdata$orig.ident=factor(Sdata$orig.ident,levels=c("Mock","Alpha","Lambda"))
Sdata$tissueColor=plyr::mapvalues(Sdata$orig.ident,from = levels(Sdata$orig.ident),to = tissueCols)

UMAPPlot(Sdata,group.by="newcluster",label=T,cols=mycol_49)+NoLegend()

UMAPPlot(Sdata,label=F,repel=T,cols=levels(Sdata$celltypeColor),group.by="Tcelltype")+NoLegend()
TSNEPlot(Sdata,label=T,repel=T,cols=mycol_41)+NoLegend()

celltypeFile$celltypeColor=mycol_49
write.csv(celltypeFile,file = paste0(sampleoutpath,"reorder_celltype_new.csv"),row.names = F)

save(Sdata,file = paste0(sampleoutpath,sample,".Rda"))
#