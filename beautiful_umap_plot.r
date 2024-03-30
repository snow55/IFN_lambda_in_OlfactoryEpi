library(ggplot2)
beauty2Dplot=function(Sdata=NULL,labelCelltype="celltype",labelID="newcluster",
                      splitPath=NULL,cols=mycol_41){
  #parapared data
  df=as.data.frame(Sdata@reductions$umap@cell.embeddings)
  df$labelID=Sdata@meta.data[rownames(df),labelID]
  df$celltype=Sdata@meta.data[rownames(df),labelCelltype]
  # get coord for labels
  dat.plot.label=as.data.frame(t(sapply(
    split(df[,1:2], df[,3]),
    function(x){
      apply(x, 2, median)
    }
  )))
  dat.plot.label$labelID=rownames(dat.plot.label)
  #plot UMAP
  pdf(paste0(splitPath,"UMAP_celltype_labeled_sample.pdf"),width = 14,height = 6)
  p=ggplot(df,aes(x=UMAP_1,y=UMAP_2, fill=labelID))+
    geom_point(shape=21,stroke=0.005,size=2,color="white")+#"grey20",alpha=0.5
    scale_fill_manual(values = cols,
                      labels = paste0(levels(Sdata@meta.data[,labelID]),". ",levels(Sdata@meta.data[,labelCelltype])))+
    geom_point(aes(x=UMAP_1,y=UMAP_2, label= labelID ), data=dat.plot.label,
               size=6, #点的半径
               show.legend = F, #不要这一层的图例
               stroke=1, #设置边的宽度
               #color 设置线的颜色
               color=cols,
               #fill="white", alpha=0.6,# alpha会影响到边框线，不如直接给fill设置不透明颜色
               fill="#FFFFFF99", #填充色，后两位是不透明度(0.6->#99)
               shape=21)+ # shape=21 就是圆圈
    geom_text(aes(x=UMAP_1,y=UMAP_2, label= labelID ), data=dat.plot.label,
              size=5,color="black",show.legend = F)+
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_blank(),
          panel.border = element_blank(),
          axis.title.x = element_text(size = 15),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.height=unit(0.7,"cm"),
          legend.key.width=unit(2,'mm'),
          legend.margin=margin(b = -0.3, unit='cm'),
          legend.title = element_text(size = 17)) +
    guides(fill = guide_legend(override.aes = list(size=4),title = '',ncol = 2))+ #可使得legend中的圆点变大
    labs(x="UMAP 1",y="UMAP 2")
  print(p)
  dev.off()
}
  
beauty2Dplot(Sdata=Sdata,splitPath = splitPath)