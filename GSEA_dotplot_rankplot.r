##--------- GSEA by clusterProfiler ------------------
library(clusterProfiler)
library(org.Mm.eg.db)
sampleoutpath = "./mergedSamples_pure/mergedSamples_pure_3/mergedSamples_rmImmuEry/"
sampleoutpath_GSEA = paste0(sampleoutpath,'GSEA/')
if(!dir.exists(sampleoutpath_GSEA)){dir.create(sampleoutpath_GSEA)}

df_all = read.csv(paste0(sampleoutpath,'DFgenes/','DFgenes_Lambda_vs_Alpha_25ct.csv'),row.names = 1)
head(df_all)
# run gsea
for(ct in unique(df_all$Tcelltype)){
  # ct='BGC'
  df = df_all[df_all$Tcelltype==ct,]
  df = df[order(df$avg_log2FC,decreasing = T),]
  geneList = df$avg_log2FC
  names(geneList)=df$gene #eg$ENTREZID
  res <- gseGO(
    geneList,    # 根据logFC排序的基因集
    ont = "BP",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
    OrgDb = org.Mm.eg.db,    # 使用人的OrgDb
    keyType = 'SYMBOL',#"ENTREZID",    # 基因id类型
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",    # p值校正方法
  )
  if(dim(res)[1]==0){next}
  save(res, file = paste0(sampleoutpath_GSEA,"GOBP_",ct,"_GSEA_result_SYMBOL.RData"))
}
# extract single signaling pathway and combined together
tmp_res = NULL
for(ct in unique(df_all$Tcelltype)){
  # ct='M cell' # 'Ionocyte'   #'Tuft cell' #
  resfile=paste0(sampleoutpath_GSEA,"GOBP_",ct,"_GSEA_result_SYMBOL.RData")
  if(file.exists(resfile)){
    load(resfile)
    res@result$celltype = ct
    tmp = res[res$Description=='response to virus',]#'defense response to virus'
    if(dim(tmp)[1]==0){
      next
    }
    if(is.null(tmp_res)){
      tmp_res=tmp
    }else{
      tmp_res=rbind(tmp_res,tmp)
    }
  }
}
head(tmp_res,2)
table(tmp_res$celltype)

tmp_res = tmp_res[-which(tmp_res$celltype %in% c('Endothelial cell','Fibroblast')),]
tmp_res[,c(1:9,ncol(tmp_res))]

## Dotplot of GSEA
res_new = tmp_res
res_new$Description=res_new$celltype
res_new[,c(1:9,ncol(res_new))]
res_new$celltype=NULL
head(res_new,2)
# res_new$Description = factor(res_new$Description,
#                              levels=c("Olfactory HBC","GBC","INP","OSN","SUS","BGC","Ciliated cell","Bpifb5+ gland cell","Fezf2+ epithelial cell"))
p = dotplotGsea(data = res_new,topn = 10, order.by = 'NES',add.seg = T)
p1=p$plot+ theme(plot.title = element_text(hjust = 0.5))+
          labs(color='P.adjust',title='GO:0009615: response to virus',x=expression('NES (IFN-'*lambda*' vs IFN-'*alpha*')') )
p1
ggsave(p1,filename = paste0(sampleoutpath_GSEA,'Dotpot_GOBP_GSEAres_all_celltype_Lambda_vs_Alpha.pdf'),w=7,h=3.5)
ggsave(p1,filename = paste0(sampleoutpath_GSEA,'Dotpot_GOBP_GSEAres_all_celltype_Lambda_vs_Alpha.pdf'),w=5.5,h=3.7)

## classic rank plot GSEA
for(ct in unique(df_all$Tcelltype)){
  # ct= 'OSN' #'SUS' #'Ciliated cell' #
  resfile=paste0(sampleoutpath_GSEA,"GOBP_",ct,"_GSEA_result_SYMBOL.RData")
  if(file.exists(resfile)){
    load(resfile)
    res@result$celltype = ct
    tmp = res[res$Description=='response to virus',]
    if(dim(tmp)[1]==0){
      next
    }
    library(GseaVis)
    # # single signaling pathway
    mygene <- c('Ifit2','Gbp4','Oas1a','Irf7',"Mx1","Isg15",'Stat1')#add gene in specific pathway
    if(tmp$NES<0){
      pval_position = c(0.05,0.2)
    }else(
      pval_position = c(0.45,0.8)
    )
    p=gseaNb(object = res,
             geneSetID = 'GO:0009615',#'response to virus',
             subPlot = 2,
             termWidth = 30,  #termWidth限制图片中 term name 的长度：
             # newGsea = T,
             addPval = T,
             pvalX = pval_position[1],
             pvalY = pval_position[2],
             pCol = 'black',
             pHjust = 0,
             htHeight=0.5,
             addGene = mygene
    )
    p1=p+labs(x=paste0('Rank in Ordered Dataset\n',ct))#;p1
    ggsave(p1,filename = paste0(sampleoutpath_GSEA,'Rankplot_',ct,'_GOBP_GSEAres_Lambda_vs_Alpha.pdf'),w=4,h=3.0)
  }
}

## 单独可视化每个细胞类型的GSEA结果
ct='Cyp4b1+ gland cell'  #'SUS'
resfile=paste0(sampleoutpath_GSEA,"GOBP_",ct,"_GSEA_result_SYMBOL.RData")
load(resfile)
dim(res)

#指定特定信号通路
useID = grep('response to virus',res$Description)
res[useID,] #[c(1,2)]
# gseaplot2(res, title = res$Description[useID[1]], geneSetID = useID[1],pvalue_table = F)
# dotplotGsea(data = res,topn = 10, order.by = 'NES',add.seg = T)

library(GseaVis)
# # single signaling pathway
# # add gene in specific pathway
mygene <- c('Ifit2','Gbp4','Oas1a','Oas2','Irf7',"Mx1","Isg15")## mygene = c(15958) #ENTREZID
p=gseaNb(object = res,
       geneSetID = 'GO:0009615',#'response to virus',
       subPlot = 2,
       termWidth = 30,  #termWidth限制图片中 term name 的长度：
       # newGsea = T,
       addPval = T,
       pvalX = 0.65,pvalY = 0.8,
       pCol = 'black',
       pHjust = 0,
       htHeight=0.5,
       addGene = mygene
       )
p+labs(x=paste0('Rank in Ordered Dataset\n',ct))


# multiple signaling pathway
useID = grep('response to virus',res$Description)
geneSetID = res[useID[c(1)],'ID']
# geneSetID=res$ID[1]
gseaNb(object = res,
       geneSetID = geneSetID,
       subPlot =2,
       # newGsea = T,
       termWidth = 25,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 0.05,pvalY = 0.05,
       curveCol = TcelltypeCols[7:8]#,
       # markTopgene=T,topGeneN = 2#,
       # ght.relHight=0.4
       # htHeight=0.5
       # addGene = mygene
       )
