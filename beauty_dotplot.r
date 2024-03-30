##-------------------- beauty dotplot of specific gene expression -----------------

# Genes = c(#"Krt18","Cyp2g1","Muc2","Sox9","Chil6","Scgb1c1",
#           'Mx1','Isg15','Stat1','Ifit1','Ifit3b',
#           # 'Ddx58','Ifih1',#RIG-I,MDA-5
#           'H2-D1','H2-K1','H2-T23','Tap1','Tap2','H2-Q4','Ifi47'
# )

Genes = c("Mx1", "Isg15", "Ifit1","Ifit3", "Rsad2", "Oasl2", "Stat1","Stat2",
          "Ifi44", "Ifit3b", "Irgm1", "Gbp7")

# Genes = rownames(Sdata[['RNA']]@data)[grep('H2-',rownames(Sdata[['RNA']]@data))]#all MHC moleculars
subData = subset(Sdata,cells=rownames(Sdata@meta.data)[Sdata$Tcelltype 
                                                       %in% c('BGC')] )#'Ciliated cell','Fezf2+ epithelial cell','Olfactory HBC'
# subData$group = paste0(subData$Tcelltype,': ',subData$orig.ident)
subData$group = subData$orig.ident
#dotplot的改造  
DotPlot(subData, features = Genes,group.by = 'group',dot.scale = 7)+ 
  # coord_flip()+
  theme_bw()+#去除背景，旋转图片  3.  
  theme(panel.grid = element_blank(), # 4.        
        axis.text = element_text(size=14,color='black'),
        axis.text.x=element_text(angle=45,hjust = 1,color='black'))+#文字90度呈现  5.  
  # scale_y_discrete(limits=c("Fezf2+ epithelial cell: Mock","Fezf2+ epithelial cell: Alpha","Fezf2+ epithelial cell: Lambda",
  #                           "Ciliated cell: Mock","Ciliated cell: Alpha","Ciliated cell: Lambda")#,
  #                   # labels=rep(c('Mock',expression("IFN-"*alpha),expression("IFN-"*lambda)),2)
  #                   )+
  scale_y_discrete(labels=rep(c('Mock',expression("IFN-"*alpha),expression("IFN-"*lambda)),2))+##limits=c("BGC: Mock","BGC: Alpha","BGC: Lambda"),,"SUS: Mock","SUS: Alpha","SUS: Lambda"
  scale_color_gradientn(colours = c('#688DCA','grey90','orange','red'))+#c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  6.
  labs(x=NULL,y=NULL,size='Percentage',colour=expression('Average\nexpression'))+#values = seq(0,1,0.2),
  guides(size=guide_legend(order=1))#+
ggsave(filename = paste0(sampleoutpath,'Dotplot_BGC_ISG-common.pdf'),h=2,w=7)
ggsave(filename = paste0(sampleoutpath,'Dotplot_BGC_ISG-common_legend.pdf'),h=4,w=7)

# guides(size=guide_legend(title='Percentage'),color=guide_legend(title='Average\nexpression'))
# ggsave(filename = paste0(sampleoutpath,'Dotplot_Ciliated-Fezf2_ISG-MHC1.pdf'),h=3,w=7.5)
# ggsave(filename = paste0(sampleoutpath,'Dotplot_SUS-BGC_ISG-MHC1.pdf'),h=3,w=7)
