##-------------------- QuSAGE 计算基因集的在两组间的差异 -------------
# transfer p value for base R plotting
transPvalue=function(p0){
  # p0=1.2e-3; p0
  if (p0<2.2e-16) {
    p0=2.2e-16
    p0=formatC(p0, format="e", digits=2);#p0 #"1.2e-52"
    # split by e
    p0=strsplit(p0, "e")[[1]]; #p0 #[1] "1.2" "-52"
    label = bquote(italic(P)~"<"~.(p0[1])~"\u00d7"~10^.(p0[2]) )
  }else if(p0<0.01){
    p0=formatC(p0, format="e", digits=2);#p0 #"1.2e-52"
    # split by e
    p0=strsplit(p0, "e")[[1]]; #p0 #[1] "1.2" "-52"
    # plot
    label = bquote(italic(P)~"="~.(p0[1])~"\u00d7"~10^.(as.numeric(p0[2]) ))
  }else{
    p0=round(p0,3); #p0 #保留3位
    label = bquote(italic(P)~"=" ~ .(p0) )
  }
  return(label)
}

# http://bioconductor.riken.jp/packages/3.0/bioc/vignettes/qusage/inst/doc/qusage.pdf
# install.packages("http://bioconductor.riken.jp/packages/3.0/bioc/src/contrib/qusage_1.6.0.tar.gz", repos=NULL, type="source")
library(qusage)
#### for single celltype
MHC_class_I = readLines('./MHC_class_I_geneset.txt')
geneSets <- list("Class I MHC"=MHC_class_I)#,"Phagocytosis"=phagosome,"Inflammatory response"=inflammatory)
geneSetsCols = c("#f01a2a", "#1abc9c","#f9a31c","#003bff")
ct="SUS"
subData = subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$orig.ident=='Mock')],invert=T)
eset=as.data.frame(subData[["RNA"]]@data[,rownames(subData@meta.data)[which(subData$Tcelltype %in% ct)]])
labels=as.character(subData$orig.ident[colnames(eset)])
contrast="Lambda-Alpha"
sf=qusage(eset,labels,contrast,geneSets)
p.vals=pdf.pVal(sf)
plot(sf, xlab="Gene Set Activation")
q.vals = p.adjust(p.vals, method="fdr")
qsTable(sf, number=10) #提取logFC, p.vals, FDR值

##LA
pdf(paste0(sampleoutpath,"genesetScore_",ct,"_QuSAGE.pdf"),w=5.5,h=4,useDingbats = F)
par(oma=c(1,1,1,1),mar=c(4,4.5,3,0))
plotDensityCurves(sf, col=geneSetsCols[1:3],cex.axis=1.2,cex.lab=1.8,cex.main=1.8,xlim = c(-0.07,0.3),#sub=ct,
                  main=expression(italic('Sirpa')^italic('-/-')*' vs WT'),
                  lwd=2, xlab=expression(Log[2](Activity)),ylab = "Density")
box(lwd=2)
legend("topright", legend=c("SUS"),lty=1,lwd=2,col=geneSetsCols[1:3],cex=1)
text(x = c(0.18), y = c(10), labels =transPvalue(q.vals[1]), col = geneSetsCols[1],cex = c(0.8))
text(x = c(-0.01), y = c(22), labels = transPvalue(q.vals[2]),col = geneSetsCols[2],cex = c(0.8))
text(x = c(0.0), y = c(50), labels =transPvalue(q.vals[3]), col = geneSetsCols[3],cex = c(0.8))
dev.off()



#### for single geneSet -- MHC class I
#LA, LM, AM
myCols = c('#E52287',  '#E6AF0C','#78A51E')
subData = Sdata
ct='Fezf2+ epithelial cell' #"Ciliated cell" #"SUS"
eset=as.data.frame(subData[["RNA"]]@data[,rownames(subData@meta.data)[which(subData$Tcelltype %in% ct)]])
labels=as.character(subData$orig.ident[colnames(eset)])
sf_LA=qusage(eset,labels,"Lambda-Alpha",geneSets)
sf_LM=qusage(eset,labels,"Lambda-Mock",geneSets)
sf_AM=qusage(eset,labels,"Alpha-Mock",geneSets)

# eset=as.data.frame(subData[["RNA"]]@data[,rownames(subData@meta.data)[which(subData$Tcelltype %in% c("SUS", "BGC"))]])
# labels=as.character(paste(Sdata$Tcelltype[colnames(eset)],Sdata$orig.ident[colnames(eset)],sep="."))
# sf_LA=qusage(eset,labels,"SUS.Lambda-SUS.Alpha",geneSets)
# sf_LM=qusage(eset,labels,"BGC.Lambda-BGC.Alpha",geneSets)

p.vals_LM=pdf.pVal(sf_LM)
q.vals_LM = p.adjust(p.vals_LM, method="fdr")
p.vals_AM=pdf.pVal(sf_AM)
q.vals_AM = p.adjust(p.vals_AM, method="fdr")
p.vals_LA=pdf.pVal(sf_LA)
q.vals_LA = p.adjust(p.vals_LA, method="fdr")

#SUS
pdf(paste0(sampleoutpath,"genesetScore_MHCI_","Lambda-Alpha-Mock_",ct,"_QuSAGE.pdf"),w=5,h=4.3,useDingbats = F)
par(oma=c(1,1,1,1),mar=c(4,4.5,3,0))
plotDensityCurves(sf_LA,col=myCols[1],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("MHC class I: "*'SUS') ,#
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.005,0.32),ylim=c(0,100))
plotDensityCurves(sf_LM,col=myCols[2],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_AM,col=myCols[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c(expression('IFN-'*lambda*' vs IFN-'*alpha),
                            expression('IFN-'*lambda*' vs Mock'),
                            expression('IFN-'*alpha*' vs Mock')),lty=1,lwd=2,col=myCols[c(1,2,3)],cex=1.2)
text(x = c(0.203), y = c(32), labels =transPvalue(q.vals_LA), col = myCols[1],cex = c(0.8))
text(x = c(0.26), y = c(51), labels =transPvalue(q.vals_LM), col = myCols[2],cex = c(0.8))
text(x = c(0.06), y = c(45), labels = transPvalue(q.vals_AM),col = myCols[3],cex = c(0.8))
dev.off()

#BGC
pdf(paste0(sampleoutpath,"genesetScore_MHCI_","Lambda-Alpha-Mock_",ct,"_QuSAGE.pdf"),w=5,h=4.3,useDingbats = F)
par(oma=c(1,1,1,1),mar=c(4,4.5,3,0))
plotDensityCurves(sf_LA,col=myCols[1],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("MHC class I: "*'BGC') ,#
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.005,0.15),ylim=c(0,100))
plotDensityCurves(sf_LM,col=myCols[2],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_AM,col=myCols[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c(expression('IFN-'*lambda*' vs IFN-'*alpha),
                            expression('IFN-'*lambda*' vs Mock'),
                            expression('IFN-'*alpha*' vs Mock')),lty=1,lwd=2,col=myCols[c(1,2,3)],cex=1.2)
text(x = c(0.038), y = c(45), labels =transPvalue(q.vals_LA), col = myCols[1],cex = c(0.8))
text(x = c(0.105), y = c(35), labels =transPvalue(q.vals_LM), col = myCols[2],cex = c(0.8))
text(x = c(0.025), y = c(60), labels = transPvalue(q.vals_AM),col = myCols[3],cex = c(0.8))
dev.off()

#Ciliated cell
pdf(paste0(sampleoutpath,"genesetScore_MHCI_","Lambda-Alpha-Mock_",ct,"_QuSAGE.pdf"),w=5,h=4.3,useDingbats = F)
par(oma=c(1,1,1,1),mar=c(4,4.5,3,0))
plotDensityCurves(sf_LA,col=myCols[1],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("MHC class I: "*'Ciliated cell') ,#
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.15,0.15),ylim=c(0,60))
plotDensityCurves(sf_LM,col=myCols[2],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_AM,col=myCols[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c(expression('IFN-'*lambda*' vs IFN-'*alpha),
                            expression('IFN-'*lambda*' vs Mock'),
                            expression('IFN-'*alpha*' vs Mock')),lty=1,lwd=2,col=myCols[c(1,2,3)],cex=1.2)
text(x = c(-0.115), y = c(20), labels =transPvalue(q.vals_LA), col = myCols[1],cex = c(0.8))
text(x = c(-0.02), y = c(32.5), labels =transPvalue(q.vals_LM), col = myCols[2],cex = c(0.8))
text(x = c(0.10), y = c(20), labels = transPvalue(q.vals_AM),col = myCols[3],cex = c(0.8))
dev.off()

#Fezf2+ epithelial cell
pdf(paste0(sampleoutpath,"genesetScore_MHCI_","Lambda-Alpha-Mock_",ct,"_QuSAGE.pdf"),w=5,h=4.3,useDingbats = F)
par(oma=c(1,1,1,1),mar=c(4,4.5,3,0))
plotDensityCurves(sf_LA,col=myCols[1],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression('Fezf2+ epithelial cell') ,#
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.2,0.3),ylim=c(0,40))
plotDensityCurves(sf_LM,col=myCols[2],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_AM,col=myCols[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c(expression('IFN-'*lambda*' vs IFN-'*alpha),
                            expression('IFN-'*lambda*' vs Mock'),
                            expression('IFN-'*alpha*' vs Mock')),lty=1,lwd=2,col=myCols[c(1,2,3)],cex=1.2)
text(x = c(-0.1), y = c(20), labels =transPvalue(q.vals_LA), col = myCols[1],cex = c(0.8))
text(x = c(0.02), y = c(18), labels =transPvalue(q.vals_LM), col = myCols[2],cex = c(0.8))
text(x = c(0.18), y = c(13), labels = transPvalue(q.vals_AM),col = myCols[3],cex = c(0.8))
dev.off()
#