library(ggplot2)
library(pheatmap)

##-------------------- heatmap of ISG genes 20240122
Genes <- c('Mx1',"Isg15",'Stat1','Ifit1','Rsad2','Gbp7')  #"Isg15" #'Mx1' #
ct=c("Olfactory HBC","GBC","INP","OSN","SUS","BGC","Oligodendrocyte",
     "RBC","Ciliated cell","Fezf2+ epithelial cell",
     "Cyp4b1+ gland cell","Dmbt1+ gland cell","Bpifb5+ gland cell")
sample_levels <- levels(Sdata$orig.ident)

# 创建一个空的列表用于存储所有的结果
result_list <- list()

for (gene in Genes) {
  df2 <- NULL
  
  for (s in sample_levels) {
    tmp <- df[df$sample == s, ]
    df1 <- tapply(tmp[, gene], tmp$celltype, mean)
    df2 <- rbind(df2, df1)
  }
  
  rownames(df2) <- sample_levels
  df2 <- t(df2)
  df2 <- df2[ct, ]
  
  # 存储结果到列表
  result_list[[gene]] <- df2
}

# 绘制每个基因的热图
heatmap_list <- lapply(result_list, function(df2) {
  pheatmap(
    df2,
    cluster_cols = F,
    cluster_rows = F,
    border_color = "white",
    scale = "row",
    show_rownames = T,
    show_colnames = T,
    labels_col = c('Mock', expression('IFN-' * alpha), expression('IFN-' * lambda)),
    gaps_row = 7,
    cellheight = 14,
    cellwidth = 16,
    angle_col = "45",
    fontsize_col = 12,
    fontsize_row = 11,
    main = colnames(df2),
    # border = 'black',
    heatmap_legend_param = list(
      title = "Expression",
      legend_height = unit(4, "cm"),
      title_position = "lefttop-rot"
    ),
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
})

# 合并所有热图到一个 grid
library(gridExtra)
library(ComplexHeatmap)
grid_combined <- grid.arrange(grobs = heatmap_list)

p_combined = p_Mx1 + p_Isg15 + p_Stat1
p_combined = draw(p_combined, ht_gap = unit(0.5, "cm"))
pdf(paste0(sampleoutpath,'Mx1-Isg15-Stat1_exp_heatmap.pdf'),w=6,h=4)
print(p_combined)
dev.off()

grid_combined <- grid.arrange(
  heatmap_list[[1]], heatmap_list[[2]], heatmap_list[[3]],
  nrow = 1
)

cowplot::plot_grid(heatmap_list[[1]]$gtable, heatmap_list[[2]]$gtable,heatmap_list[[3]]$gtable,
                   heatmap_list[[4]]$gtable, heatmap_list[[5]]$gtable,heatmap_list[[6]]$gtable,ncol= 6)

library(cowplot)

# Your existing code for creating the heatmap_list

# Function to modify gtable (remove legend and row names)
modify_gtable <- function(gtable) {
  # Check if "guide-box" exists in the layout
  if ("guide-box" %in% gtable$layout$name) {
    gtable$grobs[[which(gtable$layout$name == "guide-box")]] <- nullGrob()
  }
  # Check if "row_names" exists in the layout
  if ("row_names" %in% gtable$layout$name) {
    gtable$grobs[[which(gtable$layout$name == "row_names")]] <- nullGrob()
  }
  return(gtable)
}

# Modify gtables for the first 5 heatmaps
modified_heatmap_list <- lapply(heatmap_list[1:5], function(heatmap) {
  modify_gtable(heatmap$gtable)
})

# Use plot_grid with modified gtables
plot_grid(plotlist = modified_heatmap_list, ncol = 6)
draw(modified_heatmap_list, ht_gap = unit(0.5, "cm"))

# 将图形存储为 PDF 文件
pdf(paste0(sampleoutpath, 'Mx1-Isg15-Stat1-6Genes_exp_heatmap_20240122.pdf'), w = 6, h = 4)
print(grid_combined)
dev.off()
