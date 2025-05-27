library(reshape2)

custom_dotplot<- function(seurat_list, seurat_cluster_feature, marker_gene_list, dot.scale = 1.5, radius_size = 3, col.min = -2.5, col.max = 2.5, isMainPlot = T, cluster_subset = NULL){
  # Create dataframe
  all_markers<- unlist(marker_gene_list)
  expr_df_ls<- list()
  for(s in names(seurat_list)){
    seurat_tmp<- seurat_list[[s]]
    expr_data_tmp<- FetchData(object = seurat_tmp, vars = c(seurat_cluster_feature,all_markers))
    expr_data<- expr_data_tmp[,-1]
    clusters_tmp<- expr_data_tmp[,1]
    # expr_data<- t(apply(expr_data_tmp[,-1], 1, function(x) scale(log1p(x))))
    # colnames(expr_data)<- colnames(expr_data_tmp[,-1])
    # rownames(expr_data)<- expr_data_tmp[,1]
    pct.exp<- melt(apply(expr_data, 2, function(x) tapply(x, clusters_tmp, PercentAbove, 0)), varnames = c("cluster", "gene"), value.name = "pct_exp")
    avg.exp<- apply(expr_data, 2, function(x) tapply(x, clusters_tmp, mean, na.rm=T))
    avg.exp[]<- apply(avg.exp, 2, function(x) scale(log1p(x))) # Scale & log
    avg.exp<- melt(avg.exp, varnames = c("cluster", "gene"), value.name = "avg_exp")
    avg.exp$avg_exp<- MinMax(data = avg.exp$avg_exp, min = col.min, max = col.max) # Min/max values
    expr_df<- merge(pct.exp, avg.exp)
    expr_df$cond<- s
    expr_df_ls[[s]]<- expr_df
  }
  # expr_df<- rbind(expr_df_ls$NB1Pre, expr_df_ls$NB1Post, expr_df_ls$NB2Post)
  expr_df<- do.call("rbind", expr_df_ls)
  expr_df$gene_class<- NA
  for(clust in names(marker_gene_list)){
    expr_df$gene_class[expr_df$gene%in%marker_gene_list[[clust]]]<- clust 
  }
  expr_df$gene<- factor(expr_df$gene, levels = rev(unique(sort(all_markers))))  
  expr_df$cond<- factor(expr_df$cond, levels = names(seurat_list))
  if(isMainPlot) expr_df$cluster<- factor(expr_df$cluster, levels = names(marker_gene_list))
  expr_df$gene_class<- paste0("Markers ", expr_df$gene_class)
  expr_df$gene_class<- factor(expr_df$gene_class, levels = paste0("Markers ", names(marker_gene_list)))
  if(!is.null(cluster_subset)){
    expr_df<- expr_df[expr_df$cluster%in%cluster_subset,]
    expr_df$cluster<- factor(expr_df$cluster, levels = cluster_subset)
  } 
  
  # Plot
  ggplot(data = expr_df, mapping = aes_string(x = "cluster", y = "gene")) + 
    geom_point(mapping = aes_string(size = "pct_exp", colour = "avg_exp")) +
    scale_color_gradient2(low = "blue", mid = "gray", high = "red", breaks = c(-dot.scale,0,dot.scale)) +
    facet_grid(gene_class~cond, scales = "free", space = "free") +
    theme_cowplot() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 6, face = "italic"),
      axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.key.width = unit(0.25, "cm"),
      legend.key.height = unit(0.1, "cm"),
      legend.position = "bottom" ,
      strip.background = element_blank(),
      strip.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 5)
    ) +
    guides(size = guide_legend(title = "Percent Expressed")) + 
    scale_radius(range = c(0, radius_size))
}

