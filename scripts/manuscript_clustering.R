# Load data
###########
seurat_grp_merged_list <- read_rds("data/ST_NB_seurat.rds")
# seurat_grp_merged_list$NB2Pre<- NULL

# Color palette
#################

col_palette_main <- c(
  NE = "#A6CEE3",
  CAF = "#FA877F",
  Plasma = "#910A67", 
  Schwann = "#FF9B50",
  Endo = "#A31ACB",
  Macro = "#FF41ED",
  FZ_like = "#765827"
)

col_palette_hr <- c(
  NE1 = "#A6CEE3",
  NE2 = "#1F78B4",
  NE3 = "#864AF9",
  NE4 = "#11009E",
  CAF = "lightgrey",
  # myCAF = "#FA877F",
  # iCAF = "#F94C10",
  # imCAF = "#C70039",
  # vCAF = "#900C3F",
  Plasma = "lightgrey",
  Schwann = "lightgrey",
  Endo = "lightgrey",
  Macro = "lightgrey",
  FZ_like = "lightgrey"
)

# Plot histology
#################
p_histo_ls<- list()
for(s in c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")){
  p_tmp<- SpatialDimPlot(seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]], images = s, crop = T, alpha = 0) + 
    theme(legend.position = "none") +
    NoGrid()
  p_histo_ls[[s]]<- p_tmp
  }

# Plot main clusters
#####################
p_main_ls<- list() 
for(s in c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")){
  seurat_tmp<- seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]]
  # scale factor spots (empirically set)
  scale_f<- c(3.1, 4.1, 3.4, 3.9, 6.1, 5.9)[which(c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")==s)]
  p_tmp<- SpatialDimPlot(seurat_tmp, images = s, crop = T, group.by = "cell_type_lr", image.alpha = 0, pt.size.factor = scale_f, stroke = NA, cols = col_palette_main[levels(seurat_tmp$cell_type_lr)]) + 
    theme(legend.position = "none") +
    NoGrid()
  p_main_ls[[s]]<- p_tmp
  }
# p_main_ls$NB1Pre1
# p_main_ls$NB1Pre2
# p_main_ls$NB1Post1
# p_main_ls$NB1Post2

# Plot hr clusters
#####################
p_hr_ls<- list() 
for(s in c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")){
  seurat_tmp<- seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]]
  # scale factor spots (empirically set)
  scale_f<- c(3.1, 4.1, 3.4, 3.9, 6.1, 5.9)[which(c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")==s)]
  p_tmp<- SpatialDimPlot(seurat_tmp, images = s, crop = T, group.by = "cell_type_hr", image.alpha = 0, pt.size.factor = scale_f, stroke = NA, cols = col_palette_hr[levels(seurat_tmp$cell_type_hr)]) + 
    theme(legend.position = "none") +
    NoGrid()
  p_hr_ls[[s]]<- p_tmp
}
# p_hr_ls$NB1Pre1
# p_hr_ls$NB1Post2
# p_hr_ls$NB2Post2

# Plot MES/ADRN clusters
########################
p_ADRN_MES_ls<- list() 
for(s in c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")){
  seurat_tmp<- seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]]
  # scale factor spots (empirically set)
  scale_f<- c(3.1, 4.1, 3.4, 3.9, 6.1, 5.9)[which(c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")==s)]
  p_tmp<- SpatialFeaturePlot(seurat_tmp, images = s,  
                             features = "log2_ADRN_MES",
                             crop = T,
                             image.alpha = 0, pt.size.factor = scale_f, stroke = NA) + 
    scale_fill_gradient2(
      high= "blue", mid = "white", low = "brown", na.value = "lightgrey",  
      breaks = c(max(seurat_tmp$log2_ADRN_MES, na.rm = TRUE), min(seurat_tmp$log2_ADRN_MES, na.rm = TRUE)),
      labels = c("ADRN", "MES")
    ) +
    ggtitle("") +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.background=element_blank(),
      legend.key.size = unit(0.1, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 1))) +
    NoGrid()
  p_ADRN_MES_ls[[s]]<- p_tmp
}
# p_ADRN_MES_ls$NB1Pre1
# p_ADRN_MES_ls$NB1Post2
# p_ADRN_MES_ls$NB2Post1

# Plot UMAP
############

# Main clusters
p_umap_main_ls<- list() 
for(s in c("NB1Pre", "NB1Post", "NB2Post")){
  seurat_tmp<- seurat_grp_merged_list[[s]]
  p_tmp<- DimPlot(seurat_tmp, group.by = "cell_type_lr", pt.size = .1, cols = col_palette_main[levels(seurat_tmp$cell_type_lr)]) +
    ggtitle("") +
    theme_void() +
    theme(
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.background=element_blank(),
      legend.key.size = unit(0.1, "cm")
      ) +
    guides(color = guide_legend(override.aes = list(size = 1)))
  p_umap_main_ls[[s]]<- p_tmp
}

# hr clusters
p_umap_hr_ls<- list() 
for(s in c("NB1Pre", "NB1Post", "NB2Post")){
  seurat_tmp<- seurat_grp_merged_list[[s]]
  p_tmp<- DimPlot(seurat_tmp, group.by = "cell_type_hr", pt.size = .1, cols = col_palette_hr[levels(seurat_tmp$cell_type_hr)]) +
    ggtitle("") +
    theme_void() +
    theme(
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.background=element_blank(),
      legend.key.size = unit(0.1, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 1)))
  p_umap_hr_ls[[s]]<- p_tmp
}
# p_umap_hr_ls$NB1Pre
# p_umap_hr_ls$NB1Post
# p_umap_hr_ls$NB2Post

# CNV
p_umap_CNV_ls<- list() 
for(s in c("NB1Pre", "NB1Post", "NB2Post")){
  seurat_tmp<- seurat_grp_merged_list[[s]]
  p_tmp<- FeaturePlot(seurat_tmp, features = "inferCNV_cell_type_hr_has_cnv", pt.size = .1, max.cutoff = "q75") +
    scale_color_gradient2(
      name = "CNV score", 
      low = "lightgray", mid = "white", high= "brown", 
      breaks = c(0,5,10,15),
      midpoint = 2) +
    ggtitle("") +
    theme_void() +
    theme(
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.background=element_blank(),
      legend.key.size = unit(0.1, "cm")
    ) 
  # +
  #   guides(color = guide_legend(override.aes = list(size = 1)))
  p_umap_CNV_ls[[s]]<- p_tmp
}

# MES/ADRN
p_umap_ADRN_MES_ls<- list() 
for(s in c("NB1Pre", "NB1Post", "NB2Post")){
  seurat_tmp<- seurat_grp_merged_list[[s]]
  p_tmp<- FeaturePlot(seurat_tmp, features = "log2_ADRN_MES", pt.size = .1) +
    scale_color_gradient2(
      high= "blue", mid = "white", low = "brown", na.value = "lightgrey",  
      breaks = c(max(seurat_tmp$log2_ADRN_MES, na.rm = TRUE), min(seurat_tmp$log2_ADRN_MES, na.rm = TRUE)),
      labels = c("ADRN", "MES")
      ) +
    ggtitle("") +
    theme_void() +
    theme(
      legend.text = element_text(size = 7),
      legend.background=element_blank(),
      legend.key.size = unit(0.1, "cm")
    ) 
  # +
  #   guides(color = guide_legend(override.aes = list(size = 1)))
  p_umap_ADRN_MES_ls[[s]]<- p_tmp
}

# Corr plot MES ADRN
#####################
p_corr_ADRN_MES_ls<- list() 
for(s in c("NB1Pre", "NB1Post", "NB2Post")){
  seurat_tmp<- seurat_grp_merged_list[[s]]@meta.data
  seurat_tmp<- seurat_tmp[seurat_tmp$cell_type_lr=="NE",]
  p_tmp<- ggplot(seurat_tmp, aes(x=MES,y=ADRN)) +
    geom_point(aes(color = cell_type_hr),size = 1, alpha = 0.5) +
    scale_colour_manual(name="", values = col_palette_hr) +
    ggtitle(s) +
    theme_cowplot() +
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 6),
          plot.title = element_text(size = 8, hjust = 0.5, face = "italic"),
          legend.position = c(.95,.95),
          legend.margin = margin(0, 0, 0, 0),
          legend.justification = "right",
          # legend.justification.left = "top",
          # legend.justification.bottom = "right",
          # legend.justification.inside = c(1, 1),
          # legend.location = "plot",          
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.3, "cm")
          )
  p_corr_ADRN_MES_ls[[s]]<- p_tmp
}
# p_corr_ADRN_MES_ls$NB1Pre

# Plot
########
p_umap<- plot_grid(
  plot_grid(plotlist = p_main_ls,ncol = 2),
  plot_grid(plotlist = p_umap_main_ls,ncol = 1),
  plot_grid(plotlist = p_umap_CNV_ls,ncol = 1),
  plot_grid(plotlist = p_umap_ADRN_MES_ls,ncol = 1),
  ncol=4,
  rel_widths = c(3,2,2,2)
)

p_histo<- plot_grid(
  NA, plot_grid(plotlist = p_histo_ls,ncol = 2),
  ncol = 2,
  rel_widths = c(2,1)
)

p<- plot_grid(
  p_histo,
  p_umap,
  ncol= 1 
)

ggsave("results/figs/manuscript_fig1.pdf", width = 178, height = 0.6*265, units = "mm")

# Suppl fig with HR clusters
p_hr<- plot_grid(
  plot_grid(plotlist = p_hr_ls,ncol = 2),
  plot_grid(plotlist = p_umap_hr_ls,ncol = 1),
  ncol=2,
  rel_widths = c(3,2),
  labels = "AUTO",
  label_size = 12
)
ggsave("results/figs/manuscript_figS4_hr.pdf", p_hr, width = 5/9*178, height = 0.3*265, units = "mm")

# Suppl fig with MES/ADRN clusters spatially
p_MES_ADRN<- plot_grid(
  plotlist = p_ADRN_MES_ls,
  ncol=2
)
ggsave("results/figs/manuscript_figS4_MES_ADRN.pdf", p_MES_ADRN, width = 3/7*178, height = 0.6*265, units = "mm")

p_corr_MES_ADRN<- plot_grid(
  plotlist = p_corr_ADRN_MES_ls,
  ncol=3
)
ggsave("results/figs/manuscript_figS4_corr_ MES_ADRN.pdf", p_corr_MES_ADRN, width = 178, height = 0.2*265, units = "mm")

