# Load data
###########
seurat_grp_merged_list <- read_rds("data/ST_NB_seurat.rds")

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
      legend.position = "none"
      # legend.title = element_blank(),
      # legend.text = element_text(size = 6),
      # legend.background=element_blank(),
      # legend.key.size = unit(0.1, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 1))) +
    NoGrid()
  p_ADRN_MES_ls[[s]]<- p_tmp
}
# p_ADRN_MES_ls$NB1Pre1
# p_ADRN_MES_ls$NB1Post2
# p_ADRN_MES_ls$NB2Post1

# Plot CNV
###########
p_CNV_ls<- list()
for(s in c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")){
  seurat_tmp<- seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]]
  # scale factor spots (empirically set)
  scale_f<- c(3.1, 4.1, 3.4, 3.9, 6.1, 5.9)[which(c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")==s)]
  p_tmp<- SpatialFeaturePlot(seurat_tmp, images = s,  
                             features = "inferCNV_cell_type_hr_has_cnv",
                             crop = T,
                             image.alpha = 0, pt.size.factor = scale_f, stroke = NA, max.cutoff = "q75") + 
    scale_fill_gradient2(
      name = "CNV score", 
      low = "lightgray", mid = "white", high= "brown", 
      breaks = c(0,5, 10, 15),
      midpoint = 2) +
    ggtitle("") +
    theme(
      legend.position = "none"
    )
  p_CNV_ls[[s]]<- p_tmp
}

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
      legend.position = "none"
      # legend.text = element_text(size = 7),
      # legend.background=element_blank(),
      # legend.key.size = unit(0.1, "cm")
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
    theme(
      legend.position = "none",
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      plot.title = element_text(size = 8, hjust = 0.5, face = "italic")
      # legend.position = c(.65,.95),
      # legend.margin = margin(0, 0, 0, 0),
      # legend.justification = "left",
      # legend.text = element_text(size = 7),
      # legend.key.size = unit(0.3, "cm")
      )
  p_corr_ADRN_MES_ls[[s]]<- p_tmp
}
# p_corr_ADRN_MES_ls$NB1Pre

# Plot NBAtlas NE deconvolution 
##############################

# Seurat default colors
SpatialColors <- grDevices::colorRampPalette(
  colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
)

# Plot
p_UMAP_NE_ls<- list()
for(s in names(seurat_grp_merged_list)){
  p_UMAP_NE_ls[[s]]<- FeaturePlot(seurat_grp_merged_list[[s]], features = "RCTD_NBAtlas_full_Neuroendocrine", pt.size = 0.1) +
    scale_color_gradientn(oob = scales::squish, # Change cut-offs without graying out
                          name = NULL,
                          limits = c(0, 1), 
                          colours = SpatialColors(n = 100)
    ) +
    ggtitle("") +
    theme_void() +
    theme(
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.background=element_blank(),
      legend.key.size = unit(0.1, "cm")
    ) 
}

# Spatial plot
p_NE_ls<- list()
for(s in c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")){
  seurat_tmp<- seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]]
  # scale factor spots (empirically set)
  scale_f<- c(3.1, 4.1, 3.4, 3.9, 6.1, 5.9)[which(c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")==s)]
  p_tmp<- SpatialFeaturePlot(seurat_tmp, images = s,  
                             features = "RCTD_NBAtlas_full_Neuroendocrine",
                             crop = T,
                             image.alpha = 0, pt.size.factor = scale_f, stroke = NA) + 
    scale_fill_gradientn(
      name = NULL,
      limits = c(0, 1),
      colours = SpatialColors(n = 100)
    ) +
    ggtitle("") +
    theme(
      legend.position = "none"
    )
  p_NE_ls[[s]]<- p_tmp
}

# Rotate NB2Post1
################
library(ggplot2)
library(cowplot)
library(grid)

rotate_image <- function(p, rot_angle) {
  gt <- ggplot_gtable(ggplot_build(p))
  panel_idx <- which(gt$layout$name == "panel")
  rot_vp <- viewport(angle = rot_angle)
  gt[["grobs"]][[panel_idx]] <- editGrob(gt[["grobs"]][[panel_idx]], vp = rot_vp)
  p_rot <- ggdraw() + draw_grob(gt)
  
  return(p_rot)
}
p_NE_ls$NB2Post1<- rotate_image(p_NE_ls$NB2Post1, 90)
p_CNV_ls$NB2Post1<- rotate_image(p_CNV_ls$NB2Post1, 90)
p_ADRN_MES_ls$NB2Post1<- rotate_image(p_ADRN_MES_ls$NB2Post1, 90)

# Plot
########

# Histology
p_histo<- plot_grid(
  NA, plot_grid(plotlist = p_histo_ls,ncol = 2),
  ncol = 2,
  rel_widths = c(2,1)
)

ggsave("results/figs/manuscript_fig1B_histo.pdf", p_histo, width = 178, height = 0.6*265, units = "mm")

# UMAPs as png for Fig. 1E
p_umap2_NB1Pre<- plot_grid(
  p_umap_main_ls$NB1Pre + theme(legend.position = "none"),
  p_umap_CNV_ls$NB1Pre + theme(legend.position = "none"),
  p_UMAP_NE_ls$NB1Pre + theme(legend.position = "none"),
  ncol=3
)
p_umap2_NB1Post<- plot_grid(
  p_umap_main_ls$NB1Post + theme(legend.position = "none"),
  p_umap_CNV_ls$NB1Post + theme(legend.position = "none"),
  p_UMAP_NE_ls$NB1Post + theme(legend.position = "none"),
  ncol=3
)
p_umap2_NB2Post<- plot_grid(
  p_umap_main_ls$NB2Post + theme(legend.position = "none"),
  p_umap_CNV_ls$NB2Post + theme(legend.position = "none"),
  p_UMAP_NE_ls$NB2Post + theme(legend.position = "none"),
  ncol=3
)
ggsave("results/figs/manuscript_fig1E_umap_NB1Pre.png", plot = p_umap2_NB1Pre, width = 178/2, height = 265/9, scale = 1.5, units = "mm", dpi = 600)
ggsave("results/figs/manuscript_fig1E_umap_NB1Post.png", plot = p_umap2_NB1Post, width = 178/2, height = 265/9, scale = 1.5, units = "mm", dpi = 600)
ggsave("results/figs/manuscript_fig1E_umap_NB2Post.png", plot = p_umap2_NB2Post, width = 178/2, height = 265/9, scale = 1.5, units = "mm", dpi = 600)
ggsave("results/figs/manuscript_fig1E_umap_NE_legend.pdf", plot = p_UMAP_NE_ls$NB2Post, width = 178/6, height = 265/9, scale = 1.5, units = "mm", dpi = 600)

# Suppl fig with NE
p_NE<- plot_grid(
  plotlist = p_NE_ls,
  ncol=2
)
ggsave("results/figs/manuscript_figS3_NE.png", p_NE, width = 178/3, height = 0.3*265, units = "mm", scale = 2.5, dpi = 600, bg="white")

# Suppl fig with CNV
p_CNV<- plot_grid(
  plotlist = p_CNV_ls,
  ncol=2
)
ggsave("results/figs/manuscript_figS3_CNV.png", p_CNV, width = 178/3, height = 0.3*265, units = "mm", scale = 2.5, dpi = 600, bg="white")

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

# Suppl fig with MES/ADRN clusters 
p<- plot_grid(
  plot_grid(plotlist = p_ADRN_MES_ls,ncol = 2), NA, plot_grid(plotlist = p_umap_ADRN_MES_ls,ncol = 1), NA, plot_grid(plotlist = p_corr_ADRN_MES_ls, ncol=1), 
  ncol =5, 
  rel_widths = c(3,0.2, 2, 0.2, 2))
ggsave("results/figs/manuscript_figS5_MES_ADRN.png", p, width = 178, height = 0.5*265, units = "mm")



