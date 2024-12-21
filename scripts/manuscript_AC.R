# FZ-like analysis
###################

# Load data
###########
library(Seurat)
seurat_grp_merged_list <- read_rds("data/ST_NB_seurat.rds")
# seurat_grp_merged_list$NB2Pre<- NULL

# Plot Adrenocortical genes
#################################
p_AC_genes<- list()
for(s in c("NB2Post1", "NB2Post2")){
  for(g in c("NR5A1", "CYP11B1", "STAR")){
    seurat_tmp<- seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]]
    scale_f<- c(3.1, 4.1, 3.4, 3.9, 6.1, 5.9)[which(c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")==s)]
    p_tmp<- SpatialFeaturePlot(seurat_tmp, images = s,  
                               features = g,
                               min.cutoff = 1,
                               max.cutoff = 4,
                               crop = T,
                               image.alpha = 0, pt.size.factor = scale_f, stroke = NA) +
      theme(
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.background=element_blank(),
        legend.key.size = unit(0.1, "cm")
      ) +
      guides(color = guide_legend(override.aes = list(size = 1))) +
      NoGrid() 
    p_AC_genes[[s]][[g]]<- p_tmp
  }
}
# p_AC_genes$NB2Post1$STAR  
# p_AC_genes$NB2Post1$STAR  

# Plot umap steroid pathways
############################

p_AC_steroids<- list()
for(s in c("NB2Post1", "NB2Post2")){
  for(g in c("STAR","REACTOME_METABOLISM_OF_STEROIDS", "REACTOME_METABOLISM_OF_STEROID_HORMONES")){
    seurat_tmp<- seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]]
    scale_f<- c(3.1, 4.1, 3.4, 3.9, 6.1, 5.9)[which(c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")==s)]
    p_tmp<- SpatialFeaturePlot(seurat_tmp, images = s,  
                               features = g,
                               min.cutoff = "q05",
                               max.cutoff = "q95",
                               crop = T,
                               image.alpha = 0, pt.size.factor = scale_f, stroke = NA) +
      theme(
        title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.background=element_blank(),
        legend.key.size = unit(0.1, "cm")
      ) +
      guides(color = guide_legend(override.aes = list(size = 1))) +
      NoGrid() 
    p_AC_steroids[[s]][[g]]<- p_tmp
  }
}
# p_AC_steroids$NB2Post1$REACTOME_METABOLISM_OF_STEROIDS
# p_AC_steroids$NB2Post2$REACTOME_METABOLISM_OF_STEROIDS
# p_AC_steroids$NB2Post1$REACTOME_METABOLISM_OF_STEROID_HORMONES
# p_AC_steroids$NB2Post2$REACTOME_METABOLISM_OF_STEROID_HORMONES


# Save
#######
p_steroid<- plot_grid(
  p_AC_steroids$NB2Post2$REACTOME_METABOLISM_OF_STEROIDS,
  p_AC_steroids$NB2Post2$REACTOME_METABOLISM_OF_STEROID_HORMONES,
  ncol=1
)
ggsave("results/figs/manuscript_fig2_steroids.pdf", p_steroid, width = 0.6*178, height = 0.375*265, units = "mm")

p<- plot_grid(
  plot_grid(plotlist = p_AC_genes$NB2Post1, ncol = 3),
  plot_grid(plotlist = p_AC_genes$NB2Post2, ncol = 3),
  ncol=1
)
ggsave("results/figs/manuscript_fig2_AC_genes.pdf", p, width = 0.4*178, height = 0.25*265, units = "mm")

