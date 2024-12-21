# Show clustering tree
####################

library(clustree)
seurat_grp_merged_list <- read_rds("data/ST_NB_seurat.rds")
# seurat_grp_merged_list$NB2Pre<- NULL

resolutions <- seq(0.1, 1.5, 0.1)
seurat_obj_ls <- list()

for(s in names(seurat_grp_merged_list)){
  seurat_obj <- seurat_grp_merged_list[[s]]
  for (r in resolutions){
    seurat_obj <- FindClusters(seurat_obj, resolution = r, 
                               cluster.name = paste0("clust_res.", r),
                               dims = 1:15, min.cells = 100)
    }
  seurat_obj_ls[[s]] <- seurat_obj  
}

label_position <- function(labels) {
  if (length(unique(labels)) == 1) {
    position <- as.character(unique(labels))
  } else {
    position <- "mixed"
  }
  return(position)
}

s_res_cu<- c(0.6, 1, 0.7)
clustree_ls<- list()
for(i in 1:length(seurat_obj_ls)){
  cat(i, " ")
  s<- names(seurat_obj_ls)[i]
  seurat_obj<- seurat_obj_ls[[s]]
  seurat_obj[[paste0("clust_res.",seq(s_res_cu[i] + 0.1,1.5,0.1))]]<- NULL
  clustree_ls[[s]][["noLabels"]]<- clustree(seurat_obj, prefix = "clust_res.", node_size_range = c(1, 5), edge_width = 0.5, node_text_size = 2) + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), , legend.key.size = unit(6, 'pt'))
  clustree_ls[[s]][["withLabels"]]<- clustree(seurat_obj, prefix = "clust_res.", node_label = "cell_type_hr",
           node_label_aggr = "label_position", , node_size_range = c(1, 5), edge_width = 0.5, node_text_size = 2)
  
  clustree_ls[[s]][["umap"]]<- DimPlot(seurat_grp_merged_list[[s]], group.by = "seurat_clusters") +
    ggtitle("") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(
      # plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size=7),
      axis.text = element_text(size=6),
      legend.text = element_text(size=6), 
      legend.title = element_text(size=6), 
      legend.key.size = unit(6, 'pt'),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 2)))
}

# clustree_ls$NB1Pre$umap
# clustree_ls$NB1Pre$noLabels
# clustree_ls$NB1Pre$withLabels

p<- plot_grid(
  clustree_ls$NB1Pre$noLabels, clustree_ls$NB1Post$noLabels, clustree_ls$NB2Post$noLabels,
  clustree_ls$NB1Pre$umap, clustree_ls$NB1Post$umap, clustree_ls$NB2Post$umap,
  ncol = 3,
  rel_heights = c(1,1)
  )

ggsave("results/figs/manuscript_basic_clustering.pdf", width = 178, height = 0.6*265, units = "mm")

# clustree_ls$NB1Pre$withLabels
# clustree_ls$NB1Post$withLabels
# clustree_ls$NB2Post$withLabels

