############
# jansky_NB1Post_GSEA
###############

# Libraries
library(dplyr)
library(tidyverse)
library(patchwork)

# Merge seurat with UCell scores
#############################################
seurat<- readRDS("data/ST_NB_seurat.rds")
UC_Scores<- readRDS("temp/jansky_UCell_scores_list.rds")

for(s in names(seurat)){
  seurat[[s]]<- AddMetaData(seurat[[s]], UC_Scores[[s]], )  
}

s<- "NB1Post2"
p_ls<- list()
for(gs in c("Neuroblasts", "late Neuroblasts")){
  seurat_tmp<- seurat$NB1Post
  DefaultAssay(seurat_tmp)<- "alra"
  scale_f<- 3.9
  p_ls[[gs]]<- SpatialFeaturePlot(seurat_tmp, images = s,  
                             features = gsub(" ","\\.",paste0(gs,"_UCell")),
                             min.cutoff = "q5",
                             # max.cutoff = 4,
                             crop = T,
                             image.alpha = 0, pt.size.factor = scale_f, stroke = NA) +
    theme(
      legend.position = "none"
      # legend.title = element_text(size = 7),
      # legend.text = element_text(size = 6),
      # legend.background=element_blank(),
      # legend.key.size = unit(0.1, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 1))) +
    NoGrid() 
}


ggsave("results/figs/fig5E_jansky_lateNB.png", p_ls$`late Neuroblasts`, width = 178/2, height = 297/4, units = "mm")
ggsave("results/figs/fig5E_jansky_NB.png",  p_ls$Neuroblasts, width = 178/2, height = 297/4, units = "mm")

