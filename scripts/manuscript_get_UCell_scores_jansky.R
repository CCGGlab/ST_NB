# Run in Conda Jochim
###############

# Libraries
library(dplyr)
library(tidyverse)
library(patchwork)
library(UCell)
library(Seurat)
library(SeuratExtend)
library(ClusterFoldSimilarity)

###############
# Jansky adrenal fetal signatures

adrn_gland_sc <- read_rds("/home/joachim/projects/NB_spatialT/temp/data/adrn_gland_sub_Jansky.rds")
seurat_grp_merged_list <- read_rds("data/ST_NB.RData")

#...2. Subset SCP and Sympathoadrenal cells----
cells_select <- c("Neuroblasts", "late Neuroblasts")

Idents(adrn_gland_sc) <- "celltype"
# Find DE marker
de_markers <- FindAllMarkers(adrn_gland_sc, only.pos = TRUE, min.pct = 0.1)

jansky_NBsigs<- list()
for(c in cells_select){
  jansky_NBsigs[[c]] <- dplyr::filter(de_markers, cluster == c) %>%
    dplyr::filter(!grepl("^RPS|^RPL|^RP|^LINC", gene)) %>%
    dplyr::arrange(.,p_val_adj) %>%
    dplyr::filter(avg_log2FC>= 1) %>%
    slice_head(n = 50) %>%
    pull(gene)
}

# Run signature score on alra assay
seurat_grp_merged_list <- map(seurat_grp_merged_list, function(so){
  DefaultAssay(so) <- "alra"
  so <- UCell::AddModuleScore_UCell(so, features = jansky_NBsigs)
  return(so)
} )

# Save as dataframe
jansky_df<- list()
for(s in names(seurat_grp_merged_list)){
  jansky_df[[s]]<- data.frame(
    sample = s,
    barcode = rownames(seurat_grp_merged_list[[s]][[paste0(names(jansky_NBsigs), "_UCell")]]),
    seurat_grp_merged_list[[s]][[paste0(names(jansky_NBsigs), "_UCell")]]
  )
}

saveRDS(jansky_df, "temp/jansky_UCell_scores_list.rds")
