# Load data
###########
library(Seurat)
library(ggplot2)
library(cowplot)
library(CellChat)

seurat_grp_merged_list <- readRDS("data/ST_NB_seurat.rds")

# Show interaction
###################
col_palette_hr <- c(
  Macro = "lightgrey",
  NE1 = "#A6CEE3",
  NE2 = "#1F78B4",
  NE3 = "lightgrey",
  NE4 = "lightgrey",
  CAF = "lightgrey",
  myCAF = "lightgrey",
  iCAF = "lightgrey",
  imCAF = "lightgrey",
  vCAF = "lightgrey",
  Plasma = "lightgrey", 
  Schwann = "lightgrey",
  Endo = "lightgrey",
  FZ_like = "#765827")

seurat_tmp<- seurat_grp_merged_list$NB2Post
p_ALKAL_ls<- list() 
for(s in c("NB2Post1", "NB2Post2")){
  scale_f<- c(3.1, 4.1, 3.4, 3.9, 5.9, 5.2)[which(c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")==s)]
  p_zoom_macro_FZ<- SpatialDimPlot(seurat_tmp, images = s, crop = T, group.by = "cell_type_hr", image.alpha = 0, pt.size.factor = scale_f, stroke = NA, cols = col_palette_hr[levels(seurat_tmp$cell_type_hr)]) + 
    theme(legend.position = "none") +
    NoGrid()
  p_ALKAL_ls[[s]]<- p_zoom_macro_FZ    
}
# p_ALKAL_ls$NB2Post1
# p_ALKAL_ls$NB2Post2

# Analyse cellchat & Plot circusplot
######################################

cellchat<- readRDS("temp/cellchat_list_ST_sct_truncateMean_0.2_cell_type_hr_secreted_sig_only_20241209.rds")

# Show all pathways
pathways.show <-  cellchat$NB2Post@netP$pathways 

# Circle plot: all pathways from FZ
pdf("results/figs/manuscript_fig4_chord.pdf")
# netVisual_aggregate(cellchat$NB2Post, signaling = pathways.show, signaling.name = "All", targets.use = c("FZ_like","NE1","NE2"), sources.use = c("FZ_like","NE1","NE2"), layout = "chord",remove.isolate = T, color.use = col_palette_hr)
netVisual_aggregate(cellchat$NB2Post, signaling = pathways.show, signaling.name = "All", sources.use = c("FZ_like"), layout = "chord",remove.isolate = T, color.use = col_palette_hr)
dev.off()

# Number?
cellchat$NB2Post@net$count["FZ_like",] 
# NE1     NE2     CAF Schwann FZ_like 
# 86     172     112     116     123 

# Intersect Cell chat top DGE
###############################

# Interactions in list
AC_NE1<- sort(names(sort(cellchat$NB2Post@net$pval["FZ_like","NE1",][cellchat$NB2Post@net$pval["FZ_like","NE1",]<0.05],decreasing = F)))
# ALKAL2_ALK
# NRTN_GFRA2
# FGF9_FGFR1

AC_NE2<- sort(names(sort(cellchat$NB2Post@net$pval["FZ_like","NE2",][cellchat$NB2Post@net$pval["FZ_like","NE2",]<0.05],decreasing = F)))
# ALKAL2_ALK
# NRTN_GFRA2
# FGF9_FGFR1
# NPY_NPY1/5R

intersect(AC_NE1, AC_NE2) # 67
# ALKAL2_ALK
# NRTN_GFRA2
# FGF9_FGFR1

# ALKAL2 expression?
#####################

DE_clust3<- as.data.frame(readxl::read_excel("../ST_NB/ST_NB_Joachim/temp/tables/de_hr_filtered.xlsx",sheet = "NB2Post"))
DE_clust3[DE_clust3$gene=="ALKAL2",]
# p_val avg_log2FC pct.1 pct.2     p_val_adj cluster
# 3.288663e-224   4.857623 0.927 0.302 4.620572e-220 FZ_like

# Plot Genes
#################################
p_genes<- list()
for(s in c("NB2Post1", "NB2Post2")){
  for(g in c("ALKAL2", "ALK", "DLK1", "NOTCH3", "NRTN", "GFRA2", "PNMT", "VTN", "PLAUR", "RET", "LTK")){
    seurat_tmp<- seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]]
    DefaultAssay(seurat_tmp)<- "alra"
    scale_f<- c(5.9, 5.2)[which(c("NB2Post1", "NB2Post2")==s)]
    p_tmp<- SpatialFeaturePlot(seurat_tmp, images = s,  
                               features = g,
                               min.cutoff = "q5",
                               max.cutoff = "q95",
                               crop = T,
                               image.alpha = 0, pt.size.factor = scale_f, stroke = NA) +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.background=element_blank(),
        legend.key.size = unit(0.1, "cm")
      ) +
      guides(color = guide_legend(override.aes = list(size = 1))) +
      NoGrid() 
    p_genes[[s]][[g]]<- p_tmp
  }
}
# p_genes$NB2Post1$ALKAL2
# p_genes$NB2Post2$ALKAL2
# 
# p_genes$NB2Post1$ALK
# p_genes$NB2Post2$ALK
# 
# plot_grid(
#   p_genes$NB2Post1$ALK, p_genes$NB2Post2$ALK,
#   p_genes$NB2Post1$LTK, p_genes$NB2Post2$LTK
# )

# Save
#######
p_AC<- plot_grid(
  p_ALKAL_ls$NB2Post1,
  p_ALKAL_ls$NB2Post2,
  ncol = 1
)
ggsave("results/figs/manuscript_fig_AC.pdf", p_AC, width = 0.3*178, height = 0.25*265, units = "mm")

p_AC_comm<- plot_grid(
  p_genes$NB2Post1$ALKAL2, p_genes$NB2Post1$NRTN, NA,
  p_genes$NB2Post1$ALK, p_genes$NB2Post1$GFRA2, p_genes$NB2Post1$RET,
  p_genes$NB2Post2$ALKAL2, p_genes$NB2Post2$NRTN, NA,
  p_genes$NB2Post2$ALK, p_genes$NB2Post2$GFRA2, p_genes$NB2Post2$RET,
  ncol = 3
)
ggsave("results/figs/manuscript_fig_AC_comm.pdf", p_AC_comm, width = 0.9*178, height = 0.50*265, units = "mm")

p_NRTN_comm<- plot_grid(
  p_genes$NB2Post2$NRTN, NA,
  p_genes$NB2Post2$GFRA2, p_genes$NB2Post2$RET,
  ncol = 2
)
ggsave("results/figs/manuscript_fig_NRTN_comm.pdf", p_NRTN_comm, width = 0.6*178, height = 0.25*265, units = "mm")

p_medulla_cortex<- plot_grid(
  p_genes$NB2Post1$PNMT + theme(legend.position = "none"),
  p_genes$NB2Post2$PNMT + theme(legend.position = "none"),
  p_genes$NB2Post1$DLK1 + theme(legend.position = "none"),
  p_genes$NB2Post2$DLK1 + theme(legend.position = "none"),
  p_genes$NB2Post2$DLK1, # for legend
  ncol = 5
)
ggsave("results/figs/NE2_medulla_cortex.pdf", p_medulla_cortex, width = 1.2*178, height = 0.2*265, units = "mm")

p_NRTN<- plot_grid(
  p_genes$NB2Post1$NRTN, p_genes$NB2Post1$GFRA2,
  p_genes$NB2Post2$NRTN, p_genes$NB2Post2$GFRA2,
  ncol = 2
)
ggsave("results/figs/manuscript_fig_NRTN.pdf", p_NRTN, width = 0.6*178, height = 0.5*265, units = "mm")

# Save cellchat results as table
####################################
# Create dataframe
cc_df_ls<- list() 
for(s in c("NB1Pre", "NB1Post", "NB2Post")) cc_df_ls[[s]]<- netVisual_bubble(cellchat[[s]], thresh = 0.05, return.data = TRUE)$communication[,c("source", "target", "ligand", "receptor", "prob", "pathway_name", "evidence")]

cc_df_ls$caption<- as.data.frame(
  cbind(
    c("Table 3. Cell-cell communication results.",
      "Cell-cell communication analysis was performed using CellChat (secreted signalling interactions). Only significant interactions (P<0.05) listed",
      "Tabnames refer to the tumor sample",
      "",
      "source",
      "target",
      "ligand",
      "receptor",
      "prob",
      "pathway_name",
      "evidence"),
    c("",
      "",
      "",
      "",
      "Cluster that expresses the source, i.e. ligand",
      "Cluster that expresses the target, i.e. ligand receptor",
      "Gene encoding the ligand",
      "Gene encoding the receptor",
      "Probability of the predicted interaction",
      "Pathway the ligand & receptor are active in",      
      "Literature evidence of the ligand - receptor interaction"
    )
  )
)

WriteXLS::WriteXLS(cc_df_ls[c("caption","NB1Pre","NB1Post","NB2Post")],"../ST_NB/results/tables/manuscript_tableS3.xlsx",row.names = F)
save.image("../ST_NB/results/data/manuscript_NB2_interactions.RData")

