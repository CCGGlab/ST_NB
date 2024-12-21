# Load data
###########
library(Seurat)
library(ggplot2)
library(cowplot)
library(CellChat)
# seurat_grp_merged_list <- readRDS("ST_NB_Joachim/temp/data/seurat_grp_merged_list.rds")
# seurat_grp_merged_list$NB2Pre<- NULL

# Color palette
###############
col_palette_hr <- c(
  NE4 = "#11009E",
  Macro = "#FF41ED",
  NE1 = "lightgrey",
  NE2 = "lightgrey",
  NE3 = "lightgrey",
  CAF = "lightgrey",
  myCAF = "lightgrey",
  iCAF = "lightgrey",
  imCAF = "lightgrey",
  vCAF = "lightgrey",
  Plasma = "lightgrey", 
  Schwann = "lightgrey",
  Endo = "lightgrey",
  FZ_like = "lightgrey")


# Analyse cellchat & Plot circusplot
######################################

cellchat<- readRDS("ST_NB_Joachim/temp/data/cellchat_list_ST_sct_truncateMean_0.2_sec_sig_only_cell_type_hr.rds")

# Show all pathways
pathways.show <-  cellchat$NB1Post@netP$pathways 

# Circle plot: all pathways to NE4
par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title

pdf("results/figs/manuscript_fig3_chord.pdf")
# netVisual_aggregate(cellchat$NB1Post, signaling = pathways.show, signaling.name = "All", targets.use = "NE4",layout = "circle",remove.isolate = T)
netVisual_aggregate(cellchat$NB1Post, signaling = pathways.show, signaling.name = "All", targets.use = "NE4",layout = "chord",remove.isolate = T, color.use = col_palette_hr)
dev.off()

# Number?
cellchat$NB1Post@net$count[, "NE4"] 
# NE1    NE2    NE3    NE4    CAF Plasma   Endo  Macro 
# 13      0     35     46     55      0      0     44 

sort(cellchat$NB1Post@net$pval["Macro","NE4",],decreasing = F)[1:44]

cellchat$NB1Post@net$prob["Macro","NE4","FGF7_FGFR1"] # 0.007
cellchat$NB1Post@net$prob["Macro","NE4","CCL18_PITPNM3"] #  0.0008685348
cellchat$NB1Post@net$prob["Macro","NE4","IGF2_IGF2R"] #  0.0007526152

# Show interaction
###################
seurat_tmp<- seurat_grp_merged_list$NB1Post
scale_f<- 3.9
p_zoom_macro_NE4<- SpatialDimPlot(seurat_tmp, images = "NB1Post2", crop = T, group.by = "cell_type_hr", image.alpha = 0, pt.size.factor = scale_f, stroke = NA, cols = col_palette_hr[levels(seurat_tmp$cell_type_hr)]) + 
  theme(legend.position = "none") +
  NoGrid()

# Barplot top 10 DGE in NE4
###########################
DE_clust2<- readxl::read_excel("../ST_NB/ST_NB_Joachim/temp/tables/de_markers_ST_seurat_clust.xlsx",sheet = "NB1Post")
DE_clust2<- DE_clust2[DE_clust2$cluster=="C7",] # C7 ~ NE4
DE_clust2<- DE_clust2[1:10,] # Top 10
DE_clust2$gene<- factor(DE_clust2$gene, levels = rev(DE_clust2$gene))
p_bp<-ggplot(data=DE_clust2, aes(x=gene, y=-log10(p_val_adj))) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  ylab("-log10(Padj") +
  xlab("") +
  theme_cowplot() + 
  theme(
    axis.text.y = element_text(size = 7, face = "italic"), 
    axis.text.x = element_text(size = 6), 
    axis.title.x =  element_text(size = 7), 
    
  )

# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# Chromosomes?
# DUSP2 2q
# CALCB 11p 
# VGF 7q  
# NKX6-1 4q
# CRH 8q   
# CALCA 11p 
# CLU 8p   
# SCG2 2q   
# ARC 8q --> cyotskeleton, migration   
# WEE1 11p   

# Barplot top 10 DGE Ligands in Macro
######################################

Lig<- readRDS(file = "downloads/HGNC/Lig.rds") # All "Receptor ligands" from HGNC
DE_clust2<- readxl::read_excel("../ST_NB/ST_NB_Joachim/temp/tables/de_hr_filtered.xlsx",sheet = "NB1Post")
DE_clust2<- DE_clust2[DE_clust2$cluster=="Macro"&DE_clust2$gene%in%Lig$Approved.symbol,]
DE_clust2<- DE_clust2[1:10,] # Top 10
DE_clust2$gene<- factor(DE_clust2$gene, levels = rev(DE_clust2$gene))
p_bp_macro<-ggplot(data=DE_clust2, aes(x=gene, y=-log10(p_val_adj))) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  ylab("-log10(Padj") +
  xlab("") +
  theme_cowplot() + 
  theme(
    axis.text.y = element_text(size = 7, face = "italic"), 
    axis.text.x = element_text(size = 6), 
    axis.title.x =  element_text(size = 7), 
    
  )
p_bp_macro

# Plot Genes
#################################
p_genes<- list()
for(s in c("NB1Post1", "NB1Post2")){
  for(g in unique(c("CCL18", "PITPNM3", "FGF1", "FGF7", "FGFR1", "CD68", "TREM2", "NTRK1", "ABCB1", "THBS1","ITGAV","SDC4","IL1B","IL1R1","IL6","IL6R","WEE1","DUSP2","CRH","ARC","CALCA","CALCB","NKX6-1","NTRK1","TH"))){
    seurat_tmp<- seurat_grp_merged_list[[substr(s, 1, nchar(s)-1)]]
    DefaultAssay(seurat_tmp)<- "alra"
    scale_f<- c(3.1, 4.1, 3.4, 3.9, 6.1, 5.9)[which(c("NB1Pre1", "NB1Pre2","NB1Post1", "NB1Post2", "NB2Post1", "NB2Post2")==s)]
    p_tmp<- SpatialFeaturePlot(seurat_tmp, images = s,  
                               features = g,
                               min.cutoff = "q5",
                               # max.cutoff = 4,
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
# p_genes$NB1Post1$CCL18
# p_genes$NB1Post1$PITPNM3
# p_genes$NB1Post1$CD68
# p_genes$NB1Post1$NTRK1
# p_genes$NB1Post1$WEE1
# p_genes$NB1Post1$ABCB1
# p_genes$NB1Post1$TH
# p_genes$NB1Post1$CALCA
# p_genes$NB1Post1$CALCB
# p_genes$NB1Post1$CRH
# p_genes$NB1Post1$ARC


# 
# p_genes$NB1Post2$CCL18
# p_genes$NB1Post2$PITPNM3
# p_genes$NB1Post2$CD68
# p_genes$NB1Post2$NTRK1
# p_genes$NB1Post2$WEE1
# p_genes$NB1Post2$ABCB1
# p_genes$NB1Post2$`NKX6-1`
# p_genes$NB1Post2$TH
# p_genes$NB1Post2$DUSP2
# p_genes$NB1Post2$VGF
# p_genes$NB1Post2$ARC
# p_genes$NB1Post2$EGR4

# Save
#######
p_chr<- plot_grid(
  p_genes$NB1Post2$CRH, p_genes$NB1Post2$ARC,
  p_genes$NB1Post2$CALCB, p_genes$NB1Post2$WEE1,
  ncol = 2
)
ggsave("results/figs/manuscript_fig3_chr_8_11_genes.pdf", p_chr, width = 0.6*178, height = 0.15*265, units = "mm")

# p_diff<- plot_grid(
#   p_genes$NB1Post2$NTRK1, 
#   p_genes$NB1Post2$`NKX6-1`,
#   ncol = 1
# )
# ggsave("results/figs/manuscript_fig3_diff_genes.pdf", p_diff, width = 0.5*178, height = 0.25*265, units = "mm")

p_diff<- plot_grid(
  p_genes$NB1Post2$WEE1, p_genes$NB1Post2$`NKX6-1`, p_genes$NB1Post2$NTRK1,
  ncol = 1
)
ggsave("results/figs/manuscript_fig5_diff_genes.pdf", p_diff, width = 0.5*178, height = 0.25*265, units = "mm")

p_CCL18<- plot_grid(
  p_genes$NB1Post2$CCL18, p_genes$NB1Post2$FGF1,
  p_genes$NB1Post2$PITPNM3, p_genes$NB1Post2$FGFR1,
  p_genes$NB1Post2$THBS1, p_genes$NB1Post2$FGF7,
  p_genes$NB1Post2$ITGAV, p_genes$NB1Post2$SDC4,
  ncol = 2
)
ggsave("results/figs/manuscript_fig3_CCL18.pdf", p_CCL18, width = 178, height = 0.5*265, units = "mm")

p<- plot_grid(
  plot_grid(plotlist = p_genes$NB1Post1, ncol = 4),
  plot_grid(plotlist = p_genes$NB1Post2, ncol = 4),
  ncol=1
)
ggsave("results/figs/manuscript_fig3_CCL18_genes.pdf", p, width = 178, height = 265, units = "mm")

p<- plot_grid(
  p_zoom_macro_NE4,
  p_bp,
  p_bp_macro,
  ncol = 1
)
ggsave("results/figs/manuscript_fig3_zoom_barplot.pdf", p, width = 0.3*178, height = 0.45*265, units = "mm")


