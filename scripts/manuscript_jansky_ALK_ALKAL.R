library(cowplot)

# Load Jansky seurat objects
#############################
jansky<- readRDS("temp/adrn_med_cortex_sc_Jansky.rds")
jansky_AC_zones<- readRDS("temp/cortex_ori_sc_Jansky.rds")

# Create barplots to show % expression
#######################################

# Get ALK & ALKAL data
ALK_ALKAL<- FetchData(object = jansky, vars = c("ALK","FAM150B","celltype"))

# Stat
tapply(ALK_ALKAL$FAM150B, ALK_ALKAL$celltype=="Adrenal cortex", "mean")
# FALSE       TRUE 
# 0.01001004 0.31180493 
tapply(ALK_ALKAL$FAM150B>0, ALK_ALKAL$celltype=="Adrenal cortex", "mean")
# FALSE       TRUE 
# 0.01545768 0.43327687 
tapply(ALK_ALKAL$ALK, ALK_ALKAL$celltype=="Adrenal cortex", "mean")
# FALSE       TRUE 
# 1.05343242 0.04434088 
tapply(ALK_ALKAL$ALK>0, ALK_ALKAL$celltype=="Adrenal cortex", "mean")
# FALSE      TRUE 
# 0.7874104 0.0608131 
fisher.test(table(ALK_ALKAL$ALK>0, ALK_ALKAL$celltype=="Adrenal cortex"))$p.value #0
fisher.test(table(ALK_ALKAL$FAM150B>0, ALK_ALKAL$celltype=="Adrenal cortex"))$p.value #0

ALK_ALKAL_df<- data.frame(
  ALKAL2 = tapply(ALK_ALKAL$FAM150B>0, ALK_ALKAL$celltype, "mean"),
  ALK = tapply(ALK_ALKAL$ALK>0, ALK_ALKAL$celltype, "mean")
)
ALK_ALKAL_df$celltype <- rownames(ALK_ALKAL_df)
ALK_ALKAL_df<- melt(ALK_ALKAL_df, variable.name = "gene", value.name = "isExpr")

# Plot
ggplot(ALK_ALKAL_df, aes(fill=gene, y=isExpr, x=celltype)) + 
  geom_bar(position="dodge", stat = "identity")

ALK_ALKAL_df<- melt(ALK_ALKAL, variable.name = "gene")
# Plot
ggplot(ALK_ALKAL_df, aes(fill=gene, y=value, x=celltype)) + 
  geom_bar(position="dodge", stat = "summary", fun.y  ="mean")

# Plot ALK & ALKAL2
#####################
DefaultAssay(jansky) <- "SCT"

p_ALK_ALKAL2<- FeaturePlot(jansky, cols = c("lightgrey", "brown", "blue"), order = T, features = c("FAM150B","ALK"), pt.size = 0.001, max.cutoff = "q95", blend=T, alpha = 0.3)
p_sign_ls<- list()
for(sign in c("AP", "FZ", "DZ")){
  p_sign_ls[[sign]]<- FeaturePlot(jansky_AC_zones, order = T, features = sign, pt.size = 0.001, max.cutoff = "q95", min.cutoff = "q05") +
    scale_fill_gradient2(
      high= "blue", mid = "lightgrey", low = "lightgrey", na.value = "lightgrey") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 7),      
      legend.text = element_text(size = 6),
      legend.background=element_blank(),
      legend.key.size = unit(0.2, "cm")
    )
}

p1<- plot_grid(
  p_ALK_ALKAL2[[3]] + theme(axis.text = element_text(size = 6), axis.title = element_text(size = 7, face = "italic"), title = element_text(size=8)),
  plot_grid(
    NA,
    p_ALK_ALKAL2[[4]] + theme(axis.text = element_text(size = 6), axis.title = element_text(size = 7, face = "italic"), title = element_text(size=8)),
    rel_heights = c(3,2),
    ncol = 1
            ), # Only use combine plot
  ncol = 2,
  rel_widths = c(2,1)
)
# Save as png to avoid issue with too many points when loading
ggsave("results/figs/manuscript_jansky_main.pdf", p1, width = 0.5*178, height = 0.25*265, units = "mm")
ggsave("results/figs/manuscript_jansky_main.png", p_ALK_ALKAL2[[3]] + theme_void() + ggtitle("") + theme(legend.position = "none"), width = 0.5*178, height = 0.25*265, units = "mm", dpi = 600)
ggsave("results/figs/manuscript_jansky_main_legend.pdf", p_ALK_ALKAL2[[4]] + theme(axis.text = element_text(size = 6), axis.title = element_text(size = 7, face = "italic"), title = element_text(size=8)), width = 50, height = 40, units = "mm")

p2<- plot_grid(
  p_sign_ls$AP + theme(legend.position = "bottom") ,
  p_sign_ls$FZ + theme(legend.position = "bottom") ,
  p_sign_ls$DZ + theme(legend.position = "bottom") ,
  ncol = 3
)

ggsave("results/figs/manuscript_jansky_main_sign.png", p2, width = 178, height = 0.20*265, units = "mm", dpi = 600)
ggsave("results/figs/manuscript_jansky_main_sign_AP.png", p2, width = 178, height = 0.20*265, units = "mm", dpi = 600)

p<- plot_grid(
  p1,
  p2,
  ncol=1,
  rel_heights = c(2, 1)
)

ggsave("results/figs/manuscript_jansky.pdf", p, width = 178, height = 0.75*265, units = "mm")

p<-  DimPlot(jansky, cells.highlight = WhichCells(jansky, idents = "Neuroblasts"), pt.size = 0.1)
ggsave("results/figs/manuscript_jansky_ref_clusters.png", p, width = 178, height = 0.5*265, units = "mm")

