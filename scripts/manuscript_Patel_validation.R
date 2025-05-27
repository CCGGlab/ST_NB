# Patel validation
##################

library(cowplot)
library(Seurat)
library(reshape)

seu_obj<- readRDS( "ST_NB_Joachim/temp/data/seu_obj_LIGER_patel.rds") # Seuraty object provided by the authors

# Visualize main clusters
##########################
dimplot <- DimPlot(seu_obj, group.by = "SingleR.cluster.labels")  +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        title = element_text(size=8),
        legend.key.height = unit(3, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.text = element_text(size = 6)
  )

dimplot2 <- DimPlot(seu_obj, group.by = "annotate_refine_fine")  +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        title = element_text(size=8),
        legend.key.height = unit(3, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.text = element_text(size = 6)
  )
# Main clusters
CNVplot <- FeaturePlot(seu_obj, features = "cnv")  +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        title = element_text(size=8),
        legend.key.height = unit(3, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.text = element_text(size = 6)
  )

# CNV
p<- plot_grid(
  dimplot,
  dimplot2,
  CNVplot,
  ncol = 2
)
ggsave("results/figs/manuscript_Patel_ref.png", p, width = 178, height = 0.5*265, units = "mm", dpi = 600)


# 1. Ligand receptors
####################

for(genes in list(
  c("ALKAL2", "rna_ALK"),
  c("NRTN", "GFRA2"),
  c("NRTN", "RET"),
  c("FGF9", "FGFR1"),
  c("CCL18", "PITPNM3"),
  c("FGF7", "FGFR1")
)){
  p_tmp<- FeaturePlot(seu_obj, cols = c("lightgrey", "brown", "blue"), order = T, features = genes, pt.size = 0.001, max.cutoff = "q95", blend=T, alpha = 0.3, raster = F)
  ggsave(paste0("results/figs/manuscript_Patel_",paste0(genes, collapse = "_"),".png"), p_tmp[[3]] + theme_void() + ggtitle("") + theme(legend.position = "none"), width = 0.5*178, height = 0.25*265, units = "mm", dpi = 600)
  ggsave(paste0("results/figs/manuscript_Patel_",paste0(genes, collapse = "_"),"_legend.pdf"), p_tmp[[4]] + theme(axis.text = element_text(size = 6), axis.title = element_text(size = 7, face = "italic"), title = element_text(size=8)), width = 50, height = 40, units = "mm")
  # Sometimes difficult to seperate
  ggsave(paste0("results/figs/manuscript_Patel_",genes[1],".png"), p_tmp[[1]] + theme_void() + ggtitle("") + theme(legend.position = "none"), width = 0.5*178, height = 0.25*265, units = "mm", dpi = 600)
  ggsave(paste0("results/figs/manuscript_Patel_",genes[2],".png"), p_tmp[[2]] + theme_void() + ggtitle("") + theme(legend.position = "none"), width = 0.5*178, height = 0.25*265, units = "mm", dpi = 600)
}

# Treatment, MYCN effects
###########################

# Treatment status
seu_obj@meta.data <- seu_obj@meta.data %>% 
  dplyr::mutate(treatment_status_coarse = case_when(Treatment_Status == "treated, chemotherapy" ~ "treated",
                                                    Treatment_Status == "treated, chemo-immunotherapy" ~ "treated",
                                                    TRUE ~ as.character(Treatment_Status)
                                                    
  ))
seu_obj$treatment_status_coarse <- factor(seu_obj$treatment_status_coarse, 
                                          levels = c("untreated", "treated", "unknown"))

CCL18_expr<-  seu_obj[["RNA"]]$counts["CCL18",]
PITPNM3_expr<-  seu_obj[["RNA"]]$counts["PITPNM3",]

CCL18_df<- data.frame(
  gene = rep(c("CCL18", "PITPNM3"),each = length(CCL18_expr)),
  expr = as.numeric(c(CCL18_expr, PITPNM3_expr)),
  # expr = as.numeric(c(CCL18_expr, PITPNM3_expr)>0),
  MYCN = rep(seu_obj$MYCN,2),
  treatment = rep(seu_obj$treatment_status_coarse,2),
  macro_cluster = rep(seu_obj$SingleR.cluster.labels=="Macrophage", 2),
  neuron_cluster = rep(seu_obj$SingleR.cluster.labels=="Neurons", 2),
  adrenergic_cluster = rep(seu_obj$annotate_refine_fine=="adrenergic", 2)
)
CCL18_df$cluster<- NA
CCL18_df$cluster[CCL18_df$macro_cluster==T]<- "Macrophage"
CCL18_df$cluster[CCL18_df$adrenergic_cluster==T]<- "Adrenergic"
CCL18_df$cluster[CCL18_df$macro_cluster==F&CCL18_df$adrenergic_cluster==F]<- "Other"
CCL18_df$cluster<- factor(CCL18_df$cluster, c("Adrenergic", "Macrophage", "Other"))

p_bar_ls<- list()
for(c in c("Adrenergic","Macrophage")){
  p_bar_ls[[c]]<- ggplot(subset(CCL18_df,cluster==c & treatment!="unknown"), aes(fill=treatment, y=expr, x=MYCN)) + 
    geom_bar(position="dodge", stat="summary", fun = "mean") +
    facet_grid(gene~treatment) +
    theme_pubr() +  # use a white background
    # coord_cartesian(ylim=c(-0.5, 5)) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size=8, face = "italic"),
      axis.text = element_text(size=7),
      axis.title = element_text(size=7),
      plot.title = element_text(size = 7)
    )
}
# p_bar_ls$Macrophage
p_bar<- plot_grid(
  plotlist = p_bar_ls,
  ncol = 1
)
ggsave(paste0("results/figs/manuscript_Patel_CCL18_PITPNM3_covar_bar.pdf"), p_bar, width = 0.25*178, height = 0.4*265, units = "mm", dpi = 600, bg = "white")

# some numbers
t.test(CCL18_df$expr[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$treatment=="treated"]~CCL18_df$MYCN[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$treatment=="treated"])$p.value #4.95e-86
t.test(CCL18_df$expr[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$treatment=="untreated"]~CCL18_df$MYCN[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$treatment=="untreated"])$p.value # 0.02
t.test(CCL18_df$expr[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$MYCN=="non-amplified"&CCL18_df$treatment!="unknown"]~CCL18_df$treatment[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$MYCN=="non-amplified"&CCL18_df$treatment!="unknown"])$p.value #5.25e-81

t.test(CCL18_df$expr[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="Adrenergic"&CCL18_df$treatment=="treated"]~CCL18_df$MYCN[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="Adrenergic"&CCL18_df$treatment=="treated"])$p.value #3.008636e-14
t.test(CCL18_df$expr[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="Adrenergic"&CCL18_df$treatment=="untreated"]~CCL18_df$MYCN[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="Adrenergic"&CCL18_df$treatment=="untreated"])$p.value # 6.71967e-101
t.test(CCL18_df$expr[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="Adrenergic"&CCL18_df$MYCN=="non-amplified"&CCL18_df$treatment!="unknown"]~CCL18_df$treatment[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="Adrenergic"&CCL18_df$MYCN=="non-amplified"&CCL18_df$treatment!="unknown"])$p.value #0.00022

# Visualize UMAP
so_treated<- subset(x = seu_obj, subset = treatment_status_coarse == "treated")
p_treated<- FeaturePlot(so_treated, cols = c("lightgrey", "brown", "blue"), order = T, features = genes, pt.size = 0.001, max.cutoff = "q95", blend=T, alpha = 0.3, raster = F, split.by = "MYCN")
so_treated<- subset(x = seu_obj, subset = treatment_status_coarse == "untreated")
p_untreated<- FeaturePlot(so_treated, cols = c("lightgrey", "brown", "blue"), order = T, features = genes, pt.size = 0.001, max.cutoff = "q95", blend=T, alpha = 0.3, raster = F, split.by = "MYCN")

p_UMAP<- plot_grid(
  p_untreated[[3]] + theme_void() + ggtitle("") + theme(legend.position = "none"), p_treated[[3]]+ theme_void() + ggtitle("") + theme(legend.position = "none"), # Amp
  p_untreated[[7]]+ theme_void() + ggtitle("") + theme(legend.position = "none"), p_treated[[7]]+ theme_void() + ggtitle("") + theme(legend.position = "none"), # Non amp
  ncol = 2
)
ggsave(paste0("results/figs/manuscript_Patel_CCL18_PITPNM3_covar.png"), p_UMAP, width = 0.5*178, height = 0.3*265, units = "mm", dpi = 600, bg = "white")

# # Explore others
# WEE1_expr<-  seu_obj[["RNA"]]$counts["WEE1",]
# tapply(WEE1_expr, seu_obj$SingleR.cluster.labels, "mean")
# tapply(WEE1_expr, list(seu_obj$SingleR.cluster.labels, seu_obj$treatment_status_coarse), "mean")
# cor.test(WEE1_expr[seu_obj$treatment_status_coarse=="treated"&seu_obj$MYCN=="non-amplified"], CCL18_expr[seu_obj$treatment_status_coarse=="treated"&seu_obj$MYCN=="non-amplified"])

# AC like
##########

p_AC<- FeaturePlot(seu_obj, features = "AC_like",order = T, pt.size = 0.001, min.cutoff = "q05", raster=F) +
  scale_fill_gradient2(
    high= "blue", mid = "lightgrey", low = "lightgrey", na.value = "lightgrey") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 7),      
    legend.text = element_text(size = 6),
    legend.background=element_blank(),
    legend.key.size = unit(0.2, "cm")
  )

ggsave("results/figs/manuscript_Patel_AC.png", p_AC + theme_void() + ggtitle("") + theme(legend.position = "none"), width = 0.5*178, height = 0.25*265, units = "mm", dpi = 600)
ggsave("results/figs/manuscript_Patel_AC_legend.pdf", p_AC + theme(axis.text = element_text(size = 6), axis.title = element_text(size = 7, face = "italic"), title = element_text(size=8)), width = 50, height = 40, units = "mm")

# Overlap ALKAL2-AC
ALKAL2_expr<-  seu_obj[["RNA"]]$counts["ALKAL2",]
NRTN_expr<-  seu_obj[["RNA"]]$counts["NRTN",]

AC_class<- rep("Im", length(ALKAL2_expr))
AC_class[seu_obj$AC_like>quantile(seu_obj$AC_like, 0.9)]<- "Hi"
AC_class[seu_obj$AC_like<=quantile(seu_obj$AC_like, 0.1)]<- "Lo"

# res_t_ALKAL2<- table(AC=AC_class, ALKAL2=ALKAL2_expr>quantile(ALKAL2_expr,.9))
# signif(chisq.test(res_t_ALKAL2)$p.value,2) # 0
# 
# res_t_NRTN<- table(AC=AC_class, NRTN=NRTN_expr>quantile(ALKAL2_expr,.9))
# signif(chisq.test(res_t_NRTN)$p.value,2) # 1.1e-83
# 
# res_df<- data.frame(
#   AC_like = c("Hi","Im","Lo"),  
#   ALKAL2 = 100*prop.table(res_t_ALKAL2,1)[,"TRUE"],
#   NRTN = 100*prop.table(res_t_NRTN,1)[,"TRUE"]
# )  
# res_df$AC_like<- factor(res_df$AC_like, levels = rev(c("Lo","Im","Hi")))
# res_df
# # AC_like    ALKAL2      NRTN
# # Hi      Hi 7.6635585 0.6512137
# # Im      Im 1.2472072 0.2395834
# # Lo      Lo 0.2639787 0.1790350

ALKAL2_t<- tapply(ALKAL2_expr, AC_class, "t.test")
NRTN_t<- tapply(NRTN_expr, AC_class, "t.test")
summary(aov(ALKAL2_expr~AC_class))[[1]]["AC_class","Pr(>F)"] # P = 0
summary(aov(NRTN_expr~AC_class))[[1]]["AC_class","Pr(>F)"] # P = 2.510984e-84

res_df<- data.frame(
  AC_like = c("Hi","Im","Lo"),
  ALKAL2 = c(ALKAL2_t$Hi$estimate, ALKAL2_t$Im$estimate, ALKAL2_t$Lo$estimate),
  ALKAL2_ci = c(ALKAL2_t$Hi$conf.int[2], ALKAL2_t$Im$conf.int[2], ALKAL2_t$Lo$conf.int[2]),
  NRTN = c(NRTN_t$Hi$estimate, NRTN_t$Im$estimate, NRTN_t$Lo$estimate),
  NRTN_ci = c(NRTN_t$Hi$conf.int[2], NRTN_t$Im$conf.int[2], NRTN_t$Lo$conf.int[2])
)

# Plot
p_ALKAL2<- ggbarplot(res_df, "AC_like", "ALKAL2",
                     fill = "steelblue", color = "steelblue") +
  geom_errorbar(aes(ymin=ALKAL2, ymax=ALKAL2_ci), width=0.4, alpha=0.9, size=1.3) +
  ylab("ALKAL2 expression") +
  xlab("AC-like signature") +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "italic")
  )

p_NRTN<- ggbarplot(res_df, "AC_like", "NRTN",
                     fill = "steelblue", color = "steelblue") +
  geom_errorbar(aes(ymin=NRTN, ymax=NRTN_ci), width=0.4, alpha=0.9, size=1.3) +
  ylab("NRTN expression") +
  xlab("AC-like signature") +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "italic")
  )

p_bp<- plot_grid(
  p_ALKAL2,p_NRTN,
  ncol=2
)
ggsave("results/figs/manuscript_Patel_AC_like_cor_ALKAL2.pdf", p_bp, width = 0.375*178, height = 0.16*265, units = "mm")


