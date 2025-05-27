# NBAtlas validation
####################

library(cowplot)
library(Seurat)
library(reshape)
library(ggplot2)
library(ggpubr)
library(UCell)

NBa<- readRDS("downloads/NBAtlas/seuratObj_NBAtlas_share_v20241203.rds") # Downloaded from https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/f5969395-5f6e-4c5d-a61a-5894773d0fee/file_downloaded

# 1. Ligand receptors
####################

for(genes in list(
  c("ALKAL2", "ALK"),
  c("NRTN", "GFRA2"),
  c("NRTN", "RET"),
  c("FGF9", "FGFR1"),
  c("CCL18", "PITPNM3"),
  c("FGF7", "FGFR1")
  )){
  p_tmp<- FeaturePlot(NBa, cols = c("lightgrey", "brown", "blue"), order = T, features = genes, pt.size = 0.001, max.cutoff = "q95", blend=T, alpha = 0.3, raster = F)
  ggsave(paste0("results/figs/manuscript_NBAtlas_",paste0(genes, collapse = "_"),".png"), p_tmp[[3]] + theme_void() + ggtitle("") + theme(legend.position = "none"), width = 0.5*178, height = 0.25*265, units = "mm", dpi = 600)
  ggsave(paste0("results/figs/manuscript_NBAtlas_",paste0(genes, collapse = "_"),"_legend.pdf"), p_tmp[[4]] + theme(axis.text = element_text(size = 6), axis.title = element_text(size = 7, face = "italic"), title = element_text(size=8)), width = 50, height = 40, units = "mm")
  # Sometimes difficult to seperate
  ggsave(paste0("results/figs/manuscript_NBAtlas_",genes[1],".png"), p_tmp[[1]] + theme_void() + ggtitle("") + theme(legend.position = "none"), width = 0.5*178, height = 0.25*265, units = "mm", dpi = 600)
  ggsave(paste0("results/figs/manuscript_NBAtlas_",genes[2],".png"), p_tmp[[2]] + theme_void() + ggtitle("") + theme(legend.position = "none"), width = 0.5*178, height = 0.25*265, units = "mm", dpi = 600)
}

# Treatment, MYCN effects
###########################

table_S1<- as.data.frame(readxl::read_excel("downloads/pub/bonine_2024/1-s2.0-S2211124724011550-mmc2.xlsx", skip = 2))

# Derive treatment from table S1
NBa$treatment_status_coarse<- NA
NBa$treatment_status_coarse[NBa$Sample%in%table_S1[table_S1$Timepoint=="pre-treatment","Sample"]]<- "untreated"
NBa$treatment_status_coarse[NBa$Sample%in%table_S1[table_S1$Timepoint=="post-treatment","Sample"]]<- "treated"
NBa$treatment_status_coarse <- factor(NBa$treatment_status_coarse, levels = c("untreated", "treated"))

CCL18_expr<-  NBa[["RNA"]]$counts["CCL18",]
PITPNM3_expr<-  NBa[["RNA"]]$counts["PITPNM3",]

CCL18_df<- data.frame(
  gene = rep(c("CCL18", "PITPNM3"),each = length(CCL18_expr)),
  expr = as.numeric(c(CCL18_expr, PITPNM3_expr)),
  # expr = as.numeric(c(CCL18_expr, PITPNM3_expr)>0),
  MYCN = rep(NBa$MYCN_amplification,2),
  treatment = rep(NBa$treatment_status_coarse,2),
  macrophage_cluster = rep(NBa$Cell_type_wImmuneZoomAnnot=="Macrophage", 2),
  NE_cluster = rep(NBa$Cell_type=="Neuroendocrine", 2)
)
CCL18_df$cluster<- NA
CCL18_df$cluster[CCL18_df$macrophage_cluster==T]<- "Macrophage"
CCL18_df$cluster[CCL18_df$NE_cluster==T]<- "NE"
CCL18_df$cluster[CCL18_df$myeloid_cluster==F&CCL18_df$NE_cluster==F]<- "Other"
CCL18_df$cluster<- factor(CCL18_df$cluster, c("NE", "Macrophage", "Other"))

p_bar_ls<- list()
for(c in c("NE","Macrophage")){
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
ggsave(paste0("results/figs/manuscript_NBa_CCL18_PITPNM3_covar_bar.pdf"), p_bar, width = 0.25*178, height = 0.4*265, units = "mm", dpi = 600, bg = "white")

# some numbers
t.test(CCL18_df$expr[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$treatment=="treated"]~CCL18_df$MYCN[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$treatment=="treated"])$p.value #3.017762e-13
t.test(CCL18_df$expr[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$treatment=="untreated"]~CCL18_df$MYCN[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$treatment=="untreated"])$p.value # 7.257542e-23
t.test(CCL18_df$expr[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$MYCN=="Non-amplified"&!is.na(CCL18_df$treatment)]~CCL18_df$treatment[CCL18_df$gene=="CCL18"&CCL18_df$cluster=="Macrophage"&CCL18_df$MYCN=="Non-amplified"&!is.na(CCL18_df$treatment)])$p.value #0.00948

t.test(CCL18_df$expr[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="NE"&CCL18_df$treatment=="treated"]~CCL18_df$MYCN[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="NE"&CCL18_df$treatment=="treated"])$p.value #0.34
t.test(CCL18_df$expr[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="NE"&CCL18_df$treatment=="untreated"]~CCL18_df$MYCN[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="NE"&CCL18_df$treatment=="untreated"])$p.value # 8.989391e-07
t.test(CCL18_df$expr[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="NE"&CCL18_df$MYCN=="Non-amplified"]~CCL18_df$treatment[CCL18_df$gene=="PITPNM3"&CCL18_df$cluster=="NE"&CCL18_df$MYCN=="Non-amplified"])$p.value #0.00079

# Visualize UMAP
genes<- c("CCL18", "PITPNM3")
so_treated<- subset(x = NBa, subset = treatment_status_coarse == "treated")
p_treated<- FeaturePlot(so_treated, cols = c("lightgrey", "brown", "blue"), order = T, features = genes, pt.size = 0.001, max.cutoff = "q95", blend=T, alpha = 0.3, raster = F, split.by = "MYCN_amplification")
so_treated<- subset(x = NBa, subset = treatment_status_coarse == "untreated")
p_untreated<- FeaturePlot(so_treated, cols = c("lightgrey", "brown", "blue"), order = T, features = genes, pt.size = 0.001, max.cutoff = "q95", blend=T, alpha = 0.3, raster = F, split.by = "MYCN_amplification")

p_UMAP<- plot_grid(
  p_untreated[[3]] + theme_void() + ggtitle("") + theme(legend.position = "none"), p_treated[[3]]+ theme_void() + ggtitle("") + theme(legend.position = "none"), # Amp
  p_untreated[[7]]+ theme_void() + ggtitle("") + theme(legend.position = "none"), p_treated[[7]]+ theme_void() + ggtitle("") + theme(legend.position = "none"), # Non amp
  ncol = 2
)
ggsave(paste0("results/figs/manuscript_NBa_CCL18_PITPNM3_covar.png"), p_UMAP, width = 0.5*178, height = 0.3*265, units = "mm", dpi = 600, bg = "white")
# 
# # Explore others
# WEE1_expr<-  NBa[["RNA"]]$counts["WEE1",]
# tapply(WEE1_expr, NBa$Cell_type, "mean")
# tapply(WEE1_expr, list(NBa$Cell_type, NBa$treatment_status_coarse), "mean")
# cor.test(WEE1_expr[NBa$treatment_status_coarse=="treated"&NBa$MYCN_amplification=="Non-amplified"], CCL18_expr[NBa$treatment_status_coarse=="treated"&NBa$MYCN_amplification=="Non-amplified"])

# 2. AC-like signature
####################
AC<- list(AC = readRDS("/home/jimmy/projects/ST_NB/temp/AC_like_signature.rds"))
NBa <- AddModuleScore_UCell(
  NBa,
  AC,
  chunk.size = 5000,
  name = "",
  ncores = 1
)

p_AC<- FeaturePlot(NBa, features = "AC",order = T, pt.size = 0.001, min.cutoff = "q05", raster=F) +
  scale_fill_gradient2(
    high= "blue", mid = "lightgrey", low = "lightgrey", na.value = "lightgrey") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 7),      
    legend.text = element_text(size = 6),
    legend.background=element_blank(),
    legend.key.size = unit(0.2, "cm")
  )

ggsave("results/figs/manuscript_NBAtlas_AC.png", p_AC + theme_void() + ggtitle("") + theme(legend.position = "none"), width = 0.5*178, height = 0.25*265, units = "mm", dpi = 600)
ggsave("results/figs/manuscript_NBAtlas_AC_legend.pdf", p_AC + theme(axis.text = element_text(size = 6), axis.title = element_text(size = 7, face = "italic"), title = element_text(size=8)), width = 50, height = 40, units = "mm")

# ALK_expr<-  NBa[["RNA"]]$counts["ALK",]
# ALKAL2_expr<-  NBa[["RNA"]]$counts["ALKAL2",]
# boxplot(NBa$AC~NBa$Cell_type, las=2)
# boxplot(ALKAL2_expr~NBa$Cell_type, las=2)
# barplot(tapply(ALKAL2_expr>0, NBa$Cell_type, "mean"), las = 2)
# barplot(tapply(ALK_expr>0, NBa$Cell_type, "mean"), las = 2)
# barplot(tapply(NBa$AC>0, NBa$Cell_type, "mean"), las = 2)


# Overlap ALKAL2-AC
ALKAL2_expr<-  NBa[["RNA"]]$counts["ALKAL2",]
NRTN_expr<-  NBa[["RNA"]]$counts["NRTN",]

AC_df<- data.frame(
  AC = NBa$AC,
  clust = NBa$Cell_type_wImmuneZoomAnnot,
  ALKAL2 = ALKAL2_expr,
  NRTN = NRTN_expr
)
saveRDS(AC_df, file = "temp/NBa_AC.rds")

AC_df<- readRDS("temp/NBa_AC.rds")
AC_class<- rep("Im", nrow(AC_df))
AC_class[AC_df$AC>quantile(AC_df$AC, 0.9)]<- "Hi"
AC_class[AC_df$AC<=quantile(AC_df$AC, 0.1)]<- "Lo"

ALKAL2_t<- tapply(AC_df$ALKAL2, AC_class, "t.test")
NRTN_t<- tapply(AC_df$NRTN, AC_class, "t.test")
summary(aov(AC_df$ALKAL2~AC_class))[[1]]["AC_class","Pr(>F)"] # P = 1.413873e-120
summary(aov(AC_df$NRTN~AC_class))[[1]]["AC_class","Pr(>F)"] # P = 0.9420629

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
ggsave("results/figs/manuscript_NBa_AC_like_cor_ALKAL2.pdf", p_bp, width = 0.375*178, height = 0.16*265, units = "mm")

