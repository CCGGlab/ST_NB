# library(readr)
load("data/ST_NB.RData")
seurat <- readRDS("data/ST_NB_seurat.rds")
NBAtlas_deconv$id<- paste0(NBAtlas_deconv$orig.ident, "_",NBAtlas_deconv$barcode)

deconv_df<- NULL
for(s in c("NB1Pre", "NB1Post", "NB2Post")){
  deconv_tmp<- data.frame(
    id = names(seurat[[s]]$orig.ident),
    seurat = seurat[[s]]$seurat_clusters,
    lr = seurat[[s]]$cell_type_lr,
    hr = seurat[[s]]$cell_type_hr,
    sample = s
  )
  deconv_df<- rbind(deconv_df, deconv_tmp)
}
deconv_df<- merge(deconv_df, NBAtlas_deconv, by = "id")

# Some numbers
tapply(deconv_df$RCTD_NBAtlas_full_Neuroendocrine, list(deconv_df$lr,deconv_df$sample), "mean", na.rm=T)
# NB1Post    NB1Pre   NB2Post
# NE      0.4994327 0.4760589 0.5022684
# CAF     0.1573410 0.1148500 0.1051143
# Plasma  0.3554415 0.3452769        NA
# Schwann        NA 0.2593816 0.1122453
# Endo    0.2477847 0.2127752        NA
# Macro   0.1326388        NA        NA
# FZ_like        NA        NA 0.1445919
tapply(deconv_df$RCTD_NBAtlas_full_Fibroblast, list(deconv_df$lr,deconv_df$sample), "mean", na.rm=T)
# NB1Post     NB1Pre   NB2Post
# NE      0.1931767 0.07866458 0.2690943
# CAF     0.4956038 0.58522846 0.5450617
# Plasma  0.2767121 0.25649259        NA
# Schwann        NA 0.28763874 0.4635241
# Endo    0.3084646 0.47836452        NA
# Macro   0.3463406         NA        NA
# FZ_like        NA         NA 0.3018281
tapply(deconv_df$RCTD_NBAtlas_full_Neuroendocrine, list(deconv_df$hr,deconv_df$sample), "mean", na.rm=T)
# NB1Post    NB1Pre   NB2Post
# NE1     0.2148519 0.4445818 0.5621533
# NE2     0.5767061 0.4712177 0.0991360
# NE3     0.6421333 0.6273100        NA
# NE4     0.5203665 0.5564791        NA
# CAF     0.1573410 0.1148500 0.1051143
# Plasma  0.3554415 0.3452769        NA
# Endo    0.2477847 0.2127752        NA
# Schwann        NA 0.2593816 0.1122453
# Macro   0.1326388        NA        NA
# AC_like        NA        NA 0.1445919
tapply(deconv_df$RCTD_NBAtlas_full_Fibroblast, list(deconv_df$hr,deconv_df$sample), "mean", na.rm=T)
# NE1     0.3374005 0.04697324 0.2369286
tapply(deconv_df$RCTD_NBAtlas_full_Schwann, list(deconv_df$hr,deconv_df$sample), "mean", na.rm=T)
# NE1     0.08122158 0.08251645 0.04732303

# Color names
cols<-  c("#A6CEE3", "#FA877F", "#910A67", "#FF41ED", "#e46403", "#A31ACB", "#765827", "#a3db9a", "#145d66","#76885b", "#be6ca4", "#2e75b0")
names(cols)<- c("Neuroendocrine", "Fibroblast", "Plasma","Myeloid","Schwann","Endothelial", "Stromal_other","NK_cell","pDC", "T_cell", "B_cell","RBCs")

deconv_df<- melt(deconv_df[,-c(8:13)])
deconv_df$variable<- gsub("RCTD_NBAtlas_full_","",deconv_df$variable)
deconv_df$variable<- factor(deconv_df$variable, levels = names(cols))
deconv_df$hr<- gsub("FZ_like","AC-like", deconv_df$hr)
deconv_df$hr<- factor(deconv_df$hr, levels = c("NE1", "NE2", "NE3", "NE4", "CAF", "Plasma", "Endo", "Macro","Schwann","AC-like"))

deconv_df$sample<- factor(deconv_df$sample, levels = c("NB1Pre", "NB1Post", "NB2Post"))

# Plot
p_hr<- ggplot(deconv_df, aes(y=value, x=hr, fill=variable)) + 
  geom_bar(position="fill", stat = "summary") +
  facet_grid(.~sample, scales = "free", space = "free") +
  scale_fill_manual(values = cols, labels = names(cols),guide = guide_legend(ncol = 1)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  labs(y = "Avg. cell type proportion", x = "Subcluster",fill = "NBAtlas cell type") +
  theme_classic(12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "italic"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, 'mm')  
  )

# Plot
p_lr<- ggplot(deconv_df, aes(y=value, x=lr, fill=variable)) + 
  geom_bar(position="fill", stat = "summary") +
  facet_grid(.~sample, scales = "free", space = "free") +
  scale_fill_manual(values = cols, labels = names(cols),guide = guide_legend(ncol = 1)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  labs(y = "Avg. cell type proportion", x = "Main cluster",fill = "NBAtlas cell type") +
  theme_classic(12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "italic"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, 'mm')  
  )

# Plot
p_seurat<- list()
for(s in c("NB1Pre", "NB1Post", "NB2Post")){
  s_t<- table(deconv_df$hr[deconv_df$sample==s],deconv_df$seurat[deconv_df$sample==s])
  s_clusters<- NULL
  for(i in 1:nrow(s_t)){
    s_clusters_tmp<- names(s_t[i,][s_t[i,]!=0])
    s_clusters<- c(s_clusters,s_clusters_tmp)
  }
  s_clusters<- union(s_clusters, colnames(s_t))
  deconv_df2<- deconv_df
  deconv_df2$seurat<- factor(deconv_df2$seurat, levels = s_clusters)
  
  p_seurat[[s]]<- ggplot(deconv_df2, aes(y=value, x=seurat, fill=variable)) + 
    geom_bar(position="fill", stat = "summary") +
    facet_grid(.~sample, scales = "free", space = "free") +
    scale_fill_manual(values = cols, labels = names(cols),guide = guide_legend(ncol = 1)) +
    scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    labs(y = "Avg. cell type proportion", x = "Seurat cluster",fill = "NBAtlas cell type") +
    theme_classic(12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "italic"),
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.size = unit(3, 'mm')  
      )
}

# Merge & save
p<- plot_grid(
  p_seurat$NB1Pre,  
  p_seurat$NB1Post,  
  p_seurat$NB2Post, 
  plot_grid(p_hr + theme(legend.position = "none"), p_lr  + theme(legend.position = "none"), ncol = 2, rel_widths = c(6,4)),
  ncol=1
)
ggsave("results/figs/manuscript_figS3_deconv.pdf", width = 178, height = 0.8*265, units = "mm")
