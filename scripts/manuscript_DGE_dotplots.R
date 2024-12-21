# Libraries & functions
#######################
source("scripts/functions/custom_dotplot.R")
library(Seurat)

# Load seurat objects
#####################
seurat_grp_merged_list <- readRDS("data/ST_NB_seurat.rds")

# Load DGE data
###############
DE_clust_ls<- list()

for(s in c("NB1Pre","NB1Post","NB2Post")){
  DE_clust_ls[["seurat"]][[s]]<- readxl::read_excel("temp/de_seurat_filtered.xlsx",sheet = s)  
  DE_clust_ls[["lr"]][[s]]<- readxl::read_excel("temp/de_lr_filtered.xlsx",sheet = s)
  DE_clust_ls[["hr"]][[s]]<- readxl::read_excel("temp/de_hr_filtered.xlsx",sheet = s)
} 

# Create dotplot
################

# Select 5 genes for each cluster from leading edges (in as many as possible)
LE_gene_markers<- list(
  NE = c("NPY", "TH", "STMN2", "PRPH", "NEFM"),
  CAF = c("DCN","COL1A1","EGR1","CEBPB","MGP"),
  Plasma = c("IGKC", "IGHA1", "JCHAIN", "IGHG1","IGHM"),
  Schwann = c("SOX10", "PLP1", "S100B","MPZ", "SEMA3B"),
  Endo = c("CAV1", "MCAM", "CLDN5", "HSPG2", "CLIC4"),
  Macro = c("CD68", "TYROBP", "CD74", "SLC11A1", "FGR"),
  FZ_like = c("NR5A1","CYB5B","CYP11B1","STAR", "ALKAL2")
)

# Dotplot for lr clusters 
p<- custom_dotplot(seurat_list = seurat_grp_merged_list, seurat_cluster_feature = "cell_type_lr", marker_gene_list = LE_gene_markers, radius_size = 1, dot.scale = 1.5, col.min = -2.5, col.max = 2.5)
ggsave("results/figs/manuscript_fig1_dotplots.pdf", p, width = .3*178, height = 0.4*265, units = "mm")

# Dotplot for seurat  clusters
p<- custom_dotplot(seurat_list = seurat_grp_merged_list, seurat_cluster_feature = "seurat_clusters", marker_gene_list = LE_gene_markers, radius_size = 1, dot.scale = 1.5, col.min = -2.5, col.max = 2.5, isMainPlot = F)
ggsave("results/figs/manuscript_figS2_dotplots.pdf", p, width = .5*178, height = 0.4*265, units = "mm")

# Dotplot for hr clusters
p<- custom_dotplot(seurat_list = seurat_grp_merged_list, seurat_cluster_feature = "cell_type_hr", marker_gene_list = LE_gene_markers, radius_size = 1, dot.scale = 1.5, col.min = -2.5, col.max = 2.5, isMainPlot = F)
ggsave("results/figs/manuscript_figS4_dotplots_hr.pdf", p, width = .4*178, height = 0.4*265, units = "mm")

# Create dotplot for NE subclone-specific clusters
###########################################

# Get top 5 from each DE
genes_sel_ls<- list()
for(s in c("NB1Pre","NB1Post","NB2Post")){
  for(c in c("NE1", "NE2", "NE3", "NE4")){
    DE_tmp<- DE_clust_ls$hr[[s]][DE_clust_ls$hr[[s]]$cluster==c,]
    DE_tmp<- DE_tmp[order(DE_tmp$p_val, -1*DE_tmp$avg_log2FC),]
    if(is.na(DE_tmp$gene[1])) next
    genes_tmp<- DE_tmp$gene
    # unly unique ones
    genes_tmp<- genes_tmp[grep("DEPRECATED",genes_tmp,invert = T)]
    genes_tmp<- setdiff(genes_tmp, unlist(genes_sel_ls))
    
    # Select top 5
    genes_sel_ls[[paste0(s,"_",c)]]<- genes_tmp[1:5]
  }
}
p<- custom_dotplot(seurat_list = seurat_grp_merged_list, seurat_cluster_feature = "cell_type_hr", marker_gene_list = genes_sel_ls, radius_size = 1, dot.scale = 1.5, col.min = -2.5, col.max = 2.5, isMainPlot = F, cluster_subset = c("NE1", "NE2", "NE3", "NE4"))
ggsave("results/figs/manuscript_figS4_dotplots_hr_NE_1.pdf", p, width = .4*178, height = 0.6*265, units = "mm")

# Save data to suppl table
#############################
DE_df_ls<- list()
DE_df_ls$caption<- as.data.frame(
  cbind(
    c("Table S1. Differential gene expression results.",
      "Differential gene expression was determined between all unsupervised and annotated clusters that were identified in the 3 tumor samples.",
      "Tabnames refer to the tumor sample",
      "",
      "gene",
      "cluster",
      "p_val",
      "avg_log2FC",
      "pct.1",
      "pct.2",
      "p_val_adj"),
    c("",
      "",
      "",
      "",
      "HGNC symbol of the gene that was analysed",
      "cluster (either C0 - C13 for unsupervised clusters or as annotated)",
      "P value",
      "Average log2 fold change value between the expression of the spots in the analyzed cluster and the other clusters",
      "Proportion of the spots in the analysed cluster where the gene is expressed",
      "Proportion of the spots in the non-analysed cluster where the gene is expressed",
      "Adjusted P value"
    )
  )
)

for(cond in c("NB1Pre","NB1Post","NB2Post")){
  s<- do.call("cbind", DE_clust_ls[["seurat"]][[cond]])
  h<- do.call("cbind", DE_clust_ls[["hr"]][[cond]])
  df_tmp<- rbind(s,h)
  df_tmp<- df_tmp[,c("gene", "cluster", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
  DE_df_ls[[cond]]<- as.data.frame(df_tmp)
}

WriteXLS::WriteXLS(DE_df_ls,"results/tables/manuscript_tableS1.xlsx",row.names = F)
save.image("results/data/manuscript_clustering_DGE.RData")
