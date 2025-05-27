# Libraries & functions
########################
library(clusterProfiler)
library(fgsea)
library(ggplot2)
library(cowplot)
library(fgsea)

# Load data & results
######################
load("data/ST_NB.RData")
load("results/data/fGSEA_clusters.RData")

# Plot PDB dotplots
###################

gsea_df_ls<- list()

for(clust in c("seurat", "lr", "hr")){
  # Create df
  gsea_df_tmp_ls<- list()
  for(i in 1:3){
    s<- c("NB1Pre", "NB1Post", "NB2Post")[i]
    gsea_df<- melt(fGSEA_ls$PDB_df[[clust]][[s]][["padj"]], varnames = c("PDB_cells", "cluster"), value.name = "padj")
    # gsea_df$padj[gsea_df$padj>0.1]<- NA
    gsea_df<- na.omit(gsea_df)
    gsea_df$score<- -log10(gsea_df$padj)
    
    gsea_df_NES<- melt(fGSEA_ls$PDB_df[[clust]][[s]][["NES"]], varnames = c("PDB_cells", "cluster"), value.name = "NES")
    gsea_df_NES<- na.omit(gsea_df_NES)
    
    gsea_df<- merge(gsea_df,gsea_df_NES)
    gsea_df$cond<- s

    gsea_df$padj[gsea_df$padj>0.1]<- NA
    gsea_df<- na.omit(gsea_df)
    
    # make sure data available for each cluster
    missing_clusters<- setdiff(
      colnames(fGSEA_ls$PDB_df[[clust]][[s]][["padj"]]),
      unique(as.character(gsea_df$cluster))
    )
    if(length(missing_clusters)!=0){
      gsea_missing_df<- data.frame(
        PDB_cells = "Neurons",
        cluster = missing_clusters,
        padj = NA,
        score = NA, 
        NES = NA,
        cond = s
      )
      gsea_df_tmp_ls[[i]]<-  rbind(gsea_df, gsea_missing_df)
    }
    else gsea_df_tmp_ls[[i]]<- gsea_df 
  }
  gsea_df<- do.call("rbind", gsea_df_tmp_ls)
  gsea_df$cond<- factor(gsea_df$cond, levels = c("NB1Pre", "NB1Post", "NB2Post"))
  gsea_df$PDB_cells<- factor(gsea_df$PDB_cells, levels = rev(levels(gsea_df$PDB_cells)))
  
  # Set borders
  gsea_df$highlight<- F
  for(cond in c("NB1Pre", "NB1Post", "NB2Post")){
    if(clust!="seurat"){
      if(clust=="lr"){
        id_matrix<- matrix(F, 6, 7, dimnames = list(c("Neurons", "Fibroblasts", "Schwann_cells", "Endothelial_cells", "Plasma_cells", "Macrophages" ), c("NE", "CAF", "Plasma", "Schwann", "Endo", "Macro","FZ_like")))
        id_matrix["Neurons","NE"]<- T
      }
      if(clust=="hr"){
        id_matrix<- matrix(F, 6, 10, dimnames = list(c("Neurons", "Fibroblasts", "Schwann_cells", "Endothelial_cells", "Plasma_cells", "Macrophages" ), c("NE1", "NE2", "NE3", "NE4", "CAF", "Plasma", "Schwann", "Endo", "Macro","FZ_like")))
        id_matrix["Neurons",c("NE1", "NE2", "NE3", "NE4")]<- T
      }
      id_matrix["Fibroblasts","CAF"]<- T
      id_matrix["Macrophages","Macro"]<- T
      id_matrix["Schwann_cells","Schwann"]<- T
      id_matrix["Endothelial_cells","Endo"]<- T
      id_matrix["Plasma_cells","Plasma"]<- T
    }
    if(clust=="seurat"){
      if(cond== "NB1Pre"){
        id_matrix<- matrix(F, 6, 11, dimnames = list(c("Neurons", "Fibroblasts", "Schwann_cells", "Endothelial_cells", "Plasma_cells", "Macrophages" ), paste0("C",0:10)))
        id_matrix["Neurons",paste0("C",c(0,3,4,6,9,2,8))]<- T
        id_matrix["Fibroblasts",paste0("C",7)]<- T
        id_matrix["Schwann_cells",paste0("C",5)]<- T
        id_matrix["Endothelial_cells",paste0("C",10)]<- T
        id_matrix["Plasma_cells",paste0("C",1)]<- T
      }
      if(cond== "NB1Post"){
        id_matrix<- matrix(F, 6, 14, dimnames = list(c("Neurons", "Fibroblasts", "Schwann_cells", "Endothelial_cells", "Plasma_cells", "Macrophages" ), paste0("C",0:13)))
        id_matrix["Neurons",paste0("C",c(0,1,3,6,7,8,10))]<- T
        id_matrix["Fibroblasts",paste0("C",c(2,9,12,13))]<- T
        id_matrix["Endothelial_cells",paste0("C",5)]<- T
        id_matrix["Plasma_cells",paste0("C",11)]<- T
        id_matrix["Macrophages",paste0("C",4)]<- T
      }
      if(cond== "NB2Post"){
        id_matrix<- matrix(F, 6, 13, dimnames = list(c("Neurons", "Fibroblasts", "Schwann_cells", "Endothelial_cells", "Plasma_cells", "Macrophages" ), paste0("C",0:12)))
        id_matrix["Neurons",paste0("C",c(0,3,5,7,9,10,11))]<- T
        id_matrix["Fibroblasts",paste0("C",c(1,2,6))]<- T
        id_matrix["Schwann_cells",paste0("C",4)]<- T
      }
    }
    for(c in colnames(id_matrix)){
      PDB_cell<- names(id_matrix[,c][id_matrix[,c]])
      gsea_df$highlight[gsea_df$PDB_cells==PDB_cell & gsea_df$cluster==c & gsea_df$cond==cond]<- T
    }
  }
  gsea_df$highlight_stroke<- 0
  gsea_df$highlight_stroke[gsea_df$highlight]<- 1
  
  if(clust == "seurat") gsea_df$cluster<- factor(gsea_df$cluster, levels = paste0("C",0:13))
  if(clust == "lr") gsea_df$cluster<- factor(gsea_df$cluster, levels = c("NE", "CAF", "Endo", "FZ_like","Macro", "Plasma", "Schwann"))
  if(clust == "hr") gsea_df$cluster<- factor(gsea_df$cluster, levels = c("NE1", "NE2", "NE3", "NE4", "CAF", "Endo", "FZ_like","Macro", "Plasma", "Schwann"))
  
  gsea_df_ls[[clust]]<- gsea_df
}

# Plot
p_fGSEA_dotplot_ls<- list()
for(clust in c("seurat", "lr", "hr")){
  p<- ggplot(data = gsea_df_ls[[clust]], mapping = aes_string(x = "cluster", y = "PDB_cells")) + 
    geom_point(mapping = aes_string(size = "score", colour = "highlight", stroke = "highlight_stroke", fill = "NES"), pch = 21) +
    scale_fill_gradient(low = "gray", high = "red") +
    scale_color_manual(breaks=c(T, F), values=c("green", rgb(0,0,0,0))) +
    facet_grid(~cond, scales = "free", space = "free") +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 6, face = "italic"),
      axis.text.x = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.key.width = unit(0.25, "cm"),
      legend.key.height = unit(0.1, "cm"),
      legend.position = "bottom" ,
      strip.background = element_blank(),
      strip.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 5)
    ) +
    guides(size = guide_legend(title = "-log(Padj)"))
  p_fGSEA_dotplot_ls[[clust]]<- p
}

ggsave("results/figs/manuscript_figS2D_fGSEA.pdf", p_fGSEA_dotplot_ls$seurat, width = 178, height = 0.45*265, units = "mm")

# Get leading edges for each cluster (to be used in dotplots)

# NE
NE_1<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Pre$LE["Neurons",c("C0","C3","C4","C6","C9","C2","C8")]), collapse = ","),","))),"NA"))
NE_2<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Post$LE["Neurons",c("C0","C1","C8","C10","C3","C6","C7")]), collapse = ","),","))),"NA"))
NE_3<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB2Post$LE["Neurons",c("C0","C3","C7","C9","C10","C11")]), collapse = ","),","))),"NA"))

intersect(NE_1, intersect(NE_2, NE_3)) # NPY, PRPH, STMN2, TH, NEFM
# [1] "ARHGDIG"  "ATP6V0A1" "CACNA2D1" "CLSTN3"   "DGKB"     "EML5"    
# [7] "ENO2"     "FRRS1L"   "ISL1"     "NEFM"     "NPY"      "PCLO"    
# [13] "PRPH"     "SLC30A9"  "STMN2"    "SYT5"     "TH"       "TMEM130" 
# [19] "TMEM59L"  "TUBB3"    "VSNL1"   

# CAF
CAF_1<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Pre$LE["Fibroblasts",c("C7")]), collapse = ","),","))),"NA"))
CAF_2<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Post$LE["Fibroblasts",c("C13","C12","C2","C9")]), collapse = ","),","))),"NA"))
CAF_3<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB2Post$LE["Fibroblasts",c("C1","C6")]), collapse = ","),","))),"NA"))

intersect(CAF_1, intersect(CAF_2, CAF_3)) # DCN, COL1A1, MGP (CEBPB, EGR1: only in 2 ) 
# [1] "COL1A1"   "DCN"      "MGP"      "PDGFRA"   "RUNX1"    "SERPINH1"

# Schwann
Schwann_1<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Pre$LE["Schwann_cells",c("C5")]), collapse = ","),","))),"NA"))
Schwann_3<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB2Post$LE["Schwann_cells",c("C4")]), collapse = ","),","))),"NA"))

intersect(Schwann_1, Schwann_3) # SEMA3B, SOX10, PLP1 (MPZ, S100B: only 1) 
# [1] "CRYAB"  "EGFL8"  "ITGB4"  "PLP1"   "S100B"  "SCN7A"  "SEMA3B" "SOX10"  

# Endothelial
Endo_1<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Pre$LE["Endothelial_cells",c("C10")]), collapse = ","),","))),"NA"))
Endo_2<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Post$LE["Endothelial_cells",c("C5")]), collapse = ","),","))),"NA"))

intersect(Endo_1, Endo_2) # Only 1: CAV1, MCAM, CLDN5, HSPG2, CLIC4  
# [1] "CAV1"    "CD93"    "CLDN5"   "CLIC4"   "EPAS1"   "HSPG2"   "ITGB3"  
# [8] "MCAM"    "PDE3A"   "PODXL"   "PTPRB"   "SPARCL1" "THBD"    "TM4SF1"

# Plasma
Plasma_1<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Pre$LE["Plasma_cells",c("C1")]), collapse = ","),","))),"NA"))
Plasma_2<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Pos$LE["Plasma_cells",c("C11")]), collapse = ","),","))),"NA"))

intersect(Plasma_1, Plasma_2) # IGHA1  IGHG1  IGHM (added from 1 only)  IGKC JCHAIN
# ""IGHA1"  "IGHG1"  "IGKC"   "JCHAIN"

# Macrophage
Macro_2<- sort(setdiff(unique(unlist(strsplit(paste0((fGSEA_ls$PDB_df$seurat$NB1Post$LE["Macrophages",c("C4")]), collapse = ","),","))),"NA"))

Macro_2 # CD68, TYROBP, CD74, SLC11A1, FGR
# [1] "ADGRE5"  "AIF1"    "C5AR1"   "CCL2"    "CD14"    "CD68"    "CD74"   
# [8] "CSF1R"   "CXCL2"   "CYBB"    "CYP27A1" "FGL2"    "FGR"     "FPR1"   
# [15] "GPNMB"   "ITGAX"   "LGALS3"  "LILRA5"  "MS4A6A"  "NR4A3"   "RAB20"  
# [22] "S100A4"  "S100A8"  "SAMSN1"  "SLC11A1" "TYROBP"  "UCP2"   

# AC-like

# Plot specific RS plots
#########################

DE_clust_ls<- readRDS(file = "temp/DE_clust.rds")
s<- "NB2Post"

# Plot running score AC like
# DE_res<- readxl::read_excel("temp/de_hr_filtered.xlsx",sheet = s)
DE_res<- DE_clust_ls$hr[[s]]
DE_res$p_val[DE_res$p_val==0]<- min(DE_res$p_val[DE_res$p_val!=0])/10
DE_res<- DE_res[order(DE_res$p_val, -1*DE_res$avg_log2FC),]
stat<- -log10(DE_res$p_val[DE_res$cluster=="AC_like"])
names(stat)<- DE_res$gene[DE_res$cluster=="AC_like"]

gs_name<- "REACTOME_METABOLISM_OF_STEROIDS"
p_pw<- signif(fGSEA_all_ls$Rea_df$hr$NB2Post$AC_like$p.adjust[fGSEA_all_ls$Rea_df$hr$NB2Post$AC_like$Description==gs_name],3)
p1<- plotEnrichment(Rea_ls[[gs_name]],stat) +
  ggtitle(paste0(gs_name,"\n(Padj=",p_pw,")")) + 
  # geom_line(size=0.5, col="green") +
  theme(
    plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
  ) +
  scale_x_continuous(name = "Rank") +
  scale_y_continuous(name = "Enrichment Score", limits = c(-0.1,0.8))

gs_name<- "REACTOME_METABOLISM_OF_STEROID_HORMONES"
p_pw<- signif(fGSEA_all_ls$Rea_df$hr$NB2Post$AC_like$p.adjust[fGSEA_all_ls$Rea_df$hr$NB2Post$AC_like$Description==gs_name],3)
p2<- plotEnrichment(Rea_ls[[gs_name]],stat) +
  ggtitle(paste0(gs_name,"\n(Padj=",p_pw,")")) + 
  # geom_line(size=0.5, col="green") +
  theme(
    plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
  ) +
  scale_x_continuous(name = "Rank") +
  scale_y_continuous(name = "Enrichment Score", limits = c(-0.1,1))

p<- plot_grid(
  p1,
  p2,
  ncol= 1
)
ggsave("results/figs/manuscript_fig2A_RS.pdf", p, width = 0.2*178, height = 0.3*265, units = "mm")

# Save in xls
#############

# Merge by clustering group
fgsea_df_ls<- list()
fgsea_df_ls$caption<- as.data.frame(
  cbind(
    c("Table S2. GSEA results of all unsupervised and annotated clusters that were identified in the 3 tumor samples.",
      "GSEA was performed after sorting the genes based on their p values, as obtained from a differential gene expression analysis between each cluster and the rest of the sample. Two gene set databases were used: cell signatures from PanglaoDB (PDB) and Reactome (Rea).",
      "Tabnames refer to the gene set database and the tumor sample (seperated by an underscore).",
      "",
      "Cluster", 
      "Description",
      "setSize",
      "enrichmentScore",
      "NES",
      "pvalue",
      "p.adjust",
      "leading_edge"),
    c("",
      "",
      "",
      "",
      "cluster (either C0 - C13 for unsupervised clusters or as annotated)",
      "Gene set",
      "Number of genes retrieved in the gene set",
      "Enrichment score",
      "Normalized enrichment score",
      "P value",
      "Adjusted P value, as determined using the Benjamini Hochberg method",
      "Genes in the leading edge (i.e., highest ranked)"
      )
  )
)

for(gs_name in c("PDB_df", "Rea_df")){
  for(cond in c("NB1Pre","NB1Post","NB2Post")){
    s<- do.call("rbind", fGSEA_all_ls[[gs_name]][["seurat"]][[cond]])
    h<- do.call("rbind", fGSEA_all_ls[[gs_name]][["hr"]][[cond]])
    df_tmp<- rbind(s,h)
    df_tmp[,1]<- df_tmp$cluster
    df_tmp$cluster<- NULL
    colnames(df_tmp)[1]<- "cluster"
    # Remove less meaningful columns
    df_tmp$leading_edge<- df_tmp$core_enrichment
    df_tmp$qvalue<- NULL
    df_tmp$rank<- NULL
    df_tmp$core_enrichment<- NULL
    fgsea_df_ls[[paste0(gsub("_df","",gs_name),"_",cond)]]<- df_tmp
  }
}

WriteXLS::WriteXLS(fgsea_df_ls,"results/tables/manuscript_tableS2_fGSEA.xlsx",row.names = F)
