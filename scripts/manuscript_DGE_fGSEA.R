# Libraries & functions
#######################
load("data/ST_NB.RData")
library(reshape2)
library(cowplot)
library(clusterProfiler)

# Turn gene set lists in df (to be used by GSEA)
Rea_df<- NULL
for(gs in names(Rea_ls)){
  gs_tmp<- cbind(gs, Rea_ls[[gs]])
  Rea_df<- rbind(Rea_df, gs_tmp)
}

PDB_df<- NULL
for(gs in names(PDB_Hs_ls)){
  gs_tmp<- cbind(gs, PDB_Hs_ls[[gs]])
  if(length(gs_tmp)<=1) next
  PDB_df<- rbind(PDB_df, gs_tmp)
}

# Load DGE data
###############

DE_clust_ls<- readRDS(file = "temp/DE_clust.rds")

fGSEA_all_ls<- list()
fGSEA_ls<- list()

for(clust in c("seurat", "lr", "hr")){
  for(gs_name in c("PDB_df", "Rea_df")){
    gs<- get(gs_name)
    # fGSEA
    DE_res_all<- DE_clust_ls[[clust]]
    for(cond in c("NB1Pre","NB1Post","NB2Post")){
      DE_res<- DE_res_all[[cond]]
      DE_res$p_val[DE_res$p_val==0]<- min(DE_res$p_val[DE_res$p_val!=0])/10
      DE_res<- DE_res[order(DE_res$p_val, -1*DE_res$avg_log2FC),]
      clusters<- unique(DE_res$cluster)
      fGSEA_t<- matrix(NA, length(unique(gs[,1])), length(clusters), dimnames = list(unique(gs[,1]),clusters))
      fGSEA_t_NES<- fGSEA_t
      fGSEA_t_LE<- fGSEA_t
      for(cluster in clusters){
        stat<- -log10(DE_res$p_val[DE_res$cluster==cluster])
        names(stat)<- DE_res$gene[DE_res$cluster==cluster]
        if(length(intersect(names(stat), gs[,2]))==0) next
        # stat[stat>300]<- 300
        fGSEA_res<- GSEA(geneList = stat, TERM2GENE = gs, scoreType = "pos", pvalueCutoff=1)@result
        fGSEA_res<- fGSEA_res[order(fGSEA_res$pval),]
        if(nrow(fGSEA_res)==0) next
        
        fGSEA_t[rownames(fGSEA_res),cluster]<- fGSEA_res$p.adjust
        fGSEA_t_NES[rownames(fGSEA_res),cluster]<- fGSEA_res$NES
        fGSEA_t_LE[rownames(fGSEA_res),cluster]<- apply(fGSEA_res, 1, function(x) paste0(unlist(strsplit(x["core_enrichment"], "/")),collapse = ","))        
        fGSEA_res$cluster<- cluster
        fGSEA_all_ls[[gs_name]][[clust]][[cond]][[cluster]]<- fGSEA_res
      }
      fGSEA_ls[[gs_name]][[clust]][[cond]][["padj"]]<- fGSEA_t
      fGSEA_ls[[gs_name]][[clust]][[cond]][["NES"]]<- fGSEA_t_NES
      fGSEA_ls[[gs_name]][[clust]][[cond]][["LE"]]<- fGSEA_t_LE
    }
  }
}
save(fGSEA_all_ls, fGSEA_ls, file = "results/data/fGSEA_clusters.RData")

# View(fGSEA_all_ls$PDB_df$hr$NB1Post$NE4)  
# View(fGSEA_all_ls$PDB_df$hr$NB1Post$CAF)  
# View(fGSEA_all_ls$PDB_df$hr$NB2Post$NE2)   
# 
# View(fGSEA_all_ls$Rea_df$hr$NB1Post$NE4)  
# View(fGSEA_all_ls$Rea_df$hr$NB2Post$FZ_like)   
# 
# View(fGSEA_all_ls$Ha_df$hr$NB1Post$NE4)   


# names(fGSEA_all_ls$PDB_df$lr$NB2Post)<- gsub("FZ_like", "AC_like", names(fGSEA_all_ls$PDB_df$lr$NB2Post))
# names(fGSEA_all_ls$PDB_df$hr$NB2Post)<- gsub("FZ_like", "AC_like", names(fGSEA_all_ls$PDB_df$hr$NB2Post))
# names(fGSEA_all_ls$Rea_df$lr$NB2Post)<- gsub("FZ_like", "AC_like", names(fGSEA_all_ls$Rea_df$lr$NB2Post))
# names(fGSEA_all_ls$Rea_df$hr$NB2Post)<- gsub("FZ_like", "AC_like", names(fGSEA_all_ls$Rea_df$hr$NB2Post))
# names(fGSEA_all_ls$Ha_df$lr$NB2Post)<- gsub("FZ_like", "AC_like", names(fGSEA_all_ls$Ha_df$lr$NB2Post))
# names(fGSEA_all_ls$Ha_df$hr$NB2Post)<- gsub("FZ_like", "AC_like", names(fGSEA_all_ls$Ha_df$hr$NB2Post))
# 
# colnames(fGSEA_ls$PDB_df$lr$NB2Post$padj)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$PDB_df$lr$NB2Post$padj))
# colnames(fGSEA_ls$PDB_df$lr$NB2Post$NES)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$PDB_df$lr$NB2Post$NES))
# colnames(fGSEA_ls$PDB_df$lr$NB2Post$LE)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$PDB_df$lr$NB2Post$LE))
# colnames(fGSEA_ls$Rea_df$lr$NB2Post$padj)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Rea_df$lr$NB2Post$padj))
# colnames(fGSEA_ls$Rea_df$lr$NB2Post$NES)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Rea_df$lr$NB2Post$NES))
# colnames(fGSEA_ls$Rea_df$lr$NB2Post$LE)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Rea_df$lr$NB2Post$LE))
# colnames(fGSEA_ls$Ha_df$lr$NB2Post$padj)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Ha_df$lr$NB2Post$padj))
# colnames(fGSEA_ls$Ha_df$lr$NB2Post$NES)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Ha_df$lr$NB2Post$NES))
# colnames(fGSEA_ls$Ha_df$lr$NB2Post$LE)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Ha_df$lr$NB2Post$LE))
# colnames(fGSEA_ls$PDB_df$hr$NB2Post$padj)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$PDB_df$hr$NB2Post$padj))
# colnames(fGSEA_ls$PDB_df$hr$NB2Post$NES)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$PDB_df$hr$NB2Post$NES))
# colnames(fGSEA_ls$PDB_df$hr$NB2Post$LE)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$PDB_df$hr$NB2Post$LE))
# colnames(fGSEA_ls$Rea_df$hr$NB2Post$padj)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Rea_df$hr$NB2Post$padj))
# colnames(fGSEA_ls$Rea_df$hr$NB2Post$NES)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Rea_df$hr$NB2Post$NES))
# colnames(fGSEA_ls$Rea_df$hr$NB2Post$LE)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Rea_df$hr$NB2Post$LE))
# colnames(fGSEA_ls$Ha_df$hr$NB2Post$padj)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Ha_df$hr$NB2Post$padj))
# colnames(fGSEA_ls$Ha_df$hr$NB2Post$NES)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Ha_df$hr$NB2Post$NES))
# colnames(fGSEA_ls$Ha_df$hr$NB2Post$LE)<- gsub("FZ_like", "AC_like", colnames(fGSEA_ls$Ha_df$hr$NB2Post$LE))



