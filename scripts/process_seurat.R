### spatial transcriptomics analysis of ST data corrected for spot swapping
## Data: merge each sample group separately

library(Seurat)
library(clustree)
library(patchwork)
library(tidyverse)
library(future)
library(furrr)
library(UCell)
library(ggpubr)
library(clusterProfiler)
library(ClusterFoldSimilarity)
library(clustree)

# SEURATE PREPROCESSING
########################

# Load spotclean seurat object
seurat_decont_obj_list_NB1 <-readr::read_rds("temp/data/seurat_decont_obj_list.rds")
seurat_decont_obj_list_NB2 <-readr::read_rds("NB_spatialT/temp/data/seurat_decont_obj_list_NBLU02.rds")
seurat_decont_obj_list <-c(seurat_decont_obj_list_NB1, seurat_decont_obj_list_NB2 ) 

# Run SCtransform before merging 
seurat_decont_obj_sct_list <- seurat_decont_obj_list %>% 
            map(function(seurat_obj)
                        SCTransform(seurat_obj, assay = "Spatial", verbose = TRUE)
            )

## Merge samples based on category (pre- and post-treament) 
# Create categories 
NB1Pre <- c("NB1Pre1","NB1Pre2")
NB1Post <- c("NB1Post1","NB1Post2")
NB2Post <- c("NB2Post1","NB2Post2")

sample_names <- c("NB1Pre", "NB1Post", "NB2Post")

# Merge 
seurat_grp_merged_list <- list(NB1Pre = merge(seurat_decont_obj_sct_list$NB1Pre1, y=seurat_decont_obj_sct_list$NB1Pre2, add.cell.ids =NB1Pre),
                               NB1Post = merge(seurat_decont_obj_sct_list$NB1Post1, y=seurat_decont_obj_sct_list$NB1Post2, add.cell.ids =NB1Post),
                               NB2Post = merge(seurat_decont_obj_sct_list$NB2Post1, y=seurat_decont_obj_sct_list$NB2Post2, add.cell.ids =NB2Post)
                               )

# Rejoin layers in spatial assay
seurat_grp_merged_list <- seurat_grp_merged_list %>% 
            map(function(seurat_obj){
                        seurat_obj@assays$Spatial <- JoinLayers(seurat_obj@assays$Spatial)
                        return(seurat_obj)
            })

# Impute missing values with ALRA 
seurat_grp_merged_list <- seurat_grp_merged_list %>% map(function(seurat_obj)
            SeuratWrappers::RunALRA(seurat_obj) )

# set default assay to SCT
seurat_grp_merged_list <- seurat_grp_merged_list %>% 
            map(function(seurat_obj){
                        DefaultAssay(seurat_obj) <- "SCT"
                        return(seurat_obj)
            })

### Dimensionality reduction 

# Categorize samples
seurat_group_list <- list(NB1Pre = list(seurat_decont_obj_sct_list$NB1Pre1,seurat_decont_obj_sct_list$NB1Pre2),
                           NB1Post = list(seurat_decont_obj_sct_list$NB1Post1,seurat_decont_obj_sct_list$NB1Post2),
                           NB2Post = list(seurat_decont_obj_sct_list$NB2Post1,seurat_decont_obj_sct_list$NB2Post2))

# run PCA
seurat_grp_merged_list <- purrr::map2(seurat_grp_merged_list,seurat_group_list,function(seurat_obj, seurat_group){
            # find variable features to use for pca 
            VariableFeatures(seurat_obj) <- c(VariableFeatures(seurat_group[[1]]), VariableFeatures(seurat_group[[2]]))
            # PCA
            seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = TRUE)
            seurat_obj
            })
                                                    
#  run UMAP 
seurat_grp_merged_list <- seurat_grp_merged_list %>%  
            map(function(seurat_obj){
                        seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:15) 
                        seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:15) 
                        seurat_obj
                        
                        })

# Test different cluster resolutions
resolutions <- c(0.1,0.2,0.3, 0.4, 0.5, 0.6,0.7,0.8, 0.9, 1, 1.1, 1.2, 1.5 )
seurat_obj_test <- list()

for(i in 1:length(seurat_grp_merged_list)){
            seurat_obj <- seurat_grp_merged_list[[i]]

            for (j in 1:length(resolutions)) {
                        seurat_obj <- FindClusters(seurat_obj, resolution = resolutions[j],
                                                   cluster.name = paste0("clust_res.", resolutions[j]),
                                                   dims = 1:15, min.cells = 100)
            }

            seurat_obj_test[[i]] <- seurat_obj
}

# Name obj
seurat_obj_test <- seurat_obj_test %>%  setNames(nm = sample_names)


# Selected resolutions
# NB1Pre == 0.6,
# NB1Post == 1,
# NB2Post == 0.7

# Findclusters at selected resolutions
sample_names <- c("NB1Pre", "NB1Post", "NB2Post")
res_by_sample <- c(0.6, 1, 0.7) %>%  setNames(nm = sample_names)

seurat_grp_merged_list <- map2(seurat_grp_merged_list,res_by_sample,
                                 function(seurat_obj, res){
                                             seurat_obj <- FindClusters(seurat_obj,resolution = res)
                                             # rename seurat clusters
                                             meta <- seurat_obj@meta.data %>% 
                                                         dplyr::mutate(seurat_clusters = paste0("C", seurat_clusters)) %>% 
                                                         dplyr::select(seurat_clusters)
                                             meta$seurat_clusters <- factor( meta$seurat_clusters, levels = paste0("C", 0:length(unique(meta$seurat_clusters))))
                                             seurat_obj <- AddMetaData(seurat_obj,metadata = meta)
                                             Idents(seurat_obj) <- "seurat_clusters"
                                             return(seurat_obj)
                                         
                                 })

# Prep merged data for differential expression analysis 
seurat_grp_merged_list <- map(seurat_grp_merged_list, function(seurat_obj)
            PrepSCTFindMarkers(seurat_obj, assay = "SCT")
)

# Save seurat obj 
# write_rds(seurat_grp_merged_list,"temp/data/seurat_grp_merged_list.rds")

# FIND DE MARKERS(based on seurat_clust)
de_markers_list <- map(seurat_grp_merged_list, function(seurat_obj){
            de_markers <- FindAllMarkers(seurat_obj, assay = "SCT", min.pct = 0.1, only.pos = T)
            de_markers 
            
            })

# Save DE data/result
de_markers_list %>% 
            openxlsx::write.xlsx(., file = 'results/tables/de_markers_ST_seurat_clust.xlsx')

# Regrouping of Seurat clusters in annotated clusters
#####################################################

# Get metadata
meta_NB1Pre <- seurat_grp_merged_list$NB1Pre@meta.data
meta_NB1Post<- seurat_grp_merged_list$NB1Post@meta.data
meta_NB2Post <- seurat_grp_merged_list$NB2Post@meta.data

# NBPre
# low resolution
rename_cluster_NB1Pre_lr <- dplyr::select(meta_NB1Pre, seurat_clusters) %>% 
            dplyr::mutate(cell_type_lr=case_when(seurat_clusters %in% c("C0", "C2", "C3", "C4","C6", "C8", "C9") ~ "NE",
                                                    seurat_clusters == "C1" ~ "Plasma",
                                                    seurat_clusters == "C5" ~ "Schwann",
                                                    seurat_clusters == "C7" ~ "CAF",
                                                    seurat_clusters == "C10" ~ "Endo",
                                                 TRUE ~ as.character(seurat_clusters))) %>% 
            dplyr::mutate(cell_type_lr_cond = paste0(cell_type_lr,"_pre")) %>% 
            dplyr::mutate(condition = "Pre") %>% 
            dplyr::select(condition, cell_type_lr, cell_type_lr_cond) %>% 
            rownames_to_column(var = "barcodes")


# high resolution    
rename_cluster_NB1Pre_hr <- dplyr::select(meta_NB1Pre, seurat_clusters) %>% 
            dplyr::mutate(cell_type_hr=case_when(seurat_clusters %in% c("C0", "C3","C4" ) ~ "NE1",
                                                    seurat_clusters %in% c("C2", "C9") ~ "NE2",
                                                    seurat_clusters == "C8" ~ "NE3",
                                                    seurat_clusters == "C6" ~ "NE4",
                                                    seurat_clusters == "C1" ~ "Plasma",
                                                    seurat_clusters == "C5" ~ "Schwann",
                                                    seurat_clusters == "C7" ~ "CAF",
                                                    seurat_clusters == "C10" ~ "Endo",
                                                 TRUE ~ as.character(seurat_clusters))) %>% 
            dplyr::mutate(cell_type_hr_cond = paste0(cell_type_hr,"_pre")) %>% 
            dplyr::select(cell_type_hr, cell_type_hr_cond) %>% 
            rownames_to_column(var = "barcodes")

# NB1Post
# low resolution
rename_cluster_NB1Post_lr <-  dplyr::select(meta_NB1Post, seurat_clusters) %>% 
            dplyr::mutate(cell_type_lr=case_when(seurat_clusters %in% c("C0", "C1", "C3", "C6", "C7", "C8", "C10") ~ "NE",
                                                 seurat_clusters %in% c("C4") ~ "Macro",
                                                 seurat_clusters %in% c("C5") ~ "Endo",
                                                 seurat_clusters %in% c("C2", "C9", "C12", "C13") ~ "CAF",
                                                 seurat_clusters  %in% c("C11") ~ "Plasma",
                                                 TRUE ~ as.character(seurat_clusters))) %>% 
            dplyr::mutate(cell_type_lr_cond = paste0(cell_type_lr,"_post")) %>% 
            dplyr::mutate(condition = "Post") %>% 
            dplyr::select(condition, cell_type_lr, cell_type_lr_cond) %>% 
            rownames_to_column(var = "barcodes")

# high resolution (hr)
rename_cluster_NB1Post_hr <-  dplyr::select(meta_NB1Post, seurat_clusters) %>% 
            dplyr::mutate(cell_type_hr=case_when(seurat_clusters %in% c("C0") ~ "NE1",
                                                 seurat_clusters %in% c("C1", "C8", "C10") ~ "NE2",
                                                 seurat_clusters %in% c( "C3", "C6") ~ "NE3",
                                                 seurat_clusters %in% c("C7") ~ "NE4",
                                                 seurat_clusters %in% c("C4") ~ "Macro",
                                                 seurat_clusters %in% c("C5") ~ "Endo",
                                                 seurat_clusters %in% c("C2","C9", "C12", "C13") ~ "CAF",
                                                 seurat_clusters  %in% c("C11") ~ "Plasma",
                                                 TRUE ~ as.character(seurat_clusters))) %>% 
            dplyr::mutate(cell_type_hr_cond = paste0(cell_type_hr,"_post")) %>% 
            dplyr::select(cell_type_hr, cell_type_hr_cond) %>% 
            rownames_to_column(var = "barcodes")

# NB2Post

# low resolution
rename_cluster_NB2Post_lr <-  dplyr::select(meta_NB2Post, seurat_clusters) %>% 
            dplyr::mutate(cell_type_lr=case_when(seurat_clusters %in% c("C0", "C3", "C5", "C7", "C9", "C10", "C11" ) ~ "NE",
                                                 seurat_clusters %in% c("C1", "C2", "C6" ) ~ "CAF",
                                                 seurat_clusters %in% c("C4") ~ "Schwann",
                                                 seurat_clusters  %in% c("C8", "C12") ~ "FZ_like",
                                                 TRUE ~ as.character(seurat_clusters))) %>% 
            dplyr::mutate(cell_type_lr_cond = paste0(cell_type_lr,"_post")) %>% 
            dplyr::mutate(condition = "Post") %>% 
            dplyr::select(condition, cell_type_lr, cell_type_lr_cond) %>% 
            rownames_to_column(var = "barcodes")

# high resolution (hr)
rename_cluster_NB2Post_hr <-  dplyr::select(meta_NB2Post, seurat_clusters) %>% 
            dplyr::mutate(cell_type_hr=case_when(seurat_clusters %in% c("C0", "C3", "C5", "C9","C10", "C11") ~ "NE1",
                                                 seurat_clusters %in% c("C7") ~ "NE2",
                                                 seurat_clusters %in% c( "C1", "C2", "C6") ~ "CAF", 
                                                 seurat_clusters %in% c("C4") ~ "Schwann",
                                                 seurat_clusters %in% c("C8", "C12") ~ "FZ_like",
                                                 TRUE ~ as.character(seurat_clusters))) %>% 
            dplyr::mutate(cell_type_hr_cond = paste0(cell_type_hr,"_post")) %>% 
            dplyr::select(cell_type_hr, cell_type_hr_cond) %>% 
            rownames_to_column(var = "barcodes")


# Relevel factors
# NB1Pre
rename_cluster_NB1Pre_lr$cell_type_lr <- factor(rename_cluster_NB1Pre_lr$cell_type_lr, 
                                                   levels = c("NE", "CAF", "Plasma", "Schwann",  "Endo")) 
            
cells_NB1pre <- c("NE1","NE2", "NE3", "NE4", "CAF", "Plasma", "Endo","Schwann"  )

rename_cluster_NB1Pre_hr$cell_type_hr <- factor(rename_cluster_NB1Pre_hr$cell_type_hr, 
                                                levels = cells_NB1pre )
rename_cluster_NB1Pre_hr$cell_type_hr_cond <- factor(rename_cluster_NB1Pre_hr$cell_type_hr_cond, 
                                                     levels = paste0(cells_NB1pre, "_pre") )


annot_NB1Pre <- left_join(rename_cluster_NB1Pre_lr, rename_cluster_NB1Pre_hr, by = "barcodes" ) %>% 
            column_to_rownames(var = "barcodes")

# NB1Post
rename_cluster_NB1Post_lr$cell_type_lr <- factor(rename_cluster_NB1Post_lr$cell_type_lr, 
                                                   levels = c("NE", "CAF", "Plasma","Endo", "Macro"))
cells_NB1post <- c("NE1","NE2", "NE3", "NE4", "CAF", "Plasma","Endo", "Macro")
rename_cluster_NB1Post_hr$cell_type_hr <- factor(rename_cluster_NB1Post_hr$cell_type_hr, 
                                                 levels = cells_NB1post)
rename_cluster_NB1Post_hr$cell_type_hr_cond <- factor(rename_cluster_NB1Post_hr$cell_type_hr_cond,  
                                                      levels = paste0(cells_NB1post, "_post") )
annot_NB1Post <- left_join(rename_cluster_NB1Post_lr, rename_cluster_NB1Post_hr, by = "barcodes" ) %>% 
            column_to_rownames(var = "barcodes")

# NB2Post
rename_cluster_NB2Post_lr$cell_type_lr <- factor(rename_cluster_NB2Post_lr$cell_type_lr, 
                                                 levels = c("NE", "CAF", "Schwann","FZ_like"))
cells_NB2post <- c("NE1","NE2", "CAF", "Schwann", "FZ_like")
rename_cluster_NB2Post_hr$cell_type_hr <- factor(rename_cluster_NB2Post_hr$cell_type_hr,
                                                 levels = cells_NB2post)

rename_cluster_NB2Post_hr$cell_type_hr_cond <- factor(rename_cluster_NB2Post_hr$cell_type_hr_cond,  
                                                      levels = paste0(cells_NB2post, "_post") )

annot_NB2Post <- left_join(rename_cluster_NB2Post_lr, rename_cluster_NB2Post_hr, by = "barcodes" ) %>% 
            column_to_rownames(var = "barcodes")

# Add cell annot metadata to seurat obj 
annot_list <- list(annot_NB1Pre, annot_NB1Post,annot_NB2Pre, annot_NB2Post)

seurat_grp_merged_list <- map2(seurat_grp_merged_list,annot_list, 
                               function(seurat_obj,annot){
                                           seurat_obj <- AddMetaData( object = seurat_obj,
                                                                      metadata = annot)
                                           seurat_obj
                               })

# ADRENOCORTICAL CELL SIG. SCORES-
##################################

seurat_grp_merged_list <- read_rds("temp/data/seurat_grp_merged_list.rds")

# Load cell markers 
adrenocortical_markers_ls <- read_rds("temp/data/adrenocortical_markers_ls.rds")

# Signature score with UCell
seurat_grp_merged_list <- seurat_grp_merged_list %>% 
            map(function(seurat_obj)
                        UCell::AddModuleScore_UCell(seurat_obj, features = adrenocortical_markers_ls, name = "")
            )

# REACTOME_STEROID PATHWAYS SIGNATURE SCORE

# Reactome  steroid genesets
reactome <- msigdbr::msigdbr(species = "Homo sapiens",category = "C2",  subcategory = "CP:REACTOME") %>% 
            dplyr::select(gs_name, human_gene_symbol )

steroid_gs <- c("REACTOME_METABOLISM_OF_STEROIDS", 
                "REACTOME_METABOLISM_OF_STEROID_HORMONES") %>% 
            setNames(nm = .)

reactome_steroid_pathways_list <- steroid_gs %>% 
            map(function(gs) reactome %>% 
                            group_by(gs_name) %>% 
                            filter(gs_name == gs) %>% 
                            pull(human_gene_symbol)
                )

# UCell score
seurat_grp_merged_list <- seurat_grp_merged_list %>% 
            map(function(seurat_obj)
                        AddModuleScore_UCell(seurat_obj, 
                                             features =reactome_steroid_pathways_list,   
                                             name = "",
                                             maxRank = 3000)
                )


### Adrenergic (ADRN) and Mesenchymal (MES) signatures score
#############################################################

# ADRN and MES signatures from van Groningen et al., 2017
# Load cell makers ls 
cell_markers_list <- read_rds("temp/data/cellMarkers_merge_ls.rds")
cell_markers_list <- cell_markers_list[c("ADRN", "MES")]

# UCell score
seurat_grp_merged_list <- seurat_grp_merged_list %>% 
            map(function(seurat_obj)
                        UCell::AddModuleScore_UCell(seurat_obj, 
                                                    features = cell_markers_list, 
                                                    name = "")
                )

# Calculate log2(ADRN_norm/MES_norm) () score (VADMES score) of NE cells.
seurat_grp_merged_list <- seurat_grp_merged_list %>% map(function(seurat_obj){
            meta <- seurat_obj@meta.data %>%
                        dplyr::select(cell_type_lr, ADRN, MES)
            meta <- meta %>%
                        mutate(
                                    ADRN_sub = if_else(cell_type_lr == "NE", ADRN, NA_real_),
                                    MES_sub = if_else(cell_type_lr == "NE", MES, NA_real_)
                        )
            meta <- meta %>%
                        dplyr::mutate(ADRN_norm = ADRN_sub/mean(ADRN_sub,  na.rm = TRUE),
                                      MES_norm= MES_sub/mean(MES_sub,  na.rm = TRUE),
                                      log2_ADRN_MES = log2(ADRN_norm/MES_norm)

                        )
            seurat_obj <- AddMetaData(seurat_obj, metadata = meta)
            seurat_obj

})

saveRDS("data/ST_NB_seurat.rds")





