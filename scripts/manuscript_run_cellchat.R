# Analysis of cell-cell communications using CellChat
library(CellChat)
library(NMF)
library(presto)
library(Seurat)
library(tidyverse)

# Prepare input data for CelChat analysis
seurat_grp_merged_list <- read_rds("data/ST_NB_seurat.rds")

# Remove NB1Pre sample
seurat_grp_merged_list <- seurat_grp_merged_list[-3]

# set default assay to SCT
seurat_grp_merged_list <- seurat_grp_merged_list %>% 
            map(function(seurat_obj){
                        DefaultAssay(seurat_obj) <- "SCT"
                        return(seurat_obj)
            })


sample_grps <- list(NB1Pre= c("NB1Pre1", "NB1Pre2") %>% setNames(., .),
                    NB1Post=  c("NB1Post1", "NB1Post2") %>%  setNames(., .),
                    NB2Post=  c("NB2Post1", "NB2Post2") %>%  setNames(., .))

# Extract normalized count matrix
data_input_list  <- seurat_grp_merged_list %>% 
            map(function(seurat_obj)
                        GetAssayData(seurat_obj, layer = "data", assay = "SCT") 
            )

# Cell Metadata
meta_list  <- seurat_grp_merged_list %>% 
            map(function(seurat_obj){
                        meta <- seurat_obj@meta.data
                        # create col called slices
                        meta <- meta %>% dplyr::mutate(slices = orig.ident)
                        meta$slices <- factor(meta$slices, levels = unique(meta$slices))
                        return(meta)
            } ) 


# Spatial coordinates
spatial_coords_list <-  map2(seurat_grp_merged_list,sample_grps , 
                             function(seurat_obj, sample_grp) {
                                         df_coords <- map(sample_grp, function(sample_name){
                                                     df <- seurat_obj@images[[sample_name]]@coordinates %>% 
                                                                 dplyr::select(imagerow,imagecol)
                                                     df
                                                     
                                         })
                                         df_coords <- bind_rows(df_coords)
                                         return(df_coords)
                             })

# Scale factors and spot diameters of the full resolution images

# spaceranger folder
sample_folder_NB1 <- c("NB_LU_01_Pre_1","NB_LU_01_Pre_2","NB_LU_01_Post_1","NB_LU_01_Post_2") %>% 
            setNames(., nm= str_remove_all(.,"_"))

sample_folder_NB2 <- c("NBLU02Post1","NBLU02Post2")

# the theoretical spot size (um) in 10X Visium
spot_size = 65

#### NB1

spatial_factors_NB1_list <- sample_folder_NB1 %>% 
            map(function(folder){
                        path <- paste0("temp/data/spaceranger_out_NBLU01/",folder,"/outs/spatial")
                        spatial_factors <- jsonlite::fromJSON(txt = file.path(path,"scalefactors_json.json"))
                        conversion_factor = spot_size/spatial_factors$spot_diameter_fullres
                        spatial_factors = data.frame(ratio = conversion_factor, tol = spot_size/2)
                        return(spatial_factors)
            })
# Group
sample_grp_names <- c("NB1Pre", "NB1Post")
spatial_factors_NB1_list <- list(spatial_factors_NB1_list[1:2], spatial_factors_NB1_list[3:4]) %>%  
            setNames(nm = sample_grp_names)
# Merge
sample_grp_names_2 <- list(NB1Pre = c("NB1Pre1", "NB1Pre2"), NB1Post = c("NB1Post1", "NB1Post2"))
spatial_factors_NB1_list <- map2(spatial_factors_NB1_list, sample_grp_names_2, 
                                 function(spatial_factors, sample_grp){
                                             spatial_factors_merged <- bind_rows(spatial_factors)
                                             rownames(spatial_factors_merged) <- sample_grp
                                             return(spatial_factors_merged)
                                 })

#### NB2
spatial_factors_NB2_list <- sample_folder_NB2 %>% 
            map(function(folder){
                        path <- paste0("/home/joachim/projects/NB_spatialT/temp/data/spaceranger_out_NBLU02/",folder,"/outs/spatial")
                        spatial_factors <- jsonlite::fromJSON(txt = file.path(path,"scalefactors_json.json"))
                        conversion_factor = spot_size/spatial_factors$spot_diameter_fullres
                        spatial_factors = data.frame(ratio = conversion_factor, tol = spot_size/2)
                        return(spatial_factors)
            })
# Group
sample_grp_names <- c("NB2Post")
spatial_factors_NB2_list <- list(spatial_factors_NB2_list[1:2]) %>%  
            setNames(nm = sample_grp_names)


# Merge
sample_grp_names_2 <- list(NB2Post = c("NB2Post1", "NB2Post2"))
spatial_factors_NB2_list <- map2(spatial_factors_NB2_list, sample_grp_names_2, 
                                 function(spatial_factors, sample_grp){
                                             spatial_factors_merged <- bind_rows(spatial_factors)
                                             rownames(spatial_factors_merged) <- sample_grp
                                             return(spatial_factors_merged)
                                 })

# list of all samples
spatial_factors_list <- purrr::flatten(list(spatial_factors_NB1_list, spatial_factors_NB2_list))


# Create CellChat object
cellchat_list <- pmap(list(data_input_list, meta_list, spatial_coords_list, spatial_factors_list),
                      function(data_input,meta,spatial_coords, spatial_factors)
                                  createCellChat(object = data_input, meta = meta, 
                                                 group.by = "cell_type_hr",
                                                 datatype = "spatial", 
                                                 coordinates = spatial_coords, 
                                                 spatial.factors = spatial_factors)
)


# Set the ligand-receptor interaction database----
load( "temp/custome_LR_DB_.RDATA")
# add my custom db to CellChatDB
CellChatDB.human <- updateCellChatDB(db = db.user, merged = TRUE, species_target = "human")

# Subset CellChatDB for cell-cell communication analysis: Only SECRETED SIGNALING
CellChatDB.use <- subsetDB(CellChatDB.human,c("Secreted Signaling"), key = "annotation") 


# Set the used database in the object
cellchat_list <- cellchat_list %>% 
            map(function(cellchat){
                        cellchat@DB <- CellChatDB.use
                        return(cellchat)
                        
            })


# Subset the expression data of signaling genes in CellChatDB.use
cellchat_list <- cellchat_list %>% 
            map(function(cellchat)
                        subsetData(cellchat)  
            )



# Overexpressed genes and overexpr interactions
cellchat_list <- cellchat_list %>% 
            map(function(cellchat){
                        cellchat <- identifyOverExpressedGenes(cellchat)
                        cellchat <- identifyOverExpressedInteractions(cellchat)
                        return(cellchat)
            })





### Inference of cell-cell communication network

# Compute the communication probability and infer cellular communication network
cellchat_list <- cellchat_list %>% 
            map(function(cellchat){
                        cellchat <-  computeCommunProb(cellchat, 
                                                       type =  "truncatedMean", 
                                                       trim = 0.2,
                                                       distance.use = TRUE, 
                                                       interaction.range = 250,
                                                       scale.distance = 0.01,
                                                       contact.dependent = FALSE,
                                                       nboot = 20)
                        return(cellchat)
            })


# Filter the cell-cell communication 
cellchat_list <- cellchat_list %>% 
            map(function(cellchat)
                        filterCommunication(cellchat, min.cells = 10)
            )


# Infer the cell-cell communication at a signaling pathway level
cellchat_list <- cellchat_list %>% 
            map(function(cellchat)
                        computeCommunProbPathway(cellchat)
            )

# Calculate the aggregated cell-cell communication network----
cellchat_list <- cellchat_list %>% 
            map(function(cellchat)
                        aggregateNet(cellchat)
            )

# Compute the network centrality scores 
# identify the dominant cell senders and receivers of signaling networks
cellchat_list <- cellchat_list %>% 
            map(function(cellchat)
                        netAnalysis_computeCentrality(cellchat, slot.name = "netP")
            )


# Save cellchat obj----
write_rds(cellchat_list, "temp/data/cellchat_list_ST_sct_truncateMean_0.2_cell_type_hr_secreted_sig_only_20241209.rds") 



# Visualize

# Cell-cell communication probalities

# ST data
seurat_grp_merged_list <- read_rds("temp/data/seurat_grp_merged_list.rds")

# Cell ids
cell_idents_list <- seurat_grp_merged_list[-3] %>% map(function(seurat_obj){
            Idents(seurat_obj) <- "cell_type_hr"
            ident <- Idents(seurat_obj)
            return(ident)
})


# Get ccc probabilities as df
# CCC probalities as df
ccc_prob_df_list <-  map2(cellchat_list, cell_idents_list, function(cellchat, idents){
            inter_prob_df <- netVisual_bubble(cellchat, sources.use = idents,  return.data = TRUE)
            inter_prob_df <- inter_prob_df$communication
            return(inter_prob_df)
            # note: in the returned df: pval > 0.05 = 1, pval > 0.01 but <= 0.05 = 2, pval <= 0.01 = 3.
})

# Save data
ccc_prob_df_list %>% 
            openxlsx::write.xlsx(., file = '/temp/tables/ccc_df_secreted_sig_only.xlsx')

