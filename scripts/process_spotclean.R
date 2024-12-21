### This script removes spot swapping effects

library(SpotClean)
library(magrittr)
library(purrr)
library(dplyr)
library(S4Vectors)
library(future)
library(furrr)
library(SummarizedExperiment)

### Import spaceranger output 

# sample names
sample_names <- c("NB_LU_01_Pre_1","NB_LU_01_Pre_2", "NB_LU_01_Post_1","NB_LU_01_Post_2")
new_names <- c("NB1Pre1", "NB1Pre2", "NB1Post1", "NB1Post2")
image_slice_names <- new_names

# data directories
count_dir  <- purrr::map(sample_names, function(name) 
      paste("temp/data/spaceranger_out_NBLU01/",name,"/outs/raw_feature_bc_matrix", sep = "")) %>% 
      setNames(new_names)

# dir for spatial info
spatial_dir  <- purrr::map(sample_names, function(name) 
      paste("temp/data/spaceranger_out_NBLU01/",name,"/outs/spatial", sep = "")) %>% 
      setNames(new_names)


# load data
raw_list <- purrr::map(count_dir, function(dir)
      read10xRaw(dir))
  
# read spatial metadata    
slide_info_list <- purrr::map(spatial_dir, function(dir)
      read10xSlide(tissue_csv_file=file.path(dir, "tissue_positions.csv"),
                   tissue_img_file = file.path(dir,"tissue_lowres_image.png"),
                   scale_factor_file = file.path(dir,"scalefactors_json.json")))

#### filter outs missing barcodes
# some barcodes in raw matrix are filtered out.subset barcodes in slide_info to match those in raw matrix
slide_info_list <- purrr::map2(raw_list, slide_info_list, 
                               function(raw, slide_info){
                                           raw <- as.data.frame(raw)
                                           barcodes <- colnames(raw)
                                           slide <- slide_info$slide %>% 
                                                       dplyr::filter(barcode %in% barcodes)
                                           slide_info$slide <- slide
                                           slide_info
                               })

## create the slide object
slide_obj_list <- purrr::map2(raw_list,slide_info_list, 
                              function(raw,slide_info)
                                          createSlide(raw,slide_info))
 

## decontaminate the data

decont_obj_list  <- map(slide_obj_list,function(slide_obj)
            spotclean(slide_obj))

# save obj
readr::write_rds(decont_obj_list, "/temp/data/spotclean_decont_obj_list.rds")

## convert to Seurat object for downstream analyses
seurat_decont_obj_list <- purrr::pmap(list(decont_obj_list,spatial_dir,image_slice_names), 
                                      function(decont_obj, dir,slice_name){
                                            obj <- convertToSeurat(decont_obj,image_dir =dir, ,slice =slice_name  )
                                            obj$orig.ident <- slice_name
                                            obj
                                            })
readr::write_rds(seurat_decont_obj_list, "/temp/data/seurat_decont_obj_list.rds")
