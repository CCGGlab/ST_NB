# Get Cangelosi data from supplements
# wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7563184/bin/cancers-12-02343-s001.zip -P temp/
# unzip temp/cancers-12-02343-s001.zip -d temp/
#####################################################

# Libraries
library(singscore)
library(SummarizedExperiment)
library(tidyverse)

# Get expression data
NB_expr<- readxl::read_excel("temp/cancers-852053 supplementray layout/Table S1.xlsx", skip=2)
NB_expr<- as.data.frame(NB_expr)
rownames(NB_expr)<- NB_expr$Genes
NB_expr$Genes<- NULL

# Get clinical data Cangelosi et al
NB_clin<- readxl::read_excel("temp/cancers-852053 supplementray layout/Table S3.xlsx", skip = 2)

# Only RNA-Seq
NB_clin<- as.data.frame(NB_clin[NB_clin$Platform=="RNAseq",])
rownames(NB_clin)<- NB_clin$`Patient ID`

# Get AC signatures for Cangelosi data
######################################

# Get top50 DE genes in AC_like
DE_clust<- readxl::read_excel("../ST_NB/ST_NB_Joachim/temp/tables/de_hr_filtered.xlsx",sheet = "NB2Post")
DE_clust<- DE_clust[DE_clust$cluster == "FZ_like",] # FZ_like = AC_like (previous name)

top_markers <- DE_clust %>% 
  dplyr::filter(!str_detect(gene,"DEPRECATED")) %>% 
  dplyr::filter(avg_log2FC > 0) %>% 
  dplyr::filter(p_val_adj <= 0.01) %>%
  dplyr::filter(pct.1 >= 0.20) %>% 
  dplyr::arrange(p_val_adj) %>% # top based on padj
  dplyr::slice_head(.,n = 50) %>% 
  pull(gene)

# Calculate score
NB_expr<- as.data.frame(NB_expr[,NB_clin$`Patient ID`])
counts.se <- SummarizedExperiment(NB_expr, colData = NB_clin)
rankData <- rankGenes(counts.se)
AC_score<- simpleScore(rankData,upSet = top_markers,knownDirection=TRUE)
AC_score$AC_like<- AC_score$TotalScore

# Merge scores with clinical data
NB_clin<- merge(NB_clin, AC_score[,"AC_like", drop=F], by.x = "Patient ID", by.y = 0)
rownames(NB_clin)<- NB_clin$`Patient ID`

# Add expression to clinical data
NB_clin$ALKAL2<- as.numeric(t(NB_expr)[NB_clin$`Patient ID`,"FAM150B"]) # Only Hi & Lo?
NB_clin$NRTN<- as.numeric(t(NB_expr)[NB_clin$`Patient ID`,"NRTN"]) # Only Hi & Lo?
NB_clin$CCL18<- as.numeric(t(NB_expr)[NB_clin$`Patient ID`,"CCL18"]) # Only Hi & Lo?
NB_clin$PITPNM3<- as.numeric(t(NB_expr)[NB_clin$`Patient ID`,"PITPNM3"]) # Only Hi & Lo?

# Save
#######
saveRDS(NB_clin, "temp/cangelosi.rds")
