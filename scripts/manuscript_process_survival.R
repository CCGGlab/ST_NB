##################################################################
# Process NB survival data
##################################################################

library(singscore)
library(SummarizedExperiment)
library(tidyverse)
DE_clust_ls<- readRDS(file = "temp/DE_clust.rds")

# AC signature: Get top50 DE genes in AC_like
DE_clust<- DE_clust_ls$hr$NB2Post
DE_clust<- DE_clust[DE_clust$cluster == "AC_like",]
top_markers <- DE_clust %>%
  dplyr::filter(!str_detect(gene,"DEPRECATED")) %>%
  dplyr::filter(avg_log2FC > 0) %>%data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
  dplyr::filter(p_val_adj <= 0.01) %>%
  dplyr::filter(pct.1 >= 0.20) %>%
  dplyr::arrange(p_val_adj) %>% # top based on padj
  dplyr::slice_head(.,n = 50) %>%
  pull(gene)

surv_data_ls<- list()

# Versteeg
###########

NB_expr<- readRDS("downloads/R2/versteeg/versteeg_expr.rds")
NB_clin<- readRDS("downloads/R2/versteeg/versteeg_clin.rds")

# AC signature
NB_expr<- as.data.frame(NB_expr[,NB_clin$sample_id])
counts.se <- SummarizedExperiment(NB_expr, colData = NB_clin)
rankData <- rankGenes(counts.se)
AC_score<- simpleScore(rankData,upSet = top_markers,knownDirection=TRUE) # 3 not found, ALKAL2
AC_score$AC_like<- AC_score$TotalScore

# Merge scores with clinical data
NB_clin<- merge(NB_clin, AC_score[,"AC_like", drop=F], by.x = "sample_id", by.y = 0)
rownames(NB_clin)<- NB_clin$sample_id

# Add expression to clinical data
NB_clin$ALKAL2<- as.numeric(t(NB_expr)[NB_clin$sample_id,"FAM150B"]) # Only Hi & Lo?
NB_clin$ALK<- as.numeric(t(NB_expr)[NB_clin$sample_id,"ALK"]) # Only Hi & Lo?
NB_clin$NRTN<- as.numeric(t(NB_expr)[NB_clin$sample_id,"NRTN"]) # Only Hi & Lo?
NB_clin$CCL18<- as.numeric(t(NB_expr)[NB_clin$sample_id,"CCL18"]) # Only Hi & Lo?
NB_clin$PITPNM3<- as.numeric(t(NB_expr)[NB_clin$sample_id,"PITPNM3"]) # Only Hi & Lo?
NB_clin$RET<- as.numeric(t(NB_expr)[NB_clin$sample_id,"RET"]) # Only Hi & Lo?
NB_clin$GFRA2<- as.numeric(t(NB_expr)[NB_clin$sample_id,"GFRA2"]) # Only Hi & Lo?
NB_clin$FGF9<- as.numeric(t(NB_expr)[NB_clin$sample_id,"FGF9"]) # Only Hi & Lo?
NB_clin$FGF1<- as.numeric(t(NB_expr)[NB_clin$sample_id,"FGF1"]) # Only Hi & Lo?
NB_clin$FGFR1<- as.numeric(t(NB_expr)[NB_clin$sample_id,"FGFR1"]) # Only Hi & Lo?

# Default format
NB_clin$samplenames<-NB_clin$sample_id
NB_clin$isMYCNA<- NA
NB_clin$isMYCNA[NB_clin$MYCN_amp=="yes"]<- T 
NB_clin$isMYCNA[NB_clin$MYCN_amp=="no"]<- F 

NB_clin$dead<- NA 
NB_clin$dead[NB_clin$alive=="no"]<- T
NB_clin$dead[NB_clin$alive=="yes"]<- F
NB_clin$overall_survival<- as.numeric(NB_clin$nti_surv_overall)/365

surv_data_ls$versteeg<- NB_clin

# TARGET
#########

NB_clin<- rbind(
  readxl::read_excel("downloads/TARGET/clin/TARGET_NBL_ClinicalData_Discovery_20170525.xlsx"),
  readxl::read_excel("downloads/TARGET/clin/TARGET_NBL_ClinicalData_Validation_20190129.xlsx")
)
NB_clin<- as.data.frame(NB_clin)
NB_clin$samplenames<- NB_clin$`TARGET USI`

# Expr: Use MA data
NB_expr<- readRDS("downloads/TARGET/mRNA/microArray/expr_df.rds")

# Calculate score
common_samples<- intersect(colnames(NB_expr),NB_clin$`TARGET USI`)
NB_expr<- as.data.frame(NB_expr[,common_samples])
NB_clin<- NB_clin[NB_clin$`TARGET USI`%in%common_samples,]
counts.se <- SummarizedExperiment(NB_expr, colData = NB_clin)
rankData <- rankGenes(counts.se)
AC_score<- simpleScore(rankData,upSet = top_markers,knownDirection=TRUE) # 2 not found, ALKAL2
AC_score$AC_like<- AC_score$TotalScore

# Merge scores with clinical data
NB_clin<- merge(NB_clin, AC_score[,"AC_like", drop=F], by.x = "samplenames", by.y = 0)
rownames(NB_clin)<- NB_clin$samplenames

# Add expression to clinical data
NB_clin$ALKAL2<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FAM150B"]) # Only Hi & Lo?
NB_clin$ALK<- as.numeric(t(NB_expr)[NB_clin$samplenames,"ALK"]) # Only Hi & Lo?
NB_clin$NRTN<- as.numeric(t(NB_expr)[NB_clin$samplenames,"NRTN"]) # Only Hi & Lo?
NB_clin$CCL18<- as.numeric(t(NB_expr)[NB_clin$samplenames,"CCL18"]) # Only Hi & Lo?
NB_clin$PITPNM3<- as.numeric(t(NB_expr)[NB_clin$samplenames,"PITPNM3"]) # Only Hi & Lo?
NB_clin$RET<- as.numeric(t(NB_expr)[NB_clin$samplenames,"RET"]) # Only Hi & Lo?
NB_clin$GFRA2<- as.numeric(t(NB_expr)[NB_clin$samplenames,"GFRA2"]) # Only Hi & Lo?
NB_clin$FGF9<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FGF9"]) # Only Hi & Lo?
NB_clin$FGF1<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FGF1"]) # Only Hi & Lo?
NB_clin$FGFR1<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FGFR1"]) # Only Hi & Lo?

# Default format
NB_clin$isMYCNA<- NA
NB_clin$isMYCNA[NB_clin$`MYCN status`=="Amplified"]<- T 
NB_clin$isMYCNA[NB_clin$MYCN_amp=="Not Amplified"]<- F 

NB_clin$dead<- NA 
NB_clin$dead[NB_clin$`Vital Status`=="Dead"]<- T
NB_clin$dead[NB_clin$`Vital Status`=="Alive"]<- F
NB_clin$overall_survival<- as.numeric(NB_clin$`Overall Survival Time in Days`)/365

surv_data_ls$TARGET<- NB_clin

# Kocak
#######

NB_expr<- readRDS("downloads/R2/kocak/kocak_expr.rds")
NB_expr<- log10(NB_expr) # log scale
NB_clin<- readRDS("downloads/R2/kocak/kocak_clin.rds")

# Calculate score
NB_expr<- as.data.frame(NB_expr[,NB_clin$samplenames])
counts.se <- SummarizedExperiment(NB_expr, colData = NB_clin)
rankData <- rankGenes(counts.se)
AC_score<- simpleScore(rankData,upSet = top_markers,knownDirection=TRUE) # 7 not found, ALKAL2
AC_score$AC_like<- AC_score$TotalScore

# Merge scores with clinical data
NB_clin<- merge(NB_clin, AC_score[,"AC_like", drop=F], by.x = "samplenames", by.y = 0)
rownames(NB_clin)<- NB_clin$samplenames

# Add expression to clinical data
NB_clin$ALKAL2<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FAM150B"]) # Only Hi & Lo?
NB_clin$ALK<- as.numeric(t(NB_expr)[NB_clin$samplenames,"ALK"]) # Only Hi & Lo?
NB_clin$NRTN<- as.numeric(t(NB_expr)[NB_clin$samplenames,"NRTN"]) # Only Hi & Lo?
NB_clin$CCL18<- as.numeric(t(NB_expr)[NB_clin$samplenames,"CCL18"]) # Only Hi & Lo?
NB_clin$PITPNM3<- as.numeric(t(NB_expr)[NB_clin$samplenames,"PITPNM3"]) # Only Hi & Lo?
NB_clin$RET<- as.numeric(t(NB_expr)[NB_clin$samplenames,"RET"]) # Only Hi & Lo?
NB_clin$GFRA2<- as.numeric(t(NB_expr)[NB_clin$samplenames,"GFRA2"]) # Only Hi & Lo?
NB_clin$FGF9<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FGF9"]) # Only Hi & Lo?
NB_clin$FGF1<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FGF1"]) # Only Hi & Lo?
NB_clin$FGFR1<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FGFR1"]) # Only Hi & Lo?

# Default format
NB_clin$isMYCNA<- NA
NB_clin$isMYCNA[NB_clin$mycn=="amp"]<- T 
NB_clin$isMYCNA[NB_clin$mycn=="non_amp"]<- F 

NB_clin$dead<- NA 
NB_clin$dead[NB_clin$vital_status=="diseased"]<- T
NB_clin$dead[NB_clin$vital_status=="alive"]<- F
NB_clin$overall_survival<- as.numeric(NB_clin$Factor.Value.overall.survival.)/365

surv_data_ls$kocak<- NB_clin

# SEQC
######

# Get data
NB_expr<- readRDS("downloads/R2/SEQC/SEQC_expr.rds")
NB_expr<- log10(NB_expr) # log scale
NB_clin<- readRDS("downloads/R2/SEQC/SEQC_clin.rds")

# Calculate score
NB_expr<- as.data.frame(NB_expr[,NB_clin$samplenames])
counts.se <- SummarizedExperiment(NB_expr, colData = NB_clin)
rankData <- rankGenes(counts.se)
AC_score<- simpleScore(rankData,upSet = top_markers,knownDirection=TRUE) # 7 not found, ALKAL2
AC_score$AC_like<- AC_score$TotalScore

# Merge scores with clinical data
NB_clin<- merge(NB_clin, AC_score[,"AC_like", drop=F], by.x = "samplenames", by.y = 0)
rownames(NB_clin)<- NB_clin$samplenames

# Add expression to clinical data
NB_clin$ALKAL2<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FAM150B"]) # Only Hi & Lo?
NB_clin$ALK<- as.numeric(t(NB_expr)[NB_clin$samplenames,"ALK"]) # Only Hi & Lo?
NB_clin$NRTN<- as.numeric(t(NB_expr)[NB_clin$samplenames,"NRTN"]) # Only Hi & Lo?
NB_clin$CCL18<- as.numeric(t(NB_expr)[NB_clin$samplenames,"CCL18"]) # Only Hi & Lo?
NB_clin$PITPNM3<- as.numeric(t(NB_expr)[NB_clin$samplenames,"PITPNM3"]) # Only Hi & Lo?
NB_clin$RET<- as.numeric(t(NB_expr)[NB_clin$samplenames,"RET"]) # Only Hi & Lo?
NB_clin$GFRA2<- as.numeric(t(NB_expr)[NB_clin$samplenames,"GFRA2"]) # Only Hi & Lo?
NB_clin$FGF9<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FGF9"]) # Only Hi & Lo?
NB_clin$FGF1<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FGF1"]) # Only Hi & Lo?
NB_clin$FGFR1<- as.numeric(t(NB_expr)[NB_clin$samplenames,"FGFR1"]) # Only Hi & Lo?

# Default format
NB_clin$isMYCNA<- NA
NB_clin$isMYCNA[NB_clin$mycn_status=="mycn_amp"]<- T 
NB_clin$isMYCNA[NB_clin$mycn_status=="mycn_nonamp"]<- F 

NB_clin$dead<- NA
NB_clin$dead[NB_clin$death_from_disease=="yes"]<- T
NB_clin$dead[NB_clin$death_from_disease=="no"]<- F
NB_clin$isMYCNA<- NB_clin$mycn_status=="mycn_amp"
NB_clin$overall_survival<- NB_clin$os_day/365

surv_data_ls$SEQC<- NB_clin
saveRDS(surv_data_ls, "temp/surv_data.rds")
