# Load & process adrenal development data
############################################
library(tidyverse)
library(cowplot)

eset<- readRDS("temp/mtab_5525.rds") # Data derived from Array Express E-MTAB-5525

genes = fData(eset)$gene
.expression = exprs(eset)
.pheno = pData(eset)

# Carnegie stages
stages_dpc = c(
  CS17 = 42.5,
  CS18 = 45.75,
  CS19 = 49,
  CS20 = 51.25,
  CS21 = 53,
  CS22 = 55,
  CS23 = 57,
  F1 = 60.5,
  F2 = 66.5,
  F3 = 73.5
)

carnegie_stage_fmt <- function(x) {
  as_factor(x) %>%
    fct_relabel(~ str_replace(., "carnegie stage ", "")) %>%
    lvls_expand(names(stages_dpc)) %>%
    fct_relevel(names(stages_dpc))
}

# Pheno minimal
pheno = .pheno %>%
  as_tibble(rownames = "key") %>%
  # Keep factor variables and cleanup names
  rename_all(~ str_replace(., "\\.$", "")) %>%
  select(key, starts_with("Factor.Value.")) %>%
  rename_all(~ str_replace(., "Factor.Value.", "")) %>%
  # Order stages
  mutate_at("developmental.stage", carnegie_stage_fmt) %>%
  mutate(adrenal = organism.part == "adrenal gland") %>%
  mutate_at("adrenal",
            ~ fct_recode(as.factor(.), adrenal = "TRUE", control = "FALSE")) %>%
  mutate(dpc = map_dbl(developmental.stage, ~ stages_dpc[[.]]))

# Plots
stacked_data = .expression %>%
  as_tibble %>%
  bind_cols("gene" = genes) %>%
  gather("key", "expression",-gene) %>%
  inner_join(pheno, by = "key") %>%
  drop_na()


# Extract AC like signature
######################################

# BiocManager::install("GSVA")
library(GSVA)  
  
DE_clust<- readxl::read_excel("temp/de_hr_filtered.xlsx",sheet = "NB2Post")
DE_clust<- DE_clust[DE_clust$cluster=="FZ_like",]

geneSets<- list(
  AC_like = DE_clust$gene[1:50],
  ALK = "ALK",
  ALKAL2 = "FAM150B"
)

loc<- unique(stacked_data$key)  
genes<- unique(stacked_data$gene)  
expr_matrix<- matrix(NA, length(loc), length(genes), dimnames = list(loc, genes))
for(i in 1:nrow(expr_matrix)){
  data_tmp<- as.data.frame(stacked_data[stacked_data$key==rownames(expr_matrix)[i],])
  data_tmp<- data_tmp[!duplicated(data_tmp$gene),]
  rownames(data_tmp)<- data_tmp$gene
  expr_matrix[i,]<- data_tmp[colnames(expr_matrix),"expression"]
}

gsva_es<- as.data.frame(t(gsva(data.matrix(t(expr_matrix)), geneSets, method="ssgsea", ssgsea.norm=F)))

gsva_es_merged<- merge(pheno, gsva_es, by.x = "key", by.y = 0)

# Only focus on FZ & ALK/ALKAL in adrenal gland
gsva_es_merged<- gsva_es_merged[gsva_es_merged$organism.part%in%c("adrenal gland"),]

p_ls<- list()
for(gene in c("ALK", "ALKAL2", "AC_like")){
  gsva_es_merged$expr<- gsva_es_merged[,gene]
  
  df_corr = gsva_es_merged %>%
    group_by(organism.part) %>%
    do(
      cor.test(.$expr, .$dpc, method = "pearson") %>%
        map(~ unlist(., use.names = F)) %>%
        map_if(~ length(.) > 1, list) %>%
        as_tibble
    ) %>%
    ungroup %>%
    mutate(lbl = sprintf('%s (r = %0.2f, p = %0.1e)', organism.part, estimate, p.value))
  
  p_ls[[gene]]<- ggplot(gsva_es_merged, aes(
    x = dpc,
    y = expr,
    fill = organism.part,
    col = organism.part
  )) +
    geom_point() +
    geom_smooth(method = 'lm', se = F) +
    labs(x = "days post conception",
         y = "log2 expression intensity",
         title = gene) +
    scale_color_discrete(labels = deframe(select(df_corr, organism.part, lbl))) +
    scale_fill_discrete(labels = deframe(select(df_corr, organism.part, lbl))) +
    theme_bw() +
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 6),
          plot.title = element_text(size = 8),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.height = unit(0.15, "mm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom",
    )
  
}


p_ALK<- plot_grid(
  p_ls$ALK,
  p_ls$ALKAL2,
  p_ls$AC_like,
  ncol=3
)
ggsave("results/figs/manuscript_figS6_adrenal_dev_ALK_ALKAL_AC.pdf", p_ALK, width = 0.5*210, height = 0.15*297, units = "mm")


