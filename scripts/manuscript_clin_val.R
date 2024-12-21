##################################################################
# Clinical validation of AC-like associations & survival analysis
##################################################################

# Libraries
###########
library(ggpubr)
library(cowplot)
library(reshape2)
library(survival)
library(survminer)
library(tidyverse)

# Load data
############
NB_clin<- readRDS(file="temp/cangelosi.rds")

# Correlation AC-like Vs. gene expression
#########################################

# Stratification AC-like: Hi (>P90), Im, Lo (<P90)
th<- 0.9
NB_clin$goo<- as.numeric(NB_clin$AC_like)

NB_clin$AC_class<- NA
NB_clin$AC_class[NB_clin$goo < quantile(NB_clin$goo, 1-th, na.rm=T)]<- "Lo"
NB_clin$AC_class[NB_clin$goo >= quantile(NB_clin$goo, 1-th, na.rm=T)]<- "Im"
NB_clin$AC_class[NB_clin$goo >= quantile(NB_clin$goo, th, na.rm=T)]<- "Hi"

# Merge data in dataframe
res_df<- data.frame(
  AC_class = c(
    names(tapply(NB_clin$ALKAL2>1,  NB_clin$AC_class, "mean")),
    names(tapply(NB_clin$NRTN>1,  NB_clin$AC_class, "mean"))),  
  gene = rep(c("ALKAL2", "NRTN"), each = 3),
  expression = c(
    tapply(NB_clin$ALKAL2>1,  NB_clin$AC_class, "mean"),
    tapply(NB_clin$NRTN>1,  NB_clin$AC_class, "mean"))
)

# add p values
p_labels = c(
  NA, paste0("P = ", signif(chisq.test(table(FZHi = NB_clin$AC_class, ALKAL2_expr = NB_clin$ALKAL2>1))$p.value,2)),NA , # 5e-17
  NA, paste0("P = ", signif(chisq.test(table(FZHi = NB_clin$AC_class, NRTN_expr = NB_clin$NRTN>1))$p.value,2)),NA # p-value = 7.6e-09
)
res_df$p_labels<- p_labels

# Plot
res_df$expr_group_max<- tapply(res_df$expression, res_df$gene, "max")[res_df$gene]

p<- ggbarplot(res_df, "AC_class", "expression",
              fill = "steelblue", color = "steelblue") +
  facet_wrap(.~gene, scales = "free_y") +
  geom_text(aes(x = 2, y = expr_group_max, label=p_labels), vjust = 1, size = 6, size.unit = "pt") +
  ylab("% expressed") +
  xlab("AC-like signature") +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "italic")
  )
# p

res_df
#         AC_class   gene expression    p_labels expr_group_max
# 1       Hi ALKAL2 0.97297297        <NA>      0.9729730
# 2       Im ALKAL2 0.22184300 P = 1.4e-20      0.9729730
# 3       Lo ALKAL2 0.18918919        <NA>      0.9729730
# 4       Hi   NRTN 0.32432432        <NA>      0.3243243
# 5       Im   NRTN 0.04095563   P = 1e-10      0.3243243
# 6       Lo   NRTN 0.00000000        <NA>      0.3243243

ggsave("results/figs/manuscript_clin_val_AC_like.pdf", p, width = 0.375*178, height = 0.16*265, units = "mm")

# Survival CCL18-PITPNM3
#########################
p_res_ls<- list()

# Select data
genes<- c("CCL18", "PITPNM3")

# Create dataframe
cor_df<- data.frame(
  pt = NB_clin$`Patient ID`,
  ligand = NB_clin[,genes[1]],
  receptor = NB_clin[,genes[2]]
)
cor_df$MYCN_status<- factor(NB_clin[as.character(cor_df$pt),"MYCN status"], levels = c("no", "yes"))
cor_df<- na.omit(cor_df)

# Correlation to MYCN amp
colnames(cor_df)[2:3]<- genes
cor_df<- melt(cor_df, value.name = "expression")
bplot<- ggboxplot(cor_df, x = "MYCN_status", y = "expression",fill = "MYCN_status", palette = "jco", outlier.shape = NA, size = 0.5, facet.by = "variable") + 
  scale_x_discrete(name ="MYCN", labels=c("no" = "WT", "yes" = "Amp")) +
  stat_compare_means(label = "p.format", label.x.npc = "center", hjust = 0.5, vjust = 1, size = 6, size.unit = "pt") +
  theme(
    legend.position = "none",
    axis.text = element_text(size=6),
    axis.title = element_text(size=7, face = "italic"),
    strip.text = element_text(size=7, face = "italic"),
    strip.background = element_blank()
  )
p_res_ls[["bplot"]]<- bplot 

# Default var names
NB_clin$dead<- NB_clin$`Event overall`=="yes"
NB_clin$overall_survival<- as.numeric(NB_clin$`Overall survival`)
NB_clin$isMYCNA<- NB_clin$`MYCN status`=="yes"

# KM
surv_cols<- c(rgb(231/256, 184/256, 0/256),"black",rgb(46/256, 159/256, 223/256))
th<- 0.9
for(g in genes){
  NB_clin[[paste0("high",g)]]<- NB_clin[[g]] > median(NB_clin[[g]], na.rm=T)
  NB_clin[[paste0("survClass",g)]]<- NA
  NB_clin[[paste0("survClass",g)]][NB_clin[[g]] < quantile(NB_clin[[g]], 1-th, na.rm=T)]<- "Lo" 
  NB_clin[[paste0("survClass",g)]][NB_clin[[g]] >= quantile(NB_clin[[g]], 1-th, na.rm=T)]<- "Im" 
  NB_clin[[paste0("survClass",g)]][NB_clin[[g]] >= quantile(NB_clin[[g]], th, na.rm=T)]<- "Hi" 
  NB_clin[[paste0("survClass",g)]]<- factor(NB_clin[[paste0("survClass",g)]], levels = c("Lo","Im","Hi"))
  # Survival: univariate KM
  fit <- survfit(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", "survClass",g)), data = NB_clin)
  legend_text<- paste0(
    c(paste0("<P",100*(1-th)), paste0("P",100*(1-th),"-","P",100*(th)), paste0(">P",100*th)),
    paste0(" (n=",fit$n,")")
  )
  p_res_ls[[paste0("p_surv_",g)]]<- ggsurvplot(fit, pval = TRUE, conf.int = T, xlab = "Time (years)", legend.labs=legend_text, legend.title="", palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=6, font.title=8, legend=c(0.8,0.9), title=g, censor=F)$plot
}

# Cox multivariate
# Relevel to Im ref
for(g in genes) NB_clin[[paste0("survClass",g)]]<- factor(NB_clin[[paste0("survClass",g)]], levels = c("Im","Lo", "Hi"))
fit.coxph <- coxph(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ survClass", genes[1], " + survClass", genes[2]," + isMYCNA")), data = NB_clin)

HR_multi_df<- data.frame(
  cond = rownames(summary(fit.coxph)$coefficients),
  HR = summary(fit.coxph)$conf.int[, 1],
  ci_l = summary(fit.coxph)$conf.int[, 3],
  ci_h = summary(fit.coxph)$conf.int[, 4],
  p = summary(fit.coxph)$coefficients[,5]
)
HR_multi_df$var1<- gsub(c("Lo|Hi|TRUE"),"",HR_multi_df$cond)
HR_multi_df$var2<- gsub(paste0("survClass|",genes[1],"|",genes[2],"|isMYCNA"),"",HR_multi_df$cond)
HR_multi_df$n<- apply(HR_multi_df[c("var1","var2")], 1, function(x) table(NB_clin[[x[1]]])[x[2]])

HR_multi_df$cond_name<- paste0(gsub("survClass", "", HR_multi_df$var1), " ",gsub("TRUE","",HR_multi_df$var2), " (n=",HR_multi_df$n,")")
HR_multi_df$cond_name<- factor(HR_multi_df$cond_name, levels = rev(HR_multi_df$cond_name))

fp_multi <- ggplot(data=HR_multi_df, aes(x=cond_name, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)))) +
  geom_pointrange(size=0.25, fatten = 10, shape = 18) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  geom_text(aes(y=HR, label=signif(HR,2)), vjust=-0.8, size=2) +
  geom_text(aes(y=20, label=paste0("P=",signif(p,2))), fontface="italic", hjust=1, vjust=1.2, size=6/3) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  # facet_grid(.~gene) +
  xlab("") +
  ylab("Hazard Ratio (Mean +/- 95% CI)") +
  scale_y_continuous(limits = c(0.1,20), trans='log10') +
  theme_bw() +  # use a white background
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size=7),
    axis.title = element_text(size=7)
  )
p_res_ls[["fp_multi"]]<- fp_multi

p_surv<- plot_grid(
  p_res_ls$bplot,
  p_res_ls$fp_multi,
  ncol=1,
  align = "hv",  
  axis = "lrtb",    
  labels = c("B","C"),
  label_size = 12
)

p_surv<- plot_grid(
  plot_grid(p_res_ls$p_surv_CCL18 + theme(legend.position = "none"), p_res_ls$p_surv_PITPNM3, ncol=2),
  p_surv,
  ncol = 1,
  rel_heights = c(1,2),
  labels = c("A",NA,NA),
  label_size = 12
)

ggsave(paste0("results/figs/manuscript_clin_val_surv.pdf"),  p_surv, width = 178/2, height = 265/2, units = "mm")
  
