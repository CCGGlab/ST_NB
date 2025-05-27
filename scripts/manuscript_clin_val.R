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
load("data/ST_NB.RData")

# Perform meta-analysis on all survival data: get hi/lo expression groups per study
####################################################################################

genes<- c("AC_like","ALKAL2", "ALK", "CCL18", "PITPNM3", "NRTN", "RET", "GFRA2", "FGF1", "FGFR1", "FGF9", "FGFR1")
th<- 0.9

z_score_normalization <- function(x) {return ((x - mean(x)) / sd(x))}
NB_clin_all<- NULL
HR_multi_df_all<- NULL

# Merge studies & create per study Cox analysis
NB_clin_ls$cangelosi<- NULL # Exclude because of 100% overlap with SEQC
for(s in names(NB_clin_ls)){
  
  # Data per study
  NB_clin<- NB_clin_ls[[s]]
  for(g in genes){
    # Define expression classes
    NB_clin[[paste0("survClass",g)]]<- NA
    NB_clin[[paste0("survClass",g)]][NB_clin[[g]] < quantile(NB_clin[[g]], 1-th, na.rm=T)]<- "Lo" 
    NB_clin[[paste0("survClass",g)]][NB_clin[[g]] >= quantile(NB_clin[[g]], 1-th, na.rm=T)]<- "Im" 
    NB_clin[[paste0("survClass",g)]][NB_clin[[g]] >= quantile(NB_clin[[g]], th, na.rm=T)]<- "Hi" 
    NB_clin[[paste0("survClass",g)]]<- factor(NB_clin[[paste0("survClass",g)]], levels = c("Lo","Im","Hi"))
    
    # Z scores
    NB_clin[[paste0(g,"_z")]]<- z_score_normalization(NB_clin[[g]])
    
    # Cox: reference Im
    NB_clin[[paste0("survClass",g)]]<- factor(NB_clin[[paste0("survClass",g)]], levels = c("Im","Lo", "Hi"))
    fit.coxph <- coxph(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ survClass", g," + isMYCNA")), data = NB_clin)
    
    HR_multi_df<- data.frame(
      cond = rownames(summary(fit.coxph)$coefficients),
      HR = summary(fit.coxph)$conf.int[, 1],
      ci_l = summary(fit.coxph)$conf.int[, 3],
      ci_h = summary(fit.coxph)$conf.int[, 4],
      p = summary(fit.coxph)$coefficients[,5]
    )
    HR_multi_df$var1<- gsub(c("Lo|Hi|TRUE"),"",HR_multi_df$cond)
    HR_multi_df$var2<- gsub(paste0("survClass|",g,"|isMYCNA"),"",HR_multi_df$cond)
    HR_multi_df$n_tot<- nrow(NB_clin)
    HR_multi_df$n<- apply(HR_multi_df[c("var1","var2")], 1, function(x) table(NB_clin[[x[1]]])[x[2]])
    HR_multi_df$study<- s
    HR_multi_df$gene<- g
    HR_multi_df_all<- rbind(HR_multi_df_all, HR_multi_df)
  }
  NB_clin<- NB_clin[,c(grep("survClass|overall_survival|dead|isMYCNA|_z",colnames(NB_clin)))]
  NB_clin$study<- s
  
  if(is.null(NB_clin_all)) NB_clin_all<- NB_clin
  else NB_clin_all<- rbind(NB_clin_all, NB_clin)
}

# Cox on metadata
for(g in genes){
  
  # Cox multivariate, relevel to Im ref
  NB_clin_all[[paste0("survClass",g)]]<- factor(NB_clin_all[[paste0("survClass",g)]], levels = c("Im","Lo", "Hi"))
  fit.coxph <- coxph(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ survClass", g, "+ isMYCNA")), data = NB_clin_all)
  
  HR_multi_df<- data.frame(
    cond = rownames(summary(fit.coxph)$coefficients),
    HR = summary(fit.coxph)$conf.int[, 1],
    ci_l = summary(fit.coxph)$conf.int[, 3],
    ci_h = summary(fit.coxph)$conf.int[, 4],
    p = summary(fit.coxph)$coefficients[,5]
  )
  
  HR_multi_df$var1<- gsub(c("Lo|Hi|TRUE"),"",HR_multi_df$cond)
  HR_multi_df$var2<- gsub(paste0("survClass|",g,"|isMYCNA"),"",HR_multi_df$cond)
  HR_multi_df$n<- apply(HR_multi_df[c("var1","var2")], 1, function(x) table(NB_clin_all[[x[1]]])[x[2]])
  HR_multi_df$n_tot<- nrow(NB_clin_all)
  HR_multi_df$gene<- g
  HR_multi_df$study<- "meta"
  HR_multi_df_all<- rbind(HR_multi_df_all, HR_multi_df)
}

# Plot forest plots for all genes & studies (including meta)
############################################################

surv_cols<- c(rgb(46/256, 159/256, 223/256), rgb(231/256, 184/256, 0/256))
HR_multi_df_all$study_label<- paste0(HR_multi_df_all$study, " (n=",HR_multi_df_all$n_tot,")")
HR_multi_df_all$study_label<- factor(HR_multi_df_all$study_label, levels = rev(c("kocak (n=649)","SEQC (n=498)","TARGET (n=249)","versteeg (n=88)","meta (n=1484)")))

fp_multi <- ggplot(data=subset(HR_multi_df_all, var1!="isMYCNA") , aes(x=study_label, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)), col = var2)) +
  geom_pointrange(size=0.25, fatten = 5, shape = 18) + 
  scale_color_manual(values = surv_cols) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  geom_text(aes(y=HR, label=signif(HR,2)), vjust=-0.8, size=2) +
  geom_text(aes(y=20, label=paste0("P=",signif(p,2))), fontface="italic", hjust=1, vjust=1.2, size=6/3) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  facet_wrap(gene~., scales = "free", ncol = 2) +
  xlab("") +
  ylab("Hazard Ratio (Mean +/- 95% CI)") +
  scale_y_continuous(limits = c(0.01,20), trans='log10') +
  theme_bw() +  # use a white background
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size=8, face = "italic"),
    axis.text = element_text(size=7),
    axis.title = element_text(size=7)
  )
# fp_multi

ggsave(paste0("results/figs/manuscript_figS7_clin_val_forest_plots.pdf"),  fp_multi, width = 178, height = 265, units = "mm")

# AC-like survival analysis
############################

# Survival analysis
surv_ls<- list()
for(g in genes){
  
  surv_cols<- c(rgb(231/256, 184/256, 0/256),"black",rgb(46/256, 159/256, 223/256))
  NB_clin_all[[paste0("survClass",g)]]<- factor(NB_clin_all[[paste0("survClass",g)]], levels = c("Lo","Im","Hi"))
  fit <- survfit(as.formula(paste0("Surv(time = overall_survival, event = dead) ~ ", "survClass",g)), data = NB_clin_all)
  legend_text<- paste0(
    c(paste0("<P",100*(1-th)), paste0("P",100*(1-th),"-","P",100*(th)), paste0(">P",100*th)),
    paste0(" (n=",fit$n,")")
  )
  p_surv<- ggsurvplot(fit, pval = TRUE, conf.int = T, xlab = "Time (years)", legend.labs=legend_text, legend.title=g, palette =  surv_cols, font.x=7, font.y=7, font.tickslab=6, pval.size=7/3, font.legend=6, font.title=8, legend=c(0.8,0.9), title=g, censor=F)$plot
  surv_ls[[g]]<- p_surv
}

# Survival AC-like
ggsave("results/figs/manuscript_fig3F_clin_val_AC_like.pdf", surv_ls$AC_like, width = 0.3*178, height = 0.2*265, units = "mm")

# Survival NRTN
p_NRTN<- plot_grid(
  surv_ls$NRTN,
  surv_ls$RET,
  surv_ls$GFRA2,
  ncol=1
)
ggsave("results/figs/manuscript_fig4D_clin_val_NRTN.pdf", p_NRTN, width = 0.2*178, height = 0.4*265, units = "mm")

# AC-like correlations
#########################
res_t_ALKAL2<- table(AC=NB_clin_all$survClassAC_like, ALKAL2=NB_clin_all$survClassALKAL2=="Hi")
signif(chisq.test(res_t_ALKAL2)$p.value,2) # 1.7e-106

res_t_NRTN<- table(AC=NB_clin_all$survClassAC_like, NRTN=NB_clin_all$survClassNRTN=="Hi")
signif(chisq.test(res_t_NRTN)$p.value,2) #9.3e-33

res_t_MYCNA_AC<- table(AC=NB_clin_all$survClassAC_like, MYCNA=NB_clin_all$isMYCNA)
signif(chisq.test(res_t_MYCNA_AC)$p.value,2) # 3.0e-04

res_df<- data.frame(
  AC_like = c("Lo","Im","Hi"),  
  ALKAL2 = 100*prop.table(res_t_ALKAL2,1)[,"TRUE"],  
  NRTN = 100*prop.table(res_t_NRTN,1)[,"TRUE"], 
  MYCNA = 100*prop.table(res_t_MYCNA_AC,1)[,"TRUE"]
)  
res_df$AC_like<- factor(res_df$AC_like, levels = rev(c("Lo","Im","Hi")))
res_df
# AC_like    ALKAL2      NRTN    MYCNA
# Lo      Lo  2.013423  4.026846 22.40000
# Im      Im  4.637437  7.251265 22.22222
# Hi      Hi 61.744966 38.255034  7.03125

# Plot
p_ALKAL2<- ggbarplot(res_df, "AC_like", "ALKAL2",
                     fill = "steelblue", color = "steelblue") +
  ylab("High ALKAL2 expression (%)") +
  xlab("AC-like signature") +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "italic")
  )

p_NRTN<- ggbarplot(res_df, "AC_like", "NRTN",
                   fill = "steelblue", color = "steelblue") +
  ylab("High NRTN expression (%)") +
  xlab("AC-like signature") +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "italic")
  )

p_MYCNA<- ggbarplot(res_df, "AC_like", "MYCNA",
                    fill = "steelblue", color = "steelblue") +
  ylab("MYCN Amplified (%)") +
  xlab("AC-like signature") +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "italic")
  )

p_bp<- plot_grid(
  p_ALKAL2, p_MYCNA,
  ncol = 2
)
ggsave("results/figs/manuscript_fig3G_clin_val_AC_like_cor_ALKAL2_MYCN.pdf", p_bp, width = 0.375*178, height = 0.16*265, units = "mm")
ggsave("results/figs/manuscript_fig4C_clin_val_AC_like_cor_NRTN.pdf", p_NRTN, width = 0.18*178, height = 0.16*265, units = "mm")

# CCL18 meta analysis
#######################

# Create dataframe
cor_df<- data.frame(
  MYCN_status = NB_clin_all$isMYCNA,
  CCL18 = NB_clin_all$CCL18_z,
  PITPNM3 = NB_clin_all$PITPNM3_z
)
cor_df<- na.omit(cor_df)

# Plot
cor_df<- melt(cor_df, value.name = "expression")
p_box<- ggboxplot(subset(cor_df,variable=="CCL18"), x = "MYCN_status", y = "expression",fill = "MYCN_status", palette = "jco", outlier.shape = NA, size = 0.5) + 
  scale_x_discrete(name ="MYCN", labels=c("FALSE" = "WT", "TRUE" = "Amp")) +
  scale_y_continuous(limits = c(-3,3), name = "CCL18 Expression (z score)") +
  stat_compare_means(label = "p.format", label.x.npc = "center", hjust = 0.5, vjust = 1, size = 6, size.unit = "pt") +
  theme(
    legend.position = "none",
    axis.text = element_text(size=6),
    axis.title = element_text(size=7, face = "italic"),
    strip.text = element_text(size=7, face = "italic"),
    strip.background = element_blank()
  )

p_surv<- plot_grid(
  surv_ls$CCL18, p_box,
  ncol = 2,
  rel_widths = c(5,3)
)
ggsave(paste0("results/figs/manuscript_fig6DE_clin_CCL18.pdf"),  p_surv, width = 0.4*178, height = 265/6, units = "mm")

# Forest plot
HR_multi_df_CCL18<- subset(HR_multi_df_all, gene=="CCL18"&study=="meta")
HR_multi_df_CCL18$var_label<- c("CCL18 low (n=149)", "CCL18 high (n=149)", "MYCN Amp (n=269)")
HR_multi_df_CCL18$var_label<- factor(HR_multi_df_CCL18$var_label, levels = c("MYCN Amp (n=269)","CCL18 high (n=149)","CCL18 low (n=149)"))
fp_multi_CCL18 <- ggplot(data= HR_multi_df_CCL18, aes(x=var_label, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)))) +
  geom_pointrange(size=0.25, fatten = 10, shape = 18) + 
  scale_color_manual(values = surv_cols) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  geom_text(aes(y=HR, label=signif(HR,2)), vjust=-0.8, size=2) +
  geom_text(aes(y=20, label=paste0("P=",signif(p,2))), fontface="italic", hjust=1, vjust=1.2, size=6/3) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") +
  ylab("Hazard Ratio (Mean +/- 95% CI)") +
  scale_y_continuous(limits = c(0.3,20), trans='log10') +
  theme_bw() +  # use a white background
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size=8, face = "italic"),
    axis.text = element_text(size=7),
    axis.title = element_text(size=7)
  )
fp_multi_CCL18

ggsave(paste0("results/figs/manuscript_fig6F_clin_CCL18_fp.pdf"),  fp_multi_CCL18, width = 0.4*178, height = 0.1*265, units = "mm")

# NRTN & MYCNA?
###############

# Create dataframe
cor_df<- data.frame(
  MYCN_status = NB_clin_all$isMYCNA,
  NRTN = NB_clin_all$NRTN_z,
  RET = NB_clin_all$RET_z,
  GFRA2 = NB_clin_all$GFRA2_z
)
cor_df<- na.omit(cor_df)

# Plot
cor_df<- melt(cor_df, value.name = "expression")
p_box<- ggboxplot(cor_df, x = "MYCN_status", y = "expression",fill = "MYCN_status", palette = "jco", outlier.shape = NA, size = 0.5, facet.by = "variable") + 
  scale_x_discrete(name ="MYCN", labels=c("FALSE" = "WT", "TRUE" = "Amp")) +
  scale_y_continuous(limits = c(-3,3), name = "Expression (z score)") +
  stat_compare_means(label = "p.format", label.x.npc = "center", hjust = 0.5, vjust = 1, size = 6, size.unit = "pt") +
  theme(
    legend.position = "none",
    axis.text = element_text(size=6),
    axis.title = element_text(size=7, face = "italic"),
    strip.text = element_text(size=7, face = "italic"),
    strip.background = element_blank()
  )

ggsave(paste0("results/figs/manuscript_figx_clin_NRTN.pdf"),  p_box, width = 0.5*178, height = 265/6, units = "mm")

