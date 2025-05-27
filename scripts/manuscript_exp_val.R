# Load data
###########
load("data/ST_NB.RData")
library(BSDA)

# Plot NRTN
############

# Timelines
p_mig<- ggplot(exp_val_ls$migr$NRTN, aes(x = time_hrs, y = RWD, color = group_name, group = group_name)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = RWD - SE, ymax = RWD + SE), width = 0.25) + # Error bars
  facet_wrap(.~CL, ncol = 4, scales = "free") +
  xlab("Time (h)") +
  ylab("RWD (%)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48)) +
  # scale_linetype_manual(values = c("ctrl" = "solid", "NRTN (100 ng/ml)" = "solid", "selpercatinib (500 nM)" = "dashed", "NRTN (100 ng/ml) + selpercatinib (500 nM)" = "dashed")) + 
  scale_color_manual(
    values = c("ctrl" = "blue", "NRTN (100 ng/ml)" = "red", "selpercatinib (500 nM)" = "black", "NRTN (100 ng/ml) + selpercatinib (500 nM)" = "green")) +
  theme_bw() +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.15, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 7),
        strip.background = element_blank()
  )

p_prol<- ggplot(exp_val_ls$prol$NRTN, aes(x = time_hrs, y = relative_cell_growth, color = group_name, group = group_name)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = relative_cell_growth - SE, ymax = relative_cell_growth + SE), width = 0.25) + # Error bars
  facet_wrap(.~CL, ncol = 4, scales = "free")+
  xlab("Time (h)") +
  ylab("Relative cell growth (%)") +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96, 120, 144)) +
  # scale_linetype_manual(values = c("ctrl" = "solid", "NRTN (100 ng/ml)" = "solid", "selpercatinib (500 nM)" = "dashed", "NRTN (100 ng/ml) + selpercatinib (500 nM)" = "dashed")) + 
  scale_color_manual(
    values = c("ctrl" = "blue", "NRTN (100 ng/ml)" = "red", "selpercatinib (500 nM)" = "black", "NRTN (100 ng/ml) + selpercatinib (500 nM)" = "green")) +
  theme_bw() +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.15, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 7),
        strip.background = element_blank()
  )

p_time<- plot_grid(
  p_mig,
  p_prol,
  ncol = 1,
  labels = c("F","G"),
  label_size = 12
)

ggsave("results/figs/manuscript_NRTN_validation.pdf", p_time, width = 178, height = 0.35*297, units = "mm")

# Bar plot to compare RETi
p_ls<- list()
for(exp in c("migr", "prol")){
  for(CL in c("CLB-GA", "NB1", "SH-SY5Y", "SK-N-AS")){
    exp_sel<- exp_val_ls[[exp]]$NRTN
    if(exp=="migr"){
      t<- 12
      exp_sel$read<- exp_sel$RWD
    }
    if(exp=="prol"){
      t<- 144
      exp_sel$read<- exp_sel$relative_cell_growth
    }
    exp_sel<- exp_sel[exp_sel$CL==CL&exp_sel$time_hrs==t,]
    p_ls[[exp]][[CL]]<- ggbarplot(exp_sel, "group_name", "read",
                                  fill = "steelblue", color = "steelblue") +
      geom_errorbar(aes(ymin=read, ymax=read+SE), width=0.4, alpha=0.9, size=0.5) +
      theme_cowplot() +
      theme(
        axis.title = element_blank(),
        axis.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "italic")
      )
    for(cond in list(ctrl_NRTN=c("ctrl", "NRTN (100 ng/ml)"), selp_NRTNselp=c("selpercatinib (500 nM)", "NRTN (100 ng/ml) + selpercatinib (500 nM)"), NRTN_NRTNselp=c("NRTN (100 ng/ml)", "NRTN (100 ng/ml) + selpercatinib (500 nM)"))){
      if(!cond[2]%in%exp_sel$group_name) next
      ttemp<- tsum.test(var.equal = T,
        mean.x = exp_sel$read[exp_sel$group_name==cond[1]],
        s.x = sqrt(3)*exp_sel$SE[exp_sel$group_name==cond[1]],
        n.x = 3,
        mean.y = exp_sel$read[exp_sel$group_name==cond[2]],
        s.y = sqrt(3)*exp_sel$SE[exp_sel$group_name==cond[2]],
        n.y =3
      )
      cat(CL, exp, paste(cond, collapse = " vs "), "P=", ttemp$p.value, "Diff=", ttemp$estimate[2]-ttemp$estimate[1],"\n")
    }
  }
}
# CLB-GA migr ctrl vs NRTN (100 ng/ml) P= 0.0001277843 Diff= 22.91194 
# CLB-GA migr selpercatinib (500 nM) vs NRTN (100 ng/ml) + selpercatinib (500 nM) P= 0.04778467 Diff= -2.41588 
# CLB-GA migr NRTN (100 ng/ml) vs NRTN (100 ng/ml) + selpercatinib (500 nM) P= 6.548132e-05 Diff= -25.988 
# NB1 migr ctrl vs NRTN (100 ng/ml) P= 0.001819019 Diff= 14.63515 
# NB1 migr selpercatinib (500 nM) vs NRTN (100 ng/ml) + selpercatinib (500 nM) P= 0.8533568 Diff= -0.48032 
# NB1 migr NRTN (100 ng/ml) vs NRTN (100 ng/ml) + selpercatinib (500 nM) P= 0.002040772 Diff= -14.41901 
# SH-SY5Y migr ctrl vs NRTN (100 ng/ml) P= 0.01006616 Diff= 14.95602 
# SK-N-AS migr ctrl vs NRTN (100 ng/ml) P= 0.5403028 Diff= 0.663609 
# CLB-GA prol ctrl vs NRTN (100 ng/ml) P= 0.01641948 Diff= 3.824397 
# CLB-GA prol selpercatinib (500 nM) vs NRTN (100 ng/ml) + selpercatinib (500 nM) P= 0.500499 Diff= 1.226883 
# CLB-GA prol NRTN (100 ng/ml) vs NRTN (100 ng/ml) + selpercatinib (500 nM) P= 0.554559 Diff= -1.144782 
# NB1 prol ctrl vs NRTN (100 ng/ml) P= 0.03312831 Diff= 16.04208 
# NB1 prol selpercatinib (500 nM) vs NRTN (100 ng/ml) + selpercatinib (500 nM) P= 0.4225751 Diff= 5.60811 
# NB1 prol NRTN (100 ng/ml) vs NRTN (100 ng/ml) + selpercatinib (500 nM) P= 0.2494534 Diff= -8.19103 
# SH-SY5Y prol ctrl vs NRTN (100 ng/ml) P= 0.6654153 Diff= -0.238699 
# SK-N-AS prol ctrl vs NRTN (100 ng/ml) P= 0.1287176 Diff= 1.146124 

# p_bar<- plot_grid(
#   p_ls$migr$`CLB-GA`, p_ls$migr$NB1,
#   p_ls$prol$`CLB-GA`, p_ls$prol$NB1,
#   ncol = 1
# )
# 
# p<- plot_grid(
#   p_time, p_bar,
#   ncol=2,
#   rel_widths = c(5,1)
# )
# ggsave("results/figs/manuscript_NRTN_validation.pdf", p, width = 178, height = 0.3*297, units = "mm")

# Plot CCL18
############
p_mig<- ggplot(exp_val_ls$migr$CCL18, aes(x = time_hrs, y = RWD, color = group_name, group = group_name)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = RWD - SE, ymax = RWD + SE), width = 0.25) + # Error bars
  facet_wrap(.~CL, ncol = 4, scales = "free") +
  xlab("Time (h)") +
  ylab("RWD (%)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48)) +
  scale_color_manual(
    values = c("ctrl" = "blue", "CCL18 (100 ng/ml)" = "red")) +
  theme_bw() +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.15, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 7),
        strip.background = element_blank()
  )

p_prol<- ggplot(exp_val_ls$prol$CCL18, aes(x = time_hrs, y = relative_cell_growth, color = group_name, group = group_name)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = relative_cell_growth - SE, ymax = relative_cell_growth + SE), width = 0.25) + # Error bars
  facet_wrap(.~CL, ncol = 4, scales = "free")+
  xlab("Time (h)") +
  ylab("Relative cell growth (%)") +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96, 120, 144)) +
  scale_color_manual(
    values = c("ctrl" = "blue", "CCL18 (100 ng/ml)" = "red")) +
  theme_bw() +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.15, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 7),
        strip.background = element_blank()
  )

# P values?
for(exp in c("migr", "prol")){
  for(CL in c("CLB-GA", "NB1", "SH-SY5Y","SK-N-AS")){
    exp_sel<- exp_val_ls[[exp]]$CCL18
    if(exp=="migr"){
      t<- 48
      exp_sel$read<- exp_sel$RWD
    }
    if(exp=="prol"){
      t<- 144
      exp_sel$read<- exp_sel$relative_cell_growth
    }
    exp_sel<- exp_sel[exp_sel$CL==CL&exp_sel$time_hrs==t,]
    for(cond in list(c("ctrl", "CCL18 (100 ng/ml)"))){
      ttemp<- tsum.test(var.equal = T,
        mean.x = exp_sel$read[exp_sel$group_name==cond[1]],
        s.x = sqrt(3)*exp_sel$SE[exp_sel$group_name==cond[1]],
        n.x = 3,
        mean.y = exp_sel$read[exp_sel$group_name==cond[2]],
        s.y = sqrt(3)*exp_sel$SE[exp_sel$group_name==cond[2]],
        n.y =3
      )
      cat(CL, exp, paste(cond, collapse = " vs "), "P=", ttemp$p.value, "Diff=", ttemp$estimate[2]-ttemp$estimate[1],"\n")
    }
  }
}

# CLB-GA migr ctrl vs CCL18 (100 ng/ml) P= 0.4971924 Diff= -0.83444 
# NB1 migr ctrl vs CCL18 (100 ng/ml) P= 0.03241298 Diff= 19.61586 
# SH-SY5Y migr ctrl vs CCL18 (100 ng/ml) P= 0.9259815 Diff= -0.19525 
# SK-N-AS migr ctrl vs CCL18 (100 ng/ml) P= 0.1169082 Diff= 6.853 
# CLB-GA prol ctrl vs CCL18 (100 ng/ml) P= 0.4390459 Diff= 0.924835 
# NB1 prol ctrl vs CCL18 (100 ng/ml) P= 0.09699363 Diff= 4.11962 
# SH-SY5Y prol ctrl vs CCL18 (100 ng/ml) P= 0.74498 Diff= 0.350347 
# SK-N-AS prol ctrl vs CCL18 (100 ng/ml) P= 0.2605131 Diff= -0.636589 

p<- plot_grid(
  p_mig,
  p_prol,
  ncol = 1,
  labels = c("B","C"),
  label_size = 12
)

ggsave("results/figs/manuscript_CCL18_validation.pdf", width = 178, height = 0.35*297, units = "mm")

