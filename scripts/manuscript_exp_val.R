# Load data
###########
load("data/ST_NB.RData")

# Plot NRTN
############
p_mig<- ggplot(exp_val_ls$migr$NRTN, aes(x = time_hrs, y = RWD, color = group, group = group)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = RWD - SE, ymax = RWD + SE), width = 0.25) + # Error bars
  facet_wrap(.~CL) +
  xlab("Time (hrs)") +
  ylab("RWD (%)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48)) +
  scale_color_manual(
    values = c("ctrl" = "blue", "NRTN (100 ng/ml)" = "red")) +
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
        legend.position = c(0.1, 0.9),
        strip.text = element_text(size = 7)
  )
 
p_prol<- ggplot(exp_val_ls$prol$NRTN, aes(x = time_hrs, y = relative_cell_growth, color = group, group = group)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = relative_cell_growth - SE, ymax = relative_cell_growth + SE), width = 0.25) + # Error bars
  facet_wrap(.~CL)+
  xlab("Time (hrs)") +
  ylab("Relative cell growth (%)") +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96, 120, 144)) +
  scale_color_manual(
    values = c("ctrl" = "blue", "NRTN (100 ng/ml)" = "red")) +
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
        legend.position = c(0.1, 0.9),
        strip.text = element_text(size = 7)
  )

p<- plot_grid(
  p_prol,
  p_mig,
  ncol = 1,
  rel_heights = c(1,2),
  labels = c("B","C"),
  label_size = 12
)

ggsave("results/figs/manuscript_NRTN_validation.pdf", width = 85, height = 0.4*297, units = "mm")

# Plot CCL18
############
p_mig<- ggplot(exp_val_ls$migr$CCL18, aes(x = time_hrs, y = RWD, color = group, group = group)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = RWD - SE, ymax = RWD + SE), width = 0.25) + # Error bars
  facet_wrap(.~CL) +
  xlab("Time (hrs)") +
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
        legend.position = c(0.1, 0.9),
        strip.text = element_text(size = 7),
        strip.background = element_blank()
  )

p_prol<- ggplot(exp_val_ls$prol$CCL18, aes(x = time_hrs, y = relative_cell_growth, color = group, group = group)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = relative_cell_growth - SE, ymax = relative_cell_growth + SE), width = 0.25) + # Error bars
  facet_wrap(.~CL)+
  xlab("Time (hrs)") +
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
        legend.position = c(0.1, 0.9),
        strip.text = element_text(size = 7),
        strip.background = element_blank()
  )

p<- plot_grid(
  p_prol,
  p_mig,
  ncol = 1,
  labels = c("B","C"),
  label_size = 12
)

ggsave("results/figs/manuscript_CCL18_validation.pdf", width = 85, height = 0.27*297, units = "mm")




