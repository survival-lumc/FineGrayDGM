# Necessary packages
library(data.table)
library(ggplot2)
library(patchwork)
library(extrafont)

# Set general ggplot2 theme, source in helper functions
theme_set(theme_minimal(base_size = 15, base_family = "Roboto Condensed"))
source("helpers.R")

# Pick timepoints
t <- seq(from = 0.001, to = 10, length.out = 250)

# Colour palette "Hoiho" from Manu (github) package, see https://github.com/G-Thomson/Manu
cols <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70")


# Squeezing ---------------------------------------------------------------


# Set parameters
params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(0.5), # beta_1
    "p" = 0.2, # p
    "base_shape" = 1.25, # a_1
    "base_rate" = 1 # b_1
  ),
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(0.5), # beta_2*
    "base_shape" = 1.5, # a_2
    "base_rate" = 1 # b_2
  )
)

# Compute true hazards and cumincs conditional on X = 0, and X = 1
dat_p <- rbind(
  compute_true(
    t = t,
    model_type = "squeezing",
    newdat = list("X" = 0),
    params = params
  ),
  compute_true(
    t = t,
    model_type = "squeezing",
    newdat = list("X" = 1),
    params = params
  ),
  idcol = "X"
)

# Reformat
dat_p[, ':=' (
  X = factor(X, labels = c(0, 1)),
  cause = factor(cause, levels = c(1, 2))
)]

dat_p_long <- melt.data.table(
  data = dat_p,
  id.vars = c("time", "cause", "X"),
  variable.name = "what",
  value.name = "value"
)

# Panel C
p_squeeze_subdist <- dat_p[X == 0] |>
  ggplot(aes(time, subdist_haz, col = cause, linetype = cause)) +
  geom_line(linewidth = 1.25) +
  coord_cartesian(ylim = c(0, 2.5), xlim = c(0, 10)) +
  scale_color_manual(values = cols[c(1, 2)]) +
  labs(
    x = "Time",
    y = "Baseline subdistribution hazard",
    col = "Cause",
    linetype = "Cause"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5))+
  scale_linetype_manual(values = c(2, 1))

# Panel A
p_squeeze_cs <- dat_p[X == 0] |>
  ggplot(aes(time, cs_haz, col = cause, linetype = cause)) +
  geom_line(linewidth = 1.25) +
  coord_cartesian(ylim = c(0, 2.5), xlim = c(0, 10)) +
  scale_color_manual(values = cols[c(1, 2)]) +
  labs(
    x = "Time",
    y = "Baseline cause-specific hazard",
    col = "Cause",
    linetype = "Cause"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5))+
  scale_linetype_manual(values = c(2, 1))

# Compute hazard ratios
hr_dat <- dat_p_long[, .(
  ratio = value[X == 1] / value[X == 0]
), by = c("time", "cause", "what")][what != "cuminc"]

# Panel B
p_squeeze_cshr <- hr_dat[what == "cs_haz"] |>
  ggplot(aes(time, ratio, col = cause, linetype = cause)) +
  geom_line(linewidth = 1.25) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 10)) +
  scale_color_manual(values = cols[c(1, 2)]) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(x = "Time", y = expression(h[k](t~'|'~X==1)~'/'~h[k](t~'|'~X==0)),
       col = "Cause", linetype = "Cause") +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5))+
  scale_linetype_manual(values = c(2, 1))

# Panel D
p_squeeze_subdisthr <- hr_dat[what == "subdist_haz"] |>
  ggplot(aes(time, ratio, col = cause, linetype = cause)) +
  geom_line(linewidth = 1.25) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 10)) +
  scale_color_manual(values = cols[c(1, 2)]) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(x = "Time", y = expression(lambda[k](t~'|'~X==1)~'/'~lambda[k](t~'|'~X==0)),
       col = "Cause", linetype = "Cause") +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_linetype_manual(values = c(2, 1))


combined <- (p_squeeze_cs + p_squeeze_cshr) /
  (p_squeeze_subdist + p_squeeze_subdisthr) & theme(legend.position = "bottom")

# Figure 2
combined +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')

ggsave(
  here::here("squeeze_fig.pdf"), # adjust names with manuscript
  dpi = 300,
  units = "in",
  width = 8,
  height = 8,
  device = cairo_pdf
)


# Reduction ---------------------------------------------------------------


# Set parameters
params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(0.5), # beta_1
    "base_shape" = -2, # kappa_1
    "base_rate" = 0.5 # nu_1
  ),
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(0.25), # gamma_2
    "base_shape" = 0.5, # a_2
    "base_rate" = 1.25 # b_2
  )
)


# Compute true hazards and cumincs conditional on X = 0, and X = 1
dat_p <- rbind(
  compute_true(
    t = t,
    model_type = "reduction_factor",
    newdat = list("X" = 0),
    params = params
  ),
  compute_true(
    t = t,
    model_type = "reduction_factor",
    newdat = list("X" = 1),
    params = params
  ),
  idcol = "X"
)

# Reformat
dat_p[, ':=' (
  X = factor(X, labels = c(0, 1)),
  cause = factor(cause, levels = c(1, 2))
)]

dat_p_long <- melt.data.table(
  data = dat_p,
  id.vars = c("time", "cause", "X"),
  variable.name = "what",
  value.name = "value"
)

# Try similar plots
dat_p[, TFP := sum(cuminc), by = c("X", "time")]
dat_p[, cause := factor(cause, levels = c(2, 1))]

# P(D = 1 | X)
p1_x <- 1 - exp(params$cause1$base_rate / params$cause1$base_shape)^exp(params$cause1$betas)

p_reduct_1 <- dat_p[X == 1 & TFP <= 1] |>
  ggplot(aes(time, cuminc, fill = cause)) +
  geom_area(alpha = 0.75) +
  scale_fill_manual(values = cols[c(2, 1)]) +
  annotate(
    "text",
    x = 2,
    y = mean(c(p1_x, 1)),
    label = "F[2](t ~ '|' ~ X == 1)",
    family = "Roboto Condensed",
    parse = TRUE
  )+
  annotate(
    "text",
    x = 2,
    y = mean(c(0, p1_x)),
    label = "F[1](t ~ '|' ~ X == 1)",
    family = "Roboto Condensed",
    parse = TRUE
  ) +
  coord_cartesian(ylim = c(0, 1.1), xlim = c(0, 10)) +
  labs(x = "Time", y = "Cumulative incidence") +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(legend.position = "none")

# Get (approximately) time point at which hazards become negative
t_neg_haz <- min(dat_p[cs_haz < 0, "time"])
dat_p[, grps := fcase(
  time < t_neg_haz & X == 1 & cause == 1, 1,
  time >= t_neg_haz & X == 1 & cause == 1, 2,
  X == 1 & cause == 2, 0
)]
dat_p[, cause := factor(cause, levels = c(1, 2))]

p_reduct_2 <- dat_p[X == 1] |>
  ggplot(aes(time, cs_haz, col = cause, linetype = cause, group = grps)) +
  geom_line(linewidth = 1.25) +
  coord_cartesian(ylim = c(-5, 5), xlim = c(0, 10)) +
  scale_color_manual(values = cols[c(1, 2)]) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "Time",
    y = "Cause-specific hazard",
    col = "Cause",
    linetype = "Cause"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  geom_point(x = t_neg_haz, y = 0, size = 3, col = cols[1]) +
  annotate(
    "text",
    x = t_neg_haz + 0.5,
    y = 1.25,
    hjust = 0,
    label = "h[1](t ~ '|' ~ X == 1) ~ 'undefined'",
    family = "Roboto Condensed",
    parse = TRUE
  ) +
  geom_curve(
    aes(
      x = t_neg_haz + 0.4,
      y = 1.25,
      xend = t_neg_haz,
      yend = 0.1
    ),
    colour = "black",
    linewidth = 0.25,
    curvature = 0.3,
    arrow = arrow(length = unit(0.02, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  scale_linetype_manual(values = c(2, 1))

p_reduct_3 <- dat_p_long[, .(
  ratio = value[X == 1] / value[X == 0]
), by = c("time", "cause", "what")][what == "subdist_haz" & time <= t_neg_haz] |>
  ggplot(aes(time, ratio, col = cause, linetype = cause)) +
  geom_line(linewidth = 1.25) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 10)) +
  scale_color_manual(values = cols[c(1, 2)]) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(
    x = "Time",
    y = expression(lambda[k](t~'|'~X==1)~'/'~lambda[k](t~'|'~X==0)),
    col = "Cause",
    linetype = "Cause"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_linetype_manual(values = c(2, 1))

# Figure 1
p_reduct_1 / p_reduct_2 / p_reduct_3 + plot_annotation(tag_levels = 'A')

ggsave(
  here::here("reduc_fig.pdf"),
  dpi = 300,
  units = "in",
  width = 8,
  height = 11,
  device = cairo_pdf
)


# Two FGs -----------------------------------------------------------------


params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(0.5), # beta_1
    "p" = 0.25, # p_10
    "base_shape" = 0.75, # a_1
    "base_rate" = 1 # b_1
  ),
  "cause2" = list(
    "formula" = ~ X, # beta_2
    "betas" = c(0.5),
    "p" = 0.5, # p_20
    "base_shape" = 0.75, # a_2
    "base_rate" = 1 # b2
  )
)

dat_twofgs <- compute_true(
  t = t,
  model_type = "two_fgs",
  newdat = list(X = 0),
  params = params
)

dat_twofgs_x1 <- compute_true(
  t = t,
  model_type = "two_fgs",
  newdat = list(X = 1),
  params = params
)

p_x1 <- dat_twofgs[cause == 1] |>
  ggplot(aes(time, cuminc)) +
  geom_area(fill = cols[1], alpha = 0.75) +
  geom_ribbon(
    data = dat_twofgs[cause == 2],
    aes(ymax = cuminc + params$cause1$p, ymin = params$cause1$p),
    fill = cols[2],
    alpha = 0.75
  ) +
  coord_cartesian(ylim = c(0, 1.1), xlim = c(0, 11)) +
  geom_segment(
    y = params$cause1$p, yend = params$cause1$p, x = 0, xend = 10,
    linetype = "dashed"
  ) +
  geom_segment(
    y = params$cause1$p + params$cause2$p,
    yend = params$cause1$p + params$cause2$p, x = 0, xend = 10,
    linetype = "dashed"
  ) +
  geom_ribbon(
    aes(ymin = params$cause1$p + params$cause2$p,
        ymax = 1),
    alpha = 0.75,
    fill = cols[5]
  ) +
  annotate("text", x = 10.05, y = params$cause1$p, parse = TRUE, label = "p[10]",
           family = "Roboto Condensed", hjust = 0) +
  annotate("text", x = 10.05, y = params$cause1$p + params$cause2$p,
           parse = TRUE, label = "p[10] + p[20]", family = "Roboto Condensed",
           hjust = 0) +
  annotate("text", x = 5, y = mean(c(1, params$cause1$p + params$cause2$p)),
           label = "'Cured'", family = "Roboto Condensed", hjust = 0) +
  annotate(
    "text",
    x = 5,
    y = mean(c(params$cause1$p, params$cause1$p + params$cause2$p)),
    label = as.character(expression(F[20](t))),
    family = "Roboto Condensed",
    parse = TRUE
  ) +
  annotate(
    "text",
    x = 5,
    y = mean(c(params$cause1$p, 0)),
    label = as.character(expression(F[10](t))),
    family = "Roboto Condensed",
    parse = TRUE
  ) +
  labs(x = "Time", y = "Cumulative incidence") +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))

p_x2 <-dat_twofgs_x1[cause == 1] |>
  ggplot(aes(time, cuminc)) +
  geom_area(fill = cols[1], alpha = 0.75) +
  geom_ribbon(
    data = dat_twofgs_x1[cause == 2],
    aes(
      ymax = pmin(1, cuminc + 1 - (1 - params$cause1$p)^exp(params$cause1$betas)),
      ymin = 1 - (1 - params$cause1$p)^exp(params$cause1$betas)
    ),
    fill = cols[2],
    alpha = 0.75
  ) +
  geom_ribbon(
    data = dat_twofgs_x1[cause == 2],
    aes(
      ymax = pmax(1, cuminc + 1 - (1 - params$cause1$p)^exp(params$cause1$betas)),
      ymin = 1
    ),
    fill = cols[2],
    alpha = 0.5
  ) +
  coord_cartesian(ylim = c(0, 1.1), xlim = c(0, 11)) +
  geom_segment(
    y = 1 - (1 - params$cause1$p)^exp(params$cause1$betas),
    yend = 1 - (1 - params$cause1$p)^exp(params$cause1$betas), x = 0, xend = 10,
    linetype = "dashed"
  ) +
  annotate("text", x = 10.05,
           y = 1 - (1 - params$cause1$p)^exp(params$cause1$betas),
           parse = TRUE, label = "p[1](x)",
           family = "Roboto Condensed",
           hjust = 0) +
  geom_curve(
   aes(
      x = 4.5,
      y = 1.075,
      xend = 5,
      yend = 1.025
    ),
    colour = "black",
    linewidth = 0.25,
    curvature = -0.3,
    arrow = arrow(length = unit(0.02, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  annotate("text", x = 4.5,
           y = 1.075,
           hjust = 1,
           label = "TFP > 1",
           family = "Roboto Condensed") +
  labs(x = "Time", y = "Cumulative incidence") +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  annotate(
    "text",
    x = 5,
    y = mean(c(1 - (1 - params$cause1$p)^exp(params$cause1$betas), 1)),
    #label = "F[2](t*'|'*X)",
    label = "F[2](t ~ '|' ~ X == 1)",
    family = "Roboto Condensed",
    parse = TRUE
  )+
  annotate(
    "text",
    x = 5,
    y = mean(c(1 - (1 - params$cause1$p)^exp(params$cause1$betas), 0)),
    label = "F[1](t ~ '|' ~ X == 1)",
    family = "Roboto Condensed",
    parse = TRUE
  )

# Figure 3
p_x1 + p_x2 + plot_annotation(tag_levels = 'A')

ggsave(
  here::here("two_fgs_fig.pdf"), # also to eps/tiff?
  dpi = 300,
  units = "in",
  width = 11,
  height = 7,
  device = cairo_pdf
)

