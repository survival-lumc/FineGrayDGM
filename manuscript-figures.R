# Set-up ------------------------------------------------------------------

library(data.table) # Data wrangling
library(ggplot2) # Making plots
library(patchwork) # Arranging plots in a grid
library(extrafont) # To register Roboto Condensed font

# Set general ggplot2 theme, load helper functions
theme_set(theme_minimal(base_size = 15, base_family = "Roboto Condensed"))
source("helpers.R")

# Pick timepoints for plotting
t <- seq(from = 0.001, to = 10, length.out = 250)

# Colour palette "Hoiho" from {Manu} (Github) package,
# (see https://github.com/G-Thomson/Manu)
cols <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70")


# Figure 1 ----------------------------------------------------------------


# DGM based on specifying subdistribution hazard for event 1, the cause-specific
# hazard for event 2, and deriving the cause-specific hazard of event 1 using
# the reduction factor.

# Set parameters:
# - Fine-Gray model for cause 1, baseline Gompertz subdistribution hazard
# - Cause-specific Cox model cause 2, Weibull baseline hazard
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

# Compute true hazards and cumulative incidence functions, conditional on X = 0,
# and conditional on X = 1, for both causes
dat_reduc <- rbind(
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

# Label and reformat the data
dat_reduc[, X := factor(X, labels = c(0, 1))]
dat_reduc_long <- melt.data.table(
  data = dat_reduc,
  id.vars = c("time", "cause", "X"),
  variable.name = "what",
  value.name = "value"
)
dat_reduc[, TFP := sum(cuminc), by = c("X", "time")]
dat_reduc[, cause := factor(cause, levels = c(2, 1))]

# P(D = 1 | X), for labelling the cumulative incidence functions
p1_x <- 1 - exp(params$cause1$base_rate / params$cause1$base_shape)^exp(params$cause1$betas)

# Panel A
p_reduct_1 <- dat_reduc[X == 1 & TFP <= 1] |>
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
  ) +
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

# Get (approximately) time point at which hazards become negative, for later
# use in plots
t_neg_haz <- min(dat_reduc[cs_haz < 0, "time"])

# This is to avoid a line when it h_1(t | X) is undefined
dat_reduc[, grps := fcase(
  time < t_neg_haz & X == 1 & cause == 1, 1,
  time >= t_neg_haz & X == 1 & cause == 1, 2,
  X == 1 & cause == 2, 0
)]
dat_reduc[, cause := factor(cause, levels = c(1, 2))]

# Panel B
p_reduct_2 <- dat_reduc[X == 1] |>
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

# Panel C
p_reduct_3 <- dat_reduc_long[, .(
  ratio = value[X == 1] / value[X == 0]
), by = c("time", "cause", "what")][what == "subdist_haz" & time <= t_neg_haz] |>
  ggplot(aes(time, ratio, col = cause, linetype = cause)) +
  geom_line(linewidth = 1.25) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 10)) +
  scale_color_manual(values = cols[c(1, 2)]) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(
    x = "Time",
    y = expression(lambda[k](t ~ "|" ~ X == 1) ~ "/" ~ lambda[k](t ~ "|" ~ X == 0)),
    col = "Cause",
    linetype = "Cause"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_linetype_manual(values = c(2, 1))

# Figure 1
p_reduct_1 / p_reduct_2 / p_reduct_3 + plot_annotation(tag_levels = "A")

ggsave(
  here::here("reduc_fig.pdf"),
  dpi = 300,
  units = "in",
  width = 8,
  height = 11,
  device = cairo_pdf
)


# Figure 2 ----------------------------------------------------------------


# DGM based on specifying cumulative incidence function for event 1
# (a Fine-Gray model), and "squeezing" the cumulative incidence function
# for event 2 into the remaining probability space

# Set parameters:
# - Fine-Gray model for cause 1, baseline cumulative incidence function is a
# Weibull CDF
# - For cause 2, a Weibull hazard used for the cumulative incidence function
# conditional on D = 2
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

# Compute true hazards and cumulative incidence functions
dat_squeeze <- rbind(
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

# Reformat and label
dat_squeeze[, ":="(
  X = factor(X, labels = c(0, 1)),
  cause = factor(cause, levels = c(1, 2))
)]

dat_squeeze_long <- melt.data.table(
  data = dat_squeeze,
  id.vars = c("time", "cause", "X"),
  variable.name = "what",
  value.name = "value"
)

# Panel A
p_squeeze_cs <- dat_squeeze[X == 0] |>
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
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_linetype_manual(values = c(2, 1))

# Panel C
p_squeeze_subdist <- dat_squeeze[X == 0] |>
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
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_linetype_manual(values = c(2, 1))

# Compute hazard ratios (X = 1 vs X = 0) for remaining two panels
hr_dat <- dat_squeeze_long[, .(
  ratio = value[X == 1] / value[X == 0]
), by = c("time", "cause", "what")][what != "cuminc"]

# Panel B
p_squeeze_cshr <- hr_dat[what == "cs_haz"] |>
  ggplot(aes(time, ratio, col = cause, linetype = cause)) +
  geom_line(linewidth = 1.25) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 10)) +
  scale_color_manual(values = cols[c(1, 2)]) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(
    x = "Time",
    y = expression(h[k](t ~ "|" ~ X == 1) ~ "/" ~ h[k](t ~ "|" ~ X == 0)),
    col = "Cause",
    linetype = "Cause"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_linetype_manual(values = c(2, 1))

# Panel D
p_squeeze_subdisthr <- hr_dat[what == "subdist_haz"] |>
  ggplot(aes(time, ratio, col = cause, linetype = cause)) +
  geom_line(linewidth = 1.25) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 10)) +
  scale_color_manual(values = cols[c(1, 2)]) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(
    x = "Time",
    y = expression(lambda[k](t ~ "|" ~ X == 1) ~ "/" ~ lambda[k](t ~ "|" ~ X == 0)),
    col = "Cause",
    linetype = "Cause"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_linetype_manual(values = c(2, 1))

combined <- (p_squeeze_cs + p_squeeze_cshr) /
  (p_squeeze_subdist + p_squeeze_subdisthr) & theme(legend.position = "bottom")

# Figure 2
combined +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(
  here::here("squeeze_fig.pdf"), # adjust names with manuscript
  dpi = 300,
  units = "in",
  width = 8,
  height = 8,
  device = cairo_pdf
)


# Figure 3 ----------------------------------------------------------------


# DGM based on directly specifying the cumulative incidence function for
# both causes, such that two Fine-Gray models hold (until finite follow-up time)

# Set parameters (Fine-Gray models with Weibull baseline CDFs):
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

# Compute true values
dat_twofgs_x0 <- compute_true(
  t = t,
  model_type = "two_fgs",
  newdat = list("X" = 0),
  params = params
)

dat_twofgs_x1 <- compute_true(
  t = t,
  model_type = "two_fgs",
  newdat = list("X" = 1),
  params = params
)

# Panel A
p_twofgs_x0 <- dat_twofgs_x0[cause == 1] |>
  ggplot(aes(time, cuminc)) +
  geom_area(fill = cols[1], alpha = 0.75) +
  geom_ribbon(
    data = dat_twofgs_x0[cause == 2],
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
    aes(ymin = params$cause1$p + params$cause2$p, ymax = 1),
    alpha = 0.75,
    fill = cols[5]
  ) +
  annotate(
    "text",
    x = 10.05,
    y = params$cause1$p,
    parse = TRUE,
    label = "p[10]",
    family = "Roboto Condensed",
    hjust = 0
  ) +
  annotate(
    "text",
    x = 10.05,
    y = params$cause1$p + params$cause2$p,
    parse = TRUE,
    label = "p[10] + p[20]",
    family = "Roboto Condensed",
    hjust = 0
  ) +
  annotate(
    "text",
    x = 5,
    y = mean(c(1, params$cause1$p + params$cause2$p)),
    label = "'Cured'",
    family = "Roboto Condensed",
    hjust = 0
  ) +
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

# P(D = 1 | X)
p1_x <- 1 - (1 - params$cause1$p)^exp(params$cause1$betas)

# Panel B
p_twofgs_x1 <- dat_twofgs_x1[cause == 1] |>
  ggplot(aes(time, cuminc)) +
  geom_area(fill = cols[1], alpha = 0.75) +
  geom_ribbon(
    data = dat_twofgs_x1[cause == 2],
    aes(ymax = pmin(1, cuminc + p1_x), ymin = p1_x),
    fill = cols[2],
    alpha = 0.75
  ) +
  geom_ribbon(
    data = dat_twofgs_x1[cause == 2],
    aes(ymax = pmax(1, cuminc + p1_x), ymin = 1),
    fill = cols[2],
    alpha = 0.5
  ) +
  coord_cartesian(ylim = c(0, 1.1), xlim = c(0, 11)) +
  geom_segment(
    y = p1_x,
    yend = p1_x,
    x = 0,
    xend = 10,
    linetype = "dashed"
  ) +
  annotate("text",
    x = 10.05,
    y = p1_x,
    parse = TRUE,
    label = "p[1](x)",
    family = "Roboto Condensed",
    hjust = 0
  ) +
  geom_curve(
    aes(x = 4.5, y = 1.075, xend = 5, yend = 1.025),
    colour = "black",
    linewidth = 0.25,
    curvature = -0.3,
    arrow = arrow(length = unit(0.02, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  annotate("text",
    x = 4.5,
    y = 1.075,
    hjust = 1,
    label = "TFP > 1",
    family = "Roboto Condensed"
  ) +
  labs(x = "Time", y = "Cumulative incidence") +
  scale_x_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  annotate(
    "text",
    x = 5,
    y = mean(c(p1_x, 1)),
    label = "F[2](t ~ '|' ~ X == 1)",
    family = "Roboto Condensed",
    parse = TRUE
  ) +
  annotate(
    "text",
    x = 5,
    y = mean(c(p1_x, 0)),
    label = "F[1](t ~ '|' ~ X == 1)",
    family = "Roboto Condensed",
    parse = TRUE
  )

# Figure 3
p_twofgs_x0 + p_twofgs_x1 + plot_annotation(tag_levels = "A")

ggsave(
  here::here("two_fgs_fig.pdf"), # also to eps/tiff?
  dpi = 300,
  units = "in",
  width = 11,
  height = 7,
  device = cairo_pdf
)
