library(data.table)
library(ggplot2)
theme_set(theme_minimal(base_size = 15))
source("helpers-new.R")

newdat <- data.frame("X" = 0) # later with X = 1
t <- seq(0.001, 10, by = 0.1)

# Squeezing ---------------------------------------------------------------


params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "p" = 0.2,
    "base_rate" = 1,
    "base_shape" = 1.25
  ),
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "base_rate" = 1,
    "base_shape" = 1.5
  )
)

dat_squeeze <- compute_true(
  t = t,
  model_type = "squeezing",
  newdat = newdat,
  params = params
)

dat_squeeze_x1 <- compute_true(
  t = t,
  model_type = "squeezing",
  newdat = list("X" = 1),
  params = params
)
plot(
  t,
  dat_squeeze_x1[cause == 2][["subdist_haz"]] / dat_squeeze[cause == 2][["subdist_haz"]]
)

melt(
  data = dat_squeeze,
  id.vars = c("time", "cause"),
  variable.name = "what",
  value.name = "value"
) |>
  ggplot(aes(time, value, col = cause)) +
  geom_line(linewidth = 1.5, aes(linetype = cause)) +
  facet_wrap(~ what, scales = "free")


# Reduction ---------------------------------------------------------------


# With gompertzes
params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "base_rate" = 0.5,
    "base_shape" = -2
  ),
  # This is weib
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "base_rate" = 2,
    "base_shape" = 0.5 # Could also use Weibs for cause 2
  )
)

dat_reduction <- compute_true(
  t = t,
  model_type = "reduction_factor",
  newdat = newdat,
  params = params
)
dat_reduction[, .(max(cuminc)), by = cause][[2]] |> sum()

melt(
  data = dat_reduction,
  id.vars = c("time", "cause"),
  variable.name = "what",
  value.name = "value"
) |>
  ggplot(aes(time, value, col = cause)) +
  geom_line(linewidth = 1.5, aes(linetype = cause)) +
  facet_wrap(~ what, scales = "free")


# Two FGs -----------------------------------------------------------------


params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "p" = 0.05,
    "base_rate" = 1,
    "base_shape" = 0.75
  ),
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "p" = 0.1,
    "base_rate" = 1,
    "base_shape" = 0.75
  )
)

dat_twofgs <- compute_true(
  t = t,
  model_type = "two_fgs",
  newdat = newdat,
  params = params
)
dat_twofgs[, .(max(cuminc)), by = cause]

melt(
  data = dat_twofgs,
  id.vars = c("time", "cause"),
  variable.name = "what",
  value.name = "value"
) |>
  ggplot(aes(time, value, col = cause)) +
  geom_line(linewidth = 1.5, aes(linetype = cause)) +
  facet_wrap(~ what, scales = "free")
