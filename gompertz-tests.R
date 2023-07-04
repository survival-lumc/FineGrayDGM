# Tests with gompertz
gompertz_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    rate * exp(t * params[["base_shape"]])
  } else { # type == "cumulative" for cumulative hazard
    (rate / params[["base_shape"]]) * (exp(params[["base_shape"]] * t) - 1)
  }
}

integrate_to_t <- Vectorize(FUN = function(t, fun, ...) {
  if (!is.numeric(t) | t < 0) stop("t should be positive.")
  ifelse(
    test = t == 0,
    yes = 0,
    no = integrate(f = fun, lower = 0L, upper = t, ...)$value
  )
}, vectorize.args = "t")

params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "base_rate" = 1,
    "base_shape" = -1
  ),
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "base_rate" = 0.3,
    "base_shape" = 0.25 #0 # Probs use weibull for cause 2?
  )
)
newdat <- data.frame("X" = 0)

# Prepare model matrices (for now newdat is only one row), later vectorize
predictor_formulas <- lapply(params, "[[", "formula")
modmats <- lapply(predictor_formulas, function(form) {
  x <- model.matrix(form, data = newdat)
  x[, !colnames(x) %in% "(Intercept)"]
})
x_cause1 <- modmats[["cause1"]]
x_cause2 <- modmats[["cause2"]]

t <- seq(0.001, 10, by = 0.1)
haz_cs2 <- gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
haz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard")

integral_fun_cs1 <- function(t) {
  lambda_1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
  Lambda_1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")
  H_2 <- gompertz_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
  lambda_1 * exp(-Lambda_1 + H_2)
}

fun_haz_cs1 <- function(t, x, params, type = "hazard") {

  # Closure
  integral_fun_cs1 <- function(t) {
    lambda_1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
    Lambda_1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")
    H_2 <- gompertz_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
    lambda_1 * exp(-Lambda_1 + H_2)
  }

 integral_denom <- integrate_to_t(fun = integral_fun_cs1, t = t)


}

num <- integral_fun_cs1(t)
plot(t, num)
integral_denom <- integrate_to_t(fun = integral_fun_cs1, t = t)
plot(t, integral_denom)
haz_cs1 <- num / (1 - integral_denom)
plot(t, -log((1 - integral_denom))) # cumhaz cs 1
plot(t, haz_cs1, ylim = c(-2, 2))

plot(t, haz_subdist1, type = "l")
lines(t, haz_cs2)
lines(t, haz_cs1)

library(ggplot2)
library(data.table)

# Show also cumincs? and hazard ratios?
data.table("time" = t, haz_cs1, haz_cs2, haz_subdist1) |>
  melt(id.vars = "time", variable.name = "type", value.name = "hazard") |>
  ggplot(aes(time, hazard)) +
  geom_line(
    aes(linetype = type, col = type),
    linewidth = 1.5
  ) +
  coord_cartesian(ylim = c(-10, 10))

times <- seq(0, 10, by = 0.1)
params <- list(
  "formula" = ~ X,
  "betas" = 0,
  "base_rate" = 0.1,
  "base_shape" = -1
)

plot(
  times,
  gompertz_hazard(
    times,
    x = model.matrix(params$formula, data = data.frame("X" = 0))[, -1],
    params = params,
    type = "cumulative"
  )
)

