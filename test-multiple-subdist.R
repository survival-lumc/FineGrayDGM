# Mozumbder tings
library(data.table)
library(ggplot2)
library(Manu)
library(patchwork)
library(extrafont)
theme_set(theme_minimal(base_size = 15, base_family = "Roboto Condensed"))
source("helpers-new.R")

newdat <- data.frame("X" = 1)

# Two gompz
params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "base_rate" = 0.5,
    "base_shape" = -2
  ),
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(0.25),
    "base_rate" = 0.75,
    "base_shape" = -1.5
  )
)

predictor_formulas <- lapply(params, "[[", "formula")
modmats <- lapply(predictor_formulas, function(form) {
  x <- model.matrix(form, data = newdat)
  x[, !colnames(x) %in% "(Intercept)"]
})
x_cause1 <- modmats[["cause1"]]
x_cause2 <- modmats[["cause2"]]

ts <- seq(0.01, 10, length.out = 100)


plot(
  ts,
  gompertz_hazard(t = ts, x = x_cause1, params = params[["cause2"]], type = "hazard")
)


F1 <- 1 - exp(-gompertz_hazard(ts, x_cause1, params[["cause1"]], type = "cumulative"))
F2 <- 1 - exp(-gompertz_hazard(ts, x_cause1, params[["cause2"]], type = "cumulative"))

plot(ts, F1, ylim = c(0, 1))
lines(ts, F2)

haz_cs1 <- gompertz_hazard(ts, x_cause1, params[["cause1"]], type = "haz") * (
  1 + ((F1 + F2 - F1) / (1 - F1 - F2))
)
haz_cs2 <- gompertz_hazard(ts, x_cause2, params[["cause2"]], type = "haz") * (
  1 + ((F2 + F1 - F2) / (1 - F1 - F2))
)

plot(ts, haz_cs2)

# Continue here...
haz_cs <- function(t, x, params, cause = 1) {

  # Mod mats
  predictor_formulas <- lapply(params, "[[", "formula")
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = newdat)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats[["cause1"]]
  x_cause2 <- modmats[["cause2"]]

  # Cumincs
  F1 <- 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative"))
  F2 <- 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause2"]], type = "cumulative"))

  # Subs
  haz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "haz")
  haz_subdist2 <- gompertz_hazard(t, x_cause2, params[["cause2"]], type = "haz")

  haz_cs <- switch(cause,
    "1" = haz_subdist1 * (1 + F2 / (1 - F1 - F2)),
    "2" = haz_subdist2 * (1 + F1 / (1 - F1 - F2))
  )

  return(haz_cs)
}

library(numDeriv)

ts <- seq(0.01, 10, length.out = 100)
plot(ts, haz_cs(ts, x = newdat, params = params, cause = 1))

# Other way
F1 <- 1 - exp(-gompertz_hazard(ts, x_cause1, params[["cause1"]], type = "cumulative"))
F2 <- 1 - exp(-gompertz_hazard(ts, x_cause1, params[["cause2"]], type = "cumulative"))
f1 <- grad(
  function(t) 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")),
  x = ts
)
f2 <- grad(
  function(t) 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")),
  x = ts
)
haz_cs1 <- f1 / (1 - F1 - F2)
haz_cs2 <- f2 / (1 - F1 - F2)

plot(ts, 1 - F1 - F2)
plot(ts, haz_cs1)
plot(ts, haz_cs2)
plot(ts, haz_cs1 + haz_cs2)

cumhaz_allcause <- function(t, x, params) {

  predictor_formulas <- lapply(params, "[[", "formula")
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = newdat)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats[["cause1"]]
  x_cause2 <- modmats[["cause2"]]

  F1 <- 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative"))
  F2 <- 1 - exp(-gompertz_hazard(t, x_cause2, params[["cause2"]], type = "cumulative"))
  f1 <- grad(
    function(t) 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")),
    x = t
  )
  f2 <- grad(
    function(t) 1 - exp(-gompertz_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")),
    x = t
  )
  haz_cs1 <- f1 / (1 - F1 - F2)
  haz_cs2 <- f2 / (1 - F1 - F2)

  haz_cs1 + haz_cs2
}

plot(ts, cumhaz_allcause(ts, params = params, x = newdat))
plot(
  ts,
  exp(-integrate_to_t(fun = cumhaz_allcause, t = ts, params = params, x = newdat))
)



plot(ts, integrate_to_t(fun = haz_cs, t = ts, cause = 1, x = newdat, params = params))

# Should be exp[-cumhaz] that needs to be inverted!!
uniroot(
  f = function(t) {
    integrate_to_t(fun = haz_cs, t = t, cause = 1, x = newdat, params = params) - runif(1)
  },
  lower = 0,
  upper = 1000
)$root

library(riskRegression)
library(prodlim)
n <- 2500
X <- rbinom(n, size = 1, prob = 0.5)
X

T1 <- mapply(
  FUN = function(x) {
    uniroot(
      f = function(t) {
        exp(-integrate_to_t(fun = haz_cs, t = t, cause = 1, x = x, params = params)) -
          runif(1)
      },
      lower = 0,
      upper = 1000
    )$root
  },
  x = X
)

T2 <- mapply(
  FUN = function(x) {
    uniroot(
      f = function(t) {
        exp(-integrate_to_t(fun = haz_cs, t = t, cause = 2, x = x, params = params)) -
          runif(1)
      },
      lower = 0,
      upper = 1000
    )$root
  },
  x = X
)

time <- pmin(T1, T2)
D <- 1 + as.numeric(T2 < T1)

dat <- cbind.data.frame(X, D, time)
params
FGR(Hist(time, D)~ X, cause = 1, data = dat)$crr$coef
FGR(Hist(time, D)~ X, cause = 2, data = dat)$crr$coef



# Tests Mozumbder ---------------------------------------------------------

source("helpers-new.R")
newdat <- data.frame("X" = 1)
ts <- seq(0.01, 10, length.out = 100)

# For mixture (see section 3 of lambert et al. 2017)
# Using params however from figure 1 Mozumder et al. (2018)
params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(-0.5),
    "pi" = 0.5,
    "lambda_1" = 0.6,
    "gamma_1" = 0.5,
    "lambda_2" = 0.01,
    "gamma_2" = 0.35
  ),
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(0.2),
    "pi" = 0.5,
    "lambda_1" = 0.01,
    "gamma_1" = 0.8,
    "lambda_2" = 0.7,
    "gamma_2" = 1.45
  )
)

# From lambert paper
# params <- list(
#   # Just cause 1 in scenario 1
#   "cause1" = list(
#     "formula" = ~ X,
#     #"betas" = c(0.5),
#     "pi" = 0.4,
#     "lambda_1" = 0.02,
#     "gamma_1" = 3,
#     "lambda_2" = 0.02,
#     "gamma_2" = 1.6
#   ),
#   # Just cause 1 in scenario 2
#   "cause2" = list(
#     "formula" = ~ X,
#     #"betas" = c(0.25),
#     "pi" = 0.2,
#     "lambda_1" = 0.2,
#     "gamma_1" = 2,
#     "lambda_2" = 0.02,
#     "gamma_2" = 1.6
#   )
# )

weibull_mixture_hazard <- function(t, x, params, type = "hazard") {

  # Read-in
  lp <- drop(x %*% params[["betas"]])
  pi <- params$pi
  lambda_1 <- params$lambda_1
  gamma_1 <- params$gamma_1
  lambda_2 <- params$lambda_2
  gamma_2 <- params$gamma_2

  # Survival
  surv_fun <- pi * exp(-lambda_1 * t^gamma_1) + (1 - pi) * exp(-lambda_2 * t^gamma_2)
  num_haz <- pi * lambda_1 * gamma_1 * t^(gamma_1 - 1) * exp(-lambda_1 * t^gamma_1) +
                (1 - pi) * lambda_2 * gamma_2 * t^(gamma_2 - 1) * exp(-lambda_2 * t^gamma_2)

  # better use simsurv vignette!!!
  res <- switch(
    type,
    "cumulative" = -log(surv_fun) * exp(lp),
    "hazard" = (num_haz / surv_fun) * exp(lp)
  )

  return(res)
}

# Get mod mats
predictor_formulas <- lapply(params, "[[", "formula")
modmats <- lapply(predictor_formulas, function(form) {
  x <- model.matrix(form, data = newdat)
  x[, !colnames(x) %in% "(Intercept)"]
})
x_cause1 <- modmats[["cause1"]]
x_cause2 <- modmats[["cause2"]]

plot(
  ts,
  weibull_mixture_hazard(t = ts, x = x_cause1, params = params[["cause2"]], type = "hazard"),
  xlim = c(0, 5),
  type = "l",
  ylim = c(0, 1.5)
)
plot(
  ts,
  weibull_mixture_hazard(t = ts, x = x_cause2, params = params[["cause2"]], type = "hazard"),
  xlim = c(0, 5),
  type = "l",
  ylim = c(0, 1.5)
)




# Try comparing the calcs
F1 <- 1 - exp(-weibull_mixture_hazard(t = ts, x = x_cause1, params = params[["cause1"]], type = "cumulative"))
F2 <- 1 - exp(-weibull_mixture_hazard(t = ts, x = x_cause2, params = params[["cause2"]], type = "cumulative"))
plot(ts, F1 + F2)
abline(a = 1, b = 0)

f1 <- grad(
  function(t) 1 - exp(-weibull_mixture_hazard(t = t, x = x_cause1, params = params[["cause1"]], type = "cumulative")),
  x = ts
)
f2 <- grad(
  function(t) 1 - exp(-weibull_mixture_hazard(t = t, x = x_cause2, params = params[["cause2"]], type = "cumulative")),
  x = ts
)
haz_cs1 <- f1 / (1 - F1 - F2)
haz_cs1_bis <- weibull_mixture_hazard(t = ts, x = x_cause1, params = params[["cause1"]], type = "hazard") * (
  1 + F2 / (1 - F1 - F2)
)
haz_cs2 <- f2 / (1 - F1 - F2)
haz_cs2_bis <- weibull_mixture_hazard(t = ts, x = x_cause2, params = params[["cause2"]], type = "hazard") * (
  1 + F1 / (1 - F1 - F2)
)

# Oof this works!!
plot(ts, haz_cs1)
lines(ts, haz_cs1_bis)

plot(ts, haz_cs2)
lines(ts, haz_cs2_bis)

# Now let's try and simulate!!
haz_all_cause <- function(t, x, params) {

  # Get mod mats
  # predictor_formulas <- lapply(params, "[[", "formula")
  # modmats <- lapply(predictor_formulas, function(form) {
  #   x <- model.matrix(form, data = x)
  #   x[, !colnames(x) %in% "(Intercept)"]
  # })
  # x_cause1 <- modmats[["cause1"]]
  # x_cause2 <- modmats[["cause2"]]

  # Now the ingredients
  F1 <- 1 - exp(-weibull_mixture_hazard(t = t, x = x, params = params[["cause1"]], type = "cumulative"))
  F2 <- 1 - exp(-weibull_mixture_hazard(t = t, x = x, params = params[["cause2"]], type = "cumulative"))

  haz_cs1_bis <- weibull_mixture_hazard(t = t, x = x, params = params[["cause1"]], type = "hazard") * (
    1 + F2 / (1 - F1 - F2)
  )
  haz_cs2_bis <- weibull_mixture_hazard(t = t, x = x, params = params[["cause2"]], type = "hazard") * (
    1 + F1 / (1 - F1 - F2)
  )

  haz_cs1_bis + haz_cs2_bis
}

haz_all_cause(ts, 1, params)
plot(ts, integrate_to_t(ts, fun = haz_all_cause, params = params, x = 0))

# They restrict to 5 years
uniroot(
  f = function(t) {
    exp(-integrate_to_t(t, fun = haz_all_cause, params = params, x = 1)) - runif(1)
  },
  lower = 0,
  upper = 5
)$root


set.seed(4965)
n <- 2500
X <- rbinom(n, size = 1, prob = 0.5)
draws <- mapply(
  function(u, x) {
    time_test <- try({
      uniroot(
        f = function(t) {
          exp(-integrate_to_t(t, fun = haz_all_cause, params = params, x = x)) - u
        },
        lower = 0,
        upper = 5
      )$root
    })
    time <- if (inherits(time_test, "try-error")) 5 else time_test

    # Now see what event
    all_cause <- haz_all_cause(time, x, params)

    # Get the rest
    F1 <- 1 - exp(-weibull_mixture_hazard(t = time, x = x, params = params[["cause1"]], type = "cumulative"))
    F2 <- 1 - exp(-weibull_mixture_hazard(t = time, x = x, params = params[["cause2"]], type = "cumulative"))

    haz_cs1_bis <- weibull_mixture_hazard(t = time, x = x, params = params[["cause1"]], type = "hazard") * (
      1 + F2 / (1 - F1 - F2)
    )

    D_tilde <- 1 + rbinom(n = 1, size = 1, prob = 1 - haz_cs1_bis / all_cause)
    D <- ifelse(time == 5, 0, D_tilde)
    cbind(time, D)
  },
  u = runif(n),
  x = X,
  SIMPLIFY = FALSE
)

dat <- cbind.data.frame(
  do.call(rbind.data.frame, draws),
  X
)
dat <- dat[!is.na(dat$D), ]
dat
table(dat$D)

params
FGR(Hist(time, D)~ X, cause = 1, data = dat)$crr$coef
FGR(Hist(time, D)~ X, cause = 2, data = dat)$crr$coef

#install.packages("mets")
library(mets)
test <- doubleFGR(Event(time, D) ~ X, data = dat)
test <- doubleFGR(Event(time, D) ~ X, data = dat, restrict = 1)
print(test)
bplotdFG(test,cause=1)


cumhaz_mao <- function(t, rho) rho * (1 - exp(-t))
plot(ts, 1 - exp(-cumhaz_mao(ts, 0.1)), ylim = c(0, 1))
plot(ts, 1 - exp(-cumhaz_mao(ts, 1.5)), ylim = c(0, 1))

plot(ts, 1 - exp(-cumhaz_mao(ts, 0.1)), ylim = c(0, 1))
(ts, 1 - exp(-cumhaz_mao(ts, 1.5)), ylim = c(0, 1))



plot(ts, 1 - exp(1.5 * exp(-ts)))
plot(ts, 1 - exp(1.5 * exp(-ts)))

# You won't find a root because since survival does not go to zero!
uniroot(
  f = function(t) {
    exp(-integrate_to_t(t, fun = haz_all_cause, params = params, x = 1)) #- runif(1)
  },
  lower = 0,
  upper = 1000,
  extendInt = TRUE
)$root


