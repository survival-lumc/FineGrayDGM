# Checkiin the intergral ahain

haz_cs1 <- num / (1 - integrate_to_t(fun = integral_fun_cs1, t = t))

plot(t, haz_cs1)

# Cumulative hazard
plot(t, integrate_to_t(t, function(t) {
  integral_fun_cs1(t) /
    (1 - integrate_to_t(fun = integral_fun_cs1, t = t))
}))

cumhaz_cs1 <- integrate_to_t(t, function(t) {
  integral_fun_cs1(t) / (1 - integrate_to_t(fun = integral_fun_cs1, t = t))
})

cumhaz_cs1
-log(-integrate_to_t(fun = integral_fun_cs1, t = t))
cbind(
  cumhaz_cs1,
-log(1 - integrate_to_t(fun = integral_fun_cs1, t = t))
)


plot(t, 1 - F1 - F2)
lines(t, 1 - integrate_to_t(fun = integral_fun_cs1, t = t))
lines(t, 1 - (haz_cs2))


# Try the other way around ------------------------------------------------


library(numDeriv)
source("helpers-new.R")

params <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "base_rate" = 0.1, # 0.5 is same spot
    "base_shape" = -2
  ),
  # This is weib
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(0.5),
    "base_rate" = 0.1,
    "base_shape" = 2 #-1
  )
)
newdat <- list("X" = 0)
t <- seq(0.001, 10, length.out = 250)


predictor_formulas <- lapply(params, "[[", "formula")
modmats <- lapply(predictor_formulas, function(form) {
  x <- model.matrix(form, data = newdat)
  x[, !colnames(x) %in% "(Intercept)"]
})
x_cause1 <- modmats[["cause1"]]
x_cause2 <- modmats[["cause2"]]

haz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
plot(t, haz_subdist1, ylim = c(0, 5), type = "l", col = "blue")
haz_cs1 <- gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
lines(t, haz_cs1)
lines(t, {
  gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard") -
    gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard") -
    grad(function(t) {
      log(
        gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard") /
          gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
      )
    }, x = t)
})

fun_cs2 <- function(t) {
  gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard") -
    gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard") -
    grad(function(t) {
      log(
        gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard") /
          gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
      )
    }, x = t)
}

fun_cs2(t)

plot(t, fun_cs2(t))

plot(t, integrate_to_t(t, fun_cs2))

prod2 <- function(t, cause) {
  #browser()
  haz <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard") -
    gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard") -
    grad(function(t) {
      log(
        gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard") /
          gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
      )
    }, x = t)
  cumhaz_cause1 <- gompertz_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
  cumhaz_cause2 <- integrate_to_t(t, fun_cs2)
  haz * exp(-cumhaz_cause1 - cumhaz_cause2)
}

# Calculate both cumulative incidences
# Could use just gompertz property to get cumulative incidence?
#F1 <- integrate_to_t(fun = prod, t = t, cause = 1)
F1 <- 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative"))
F2 <- try(integrate_to_t(fun = prod2, t = t, cause = 2))

plot(t, F1 + F2)
plot(t, F1)
plot(t, F2)





plot(t, haz_cs1, type = 'l', col = "blue", ylim = c(0, 3))
lines(t, haz_cs2, col = "red")
lines(t, haz_cs1 + haz_cs2, lwd = 2)
lines(t, haz_subdist1, lty = 2)


plot(t, haz_cs1 + haz_cs2, lwd = 2, ylim = c(0, 3), type = "l")
lines(t, haz_cs1 + haz_cs2 - haz_subdist1)

# Cumhazes
cumhaz_cs1 <- -log((1 - integrate_to_t(fun = integral_fun_cs1, t = t)))
cumhaz_cs2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
cumhaz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")

plot(t, cumhaz_cs1, type = 'l', col = "blue", ylim = c(0, 5),
     ylab = "Cumhaz")
lines(t, cumhaz_cs2, col = "red")
lines(t, cumhaz_cs1 + cumhaz_cs2, lwd = 2)
lines(t, cumhaz_subdist1, lty = 2)
lines(t, cumhaz_subdist1 - cumhaz_cs2, lty = 2)

plot(
  t,
  integrate_to_t(fun = integral_fun_cs1, t = t)

)

# exp(-H1(t)) =
1- integrate_to_t(fun = integral_fun_cs1, t = t)#

#.. so all cause is
cumhaz_cs2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
EFS <- (1 - integrate_to_t(fun = integral_fun_cs1, t = t)) * exp(-cumhaz_cs2)

plot(
  t,
  1 - EFS, # TFP!
  ylim = c(0.9, 1.1),
  type = "l"
)
plot(t, 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")), ylim = c(0, 1))
# P(D = 1 | X)
# 1 - exp(params[["cause1"]]$base_rate / params[["cause1"]]$base_shape)

F1 <- 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative"))
plot(
  t,
  F1,
  ylim = c(0, 1),
  type = "l"
)
d1_inf <- 1 - (1 - (1 - exp(params[["cause1"]]$base_rate / params[["cause1"]]$base_shape)))^exp(params[["cause1"]]$betas)
abline(
  b = 0,
  lty = 2,
  a = d1_inf
)
lines(
  t,
  d1_inf + (1 - EFS - F1)
    #((1 - EFS) - (1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative"))))
)



plot(t, exp(-cumhaz_cs1 - cumhaz_cs2), type = 'l', col = "blue", ylim = c(0, 1),
     ylab = "Cumhaz")
lines(t, exp(-cumhaz_subdist1), col = "red")
lines(t, cumhaz_cs1 + cumhaz_cs2, lwd = 2)
lines(t, cumhaz_subdist1, lty = 2)
lines(t, cumhaz_subdist1 - cumhaz_cs2, lty = 2)
