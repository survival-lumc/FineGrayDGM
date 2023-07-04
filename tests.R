#https://cran.rstudio.com/web/packages/flexsurv/flexsurv.pdf

# Use: https://github.com/survival-lumc/FineGrayCovarMI/blob/bc146cba03f2a14edf59fa12a8da759c157b6c42/analysis/testing-all-DGDs.R

haz_gompertz_KM <- function(theta, alpha, t) theta * exp(alpha * t)
t <- seq(0, 10, length.out = 100)
plot(t, haz_gompertz_KM(theta = 0.01, alpha = -1, t), ylim = c(0, 0.05))
# alpha is steepness, theta is basehaz at time = 0
points(t, haz_gompertz_KM(theta = 0.01, alpha = 0.1, t))

# More concrete
alpha_1 <- -1
alpha_2 <- -0.95 #0.05
theta_1 <- theta_2 <- 0.75


plot(
  t,
  {
    haz_gompertz(theta = theta_1, alpha = alpha_1, t) -
      haz_gompertz(theta = theta_2, alpha = alpha_2, t) -
      (alpha_1 - alpha_2)
  },
  col = "black",
  lty = 3,
  lwd = 3,
  type = "l",
  ylim = c(-3, 1.25)
)

lines(
  t,
  haz_gompertz_KM(theta = theta_1, alpha = alpha_1, t),
  ylim = c(0, theta_1 * 5),
  col = "blue",
  lty = 1,
  type = "l",
  lwd = 3
)
lines(
  t,
  haz_gompertz_KM(theta = theta_2, alpha = alpha_2, t),
  col = "red",
  lty = 2,
  lwd = 3
)


testo <- function(theta_1, theta_2, alpha_1, alpha_2, t) {
  haz_gompertz_KM(theta = theta_1, alpha = alpha_1, t) -
    haz_gompertz_KM(theta = theta_2, alpha = alpha_2, t) -
    (alpha_1 - alpha_2) - 0
}

uniroot(
  f = testo,
  interval = c(0, 1000),
  alpha_1 = alpha_1,
  alpha_2 = alpha_2,
  theta_1 = theta_1,
  theta_2 = theta_2
)$root

cumhaz_gompertz <- function(theta, alpha, t) (theta / alpha) * (exp(alpha * t) - 1)

plot(t, cumhaz_gompertz(theta = 0.75, alpha = -1, t = t))
plot(
  t,
  1 - exp(-cumhaz_gompertz(theta = 0.75, alpha = -1, t = t))
)
1 - exp(0.75/-1) # cuminc at inf



# Clean -------------------------------------------------------------------


library(survival)
library(data.table)
library(ggplot2)

haz_gompertz <- function(t, alpha, theta) theta * exp(alpha * t)
cumhaz_gompertz <- function(t, alpha, theta) (theta / alpha) * (exp(alpha * t) - 1)
# This next one to make more general
#haz_cs2_derived <- function()
#haz_cs2_gompertz <-
times <- seq(0, 10, length.out = 1000)
hr_subdist <- exp(0.5)
alpha_subdist1 <- -1
hr_cs1 <- exp(0.5)
alpha_cs1 <- -2
theta_cs1 <- theta_subdist1 <- 0.75


plot(
  times,
  haz_gompertz(t = times, alpha = alpha_subdist1, theta = theta_subdist1),
  ylim = c(-10, 10)
)
points(
  times,
  haz_gompertz(t = times, alpha = alpha_subdist1, theta = theta_subdist1 * hr_subdist)
)
points(
  times,
  haz_gompertz(t = times, alpha = alpha_cs1, theta = theta_cs1)
)
points(
  times,
  haz_gompertz(t = times, alpha = alpha_cs1, theta = theta_cs1 * hr_cs1)
)
points(
  times,
  {
    haz_gompertz(t = times, alpha = alpha_subdist1, theta = theta_subdist1) -
      haz_gompertz(t = times, alpha = alpha_cs1, theta = theta_cs1) -
      (alpha_subdist1 - alpha_cs1)
  }
)
points(
  times,
  {
    haz_gompertz(t = times, alpha = alpha_subdist1, theta = theta_subdist1 * hr_subdist) -
      haz_gompertz(t = times, alpha = alpha_cs1, theta = theta_cs1 * hr_cs1) -
      (alpha_subdist1 - alpha_cs1)
  }
)

# Try simulating
cumhaz_cs2 <- function(t,
                       alpha_subdist1,
                       alpha_cs1,
                       theta_subdist1,
                       theta_cs1) {
  cumhaz_subdist1 <- cumhaz_gompertz(t, alpha_subdist1, theta_subdist1)
  cumhaz_cs1 <- cumhaz_gompertz(t, alpha_cs1, theta_cs1)
  cumhaz_subdist1 - cumhaz_cs1 - t * (alpha_subdist1 - alpha_cs1)
}

plot(times, cumhaz_cs2(times, -1, -2, 0.75, 0.75))
n <- 2000
X <- rbinom(n, size = 1, prob = 0.5)
mapply(function(x) {
  uniroot(
    f = function(x, u, ...) cumhaz_cs2 - u
  )
}, x = X)

#rgompertz <- function(n, alpha, theta) {
#  u <- runif(n)
#  log(-alpha * log(u) / theta + 1) / alpha
#}

n <- 10000
X <- rbinom(n, size = 1, prob = 0.5)
timez <- flexsurv::rgompertz(n, shape = -1, rate = 1 * exp(X)) #rgompertz(n, alpha = 5, theta = 0.75)
hist(timez)
dat <- cbind.data.frame(timez, D = 1, X)
coxph(Surv(timez, D) ~ X, data = dat[dat$timez != Inf, ])




# Two FGs -----------------------------------------------------------------



# Squeezing ---------------------------------------------------------------


source("helpers.R")
times <- seq(1e-6, 10, length.out = 100)

params_squeeze <- list(
  "cause1" = list(
    "formula" = ~ X,
    "betas" = c(1),
    "p" = 0.3,
    "base_rate" = 1,
    "base_shape" = 0.75
  ),
  "cause2" = list(
    "formula" = ~ X,
    "betas" = c(1),
    "base_rate" = 1,
    "base_shape" = 0.75
  )
)

df_base <- compute_true(
  t = times,
  newdat = cbind.data.frame(X = 0),
  params = params_squeeze,
  model_type = "correct_FG"
)

df_x1 <- compute_true(
  t = times,
  newdat = cbind.data.frame(X = 1),
  params = params_squeeze,
  model_type = "correct_FG"
)

# Better compute true funcs?
# Should allow multiple patient, and give times directly in new dat

# Possibly try numerical derivatives in R:
# https://stackoverflow.com/questions/18494302/numerical-derivatives-of-an-arbitrarily-defined-function
