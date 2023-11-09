# General helpers ---------------------------------------------------------


# Weibull (cumulative) hazard, shape and rate parametrisation
# (see Klein and Moeschberger book)
weibull_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    rate * params[["base_shape"]] * t^(params[["base_shape"]] - 1)
  } else { # type == "cumulative" for cumulative hazard
    rate * t^params[["base_shape"]]
  }
}

# Gompertz (cumulative) hazard, shape and rate parametrisation
# (e.g. see documentation of {flexsurv} package)
gompertz_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    rate * exp(t * params[["base_shape"]])
  } else { # type == "cumulative" for cumulative hazard
    (rate / params[["base_shape"]]) * (exp(params[["base_shape"]] * t) - 1)
  }
}

# Generic function to perform integration from 0 to t e.g. for cumulative hazards,
# or for the net risk of cause 1 in the reduction factor formula (2)
integrate_to_t <- Vectorize(FUN = function(t, fun, ...) {
  if (!is.numeric(t) | t < 0) stop("t should be positive.")
  ifelse(
    test = t == 0,
    yes = 0,
    no = integrate(f = fun, lower = 0L, upper = t, ...)$value
  )
}, vectorize.args = "t")


# Workhorse function ------------------------------------------------------


# Computes (for both causes) at given timepoints the true:
# - cumulative incidence functions
# - subdistribution hazards
# - cause-specific hazards
compute_true <- function(t,
                         model_type = c("squeezing", "reduction_factor", "two_fgs"),
                         newdat,
                         params) {

  # Get model matrices for all causes
  predictor_formulas <- lapply(params, "[[", "formula")
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = newdat)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats[["cause1"]]
  x_cause2 <- modmats[["cause2"]]

  # Pick a DGM:

  # - Specifying subdistribution cause 1, cause-specific hazard cause 2
  if (model_type == "reduction_factor") {

    # This is what we specify
    haz_cs2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
    haz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard")

    # Use reduction factor to get cause-specific hazard for cause 1
    integral_fun_cs1 <- function(t) {
      haz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
      cumhaz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")
      cumhaz_cs2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
      haz_subdist1 * exp(-cumhaz_subdist1 + cumhaz_cs2)
    }
    num <- integral_fun_cs1(t)
    haz_cs1 <- num / (1 - integrate_to_t(fun = integral_fun_cs1, t = t)) # Formula (2)

    # To get point at which neg hazards occur
    # uniroot(f = function(t) {integral_fun_cs1(t) / (1 - integrate_to_t(fun = integral_fun_cs1, t = t))}, interval = c(0, 10))

    # Cumulative incidence cause 1 simple since it is a Fine-Gray model
    F1 <- 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative"))

    # For cause 2, simple to calculate event-free survival (EFS) first,
    # and then subtract it + F1 from 1
    cumhaz_cause2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
    EFS <- (1 - integrate_to_t(fun = integral_fun_cs1, t = t)) * exp(-cumhaz_cause2)
    F2 <- 1 - EFS - F1

    # For the subdistribution hazard for cause 2: compute subdensity first
    subdens_2 <- haz_cs2 * EFS
    haz_subdist2 <- subdens_2 / (1 - F2)

    # (Not used as example in article, but can use all-cause hazard too)
  } else if (model_type == "all_cause") {

    # Specify subdistribution cause 1, and proportional all-cause model
    haz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
    EFS <- exp(-weibull_hazard(t, x_cause2, params[["cause2"]], type = "cumulative"))

    # Cumulative incidence cause 1 still the same
    F1 <- 1 - exp(-gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative"))

    # haz_subdist1 * (1 - F1) is equal to the subdensity! (re-arranging expression
    # for the subdistribution hazard) - therefore can compute cause-specific hazard 1
    haz_cs1 <- haz_subdist1 * (1 - F1) / EFS

    # Subtract to obtain cause-specific cause 1
    haz_cs2 <- weibull_hazard(t, x_cause2, params[["cause2"]], type = "hazard") - haz_cs1
    F2 <- 1 - EFS - F1

    # For cause 2, again use the subdensity
    subdens_2 <- haz_cs2 * EFS
    haz_subdist2 <- subdens_2 / (1 - F2)

    # - Directly specifying two Fine-Gray models
  } else if (model_type == "two_fgs") {

    # Cause 1
    hr_subdist1 <- drop(exp(x_cause1 %*% params[["cause1"]][["betas"]]))
    p1 <- params[["cause1"]][["p"]]
    shape_1 <- params[["cause1"]][["base_shape"]]
    rate_1 <- params[["cause1"]][["base_rate"]]
    F1 <- 1 - (1 - p1 * (1 - exp(-rate_1 * t^shape_1)))^hr_subdist1

    # Cause 2
    hr_subdist2 <- drop(exp(x_cause2 %*% params[["cause2"]][["betas"]]))
    p2 <- params[["cause2"]][["p"]]
    shape_2 <- params[["cause2"]][["base_shape"]]
    rate_2 <- params[["cause2"]][["base_rate"]]
    F2 <- 1 - (1 - p2 * (1 - exp(-rate_2 * t^shape_2)))^hr_subdist2

    # Both sub densities obtained simply by taking derivatives of F1 and F2, w.r.t time
    subdens_1 <- hr_subdist1 * (1 - p1 + p1 * exp(-rate_1 * t^shape_1))^(hr_subdist1 - 1) *
      p1 * exp(-rate_1 * t^shape_1) * rate_1 * shape_1 * t^(shape_1 - 1)

    subdens_2 <- hr_subdist2 * (1 - p2 + p2 * exp(-rate_2 * t^shape_2))^(hr_subdist2 - 1) *
      p2 * exp(-rate_2 * t^shape_2) * rate_2 * shape_2 * t^(shape_2 - 1)

    # Compute all hazards
    haz_subdist1 <- subdens_1 / (1 - F1)
    haz_subdist2 <- subdens_2 / (1 - F2)
    haz_cs1 <- subdens_1 / (1 - F1 - F2)
    haz_cs2 <- subdens_2 / (1 - F1 - F2)

    # - "Squeezing" DGM, same flavour as two_fgs above
  } else if (model_type == "squeezing") {

    # Calculate cumulative incidences directly - here for cause 1
    hr_subdist <- drop(exp(x_cause1 %*% params[["cause1"]][["betas"]]))
    p <- params[["cause1"]][["p"]]
    shape_1 <- params[["cause1"]][["base_shape"]]
    rate_1 <- params[["cause1"]][["base_rate"]]
    F1 <- 1 - (1 - p * (1 - exp(-rate_1 * t^shape_1)))^hr_subdist

    # Cause 2
    p2_inf <- (1 - p)^hr_subdist
    hr_condit <- drop(exp(x_cause2 %*% params[["cause2"]][["betas"]]))
    shape_2 <- params[["cause2"]][["base_shape"]]
    rate_2 <- params[["cause2"]][["base_rate"]]
    cumhaz_condit <- rate_2 * hr_condit * t^shape_2
    F2 <- (1 - exp(-cumhaz_condit)) * p2_inf

    # Compute both subdensities
    subdens_1 <- hr_subdist * (1 - p + p * exp(-rate_1 * t^shape_1))^(hr_subdist - 1) *
      p * exp(-rate_1 * t^shape_1) * rate_1 * shape_1 * t^(shape_1 - 1)

    subdens_2 <- hr_condit * rate_2 * shape_2 * t^(shape_2 - 1) *
      exp(-hr_condit * rate_2 * t^shape_2) *
      p2_inf # important

    # Note this can be checked with e.g.:
    # library(numDeriv)
    # F1_fun <- function(t) {
    #  1 - (1 - p * (1 - exp(-rate_1 * t^shape_1)))^hr_subdist
    # }
    # plot(t, grad(func = F1_fun, t))

    # Compute all hazards
    haz_subdist1 <- subdens_1 / (1 - F1)
    haz_subdist2 <- subdens_2 / (1 - F2)
    haz_cs1 <- subdens_1 / (1 - F1 - F2)
    haz_cs2 <- subdens_2 / (1 - F1 - F2)
  }

  # Return everything in long format
  dat_ev1 <- data.table(
    "time" = t,
    "cause" = "1",
    "cuminc" = F1,
    "subdist_haz" = haz_subdist1,
    "cs_haz" = haz_cs1
  )
  dat_ev2 <- data.table(
    "time" = t,
    "cause" = "2",
    "cuminc" = F2,
    "subdist_haz" = haz_subdist2,
    "cs_haz" = haz_cs2
  )
  rbind(dat_ev1, dat_ev2)
}
