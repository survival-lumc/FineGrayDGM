weibull_hazard <- function(t, x, params, type = "hazard") {
  lp <- drop(x %*% params[["betas"]])
  rate <- params[["base_rate"]] * exp(lp)
  if (type == "hazard") {
    rate * params[["base_shape"]] * t^(params[["base_shape"]] - 1)
  } else { # type == "cumulative" for cumulative hazard
    rate * t^params[["base_shape"]]
  }
}

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


# Test ones
compute_true <- function(t,
                         model_type = c("squeezing", "reduction_factor", "two_fgs"),
                         newdat,
                         params) {

  predictor_formulas <- lapply(params, "[[", "formula")
  modmats <- lapply(predictor_formulas, function(form) {
    x <- model.matrix(form, data = newdat)
    x[, !colnames(x) %in% "(Intercept)"]
  })
  x_cause1 <- modmats[["cause1"]]
  x_cause2 <- modmats[["cause2"]]

  if (model_type == "reduction_factor") {

    haz_cs2 <- gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
    haz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard")

    integral_fun_cs1 <- function(t) {
      haz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "hazard")
      cumhaz_subdist1 <- gompertz_hazard(t, x_cause1, params[["cause1"]], type = "cumulative")
      cumhaz_cs2 <- gompertz_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
      haz_subdist1 * exp(-cumhaz_subdist1 + cumhaz_cs2)
    }
    num <- integral_fun_cs1(t)
    haz_cs1 <- num / (1 - integrate_to_t(fun = integral_fun_cs1, t = t))

    # Now for the CIs
    prod <- function(t, cause) {
      haz <- switch(
        cause,
        "1" = integral_fun_cs1(t) / (1 - integrate_to_t(fun = integral_fun_cs1, t = t)),
        "2" = gompertz_hazard(t, x_cause2, params[["cause2"]], type = "hazard")
      )
      cumhaz_cause1 <- -log(1 - integrate_to_t(fun = integral_fun_cs1, t = t))
      cumhaz_cause2 <- gompertz_hazard(t, x_cause2, params[["cause2"]], type = "cumulative")
      haz * exp(-cumhaz_cause1 - cumhaz_cause2)
    }

    # Calculate both cumulative incidences
    F1 <- integrate_to_t(fun = prod, t = t, cause = 1)
    F2 <- integrate_to_t(fun = prod, t = t, cause = 2)
    subdens_2 <- haz_cs2 * (1 - F1 - F2)
    haz_subdist2 <- subdens_2 / (1 - F2)

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

    get_subdisthaz_twofgs <- function(t, cause) {
      if (cause == 1) {
        nom <- p1 * shape_1 * rate_1 * exp(-rate_1 * t^shape_1) * t^(shape_1 - 1) #baseline subdens
        denom <- 1 - p1 * (1 - exp(-rate_1 * t^shape_1))
        hr_subdist1 * nom / denom
      } else {
        nom <- p2 * shape_2 * rate_2 * exp(-rate_2 * t^shape_2) * t^(shape_2 - 1) #baseline subdens
        denom <- 1 - p2 * (1 - exp(-rate_2 * t^shape_2)) # this I think is still ok even if in fraction
        hr_subdist2 * nom / denom
      }
    }
    haz_subdist1 <- get_subdisthaz_twofgs(t, cause = 1) # maybe these subdist hazards are wrong
    haz_subdist2 <- get_subdisthaz_twofgs(t, cause = 2)

    tfp_suscep <- p1^hr_subdist1 + p2^hr_subdist2 # this is like 1

    # Now for cs hazards.. ask Hein here...
    # p1^hr_subdist1 + p2^hr_subdist2 is max TFP??
    # I think it's CS hazards within the fraction?
    # See vertical model
    #haz_cs1 <- 0 # for now
   # haz_cs2 <- 0

    # Trying stuff
    haz_cs1 <- haz_subdist1 * (tfp_suscep - F1) / (tfp_suscep - F1 - F2)
    haz_cs2 <- haz_subdist2 * (tfp_suscep - F2) / (tfp_suscep - F1 - F2)

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

    # Closure for subdistribution hazard in this mechanism
    subdens_2 <- hr_condit * rate_2 * shape_2 * t^(shape_2 - 1) * exp(-hr_condit * rate_2 * t^shape_2)
    get_subdisthaz_squeezing <- function(t, cause) {
      if (cause == 1) {
        nom <- p * shape_1 * rate_1 * exp(-rate_1 * t^shape_1) * t^(shape_1 - 1)
        denom <- 1 - p * (1 - exp(-rate_1 * t^shape_1))
        hr_subdist * nom / denom
      } else {
        subdens_2 / (1 - F2)
      }
    }

    get_cshaz_squeezing <- function(t, cause) { # use switch??
      if (cause == 1) {
        # This is reduction factor (would it work both ways?)
        get_subdisthaz_squeezing(t, cause = 1) * (1 + F2 / (1 - F1 - F2)) # check for floating point issues?
      } else {
        subdens_2 * (1 - F1 - F2)
      }
    }

    haz_subdist1 <- get_subdisthaz_squeezing(t, cause = 1)
    haz_subdist2 <- get_subdisthaz_squeezing(t, cause = 2)
    haz_cs1 <- get_cshaz_squeezing(t, cause = 1)
    haz_cs2 <- get_cshaz_squeezing(t, cause = 2)
  }

  # Return everything; do long format?
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
