set.seed(5565)


library(flexsurv)
library(survival)
n <- 100000
X <- rbinom(n, size = 1, prob = 0.5)

t_plot <- seq(0.001, 10, length.out = 500)
plot(t_plot, hgompertz(t_plot, shape = -0.25, rate = 1), type = "l", col = "blue")
lines(t_plot, hgompertz(t_plot, shape = -0.1, rate = 1))


T1 <- rgompertz(n, shape = -0.1, rate = exp(0.5 * X))
T2 <- rgompertz(n, shape = -0.25, rate = exp(0.25 * X))
time <- pmin(T1, T2)
D <- 1 + as.numeric(T2 < T1)
table(D)

dat <- cbind.data.frame(time, D, X)

coxph(Surv(time, D == 1) ~ X, data = dat)
coxph(Surv(time, D == 2) ~ X, data = dat)


