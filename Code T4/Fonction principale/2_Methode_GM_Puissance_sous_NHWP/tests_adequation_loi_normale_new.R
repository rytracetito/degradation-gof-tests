# Les différents tests pour processus de Wiener

# Test de Anderson-Darling
gof.wiener.anderson <- function(dx, dt_eq) {
  if (!requireNamespace("nortest", quietly = TRUE)) install.packages("nortest")
  library(nortest)
  
  params <- estimer_wiener(seq_along(dx) * dt_eq, cumsum(c(0, dx))[-1])
  mu_est <- params["mu"]
  sigma_est <- params["sigma"]
  
  mean_ref <- mu_est * dt_eq
  sd_ref <- sigma_est * sqrt(dt_eq)
  
  dx_std <- (dx - mean_ref) / sd_ref
  ad.test(dx_std)
}

# Test de Cramér-von Mises
gof.wiener.cramer <- function(dx, dt_eq) {
  if (!requireNamespace("nortest", quietly = TRUE)) install.packages("nortest")
  library(nortest)
  
  params <- estimer_wiener(seq_along(dx) * dt_eq, cumsum(c(0, dx))[-1])
  mu_est <- params["mu"]
  sigma_est <- params["sigma"]
  
  mean_ref <- mu_est * dt_eq
  sd_ref <- sigma_est * sqrt(dt_eq)
  
  dx_std <- (dx - mean_ref) / sd_ref
  cvm.test(dx_std)
}

# Test de Shapiro-Wilk
gof.wiener.shapiro <- function(dx) {
  shapiro.test(dx)
}

# Test de Lilliefors
gof.wiener.lilliefors <- function(dx) {
  if (!requireNamespace("nortest", quietly = TRUE)) install.packages("nortest")
  library(nortest)
  
  lillie.test(dx)
}

# Test de D'Agostino
gof.wiener.agostino <- function(dx) {
  if (!requireNamespace("moments", quietly = TRUE)) install.packages("moments")
  library(moments)
  
  agostino.test(dx)
}