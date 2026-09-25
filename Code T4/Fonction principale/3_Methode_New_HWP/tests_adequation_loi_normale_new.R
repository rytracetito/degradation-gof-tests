# Les différents tests pour processus de Wiener

# Test d'Anderson-Darling avec estimation interne
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

# Test de Cramér-von Mises avec estimation interne
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

# Test de Jarque-Bera
gof.wiener.jarquebera <- function(dx) {
  if (!requireNamespace("tseries", quietly = TRUE)) install.packages("tseries")
  library(tseries)
  
  jarque.bera.test(dx)
}

# Test de D'Agostino-Pearson (omnibus) avec fBasics
gof.wiener.agostino <- function(dx) {
  if (!requireNamespace("fBasics", quietly = TRUE)) install.packages("fBasics")
  library(fBasics)
  
  res <- dagoTest(dx)
  p_omnibus <- res@test$p.value[1]
  stat_omnibus <- res@test$statistic[1]
  
  result <- list(
    statistic = c(`Chi-squared` = stat_omnibus),
    parameter = c(df = 2),
    p.value = p_omnibus,
    method = "D'Agostino-Pearson omnibus normality test",
    data.name = deparse(substitute(dx))
  )
  class(result) <- "htest"
  return(result)
}