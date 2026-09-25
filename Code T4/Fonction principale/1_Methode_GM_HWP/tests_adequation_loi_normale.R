
# Test d'Anderson-Darling 
gof.wiener.anderson <- function(dx) {
  if (!requireNamespace("nortest", quietly = TRUE)) install.packages("nortest")
  library(nortest)
  
  res <- ad.test(dx)   
  return(res)
}

# Test de Cramér-von Mises 
gof.wiener.cramer <- function(dx) {
  if (!requireNamespace("nortest", quietly = TRUE)) install.packages("nortest")
  library(nortest)
  
  res <- suppressWarnings(cvm.test(dx))   
  return(res)
}

## Test de Shapiro-Wilk 
gof.wiener.shapiro <- function(dx) {
  shapiro.test(dx)
}

## Test de Lilliefors 
gof.wiener.lilliefors <- function(dx) {
  if (!requireNamespace("nortest", quietly = TRUE)) install.packages("nortest")
  library(nortest)
  
  lillie.test(dx)
}

## Test de Jarque-Bera 
gof.wiener.jarquebera <- function(dx) {
  if (!requireNamespace("tseries", quietly = TRUE)) install.packages("tseries")
  library(tseries)
  
  jarque.bera.test(dx)
}


# Test D'Agostino-Pearson 
gof.wiener.dagostino_pearson <- function(dx) {
  if (!requireNamespace("fBasics", quietly = TRUE)) install.packages("fBasics")
  library(fBasics)
  
  res <- dagoTest(dx)  
  return(res)
}

# Test complet de D'Agostino-Pearson 
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



