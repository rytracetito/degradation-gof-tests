gof.wp <- function(X, test, B = 200) {
  t <- X[, 1]
  x <- X[, 2]
  dt <- diff(t)
  tol <- 1e-8
  
  if (all(abs(dt - dt[1]) < tol)) {
    dt_eq <- dt[1]
    dx_test <- diff(x)
  } else {
    params <- estimer_wiener(t, x)
    sigma_est <- params["sigma"]
    
    interp <- wiener_bridge(X = x, t = t, sigma = sigma_est)
    dt_eq <- interp$teq
    dx_test <- diff(interp$new_obs)
  }
  
  if (test == "Shapiro") return(gof.wiener.shapiro(dx_test))
  if (test == "Lilliefors") return(gof.wiener.lilliefors(dx_test))
  if (test == "Anderson") return(gof.wiener.anderson(dx_test, dt_eq))
  if (test == "Cramer") return(gof.wiener.cramer(dx_test, dt_eq))
  if (test == "Agostino") return(gof.wiener.agostino(dx_test))
  
  stop("Test inconnu. Choix: Shapiro, Lilliefors, Anderson, Cramer, Agostino")
}