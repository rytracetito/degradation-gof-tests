gof.wp <- function(X, test) {
  t <- X[, 1]
  x <- X[, 2]
  dt <- diff(t)
  tol <- 1e-8
  
  if (all(abs(dt - dt[1]) < tol)) {
    dt_eq <- dt[1]
    dx_test <- diff(x)
  } else {
    resultat_bridge <- simulate_wiener_bridge(X = X, digits = 1)
    dt_eq <- resultat_bridge$grille$pas_fine
    dx_test <- diff(resultat_bridge$X_fine)
  }
  
  if (test == "Shapiro") return(gof.wiener.shapiro(dx_test))
  if (test == "Lilliefors") return(gof.wiener.lilliefors(dx_test))
  if (test == "JarqueBera") return(gof.wiener.jarquebera(dx_test))
  if (test == "Anderson") return(gof.wiener.anderson(dx_test, dt_eq))
  if (test == "Cramer") return(gof.wiener.cramer(dx_test, dt_eq))
  if (test == "Agostino") return(gof.wiener.agostino(dx_test))
  
  stop("Test inconnu. Choix: Shapiro, Lilliefors, JarqueBera, Anderson, Cramer, Agostino")
}