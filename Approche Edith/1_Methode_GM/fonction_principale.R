gof.gp <- function(X, test) {
  t <- X[, 1]
  x <- X[, 2]
  dt <- diff(t)
  dx <- diff(x)
  tol <- 1e-8
  
  if (all(abs(dt - dt[1]) < tol)) {
    dx_test <- dx
  } else {
    params <- estimer(t, x)
    a_est <- as.numeric(params["a"])
    interp <- gamma_bridge(x, t, a_est)
    dx_test <- diff(interp$new_obs)
  }
  
  switch(test,
         "Henze" = gof.gp.henze(dx_test),
         "VGE"   = gof.gp.vge(dx_test),
         "KS"    = gof.gp.ks(dx_test),
         "BET"   = gof.gp.betsch(dx_test),
         stop("Test inconnu : '", test, "'. Choix possibles : 'Henze', 'VGE', 'KS', 'BET'.")
  )
}