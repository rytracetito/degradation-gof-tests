set.seed(123)

# ---- 1. Pont gamma sur une grille fournie ----

gamma_bridge_on_grid <- function(X, t, a, t_star) {
  n_star <- length(t_star)
  Y <- numeric(n_star)
  Y[1] <- X[1]
  Y[n_star] <- X[length(X)]
  
  if (n_star > 2) {
    for (i in 2:(n_star - 1)) {
      ts <- t_star[i]
      inter <- max(t[t <= ts])
      if (inter > t_star[i - 1]) {
        t_inf <- inter
        x_inf <- X[t == inter]
      } else {
        t_inf <- t_star[i - 1]
        x_inf <- Y[i - 1]
      }
      t_sup <- min(t[t > ts])
      x_sup <- X[t == t_sup]
      y <- rbeta(1, a * (ts - t_inf), a * (t_sup - ts))
      Y[i] <- x_inf + (x_sup - x_inf) * y
    }
  }
  Y
}

# ---- 2. Pont gamma multi-trajectoires, grille fine régulière globale ----

gamma_bridge_multi <- function(List.Df, a) {
  
  tol = 1e-8
  tps <- NULL
  for (i in seq_along(List.Df)) tps <- c(tps, List.Df[[i]]$t)
  tps <- sort(unique(tps))
  
  n   <- length(tps) - 1
  t_n <- tps[length(tps)]
  teq <- t_n / n
  t_star_global <- seq(0, t_n, teq)
  
  resultats <- vector("list", length(List.Df))
  for (i in seq_along(List.Df)) {
    t_i <- List.Df[[i]]$t
    X_i <- List.Df[[i]]$x
    t_n_i <- max(t_i)
    
    t_star_i <- t_star_global[t_star_global <= t_n_i + tol]
    if (t_n_i - max(t_star_i) > tol) {
      t_star_i <- c(t_star_i, t_n_i)
    }
    
    Y_i <- gamma_bridge_on_grid(X_i, t_i, a, t_star_i)
    resultats[[i]] <- data.frame(t = t_star_i, x = Y_i)
  }
  resultats
}

# ---- 3. Simulation de 2 trajectoires gamma à instants irréguliers ----

simuler_traj_gamma <- function(t_obs, shape_rate, scale) {
  n <- length(t_obs) - 1
  increments <- rgamma(n, shape = diff(t_obs) * shape_rate, scale = scale)
  x <- c(0, cumsum(increments))
  data.frame(t = t_obs, x = x)
}

t1 <- c(0, 2, 5, 9, 11)
traj1 <- simuler_traj_gamma(t1, shape_rate = 1, scale = 2)

t2 <- c(0, 1, 3, 4, 7, 10, 13, 15)
traj2 <- simuler_traj_gamma(t2, shape_rate = 1, scale = 2)

List.Df <- list(traj1, traj2)

# ---- 4. Application du pont ----

a <- 5
resultats <- gamma_bridge_multi(List.Df, a)

# Vérifications utiles avant de tracer
print(resultats[[1]])
print(resultats[[2]])

# ---- 5. Affichage : observé vs simulé ----

plot(NA, xlim = c(0, 15),
     ylim = c(0, max(sapply(List.Df, function(d) max(d$x))) * 1.1),
     xlab = "Temps", ylab = "Dégradation X(t)",
     main = "Pont Gamma multi-trajectoires (grille fine régulière globale)")

couleurs <- c("blue", "darkgreen")

for (i in seq_along(List.Df)) {
  points(List.Df[[i]]$t, List.Df[[i]]$x, col = couleurs[i], pch = 16, cex = 1.3)
  lines(List.Df[[i]]$t, List.Df[[i]]$x, col = couleurs[i], lty = 1)
  
  points(resultats[[i]]$t, resultats[[i]]$x, col = couleurs[i], pch = 4)
  lines(resultats[[i]]$t, resultats[[i]]$x, col = couleurs[i], lty = 3)
}


# Accroissement 

increments_list <- lapply(resultats, function(df) diff(df$x))

Z <- unlist(increments_list)
Z <- Z[Z > 1e-6] 
print(Z)

# Test d'adéquation

library(gofgamma)
library(goft)
## Test de Henze 
res_henze <- test.HME1(Z, a = 1, boot = 200, alpha = 0.05)
print(res_henze)

## Test de Villaseñor 
res_villasenor <- gamma_test(Z)
print(res_villasenor)



