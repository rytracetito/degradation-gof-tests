gamma_bridge_multi_vis <- function(trajectoires, a) {
  total_time <- sum(sapply(trajectoires, function(tr) max(tr$t)))
  total_obs <- sum(sapply(trajectoires, function(tr) length(tr$t) - 1))
  h <- total_time / total_obs
  
  resultats <- list()
  for (i in seq_along(trajectoires)) {
    tr <- trajectoires[[i]]
    t_n <- max(tr$t)
    K <- floor(t_n / h)
    t_star <- seq(0, by = h, length.out = K + 1)
    X_star <- numeric(K + 1)
    X_star[1] <- 0
    
    for (j in 2:(K + 1)) {
      if (t_star[j] == t_n) {
        X_star[j] <- tr$X[length(tr$X)]
        next
      }
      inter <- max(tr$t[tr$t <= t_star[j]])
      if (inter > t_star[j - 1]) {
        t_inf <- inter
        x_inf <- tr$X[tr$t == inter]
      } else {
        t_inf <- t_star[j - 1]
        x_inf <- X_star[j - 1]
      }
      t_sup <- min(tr$t[tr$t > t_star[j]])
      x_sup <- tr$X[tr$t == t_sup]
      y <- rbeta(1, a * (t_star[j] - t_inf), a * (t_sup - t_star[j]))
      X_star[j] <- x_inf + (x_sup - x_inf) * y
    }
    
    if (t_star[length(t_star)] < t_n) {
      t_star <- c(t_star, t_n)
      X_star <- c(X_star, tr$X[length(tr$X)])
    }
    
    diff_X_star <- diff(X_star)
    resultats[[i]] <- list(t_star = t_star, X_star = X_star, diff_X_star = diff_X_star)
  }
  return(resultats)
}


# Application 1
# Simulation de deux trajectoires

set.seed(73)
simuler_gamma <- function(a, b, temps_obs) {
  delta_t <- diff(c(0, temps_obs))
  delta_X <- rgamma(length(delta_t), shape = a * delta_t, scale = 1/b)
  X <- cumsum(c(0, delta_X))
  list(t = c(0, temps_obs), X = X)
}


t_obs1 <- c(1, 3, 6, 10)
t_obs2 <- c(2, 5, 9)
traj1 <- simuler_gamma(2, 0.5, t_obs1)
traj2 <- simuler_gamma(2, 0.5, t_obs2)
trajectoires <- list(traj1, traj2)

# Application de la fonction
res <- gamma_bridge_multi_vis(trajectoires, a = 2)

# Extraction de tous les incréments i.i.d.
all_increments <- unlist(lapply(res, function(x) x$diff_X_star))

# Test de Kolmogorov-Smirnov 
h <- sum(sapply(trajectoires, function(tr) max(tr$t))) / sum(sapply(trajectoires, function(tr) length(tr$t)-1))
ks.test(all_increments, pgamma, shape = 2*h, scale = 2)


# Visualisation de la première trajectoire
plot(traj1$t, traj1$X, type = "b", col = "red", pch = 16, 
     xlab = "Temps", ylab = "Dégradation", 
     main = "Deux trajectoires - observé (cercles) et pont gamma (croix)",
     ylim = range(c(traj1$X, traj2$X, res[[1]]$X_star, res[[2]]$X_star)))
lines(res[[1]]$t_star, res[[1]]$X_star, type = "b", col = "blue", pch = 4)
lines(traj2$t, traj2$X, type = "b", col = "darkred", pch = 17)
lines(res[[2]]$t_star, res[[2]]$X_star, type = "b", col = "cyan", pch = 5)


## Application 2 

n1 <- 15
t_obs1 <- sort(runif(n1, min = 1, max = 10))
n2 <- 9
t_obs2 <- seq(2, 10, length.out = n2)

# Simulation
traj1 <- simuler_gamma(2, 0.5, t_obs1)
traj2 <- simuler_gamma(2, 0.5, t_obs2)
trajectoires <- list(traj1, traj2)


# Application du pont gamma multi
res <- gamma_bridge_multi_vis(trajectoires, a = 2)

# Extraction de tous les incréments i.i.d.
all_increments <- unlist(lapply(res, function(x) x$diff_X_star))

# Test KS
h <- sum(sapply(trajectoires, function(tr) max(tr$t))) / sum(sapply(trajectoires, function(tr) length(tr$t)-1))
ks.test(all_increments, pgamma, shape = 2*h, scale = 2)

# Visualisation des deux trajectoires sur une même figure
plot(traj1$t, traj1$X, type = "b", col = "red", pch = 16, 
     xlab = "Temps", ylab = "Dégradation", 
     main = "Deux trajectoires (observé = cercles, pont gamma = croix)",
     ylim = range(c(traj1$X, traj2$X, res[[1]]$X_star, res[[2]]$X_star)))
lines(res[[1]]$t_star, res[[1]]$X_star, type = "b", col = "blue", pch = 4)
lines(traj2$t, traj2$X, type = "b", col = "darkred", pch = 17)
lines(res[[2]]$t_star, res[[2]]$X_star, type = "b", col = "cyan", pch = 5)
legend("topleft", 
       legend = c("Traj1 observée", "Traj1 pont gamma", "Traj2 observée", "Traj2 pont gamma"),
       col = c("red", "blue", "darkred", "cyan"), 
       pch = c(16, 4, 17, 5), lty = 1)


## Application (Données réelles)

library(readr)
# Lecture des données
fichier <- "C:/Users/ADMIN/Desktop/Stage UPPA/Code R/Code T6/Dim_1.csv"
M <- as.matrix(read.csv2(fichier, header = FALSE))

temps_trajectoires <- list(
  c(0, 1.0, 2.2, 3.5, 4.8, 6.0, 7.1, 8.3, 9.2, 10.0),   
  c(0, 0.8, 1.9, 3.0, 4.2, 5.5, 6.7, 7.9, 9.0, 10.2),   
  c(0, 1.2, 2.5, 3.8, 5.0, 6.2, 7.4, 8.6, 9.8, 11.0),   
  c(0, 0.9, 2.0, 3.2, 4.5, 5.8, 7.0, 8.2, 9.5, 10.5),   
  c(0, 1.1, 2.3, 3.6, 4.9, 6.1, 7.3, 8.5, 9.7, 10.8),   
  c(0, 0.7, 1.8, 2.9, 4.1, 5.3, 6.5, 7.7, 8.9, 10.1),   
  c(0, 1.3, 2.6, 3.9, 5.2, 6.4, 7.6, 8.8, 10.0, 11.2),  
  c(0, 0.6, 1.7, 2.8, 4.0, 5.2, 6.4, 7.6, 8.8, 10.0),   
  c(0, 1.4, 2.7, 4.0, 5.3, 6.6, 7.9, 9.2, 10.4, 11.6),  
  c(0, 0.5, 1.6, 2.7, 3.9, 5.1, 6.3, 7.5, 8.7, 9.9)     
)


trajectoires <- list()
for (i in 1:n_traj) {
  trajectoires[[i]] <- list(
    t = temps_trajectoires[[i]],
    X = M[, i]   
  )
}

# Estimation des paramètres (a, b) par méthode des moments agrégée
estimer_multi <- function(all_dt, all_dx) {
  mu_hat <- sum(all_dx) / sum(all_dt)
  denom <- sum(all_dt) - sum(all_dt^2) / sum(all_dt)
  if (denom <= 0) stop("Dénominateur nul ou négatif")
  nu_hat <- sum((all_dx - mu_hat * all_dt)^2) / denom
  a_hat <- mu_hat^2 / nu_hat
  b_hat <- mu_hat / nu_hat
  c(a = a_hat, b = b_hat)
}

# Collecte des incréments bruts de toutes les trajectoires
all_dt <- c()
all_dx <- c()
for (i in 1:length(trajectoires)) {
  dt <- diff(trajectoires[[i]]$t)
  dx <- diff(trajectoires[[i]]$X)
  all_dt <- c(all_dt, dt)
  all_dx <- c(all_dx, dx)
}
params <- estimer_multi(all_dt, all_dx)
a_hat <- params["a"]
b_hat <- params["b"]

# Application du pont gamma multi‑trajectoires
res <- gamma_bridge_multi_vis(trajectoires, a = a_hat)

# Récupération des incréments réguliers i.i.d.
all_increments <- unlist(lapply(res, function(x) x$diff_X_star))

# Calcul du pas commun h (identique à celui utilisé dans la fonction)
total_time <- sum(sapply(trajectoires, function(tr) max(tr$t)))
total_obs <- sum(sapply(trajectoires, function(tr) length(tr$t) - 1))
h <- total_time / total_obs

# Test KS
ks.test(all_increments, pgamma, shape = a_hat * h, scale = 1 / b_hat)

# Visualisation 
#traj 1 
plot(trajectoires[[6]]$t, trajectoires[[6]]$X,
     type = "b", col = "red", pch = 16,
     xlab = "Temps", ylab = "Dégradation", main = "Trajectoires")
lines(res[[6]]$t_star, res[[6]]$X_star,
      type = "b", col = "blue", pch = 16)

#traj 2 
lines(trajectoires[[2]]$t, trajectoires[[2]]$X,
     type = "b", col = "red", pch = 16,
     xlab = "Temps", ylab = "Dégradation", main = "Trajectoires")
lines(res[[2]]$t_star, res[[2]]$X_star,
      type = "b", col = "blue", pch = 16)


#traj 3 
lines(trajectoires[[3]]$t, trajectoires[[3]]$X,
      type = "b", col = "red", pch = 16,
      xlab = "Temps", ylab = "Dégradation", main = "Trajectoires")
lines(res[[3]]$t_star, res[[3]]$X_star,
      type = "b", col = "blue", pch = 16)

#traj 4 
lines(trajectoires[[1]]$t, trajectoires[[1]]$X,
      type = "b", col = "red", pch = 16,
      xlab = "Temps", ylab = "Dégradation", main = "Trajectoires")
lines(res[[1]]$t_star, res[[1]]$X_star,
      type = "b", col = "blue", pch = 16)

legend("topleft", legend = c("Observé", "Pont gamma"), col = c("red", "blue"), pch = c(16, 16))
