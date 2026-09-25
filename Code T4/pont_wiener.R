Wiener_bridge <- function(T, mu, sigma, N) {
  
  # T : Temps final ([0,T])
  # sigma : Volatilité  
  # mu : Dérive
  # N : Niveau de raffinement (résolution = 2^N intervalles)
  
  n_points <- 2^N + 1
  
  # Grille de temps
  t <- seq(0, T, length.out = n_points)
  
  X <- numeric(n_points)
  
  normals <- rnorm(2^N - 1)
  idx <- 1
  
  # Initialisation
  X[1] <- 0
  X[n_points] <- mu * T + sigma * sqrt(T) * normals[idx]
  idx <- idx + 1
  
  
  known_idx <- c(1, n_points)
  
  
  for (n in 1:N) {
    new_idx <- numeric(length(known_idx) - 1)
    pos <- 1
    
    for (k in 1:(length(known_idx) - 1)) {
      left <- known_idx[k]
      right <- known_idx[k + 1]
      
      mid <- as.integer((left + right) / 2)
      
      delta_t <- t[right] - t[left]
      
      X[mid] <- (X[left] + X[right]) / 2 +
        sigma * sqrt(delta_t) / 2 * normals[idx]
      
      idx <- idx + 1
      new_idx[pos] <- mid
      pos <- pos + 1
    }
    
    known_idx <- sort(c(known_idx, new_idx))
  }
  
  return(list(time = t, X = X))
}
# Application 
library(ggplot2)

# Simulation
result <- Wiener_bridge(T = 1, mu = 0.5, sigma = 0.3, N = 20)

# Graphique
df <- data.frame(temps = result$time, X = result$X)

ggplot(df, aes(x = temps, y = X)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(intercept = 0, slope = 0.5, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Simulation d'un Wiener avec dérive",
       x = "Temps", y = "X(t)") +
  theme_minimal()






# ============================================================================
# SIMULATION D'UN MOUVEMENT BROWNIEN PAR PONT BROWNIEN (ALGORITHM 1.9)
# ============================================================================
# Modèle : B(t) mouvement brownien standard (sans dérive)
# ============================================================================

brownian_bridge_simulation <- function(T, N) {
  # T : Temps final
  # N : Niveau de raffinement (nombre de subdivisions = 2^N)
  # Sortie : liste contenant les temps et les valeurs du brownien
  
  n_points <- 2^N + 1
  
  # Grille de temps
  t <- seq(0, T, length.out = n_points)
  
  # Vecteur des valeurs
  B <- numeric(n_points)
  
  # Génération des variables normales standards (2^N - 1 suffisent)
  normals <- rnorm(2^N - 1)
  idx <- 1
  
  # Step 0
  t[1] <- 0
  t[n_points] <- T
  B[1] <- 0
  B[n_points] <- sqrt(T) * normals[idx]
  idx <- idx + 1
  
  # Liste des indices connus
  known_idx <- c(1, n_points)
  
  # Boucle principale (Step n+1)
  for (n in 1:N) {
    new_idx <- numeric(length(known_idx) - 1)
    pos <- 1
    
    for (k in 1:(length(known_idx) - 1)) {
      left <- known_idx[k]
      right <- known_idx[k + 1]
      
      mid <- as.integer((left + right) / 2)
      
      delta_t <- t[right] - t[left]
      
      # Formule du pont brownien au milieu
      B[mid] <- (B[left] + B[right]) / 2 + sqrt(delta_t) / 2 * normals[idx]
      
      idx <- idx + 1
      new_idx[pos] <- mid
      pos <- pos + 1
    }
    
    known_idx <- sort(c(known_idx, new_idx))
  }
  
  return(list(time = t, B = B))
}

# ============================================================================
# APPLICATION
# ============================================================================

library(ggplot2)

# Simulation d'une trajectoire
result <- brownian_bridge_simulation(T = 1, N = 10)

# Graphique
df <- data.frame(temps = result$time, B = result$B)

ggplot(df, aes(x = temps, y = B)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  labs(title = "Simulation d'un mouvement brownien standard par pont brownien",
       x = "Temps", y = "B(t)") +
  theme_minimal()




wiener_bridge <- function(t_obs, x_obs, sigma, dt_eq) {
  
  
  # Vérifications des entrées
  if (any(diff(t_obs) <= 0)) stop("Les temps doivent être strictement croissants")
  if (sigma <= 0) stop("Le paramètre sigma doit être positif")
  if (dt_eq <= 0) stop("dt_eq doit être positif")
  
  t_full <- c(0, t_obs)
  x_full <- c(0, x_obs)
  
  # Nombre de points à générer
  Np <- floor(max(t_obs) / dt_eq)
  
  t_gen <- numeric(Np)
  x_gen <- numeric(Np)
  
  t_prev <- 0
  x_prev <- 0
  
  for (i in 1:Np) {
    ti <- i * dt_eq
    
    # Trouver l'intervalle contenant ti
    j <- which(t_full >= ti)[1]
    t_sup <- t_full[j]
    x_sup <- x_full[j]
    
    # Déterminer la borne inférieure (corrigé)
    if (j > 1 && t_full[j - 1] > t_prev) {
      t_inf <- t_full[j - 1]
      x_inf <- x_full[j - 1]
    } else {
      t_inf <- t_prev
      x_inf <- x_prev
    }
    
    # Paramètres du pont brownien
    delta_t <- t_sup - t_inf
    mean_cond <- ((t_sup - ti) / delta_t) * x_inf + ((ti - t_inf) / delta_t) * x_sup
    var_cond <- sigma^2 * ((ti - t_inf) * (t_sup - ti)) / delta_t
    
    # Génération de la valeur
    x_new <- mean_cond + sqrt(var_cond) * rnorm(1)
    
    t_gen[i] <- ti
    x_gen[i] <- x_new
    t_prev <- ti
    x_prev <- x_new
  }
  
  return(data.frame(time = t_gen, value = x_gen))
}


# Exemple d'utilisation
# ============================================================================
# VISUALISATION : TRAJECTOIRE INTERPOLÉE + OBSERVATIONS ORIGINALES
# ============================================================================

library(ggplot2)

# Exemple de données
t_obs <- c(1, 3, 6, 10)
x_obs <- c(2, 5, 9, 15)
sigma <- 0.5
dt_eq <- 0.5

# Génération de la trajectoire interpolée
resultat <- wiener_bridge(t_obs, x_obs, sigma, dt_eq)


ggplot() +
  # Trajectoire interpolée (pont brownien)
  geom_line(data = resultat, aes(x = time, y = value), 
            color = "blue", size = 1) +
  # Points des observations originales
  geom_point(data = data.frame(t_obs, x_obs), 
             aes(x = t_obs, y = x_obs), 
             color = "red", size = 3) +
  # Étiquettes
  labs(title = "Interpolation par pont brownien",
       subtitle = paste("sigma =", sigma, ", Delta t_eq =", dt_eq),
       x = "Temps", y = "Valeur") +
  theme_minimal() +
  # Ajouter une légende manuelle
  annotate("text", x = max(t_obs) * 0.8, y = max(x_obs) * 0.9, 
           label = "Trajectoire interpolée", color = "blue", size = 4) +
  annotate("text", x = max(t_obs) * 0.8, y = max(x_obs) * 0.85, 
           label = "Observations originales", color = "red", size = 4)







