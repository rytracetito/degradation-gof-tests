
# Pont Brownien (Grall Maes)

wiener_bridge <- function(X, t, sigma) {
  # Vérifications des entrées
  if (length(X) != length(t)) stop("X et t doivent avoir la même longueur")
  if (any(diff(t) <= 0)) stop("Les temps doivent être strictement croissants")
  if (sigma <= 0) stop("Le paramètre sigma doit être positif")
  if (X[1] != 0) warning("La première valeur de X devrait être 0 pour un processus de Wiener")
  
  n <- length(X) - 1
  t_n <- t[length(t)]
  teq <- t_n / n
  
  t_star <- seq(0, t_n, teq)
  
  Y <- numeric(length(t_star))
  Y[1] <- 0
  
  for (i in 2:(length(t_star)-1)) {
    
    # Trouver la borne inférieure
    inter <- max(t[t <= t_star[i]])
    if (inter > t_star[i - 1]) {
      t_inf <- inter
      x_inf <- X[t == inter]
    } else {
      t_inf <- t_star[i - 1]
      x_inf <- Y[i - 1]
    }
    
    # Trouver la borne supérieure
    t_sup <- min(t[t > t_star[i]])
    x_sup <- X[t == t_sup]
    
    # Paramètres du pont brownien
    mu <- x_inf + (x_sup - x_inf) * (t_star[i] - t_inf) / (t_sup - t_inf)
    var_cond <- sigma^2 * (t_star[i] - t_inf) * (t_sup - t_star[i]) / (t_sup - t_inf)
    
    # Génération de la valeur
    y <- rnorm(1, mean = mu, sd = sqrt(var_cond))
    
    Y[i] <- y
  }
  
  Y[length(t_star)] <- X[length(X)]
  
  result <- list(teq = teq, new_obs = Y, time = t_star)
  return(result)
}

# Application 
library(ggplot2)

# Paramètres
set.seed(42)
sigma <- 1.5
a_drift <- 0.5  

# Temps irréguliers
t_obs <- c(0, 1.2, 2, 3.8, 5.5, 7.1, 10.4, 12.9, 15.3, 18.7, 22.5)

# Simulation d'un Wiener avec dérive
x_obs <- a_drift * t_obs + cumsum(c(0, rnorm(length(t_obs)-1, mean = 0, sd = sigma * sqrt(diff(t_obs)))))

# Application du pont brownien
bridge <- wiener_bridge(X = x_obs, t = t_obs, sigma = sigma)

# Graphique
df_obs <- data.frame(t = t_obs, x = x_obs, type = "Observé")
t_reg <- seq(0, max(t_obs), by = bridge$teq)
df_sim <- data.frame(t = t_reg, x = bridge$new_obs, type = "Pont brownien")

ggplot() +
  geom_point(data = df_obs, aes(t, x, color = type), size = 3) +
  geom_point(data = df_sim, aes(t, x, color = type), size = 2, shape = 16) +
  scale_color_manual(values = c("Observé" = "red", "Pont brownien" = "blue")) +
  labs(title = "Pont brownien avec dérive",
       x = "Temps", y = "X(t)") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
