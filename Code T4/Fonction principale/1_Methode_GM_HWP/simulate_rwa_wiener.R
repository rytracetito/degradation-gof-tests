

# Simulation d'un processus de Wiener par RWA

simulate_wiener_rwa <- function(times, mu, sigma, nsim = 1) {
  
  # Vérification des entrées
  if (!is.numeric(times) || !is.vector(times)) {
    stop("'times' doit être un vecteur numérique.")
  }
  
  if (length(times) < 2) {
    stop("'times' doit contenir au moins deux points.")
  }
  
  if (any(diff(times) <= 0)) {
    stop("'times' doit être strictement croissant.")
  }
  
  if (times[1] != 0) {
    warning("Le premier temps n'est pas 0. Ajout de 0 au début.")
    times <- c(0, times)
  }
  
  if (!is.numeric(mu)) {
    stop("'mu' doit être numérique.")
  }
  
  if (!is.numeric(sigma) || sigma <= 0) {
    stop("'sigma' doit être un nombre positif.")
  }
  
  if (!is.numeric(nsim) || nsim < 1 || nsim != floor(nsim)) {
    stop("'nsim' doit être un entier positif.")
  }
  
  n <- length(times) - 1
  dt <- diff(times)
  
  increments <- matrix(
    data = rnorm(n = n * nsim, mean = 0, sd = 1),
    nrow = n,
    ncol = nsim
  )
  
  
  increments <- mu * dt + sigma * sqrt(dt) * increments
  
  trajectories <- rbind(
    rep(0, nsim),
    apply(increments, 2, cumsum)
  )
  
  rownames(trajectories) <- times
  colnames(trajectories) <- paste0("Traj_", 1:nsim)
  
  return(trajectories)
}




# Exemple d'utilisation
library(ggplot2)
library(tidyr)

set.seed(124)
times <- seq(0, 10, by = 0.1)
mu <- 2
sigma <- 1.5

resultat <- simulate_wiener_rwa(times, mu, sigma, nsim = 3)

df <- as.data.frame(resultat)
df$temps <- times

df_long <- pivot_longer(df, cols = starts_with("Traj_"), 
                        names_to = "Trajectoire", 
                        values_to = "X")

ggplot(df_long, aes(x = temps, y = X, color = Trajectoire)) +
  geom_step(linewidth = 1) +
  labs(title = paste("Processus de Wiener (mu =", mu, ", sigma =", sigma, ")"),
       x = "Temps t", y = "X(t)") +
  theme_minimal() +
  scale_color_manual(values = rainbow(3))
