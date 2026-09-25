# Fonction PGCD
pgcd <- function(a, b) {
  a <- abs(a)
  b <- abs(b)
  while (b != 0) {
    temp <- b
    b <- a %% b
    a <- temp
  }
  return(a)
}

# Construire grille
construire_grille <- function(t_obs, digits = 1) {
  
  if (!is.numeric(t_obs) || length(t_obs) < 2) {
    stop("t_obs doit être un vecteur numérique d'au moins 2 éléments.")
  }
  if (any(diff(t_obs) <= 0)) {
    stop("t_obs doit être strictement croissant.")
  }
  if (!is.numeric(digits) || digits < 0) {
    stop("digits doit être un nombre positif ou nul.")
  }
  
  t_obs_round <- round(t_obs, digits)
  t_obs_unique <- unique(t_obs_round)
  
  if (length(t_obs_unique) < 2) {
    stop("Après arrondi à ", digits, " décimales, moins de 2 temps distincts restent : ",
         "augmente 'digits'.")
  }
  
  facteur <- 10^digits
  t_obs_int <- round(t_obs_unique * facteur)
  differences <- diff(t_obs_int)
  
  delta_int <- differences[1]
  if (length(differences) > 1) {
    for (i in 2:length(differences)) {
      delta_int <- pgcd(delta_int, differences[i])
    }
  }
  
  pas_fine <- delta_int / facteur
  t_max <- max(t_obs_round)
  t_fine <- seq(0, t_max, by = pas_fine)
  
  list(
    t_fine   = t_fine,
    pas_fine = pas_fine
  )
}

# Pont Gamma 

simulate_gamma_grille <- function(X, t, a) {
  
  grille <- construire_grille(t)
  t_star <- grille$t_fine
  pas_fine <- grille$pas_fine
  n_star <- length(t_star) - 1
  
  Y <- numeric(length(t_star))
  Y[1] <- 0
  
  if (n_star > 1) {
    for (i in 2:n_star) {
      inter <- max(t[t <= t_star[i]])
      if (inter > t_star[i - 1]) {
        t_inf <- inter
        x_inf <- X[t == inter]
      } else {
        t_inf <- t_star[i - 1]
        x_inf <- Y[i - 1]
      }
      t_sup <- min(t[t > t_star[i]])
      x_sup <- X[t == t_sup]
      y <- rbeta(1, a * (t_star[i] - t_inf), a * (t_sup - t_star[i]))
      Y[i] <- x_inf + (x_sup - x_inf) * y
    }
  }
  
  Y[length(t_star)] <- X[length(X)]
  
  list(pas_fine = pas_fine, t_fine = t_star, new_obs = Y)
}

