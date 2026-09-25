

# ============================================================================
# Construction d'une grille fine pour temps entiers
# ============================================================================

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
# Fonction construire_grille pour temps entiers
construire_grille <- function(t_obs) {
  
  # Vérifications des entrées
  if (!is.numeric(t_obs) || length(t_obs) < 2) {
    stop("t_obs doit être un vecteur numérique d'au moins 2 éléments.")
  }
  if (any(diff(t_obs) <= 0)) {
    stop("t_obs doit être strictement croissant.")
  }
  
  # Calcul du PGCD des différences
  differences <- diff(t_obs)
  delta <- differences[1]
  for (i in 2:length(differences)) {
    delta <- pgcd(delta, differences[i])
  }
  
  # Construction de la grille fine (pas_fine = delta)
  pas_fine <- delta
  t_max <- max(t_obs)
  t_fine <- seq(0, t_max, by = pas_fine)
  idx_obs <- match(t_obs, t_fine)
  
  return(list(
    t_fine = t_fine,
    delta = delta,
    pas_fine = pas_fine,
    idx_obs = idx_obs,
    n_fine = length(t_fine),
    n_obs = length(t_obs)
  ))
}

# Fonction construire_grille_dec pour temps décimaux
construire_grille_dec <- function(t_obs, digits = 1) {
  
  if (!is.numeric(t_obs) || length(t_obs) < 2) {
    stop("t_obs doit être un vecteur numérique d'au moins 2 éléments.")
  }
  if (any(diff(t_obs) <= 0)) {
    stop("t_obs doit être strictement croissant.")
  }
  if (!is.numeric(digits) || digits < 0) {
    stop("digits doit être un nombre positif ou nul.")
  }
  
  n_dec <- function(x) {
    if (grepl("\\.", as.character(x))) {
      return(nchar(strsplit(as.character(x), "\\.")[[1]][2]))
    } else {
      return(0)
    }
  }
  
  max_dec <- max(sapply(t_obs, n_dec))
  
  if (max_dec <= digits) {
    facteur <- 10^max_dec
    t_obs_int <- round(t_obs * facteur)
    diff_int <- diff(t_obs_int)
    delta_int <- diff_int[1]
    for (i in 2:length(diff_int)) {
      delta_int <- pgcd(delta_int, diff_int[i])
    }
    delta <- delta_int / facteur
  } else {
    t_obs_round <- round(t_obs, digits)
    t_obs_round <- unique(t_obs_round)
    facteur <- 10^digits
    t_obs_int <- round(t_obs_round * facteur)
    diff_int <- diff(t_obs_int)
    delta_int <- diff_int[1]
    for (i in 2:length(diff_int)) {
      delta_int <- pgcd(delta_int, diff_int[i])
    }
    delta <- delta_int / facteur
  }
  
  pas_fine <- delta
  
  t_max <- max(t_obs)
  t_fine <- seq(0, t_max, by = pas_fine)
  
  idx_obs <- sapply(t_obs, function(x) which.min(abs(t_fine - x)))
  
  return(list(
    t_fine = t_fine,
    delta = delta,
    pas_fine = pas_fine,
    idx_obs = idx_obs,
    t_obs = t_obs,
    n_fine = length(t_fine),
    n_obs = length(t_obs)
  ))
}