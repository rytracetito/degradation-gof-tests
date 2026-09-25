

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
construire_grille_dec <- function(t_obs, digits = 3) {
  
  # Vérifications des entrées
  if (!is.numeric(t_obs) || length(t_obs) < 2) {
    stop("t_obs doit être un vecteur numérique d'au moins 2 éléments.")
  }
  if (any(diff(t_obs) <= 0)) {
    stop("t_obs doit être strictement croissant.")
  }
  if (!is.numeric(digits) || digits < 0) {
    stop("digits doit être un nombre positif ou nul.")
  }
  
  # Fonction pour compter le nombre de décimales
  n_dec <- function(x) {
    if (grepl("\\.", as.character(x))) {
      return(nchar(strsplit(as.character(x), "\\.")[[1]][2]))
    } else {
      return(0)
    }
  }
  
  # Compter le nombre maximum de décimales
  max_dec <- max(sapply(t_obs, n_dec))
  
  # Cas où le nombre de décimales est petit (<= 3)
  if (max_dec <= 3) {
    facteur <- 10^max_dec
    t_obs_int <- round(t_obs * facteur)
    diff_int <- diff(t_obs_int)
    delta_int <- diff_int[1]
    for (i in 2:length(diff_int)) {
      delta_int <- pgcd(delta_int, diff_int[i])
    }
    delta <- delta_int / facteur
    
    # Cas où il y a trop de décimales
  } else {
    t_obs <- round(t_obs, digits)
    t_obs <- unique(t_obs)
    facteur <- 10^digits
    t_obs_int <- round(t_obs * facteur)
    diff_int <- diff(t_obs_int)
    delta_int <- diff_int[1]
    for (i in 2:length(diff_int)) {
      delta_int <- pgcd(delta_int, diff_int[i])
    }
    delta <- delta_int / facteur
  }
  
  # Calcul du pas de la grille fine (pas_fine = delta)
  pas_fine <- delta
  
  # Création de la grille fine
  t_max <- max(t_obs)
  t_fine <- seq(0, t_max, by = pas_fine)
  
  # Trouver les indices des temps d'observation
  idx_obs <- match(round(t_obs / pas_fine, 10), round(t_fine / pas_fine, 10))
  
  # Retourner les résultats
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