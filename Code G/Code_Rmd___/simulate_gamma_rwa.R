# ============================================================================
# Fonctions pour la simulation du processus gamma (grille fine)
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

# Construction de la grille fine pour des temps entiers
construire_grille <- function(t_obs, rho = 1) {
  if (!is.numeric(t_obs) || length(t_obs) < 2) stop("t_obs doit être un vecteur numérique d'au moins 2 éléments.")
  if (any(diff(t_obs) <= 0)) stop("t_obs doit être strictement croissant.")
  if (!is.numeric(rho) || rho <= 0) stop("rho doit être un nombre positif.")
  differences <- diff(t_obs)
  delta <- differences[1]
  for (i in 2:length(differences)) delta <- pgcd(delta, differences[i])
  pas_fine <- delta / rho
  t_max <- max(t_obs)
  t_fine <- seq(0, t_max, by = pas_fine)
  idx_obs <- match(t_obs, t_fine)
  list(t_fine = t_fine, delta = delta, pas_fine = pas_fine,
       idx_obs = idx_obs, n_fine = length(t_fine), n_obs = length(t_obs))
}

# Construction de la grille fine pour des temps décimaux
construire_grille_dec <- function(t_obs, rho = 1, digits = 3) {
  if (!is.numeric(t_obs) || length(t_obs) < 2) stop("t_obs doit être un vecteur numérique d'au moins 2 éléments.")
  if (any(diff(t_obs) <= 0)) stop("t_obs doit être strictement croissant.")
  if (!is.numeric(rho) || rho <= 0) stop("rho doit être un nombre positif.")
  n_dec <- function(x) ifelse(grepl("\\.", as.character(x)), nchar(strsplit(as.character(x), "\\.")[[1]][2]), 0)
  max_dec <- max(sapply(t_obs, n_dec))
  if (max_dec <= 3) {
    facteur <- 10^max_dec
    t_obs_int <- round(t_obs * facteur)
    diff_int <- diff(t_obs_int)
    delta_int <- diff_int[1]
    for (i in 2:length(diff_int)) delta_int <- pgcd(delta_int, diff_int[i])
    delta <- delta_int / facteur
  } else {
    t_obs <- round(t_obs, digits)
    t_obs <- unique(t_obs)
    facteur <- 10^digits
    t_obs_int <- round(t_obs * facteur)
    diff_int <- diff(t_obs_int)
    delta_int <- diff_int[1]
    for (i in 2:length(diff_int)) delta_int <- pgcd(delta_int, diff_int[i])
    delta <- delta_int / facteur
  }
  pas_fine <- delta / rho
  t_max <- max(t_obs)
  t_fine <- seq(0, t_max, by = pas_fine)
  idx_obs <- match(round(t_obs / pas_fine, 10), round(t_fine / pas_fine, 10))
  list(t_fine = t_fine, delta = delta, pas_fine = pas_fine,
       idx_obs = idx_obs, t_obs_original = t_obs, n_fine = length(t_fine), n_obs = length(t_obs))
}

# Simulation par Random Walk Approximation (avec paramètre d'échelle)
simulate_gamma_rwa <- function(times, A = NULL, a = NULL, scale = 1, nsim = 1) {
  if (!is.numeric(times) || !is.vector(times)) stop("'times' doit être un vecteur numérique.")
  if (length(times) < 2) stop("'times' doit contenir au moins deux points.")
  if (any(diff(times) <= 0)) stop("'times' doit être strictement croissant.")
  if (!is.numeric(nsim) || nsim < 1 || nsim != floor(nsim)) stop("'nsim' doit être un entier positif.")
  if (scale <= 0) stop("'scale' doit être positif.")
  if (!is.null(A) && is.function(A)) {
    A_vals <- A(times)
  } else if (!is.null(a) && is.numeric(a) && length(a) == 1) {
    A_vals <- a * times
  } else {
    stop("Il faut fournir soit une fonction 'A', soit un scalaire 'a'.")
  }
  if (any(diff(A_vals) <= 0)) stop("La fonction A doit être strictement croissante.")
  diff_A <- diff(A_vals)
  n_intervals <- length(diff_A)
  increments <- matrix(rgamma(n_intervals * nsim, shape = rep(diff_A, each = nsim), scale = scale),
                       nrow = n_intervals, ncol = nsim, byrow = TRUE)
  X <- rbind(rep(0, nsim), apply(increments, 2, cumsum))
  rownames(X) <- paste0("t", seq_along(times) - 1)
  colnames(X) <- paste0("sim", seq_len(nsim))
  return(X)
}

# Fonction principale de simulation (avec paramètre b pour le taux)
simuler_gamma <- function(t_obs, A = NULL, a = NULL, b = 1, rho = 10, nsim = 1, digits = 3) {
  if (all(t_obs == round(t_obs))) {
    grille <- construire_grille(t_obs, rho)
    t_obs_clean <- t_obs
  } else {
    grille <- construire_grille_dec(t_obs, rho, digits)
    t_obs_clean <- grille$t_obs_original
  }
  if (!is.null(A) && is.function(A)) {
    A_func <- A
  } else if (!is.null(a) && is.numeric(a) && length(a) == 1) {
    A_func <- function(t) a * t
  } else {
    stop("Il faut fournir soit une fonction 'A', soit un scalaire 'a'.")
  }
  X_fine <- simulate_gamma_rwa(times = grille$t_fine, A = A_func, scale = 1/b, nsim = nsim)
  if (nsim == 1) {
    X_obs <- X_fine[grille$idx_obs, 1]
    obs <- data.frame(t = t_obs_clean, X = X_obs)
    traj_fine <- data.frame(t = grille$t_fine, X = X_fine[, 1])
  } else {
    X_obs <- X_fine[grille$idx_obs, ]
    obs <- data.frame(t = t_obs_clean, X_obs)
    colnames(obs)[-1] <- paste0("sim", 1:nsim)
    traj_fine <- data.frame(t = grille$t_fine, X_fine)
  }
  return(list(obs = obs, traj_fine = traj_fine, grille = grille))
}

# ============================================================================
# Scénarios de simulation (inspections régulières et irrégulières)
# ============================================================================

# Inspections régulières, processus homogène
sim_HGP <- function(n, a, b, tmax = 100) {
  t <- seq(0, tmax, length.out = n)
  X <- simuler_gamma(t_obs = t[-1], a = a, b = b, nsim = 1)
  data.frame(t = t, X = c(0, X$obs$X))
}

# Inspections régulières, processus non homogène
sim_NHGP <- function(n, a, b, alpha, tmax = 100) {
  t <- seq(0, tmax, length.out = n)
  A_func <- function(t) a * t^alpha
  X <- simuler_gamma(t_obs = t[-1], A = A_func, b = b, nsim = 1)
  data.frame(t = t, X = c(0, X$obs$X))
}

# Inspections irrégulières, processus homogène
sim_HGP_irregular <- function(t_obs, a, b) {
  X <- simuler_gamma(t_obs = t_obs, a = a, b = b, nsim = 1)
  data.frame(t = c(0, t_obs), X = c(0, X$obs$X))
}

# Inspections irrégulières, processus non homogène
sim_NHGP_irregular <- function(t_obs, a, b, alpha) {
  A_func <- function(t) a * t^alpha
  X <- simuler_gamma(t_obs = t_obs, A = A_func, b = b, nsim = 1)
  data.frame(t = c(0, t_obs), X = c(0, X$obs$X))
}
