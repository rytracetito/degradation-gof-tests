simulate_wiener_bridge <- function(X, digits = 1) {
  
  t_obs <- X[, 1]
  x_obs <- X[, 2]
  
  if (length(t_obs) != length(x_obs)) {
    stop("t_obs et x_obs doivent avoir la même longueur.")
  }
  if (length(t_obs) < 2) {
    stop("Il faut au moins 2 observations.")
  }
  
  params <- estimer_wiener(t_obs, x_obs)
  sigma <- params["sigma"]
  
  grille <- construire_grille_dec(t_obs, digits)
  
  delta <- grille$pas_fine
  t_fine <- grille$t_fine
  idx_obs <- grille$idx_obs
  
  n_fine <- length(t_fine)
  n_obs <- length(t_obs)
  
  X_fine <- numeric(n_fine)
  
  for (i in 1:n_obs) {
    X_fine[idx_obs[i]] <- x_obs[i]
  }
  
  for (i in 1:(n_obs - 1)) {
    idx_start <- idx_obs[i]
    idx_end <- idx_obs[i + 1]
    
    if (is.na(idx_start) || is.na(idx_end)) next
    
    t_start <- t_fine[idx_start]
    t_end <- t_fine[idx_end]
    
    x_start <- x_obs[i]
    x_end <- x_obs[i + 1]
    
    K <- idx_end - idx_start
    
    if (!is.na(K) && K >= 2) {
      X_current <- x_start
      
      for (k in 0:(K - 2)) {
        s_k <- t_start + k * delta
        s_kp1 <- s_k + delta
        
        mean_cond <- X_current + (delta / (t_end - s_k)) * (x_end - X_current)
        var_cond <- sigma^2 * delta * (t_end - s_kp1) / (t_end - s_k)
        
        X_next <- rnorm(1, mean = mean_cond, sd = sqrt(var_cond))
        
        X_fine[idx_start + k + 1] <- X_next
        X_current <- X_next
      }
    }
    
    X_fine[idx_end] <- x_end
  }
  
  return(list(
    t_fine = t_fine,
    X_fine = X_fine,
    idx_obs = idx_obs,
    grille = grille,
    sigma_estime = sigma
  ))
}