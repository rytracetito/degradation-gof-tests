# Processus Gamma Homogène 
sim_HGP <- function(t_obs, a, b) {
  X <- simulate_gamma_rwa(times = t_obs, a = a, b = b, nsim = 1)
  return(data.frame(t = t_obs, X = X[, 1]))
}

# Processus Gamma Non-Homogène 
sim_NHGP <- function(t_obs, a, b, alpha) {
  A_func <- function(t) a * t_obs^alpha
  X <- simulate_gamma_rwa(times = t_obs, A = A_func, b = b, nsim = 1)
  return(data.frame(t = t_obs, X = X[, 1]))
}

# Processus Gaussien Inverse Homogène
sim_IGP <- function(t_obs, mu, lambda) {
  if (mu <= 0)     stop("mu doit être strictement positif")
  if (lambda <= 0) stop("lambda doit être strictement positif")
  
  n <- length(t_obs) - 1
  delta_t <- diff(t_obs)
  
  delta_X <- rinvgauss(n, mean = mu * delta_t, shape = lambda * delta_t^2)
  
  if (any(delta_X <= 0) || any(!is.finite(delta_X))) {
    stop("Incréments invalides générés (négatifs, nuls ou non finis) : ",
         "vérifie mu, lambda et t_obs")
  }
  
  X <- c(0, cumsum(delta_X))
  return(data.frame(t = t_obs, X = X))
}



