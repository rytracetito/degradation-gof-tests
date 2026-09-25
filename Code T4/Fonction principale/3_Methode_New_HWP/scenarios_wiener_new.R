# H0 : Wiener regulier 
sim_Wiener_reg <- function(n, mu, sigma, tmax = 100) {
  t <- seq(0, tmax, length.out = n)
  dt <- diff(t)
  dx <- rnorm(length(dt), mean = mu * dt, sd = sigma * sqrt(dt))
  x <- c(0, cumsum(dx))
  return(cbind(t, x))
}

# H0 : Wiener irregulier 
sim_Wiener_irreg <- function(t_obs, mu, sigma) {
  t_full <- c(0, t_obs)
  dt <- diff(t_full)
  dx <- rnorm(length(dt), mean = mu * dt, sd = sigma * sqrt(dt))
  x <- c(0, cumsum(dx))
  return(cbind(t_full, x))
}

# H1 : Gamma regulier
sim_Gamma_reg <- function(n, a, b, tmax = 100) {
  t <- seq(0, tmax, length.out = n)
  dt <- diff(t)
  dx <- rgamma(length(dt), shape = a * dt, rate = b)
  x <- c(0, cumsum(dx))
  return(cbind(t, x))
}

# H1 : Gamma irregulier
sim_Gamma_irreg <- function(t_obs, a, b) {
  t_full <- c(0, t_obs)
  dt <- diff(t_full)
  dx <- rgamma(length(dt), shape = a * dt, rate = b)
  x <- c(0, cumsum(dx))
  return(cbind(t_full, x))
}