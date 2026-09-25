# Scénarios pour l'étude Monte Carlo

# H0 : Wiener régulier 
sim_Wiener_reg <- function(n, mu, sigma, tmax = 100) {
  t <- seq(0, tmax, length.out = n)
  dt <- diff(t)
  dx <- rnorm(length(dt), mean = mu * dt, sd = sigma * sqrt(dt))
  x <- c(0, cumsum(dx))
return(cbind(t, x))
}

# H0 : Wiener irrégulier 
sim_Wiener_irreg <- function(t_obs, mu, sigma) {
  t_full <- c(0, t_obs)
  dt <- diff(t_full)
  dx <- rnorm(length(dt), mean = mu * dt, sd = sigma * sqrt(dt))
  x <- c(0, cumsum(dx))
  return(cbind(t_full, x))
}

# H1 : Gamma régulier
sim_Gamma_reg <- function(n, a, b, tmax = 100) {
  t <- seq(0, tmax, length.out = n)
  dt <- diff(t)
  dx <- rgamma(length(dt), shape = a * dt, rate = b)
  x <- c(0, cumsum(dx))
  return(cbind(t, x))
}

# H1 : Gamma irrégulier
sim_Gamma_irreg <- function(t_obs, a, b) {
  t_full <- c(0, t_obs)
  dt <- diff(t_full)
  dx <- rgamma(length(dt), shape = a * dt, rate = b)
  x <- c(0, cumsum(dx))
  return(cbind(t_full, x))
}


# H1 : Wiener non homogène (NHWP)

sim_Wiener_nonhomog <- function(t_obs, mu_a, delta_a, sigma_true) {
  a_nhom <- rnorm(1, mean = mu_a, sd = delta_a)
  tau <- t_obs^2
  dtau <- diff(c(0, tau))
  z2 <- rnorm(length(t_obs))
  dx_nhom <- a_nhom * dtau + sigma_true * sqrt(dtau) * z2
  x_nhom <- c(0, cumsum(dx_nhom))
  return(cbind(t = c(0, t_obs), x = x_nhom))
}