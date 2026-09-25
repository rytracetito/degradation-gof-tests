estimer_wiener <- function(t, x) {
  if (length(t) != length(x)) stop("t et x doivent avoir la même longueur")
  if (length(t) < 2) stop("Il faut au moins 2 observations")
  if (any(diff(t) <= 0)) stop("Les temps doivent être strictement croissants")
  
  dt <- diff(t)
  dx <- diff(x)
  N <- length(dt)
  
  mu_hat <- sum(dx) / sum(dt)
  
  residus <- (dx - mu_hat * dt)^2 / dt
  sigma2_hat <- sum(residus) / N
  sigma_hat <- sqrt(sigma2_hat)
  
  c(mu = mu_hat, sigma2 = sigma2_hat, sigma = sigma_hat)
}


#  Fonction de tri
tri_sample <- function(dt, dx, dt_eq) {
  N <- length(dt)
  dt_c <- dt
  dx_c <- dx
  dt_s <- numeric(N)
  dx_s <- numeric(N)
  max_i <- which.max(dt_c)
  dt_s[1] <- dt_c[max_i]
  dx_s[1] <- dx_c[max_i]
  dt_c <- dt_c[-max_i]
  dx_c <- dx_c[-max_i]
  for (i in 2:N) {
    cum_time <- sum(dt_s[1:(i - 1)])
    diffs <- (cum_time + dt_c) - i * dt_eq
    valid <- which(diffs >= 0)
    if (length(valid) > 0) {
      choix <- valid[which.min(diffs[valid])]
    } else {
      choix <- which.max(diffs)
    }
    dt_s[i] <- dt_c[choix]
    dx_s[i] <- dx_c[choix]
    dt_c <- dt_c[-choix]
    dx_c <- dx_c[-choix]
  }
  list(dt = dt_s, dx = dx_s)
}