# Définition des fonctions


# Fonction : gamma_bridge (interpolation par pont gamma)
gamma_bridge <- function(X, t, a) {
  n <- length(X) - 1
  t_n <- t[length(t)]
  teq <- t_n / n
  
  t_star <- seq(0, t_n, teq)
  
  Y <- numeric(length(t_star))
  Y[1] <- 0
  
  for (i in 2:n) {
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
  Y[n + 1] <- X[n + 1]
  
  list(teq = teq, new_obs = Y)
}