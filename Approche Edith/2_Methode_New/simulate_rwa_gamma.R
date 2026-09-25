

# Simuler un processus gamma par la méthode Random Walk Approximation (RWA)
simulate_gamma_rwa <- function(times, A = NULL, a = NULL, nsim = 1, b = 1) {
  if (!is.numeric(times) || !is.vector(times)) {
    stop("'times' doit être un vecteur numérique.")
  }
  if (length(times) < 2) {
    stop("'times' doit contenir au moins deux points.")
  }
  if (any(diff(times) <= 0)) {
    stop("'times' doit être strictement croissant.")
  }
  if (!is.numeric(nsim) || nsim < 1 || nsim != floor(nsim)) {
    stop("'nsim' doit être un entier positif.")
  }
  
  if (!is.null(A) && is.function(A)) {
    A_vals <- A(times)
  } else if (!is.null(a) && is.numeric(a) && length(a) == 1) {
    A_vals <- a * times
  } else {
    stop("Il faut fournir soit une fonction 'A', soit un scalaire 'a'.")
  }
  
  if (any(diff(A_vals) <= 0)) {
    stop("La fonction A doit être strictement croissante.")
  }
  
  diff_A <- diff(A_vals)
  n_intervals <- length(diff_A)
  
  increments <- matrix(rgamma(n_intervals * nsim, shape = rep(diff_A, each = nsim), rate = b),
                       nrow = n_intervals, ncol = nsim, byrow = TRUE)
  
  X <- rbind(rep(0, nsim), apply(increments, 2, cumsum))
  
  rownames(X) <- paste0("t", seq_along(times) - 1)
  colnames(X) <- paste0("sim", seq_len(nsim))
  
  return(X)
}

