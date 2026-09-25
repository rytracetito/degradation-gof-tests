# Cas reguliers
MC_test_wiener <- function(N = 50, scenario, tmax = 100, n = 100,
                           mu = NULL, sigma = NULL, a = NULL, b = NULL) {
  
  une_iter <- function(i) {
    if (!is.null(a) && !is.null(b)) {
      X <- scenario(tmax = tmax, n = n, a = a, b = b)
    } else {
      X <- scenario(tmax = tmax, n = n, mu = mu, sigma = sigma)
    }
    
    res_s <- gof.wp(X, "Shapiro")
    rej_s <- as.integer(res_s$p.value < 0.05)
    
    res_l <- gof.wp(X, "Lilliefors")
    rej_l <- as.integer(res_l$p.value < 0.05)
    
    res_g <- gof.wp(X, "Agostino")
    rej_g <- as.integer(res_g$p.value < 0.05)
    
    res_a <- gof.wp(X, "Anderson")
    rej_a <- as.integer(res_a$p.value < 0.05)
    
    res_c <- gof.wp(X, "Cramer")
    rej_c <- as.integer(res_c$p.value < 0.05)
    
    c(Shapiro = rej_s, Lilliefors = rej_l, Agostino = rej_g,
      Anderson = rej_a, Cramer = rej_c)
  }
  
  res <- sapply(1:N, une_iter)
  taux <- rowMeans(res)
  return(taux)
}

# Cas irreguliers
MC_test_wiener_irreg <- function(N = 50, scenario, t_obs,
                                 mu = NULL, sigma = NULL, a = NULL, b = NULL) {
  
  une_iter <- function(i) {
    if (!is.null(a) && !is.null(b)) {
      X <- scenario(t_obs = t_obs, a = a, b = b)
    } else {
      X <- scenario(t_obs = t_obs, mu = mu, sigma = sigma)
    }
    
    res_s <- gof.wp(X, "Shapiro")
    rej_s <- as.integer(res_s$p.value < 0.05)
    
    res_l <- gof.wp(X, "Lilliefors")
    rej_l <- as.integer(res_l$p.value < 0.05)
    
    res_g <- gof.wp(X, "Agostino")
    rej_g <- as.integer(res_g$p.value < 0.05)
    
    res_a <- gof.wp(X, "Anderson")
    rej_a <- as.integer(res_a$p.value < 0.05)
    
    res_c <- gof.wp(X, "Cramer")
    rej_c <- as.integer(res_c$p.value < 0.05)
    
    c(Shapiro = rej_s, Lilliefors = rej_l, Agostino = rej_g,
      Anderson = rej_a, Cramer = rej_c)
  }
  
  res <- sapply(1:N, une_iter)
  taux <- rowMeans(res)
  return(taux)
}