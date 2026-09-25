library(future.apply)

MC_test <- function(N = 50, scenario, t_obs, a = 2, b = 1, alpha = NULL,
                    ncores = parallel::detectCores() - 1) {
  
  plan(multisession, workers = ncores)
  
  une_iter <- function(i, scenario, t_obs, a, b, alpha) {
    if (!is.null(alpha)) {
      X <- scenario(t_obs = t_obs, a = a, b = b, alpha = alpha)
    } else {
      X <- scenario(t_obs = t_obs, a = a, b = b)
    }
    
    res_h <- gof.gp(X, "Henze")
    rej_h <- as.integer(res_h$T.value > res_h$cv)
    
    res_v <- gof.gp(X, "VGE")
    rej_v <- as.integer(res_v$p.value < 0.05)
    
    res_k <- gof.gp(X, "KS")
    rej_k <- as.integer(res_k$T.value > res_k$cv)
    
    res_b <- gof.gp(X, "BET")
    rej_b <- as.integer(res_b$T.value > res_b$cv)
    
    c(Henze = rej_h, VGE = rej_v, KS = rej_k, BET = rej_b)
  }
  
  res <- future_sapply(1:N, une_iter,
                       scenario = scenario, t_obs = t_obs,
                       a = a, b = b, alpha = alpha,
                       future.seed = TRUE)
  taux <- rowMeans(res)
  names(taux) <- c("Henze", "VGE", "KS", "BET")
  return(taux)
}