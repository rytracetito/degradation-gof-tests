
library(gofgamma)
library(goft)

 
K = 100
# --- Test de Henze-Meintanis-Ebner (HME1) ---
gof.gp.henze <- function(dx, a = 1, boot = K, alpha = 0.05) {
  test.HME1(dx, a = a, boot = boot, alpha = alpha)
}

# --- Test de Villaseñor-González-Estrada (VGE) ---


gof.gp.vge <- function(dx, B = K) {
  if (any(dx <= 0)) stop("dx doit contenir des valeurs positives")
  
  fit   <- gamma_fit(dx)
  shape <- fit["shape", 1]   
  scale <- fit["scale", 1]   
  
  T_obs <- abs(as.numeric(suppressWarnings(gamma_test(dx))$statistic))
  
  T_boot <- replicate(B, {
    dx_star <- rgamma(length(dx), shape = shape, scale = scale)
    abs(as.numeric(suppressWarnings(gamma_test(dx_star))$statistic))
  })
  
  p_val <- mean(T_boot >= T_obs, na.rm = TRUE)
  
  list(statistic = T_obs, p.value = p_val, shape = shape, scale = scale)
}

# --- Test de Kolmogorov-Smirnov 
gof.gp.ks <- function(dx, boot = K, alpha = 0.05) {
  test.KS(dx, boot = boot, alpha = alpha)
}

# --- Test de Betsch-Ebner ---
gof.gp.betsch <- function(dx, B = K) {
  test.BE(dx, boot = B)
}