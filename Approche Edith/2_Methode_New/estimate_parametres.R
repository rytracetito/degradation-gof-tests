# Méthode des moments 
estimer <- function(t, x) {
  dt <- diff(t)
  dx <- diff(x)
  
  mu_hat <- sum(dx) / sum(dt)
  if (mu_hat <= 0) {
    stop("mu_hat <= 0 : tendance de dégradation non positive, modèle Gamma inadapté")
  }
  
  denom <- sum(dt) - sum(dt^2) / sum(dt)
  if (denom <= 0) {
    stop("Dénominateur nul ou négatif : données inadaptées (un seul intervalle ?)")
  }
  
  nu_hat <- sum((dx - mu_hat * dt)^2) / denom
  if (nu_hat <= 0) {
    stop("nu_hat <= 0 : estimation dégénérée, vérifier les données")
  }
  
  a_hat <- mu_hat^2 / nu_hat
  b_hat <- mu_hat / nu_hat
  
  c(a = a_hat, b = b_hat)
}
