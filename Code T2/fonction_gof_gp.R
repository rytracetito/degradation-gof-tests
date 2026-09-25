# Implémentation de la fonction gp.gof


############################################################
#  Fonction principale : gof.gp(X , test)
#  X = matrice contenant (temps ti , observation xi)
#  test = "Henze" , "VGE" ou "KS"
############################################################


############################################################
#  Fonction : gamma_bridge
#  Interpolation d'un processus gamma sur une grille régulière
############################################################

gamma_bridge <- function(t, x, a, dt_eq){
  
  # Vérification : les temps doivent être strictement croissants
  if(any(diff(t) <= 0))
    stop("Les temps doivent être strictement croissants")
  
  # Vérification du paramètre du processus gamma
  if(a <= 0)
    stop("Le paramètre a doit être positif")
  
  # Nombre de points dans la grille régulière
  Np <- floor(max(t) / dt_eq)
  
  # Construction de la grille temporelle régulière
  t_eq <- seq(dt_eq , Np*dt_eq , by = dt_eq)
  
  # Vecteur qui contiendra les valeurs interpolées
  x_eq <- numeric(Np)
  
  # Boucle d'interpolation
  for(i in 1:Np){
    
    ti <- t_eq[i]
    
    # Premier temps supérieur à ti
    sup_i <- which(t >= ti)[1]
    
    t_sup <- t[sup_i]
    x_sup <- x[sup_i]
    
    # Détermination du point inférieur
    if(sup_i == 1){
      t_inf <- 0
      x_inf <- 0
    }else{
      
      t_inf_candidate <- t[sup_i-1]
      
      if(t_inf_candidate > (i-1)*dt_eq){
        t_inf <- t_inf_candidate
        x_inf <- x[sup_i-1]
      }else{
        t_inf <- (i-1)*dt_eq
        
        if(i > 1){
          x_inf <- x_eq[i-1]
        }else{
          x_inf <- 0
        }
      }
    }
    
    # Paramètres de la loi Beta
    alpha1 <- a * max(ti - t_inf , 1e-6)
    alpha2 <- a * max(t_sup - ti , 1e-6)
    
    # Tirage beta pour interpolation gamma bridge
    y <- rbeta(1 , alpha1 , alpha2)
    
    # Valeur interpolée
    x_eq[i] <- x_inf + (x_sup - x_inf) * y
  }
  
  return(data.frame(time = t_eq , value = x_eq))
}


tri_sample <- function(dt , dx , dt_eq){
  
  N <- length(dt)
  
  dt_c <- dt
  dx_c <- dx
  
  dt_s <- numeric(N)
  dx_s <- numeric(N)
  
  # On sélectionne d'abord le plus grand incrément de temps
  max_i <- which.max(dt_c)
  
  dt_s[1] <- dt_c[max_i]
  dx_s[1] <- dx_c[max_i]
  
  dt_c <- dt_c[-max_i]
  dx_c <- dx_c[-max_i]
  
  # Construction progressive du nouvel ordre
  for(i in 2:N){
    
    cum_time <- sum(dt_s[1:(i-1)])
    
    diffs <- (cum_time + dt_c) - i*dt_eq
    
    valid <- which(diffs >= 0)
    
    if(length(valid) > 0){
      choix <- valid[which.min(diffs[valid])]
    }else{
      choix <- which.max(diffs)
    }
    
    dt_s[i] <- dt_c[choix]
    dx_s[i] <- dx_c[choix]
    
    dt_c <- dt_c[-choix]
    dx_c <- dx_c[-choix]
  }
  
  return(list(dt = dt_s , dx = dx_s))
}



############################################################
# Fonction : estimer_a
#  Estimation du paramètre du processus gamma
############################################################

# Métode 1: Méthode de Villaseñor-González

estimer_params_VG <- function(dx) {
  # Supprimer les NA
  dx <- dx[!is.na(dx) & is.finite(dx)]
  
  if(length(dx) < 2) {
    return(list(shape = 1, rate = 1))
  }
  
  m <- mean(dx)
  log_dx <- log(dx)
  ml <- mean(log_dx)
  b_hat <- mean((dx - m) * (log_dx - ml))
  
  if(is.na(b_hat) || b_hat <= 0) b_hat <- 1e-6
  
  a_hat <- m / b_hat
  if(is.na(a_hat) || a_hat <= 0) a_hat <- 1e-6
  
  list(shape = a_hat, rate = 1 / b_hat)
}

# Méthode 2: Méthode des moments 
estimer_a <- function(t , x){
  
  dt <- diff(t)
  dx <- diff(x)
  
  # Estimateur de la moyenne
  m_hat <- sum(dx) / sum(dt)
  
  # Estimateur de la variance
  n_hat <- sum((dx - m_hat*dt)^2) /
    (sum(dt) - sum(dt^2)/sum(dt))
  
  # Estimateur du paramètre a
  a_hat <- m_hat^2 / n_hat
  
  return(a_hat)
}


############################################################
#  Test d'adéquation de Henze
############################################################

gof.gp.henze <- function(dx){
  
  if(!requireNamespace("gofgamma", quietly = TRUE))
    install.packages("gofgamma")
  
  library(gofgamma)
  
  res <- test.HME1(dx , a = 1 , boot = 500
                   , alpha = 0.05)
  
  return(res)
}



############################################################
#  Test VGE 
############################################################
gof.gp.vge <- function(dx) {
  if (!requireNamespace("goft", quietly = TRUE))
    install.packages("goft")
  library(goft)
  gamma_test(dx)
}
############################################################
#  Test de Kolmogorov-Smirnov 
############################################################

gof.gp.ks <- function(dx, B = 100) {
  params <- estimer_params_VG(dx)
  shape_hat <- params$shape
  rate_hat <- params$rate
  
  D_obs <- suppressWarnings(
    ks.test(dx, "pgamma", shape = shape_hat, rate = rate_hat)$statistic
  )
  
  D_boot <- replicate(B, {
    dx_boot <- rgamma(length(dx), shape = shape_hat, rate = rate_hat)
    tryCatch({
      params_b <- estimer_params_VG(dx_boot)
      ks.test(dx_boot, "pgamma", shape = params_b$shape, rate = params_b$rate)$statistic
    }, error = function(e) NA)
  })
  
  p_val <- mean(D_boot >= D_obs, na.rm = TRUE)
  
  structure(list(
    statistic = D_obs,
    p.value = p_val,
    method = "Kolmogorov-Smirnov test with parametric bootstrap (VG)",
    data.name = deparse(substitute(dx))
  ), class = "htest")
}



############################################################
#  Test de Steffen Betsch & Bruno Ebner 
############################################################


gof.gp.betsch <- function(dx, B = 500) {
  
  if(!requireNamespace("gofgamma", quietly = TRUE))
    install.packages("gofgamma")
  
  library(gofgamma)
  
  # Test de Betsch & Ebner
  res <- test.BE(dx, boot = B)
  
  return(res)
}



############################################################
#  Fonction principale : gof.gp
############################################################

gof.gp <- function(X , test){
  
  # Extraction du temps et des observations
  t <- X[,1]
  x <- X[,2]
  
  dt <- diff(t)
  dx <- diff(x)
  
  tol <- 1e-8
  
  ##########################################################
  # Cas 1 : observations à temps réguliers
  ##########################################################
  
  if(all(abs(dt - dt[1]) < tol)){
    
    dx_test <- dx
    
  }else{
    
    ##########################################################
    # Cas 2 : observations à temps irréguliers
    ##########################################################
    
    dt_eq <- max(t)/length(dt)
    
    tri <- tri_sample(dt , dx , dt_eq)
    
    t_tri <- c(0 , cumsum(tri$dt))
    x_tri <- c(0 , cumsum(tri$dx))
    
    # Estimation du paramètre du processus
    a_ref <- estimer_a(t_tri , x_tri)
    
    # Interpolation gamma bridge
    interp <- gamma_bridge(t_tri , x_tri , a_ref , dt_eq)
    
    dx_test <- diff(c(0 , interp$value))
  }
  
  ##########################################################
  # Sélection du test
  ##########################################################
  
  
  if(test == "Henze")
    return(gof.gp.henze(dx_test))
  
  if(test == "VGE")
    return(gof.gp.vge(dx_test))
  
  if(test == "KS")
    return(gof.gp.ks(dx_test))
  
  if(test == "BET")
    return(gof.gp.betsch(dx_test))
  
  stop("Test inconnu")
}

