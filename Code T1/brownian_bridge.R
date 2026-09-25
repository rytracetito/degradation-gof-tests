
############################################################
#  Adaptation de la méthodologie proposée par Grall et Maës 
# au processus de Wiener avec dérive
############################################################

# Cas régulier 

############################################
# Simulation du processus de Wiener
############################################

simulate_wiener <- function(a, sigma, delta_t, T){
  
  time <- seq(0, T, by = delta_t)
  n <- length(time)
  
  X <- numeric(n)
  X[1] <- 0
  
  for(i in 2:n){
    
    X[i] <- X[i-1] +
      a * delta_t +
      sigma * sqrt(delta_t) * rnorm(1)
  }
  
  return(data.frame(time = time, value = X))
}

############################################
# Calcul des accroissements
############################################

compute_increments <- function(time, values){
  
  delta_x <- diff(values)
  delta_t <- diff(time)
  
  return(data.frame(delta_t = delta_t, delta_x = delta_x))
}

############################################
# Estimation des paramètres
############################################

estimate_parameters <- function(delta_t, delta_x){
  
  a_hat <- sum(delta_x) / sum(delta_t)
  
  sigma_hat <- sqrt(mean((delta_x - a_hat * delta_t)^2 / delta_t))
  
  return(list(a = a_hat, sigma = sigma_hat))
}

############################################
# Standardisation des accroissements
############################################

standardize_increments <- function(delta_x, delta_t, a, sigma){
  
  Z <- (delta_x - a * delta_t) / (sigma * sqrt(delta_t))
  
  return(Z)
}

############################################
# Tests de normalité
############################################

normality_tests <- function(Z){
  
  library(tseries)
  library(nortest)
  
  shapiro <- shapiro.test(Z)
  jb <- jarque.bera.test(Z)
  ad <- ad.test(Z)
  
  return(list(
    Shapiro_Wilk = shapiro,
    Jarque_Bera = jb,
    Anderson_Darling = ad
  ))
}

############################################
# Analyse complète
############################################

final_analysis <- function(time, values){
  
  inc <- compute_increments(time, values)
  
  params <- estimate_parameters(inc$delta_t, inc$delta_x)
  
  Z <- standardize_increments(
    inc$delta_x,
    inc$delta_t,
    params$a,
    params$sigma
  )
  
  tests <- normality_tests(Z)
  
  return(list(
    increments = inc,
    parameters = params,
    Z = Z,
    tests = tests
  ))
}

############################################
# APPLICATION 1
############################################

# paramètres du modèle

a <- 0.3
sigma <- 0.5
delta_t <- 0.1
T <- 2

# simulation

sim <- simulate_wiener(a, sigma, delta_t, T)

time <- sim$time
values <- sim$value

############################################
# Visualisation
############################################

plot(time, values,
     type="l",
     lwd=2,
     col="blue",
     main="Trajectoire d'un processus de Wiener",
     xlab="Temps",
     ylab="X(t)")

############################################
# Analyse statistique
############################################

result <- final_analysis(time, values)

# paramètres estimés

print(result$parameters)

# résultats des tests

print(result$tests)

############################################
# Histogramme
############################################

hist(result$Z,
     probability=TRUE,
     main="Histogramme des Z")

curve(dnorm(x), add=TRUE, lwd=2)

############################################
# QQ plot
############################################

qqnorm(result$Z)
qqline(result$Z, col="red")



############################################
# APPLICATION 2 : PROCESSUS GAMMA
############################################

alpha <- 3
beta <- 0.4
delta_t <- 0.1
T <- 2

time_g <- seq(0, T, by = delta_t)
n_g <- length(time_g)

values_g <- numeric(n_g)
values_g[1] <- 0

for(i in 2:n_g){
  
  increment <- rgamma(1,
                      shape = alpha * delta_t,
                      scale = beta)
  
  values_g[i] <- values_g[i-1] + increment
}

############################################
# Visualisation
############################################

plot(time_g, values_g,
     type="l",
     col="darkgreen",
     lwd=2,
     main="Trajectoire du processus Gamma",
     xlab="Temps",
     ylab="G(t)")

############################################
# Analyse 
############################################

result_gamma <- final_analysis(time_g, values_g)

# paramètres estimés

print(result_gamma$parameters)

# résultats des tests

print(result_gamma$tests)

############################################
# Histogramme
############################################

hist(result_gamma$Z,
     probability=TRUE,
     main="Histogramme des Z (processus Gamma)")

curve(dnorm(x), add=TRUE, lwd=2)

############################################
# QQ plot
############################################

qqnorm(result_gamma$Z)
qqline(result_gamma$Z, col="red")













# Cas irrégulier 


############################################
# SORT ALGORITHM 
############################################

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
  
  return(list(dt = dt_s, dx = dx_s))
}

############################################
# BROWNIAN BRIDGE INTERPOLATION
############################################

brownian_bridge_interpolation <- function(t_obs, x_obs, sigma, delta_t){
  
  t_new <- seq(min(t_obs), max(t_obs), by = delta_t)
  x_new <- numeric(length(t_new))
  
  # imposer les observations
  for(i in 1:length(t_obs)){
    idx <- which(abs(t_new - t_obs[i]) < 1e-10)
    if(length(idx)>0){
      x_new[idx] <- x_obs[i]
    }
  }
  
  for(i in 1:(length(t_obs)-1)){
    
    tA <- t_obs[i]
    tB <- t_obs[i+1]
    
    xA <- x_obs[i]
    xB <- x_obs[i+1]
    
    idx <- which(t_new > tA & t_new < tB)
    
    for(j in idx){
      
      t <- t_new[j]
      
      mu <- xA + (t - tA)/(tB - tA)*(xB - xA)
      
      v <- sigma^2*(t - tA)*(tB - t)/(tB - tA)
      
      x_new[j] <- mu + sqrt(v)*rnorm(1)
    }
  }
  
  return(data.frame(time=t_new,value=x_new))
}

############################################
# CALCUL DES ACCROISSEMENTS
############################################

compute_increments <- function(time, values){
  
  delta_x <- diff(values)
  delta_t <- diff(time)
  
  return(data.frame(delta_t = delta_t, delta_x = delta_x))
}

############################################
# ESTIMATION DES PARAMÈTRES
############################################

estimate_parameters <- function(delta_t, delta_x){
  
  a_hat <- sum(delta_x) / sum(delta_t)
  
  sigma_hat <- sqrt(mean((delta_x - a_hat * delta_t)^2 / delta_t))
  
  return(list(a = a_hat, sigma = sigma_hat))
}

############################################
# STANDARDISATION
############################################

standardize_increments <- function(delta_x, delta_t, a, sigma){
  
  Z <- (delta_x - a * delta_t) / (sigma * sqrt(delta_t))
  
  return(Z)
}

############################################
# TESTS STATISTIQUES
############################################

normality_tests <- function(Z){
  
  library(tseries)
  library(nortest)
  
  shapiro <- shapiro.test(Z)
  jb <- jarque.bera.test(Z)
  ad <- ad.test(Z)
  
  return(list(
    Shapiro_Wilk = shapiro,
    Jarque_Bera = jb,
    Anderson_Darling = ad
  ))
}

############################################
# ANALYSE COMPLÈTE
############################################

final_analysis <- function(time, values){
  
  inc <- compute_increments(time, values)
  
  params <- estimate_parameters(inc$delta_t, inc$delta_x)
  
  Z <- standardize_increments(
    inc$delta_x,
    inc$delta_t,
    params$a,
    params$sigma
  )
  
  tests <- normality_tests(Z)
  
  return(list(
    increments = inc,
    parameters = params,
    Z = Z,
    tests = tests
  ))
}

############################################
# APPLICATION
############################################

# données observées (temps irréguliers)

t_obs <- c(0,0.9,7.7,1.5,7)
x_obs <- c(0,0.2,0.5,0.3,0.9)

sigma <- 0.5
delta_t <- 0.1

############################################
# CALCUL DES INCRÉMENTS OBSERVÉS
############################################

dt <- diff(t_obs)
dx <- diff(x_obs)

dt_eq <- mean(dt)

############################################
# SORT ALGORITHM
############################################

tri <- tri_sample(dt, dx, dt_eq)

dt_sorted <- tri$dt
dx_sorted <- tri$dx

############################################
# RECONSTRUCTION DES DONNÉES TRIÉES
############################################

t_sorted <- c(t_obs[1], t_obs[1] + cumsum(dt_sorted))
x_sorted <- c(x_obs[1], x_obs[1] + cumsum(dx_sorted))

############################################
# BROWNIAN BRIDGE
############################################

sim <- brownian_bridge_interpolation(t_sorted, x_sorted, sigma, delta_t)

############################################
# VISUALISATION
############################################

plot(sim$time, sim$value,
     type="l",
     col="blue",
     lwd=2,
     main="Interpolation Brownian Bridge",
     xlab="t",
     ylab="X(t)")

points(t_sorted, x_sorted, col="red", pch=19)

############################################
# ANALYSE STATISTIQUE
############################################

result <- final_analysis(sim$time, sim$value)

print(result$parameters)

print(result$tests)

############################################
# HISTOGRAMME
############################################

hist(result$Z,
     probability = TRUE,
     main="Histogramme des Z")

curve(dnorm(x), add=TRUE, lwd=2)

############################################
# QQ PLOT
############################################

qqnorm(result$Z)
qqline(result$Z, col="red")



# Test de KS

ks_normal_test <- function(Z){
  
  ks <- ks.test(Z, "pnorm", mean = 0, sd = 1)
  
  return(ks)
}

result <- final_analysis(sim$time, sim$value)

ks_result <- ks_normal_test(result$Z)

print(ks_result)
