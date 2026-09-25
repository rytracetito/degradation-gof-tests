




# Estimation du paramètre a par la méthode de Villaseñor-González 
estimer_a_VG <- function(dx) {
  dx <- dx[!is.na(dx) & is.finite(dx)]
  
  if (length(dx) < 2) {
    return(c(shape = 1, rate = 1))
  }
  
  m <- mean(dx)
  log_dx <- log(dx)
  ml <- mean(log_dx)
  b_hat <- mean((dx - m) * (log_dx - ml))
  
  if (is.na(b_hat) || b_hat <= 0) b_hat <- 1e-6
  
  a_hat <- m / b_hat
  if (is.na(a_hat) || a_hat <= 0) a_hat <- 1e-6
  
  c(shape = a_hat, rate = 1 / b_hat)
}

# Fonction gamma_bridge 
gamma_bridge <- function(t, x, a, dt_eq){
  
  if(any(diff(t) <= 0))
    stop("Les temps doivent être strictement croissants")
  
  if(a <= 0)
    stop("Le paramètre a doit être positif")
  
  Np <- floor(max(t) / dt_eq)
  t_eq <- seq(dt_eq, Np*dt_eq, by = dt_eq)
  x_eq <- numeric(Np)
  
  for(i in 1:Np){
    
    ti <- t_eq[i]
    sup_i <- which(t >= ti)[1]
    t_sup <- t[sup_i]
    x_sup <- x[sup_i]
    
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
    
    alpha1 <- a * max(ti - t_inf, 1e-6)
    alpha2 <- a * max(t_sup - ti, 1e-6)
    y <- rbeta(1, alpha1, alpha2)
    x_eq[i] <- x_inf + (x_sup - x_inf) * y
  }
  
  return(data.frame(time = t_eq, value = x_eq))
}

# Les différents tests
############################################################
# Test d'adéquation de Henze
############################################################
gof.gp.henze <- function(dx) {
  if (!requireNamespace("gofgamma", quietly = TRUE)) install.packages("gofgamma")
  library(gofgamma)
  test.HME1(dx, a = 1, boot = 500, alpha = 0.05)
}

############################################################
# Test VGE
############################################################
gof.gp.vge <- function(dx) {
  if (!requireNamespace("goft", quietly = TRUE)) install.packages("goft")
  library(goft)
  gamma_test(dx)
}

############################################################
# Test de Kolmogorov-Smirnov
############################################################
gof.gp.ks <- function(dx, B = 100) {
  params <- estimer_a_VG(dx)
  shape_hat <- params["shape"]
  rate_hat <- params["rate"]
  
  D_obs <- ks.test(dx, "pgamma", shape = shape_hat, rate = rate_hat)$statistic
  
  D_boot <- replicate(B, {
    dx_boot <- rgamma(length(dx), shape = shape_hat, rate = rate_hat)
    params_b <- estimer_a_VG(dx_boot)
    ks.test(dx_boot, "pgamma", shape = params_b["shape"], rate = params_b["rate"])$statistic
  })
  
  p_val <- mean(D_boot >= D_obs)
  
  structure(list(
    statistic = D_obs,
    p.value = p_val,
    method = "Kolmogorov-Smirnov with parametric bootstrap",
    data.name = deparse(substitute(dx))
  ), class = "htest")
}

############################################################
# Test de Steffen Betsch & Bruno Ebner 
############################################################
gof.gp.betsch <- function(dx, B = 500) {
  if(!requireNamespace("gofgamma", quietly = TRUE)) install.packages("gofgamma")
  library(gofgamma)
  test.BE(dx, boot = B)
}

############################################################
# Fonction principale : gof.gp
############################################################
gof.gp <- function(X, test) {
  
  t <- X[, 1]
  x <- X[, 2]
  
  dt <- diff(t)
  dx <- diff(x)
  
  tol <- 1e-8
  
  # Cas 1 : observations à temps réguliers
  if (all(abs(dt - dt[1]) < tol)) {
    dx_test <- dx
  } 
  # Cas 2 : observations à temps irréguliers
  else {
    dt_eq <- max(t) / length(dt)
    params <- estimer_a_VG(dx)
    a_ref <- params["shape"]
    interp <- gamma_bridge(t, x, a_ref, dt_eq)
    dx_test <- diff(c(0, interp$value))
  }
  
  if (test == "Henze") return(gof.gp.henze(dx_test))
  if (test == "VGE")   return(gof.gp.vge(dx_test))
  if (test == "KS")    return(gof.gp.ks(dx_test))
  if (test == "BET")   return(gof.gp.betsch(dx_test))
  stop("Test inconnu")
}



# Simulations et Application

# Cas Homogène 

# Paramètres
set.seed(123)
a <- 2
b <- 1
tmax <- 100
n <- 100

# Génération d'un processus gamma homogène avec temps réguliers
t_reg <- seq(0, tmax, length.out = n)
dt_reg <- diff(t_reg)
dx_reg <- rgamma(n-1, shape = a * dt_reg, rate = b)
x_reg <- c(0, cumsum(dx_reg))
X_reg <- cbind(t_reg, x_reg)


# Test 

gof.gp(X_reg, "KS")
gof.gp(X_reg, "Henze")
gof.gp(X_reg, "VGE")
gof.gp(X_reg, "BET")


# Cas non Homogène 

# Cas 1: alpha = 1.5

set.seed(123)
a <- 2
b <- 1
tmax <- 100
n <- 100


alpha_nh <- 1.5
t <- seq(0.01, tmax, length.out = n)
dt <- diff(t^alpha_nh)
dx <- rgamma(n-1, shape = a * dt, rate = b)
x <- c(0, cumsum(dx))
X_nh <- cbind(t, x)

# Test
gof.gp(X_nh, "KS")
gof.gp(X_nh, "Henze")
gof.gp(X_nh, "VGE")
gof.gp(X_nh, "BET")

# Cas 2 :  alpha = 0.8

set.seed(123)
a <- 2
b <- 1
tmax <- 100
n <- 1000

# Génération d'un processus gamma NON homogène 

alpha_nh <- 0.03
t <- seq(0.01, tmax, length.out = n)
dt <- diff(t^alpha_nh)
dx <- rgamma(n-1, shape = a * dt, rate = b)
x <- c(0, cumsum(dx))
X_nh <- cbind(t, x)

# Test
gof.gp(X_nh, "KS")
gof.gp(X_nh, "Henze")
gof.gp(X_nh, "VGE")
gof.gp(X_nh, "BET")
