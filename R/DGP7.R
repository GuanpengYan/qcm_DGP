dataGeneratingProcess7 <- function(N, T0, T1){
  K <- 2
  beta <- c(1, 4)
  T <- T0 + T1
  gamma <- array(dim = c(3, N), data = rnorm(N * 3, mean = 0, sd = 1))
  rho <- array(dim = c(K, N), data = runif(N * K, min = 0.1, max = 0.9))
  c <- array(dim = c(3, N), data = runif(N * 2, min = 1, max = 2))
  eta <- array(dim = c(K, N, T), data = rnorm(N * T * K, mean = 0, sd = 1))
  zeta <- array(dim = c(3, T), data = rnorm(n = 3 * T, mean = 0, sd = 1))
  f <- array(dim = c(3, T))
  for(t in 1:T){
    if(t == 1){
      f[2, t] <- zeta[2, t]
      f[1, t] <- 0.8 * f[2, t] + zeta[1, t]
      f[3, t] <- zeta[3, t]
    }else if(t == 2){
      f[2, t] <- 0.5 * f[2, t - 1] + zeta[2, t] + 0.7 * zeta[2, t - 1]
      f[1, t] <- 0.8 * f[2, t] + zeta[1, t]
      f[3, t] <- zeta[3, t] + 0.9 * zeta[3, t - 1]
    }else{
      f[2, t] <- 0.5 * f[2, t-1] + zeta[2, t] + 0.7 * zeta[2, t - 1]
      f[1, t] <- 0.8 * f[2, t] + zeta[1, t]
      f[3, t] <- zeta[3, t] + 0.9 * zeta[3, t - 1] + 0.4 *zeta[3, t - 2]
    }
  }
  u <- array(dim = c(N, T), data = rnorm(N * T, mean = 0, sd = 1))
  
  #Calculate x
  x <- array(dim = c(K, N, T))
  for(t in 1:T) for(i in 1:N) for(k in 1:K) x[k, i, t] <- ifelse(
    t == 1,
    1 + sum(c[, i] * f[, t]) + eta[k, i, t],
    1 + rho[k, i] * x[k, i, t - 1] + c[1, i] * gamma[k, i] + c[2, i] * f[k, t] + eta[k, i, t])
  #Calculate y.ctfl
  y.ctfl <- array(dim = c(N, T))
  for(t in 1:T) for(i in 1:N) 
    y.ctfl[i, t] <- sum(x[, i, t] * beta) + sum(gamma[, i] * f[, t]) + u[i, t]
  #Calculate y.actl
  y.actl <- y.ctfl; for(t in (T0 + 1):T) y.actl[1, t] <- y.actl[1, t] + 1
  
  return(list(
    N = N,
    T0 = T0,
    T1 = T1,
    beta = beta,
    K = K,
    gamma = gamma,
    rho = rho,
    c = c,
    eta = eta,
    zeta = zeta,
    f = f,
    u = u,
    x = x,
    y.ctfl = y.ctfl,
    y.actl = y.actl
  ))
}
# source("R/.Source.R")
# 
# dataGeneratingProcess7(40, 40, 10)$y.actl
# save(DGP1_data,file = "DGP1_data.Rdata")
