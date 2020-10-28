dataGeneratingProcess2 <- function(N, T0, T1){
  K <- 2
  beta <- c(1, 4)
  T <- T0 + T1
  gamma <- array(dim = c(3, N), data = rnorm(N * 3, mean = 0, sd = 1))
  rho <- array(dim = c(K, N), data = runif(N * K, min = 0.1, max = 0.9))
  c <- array(dim = c(3, N), data = runif(N * 2, min = 1, max = 2))
  eta <- array(dim = c(K, N, T), data = rnorm(N * T * K, mean = 0, sd = 1))
  f <- array(dim = c(3, T), data = rnorm(n = 3 * T, mean = 0, sd = 1))
  u <- array(dim = c(N, T), data = rnorm(N * T, mean = 0, sd = 1))
  phi <- array(dim = c(K, N), data = runif(N * K, min = 0.1, max = 0.9))
  #Calculate x
  x <- array(dim = c(K, N, T))
  for(t in 1:T) for(i in 1:N) for(k in 1:K) x[k, i, t] <- ifelse(
    t == 1,
    1 + eta[k, i, t],
    1 + rho[k, i] * x[k, i, t - 1] + eta[k, i, t] + phi[k,i] * eta[k, i, t])
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
    f = f,
    u = u,
    phi = phi,
    x = x,
    y.ctfl = y.ctfl,
    y.actl = y.actl
  ))
}
# source("R/.Source.R")
# 
# dataGeneratingProcess1(40, 40, 10)
# save(DGP1_data,file = "DGP1_data.Rdata")
