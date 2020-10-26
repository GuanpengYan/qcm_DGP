### Confidence Intervals via Quantile Random Forest for DGP 1 in Hsiao et al. (2018,p.8)
# N = 10,20,30,40,50,60 (number of individuals), only the first individual is treated in the posttreatment period
# T0 = 10,20,30,40,50,60,70,80,90 (pretreatment period)

# T1 = 1 (posttreatment period)
# K = 2 (number of factors)
# reps = 1000 (number of simulations)

# Define a function ci.qcm of T0 to compute the coverage rate for T1 posttreatment periods
# Other parameters and DGP are all fixed 
# ci.qcm refers to "confidence interval for quantile control method"

DGP1 <- function(N, T0){
  T1 <- 10
  K <- 2
  beta <- c(1, 2)
  T <- T0 + T1
  gamma <- array(dim = c(3, N), data = rnorm(N * 3, mean = 0, sd = 1))
  rho <- array(dim = c(K, N), data = runif(N * K, min = 0.1, max = 0.9))
  c <- array(dim = c(2, N), data = runif(N * 2, min = 1, max = 2))
  eta <- array(dim = c(K, N, T), data = rchisq(N * T * K, df = 1) - 1)
  f <- array(dim = c(3, T), data = rnorm(n = 3 * T, mean = 0, sd = 1))
  #Calculate u
  {
    sigma2 <- 0.5 * (rchisq(N, df = 1) + 1)
    v <- array(dim = c(N, T)); for(i in 1:N) v[i, ] <- rnorm(n = T, mean = 0, sd = sqrt(sigma2[i]));
    u <- array(dim = c(N, T), data = v); b = 1; 
    for(i in 2:(N - 1)) u[i, ] <- (1 + b ^ 2) * v[i, ] + b * v[i + 1, ] + b * v[i - 1,]
  }
  #Calculate x
  x <- array(dim = c(K, N, T))
  for(t in 1:T) for(i in 1:N) for(k in 1:K) x[k, i, t] <- ifelse(
    t == 1,
    1 + sum(c[, i]*c(gamma[k, i],f[k, t]))+eta[k, i, t],
    1 + rho[k, i]* x[k, i, t - 1] + c[1, i] * gamma[k, i] + c[2, i] * f[k, t] + eta[k, i, t])
  #Calculate y.ctfl
  y.ctfl <- array(dim = c(N, T))
  for(t in 1:T) for(i in 1:N) 
    y.ctfl[i, t] <- sum(x[, i, t] * beta)+sum(gamma[ ,i] * f[, t]) + u[i, t]
  #Calculate y.actl
  y.actl <- y.ctfl; for(t in (T0 + 1):T) y.actl[1, t] <- y.actl[1, t] + 1
  #return data
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
    sigma2 = sigma2,
    v = v,
    f = f,
    u = u,
    x = x,
    y.ctfl = y.ctfl,
    y.actl = y.actl
  ))
}

BS <- 1:100
T0 <- seq(10, 90, 10)
N <- seq(60, 90, 10)
DGP1_data  <- list()
set.seed(1)
for(n in N){
  T_list <- list()
  for(t0 in T0){
    BS_list <- list()
    for(bs in BS){
      cat("N = ", n, ", T0 = ", t0, ", BS = ", bs, "\n")
      BS_list[[bs]] <- DGP1(n, t0)
    }
    T_list[[t0]] <- BS_list
  }
  DGP1_data[[n]] <- T_list
}
save(DGP1_data, file = paste("DGP1_data(", max(BS), ").Rdata", sep = ""))
