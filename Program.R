rm(list = ls())
library(ggplot2)
source("R/.Source.R")

# Default Setting
N <- c(40,70)
T0 <- c(40,70)
T1 <- 10
BS <- 1:10
nthreads <- 2
seed <- 1
Process <- 1:7

for(i in Process){
  file <- paste0("DGP", i, "_data.Rdata")
  if(!file %in% dir(path = "data")){
    object <- getQCM(getDGP(
      Process = i,
      nthreads = nthreads,
      seed = seed,
      N = N,
      T0 = T0,
      T1 = T1,
      BS = BS),nthreads = nthreads, seed = seed)
    save(object, file = paste0("data/", file))
  }else load(paste("data/", file, sep = ""))
  
  result <- getIndicator(
    data = object, 
    Process = i, 
    N = N, 
    T0 = T0, 
    T1 = T1,
    BS = BS)
  openxlsx::write.xlsx(result, file = paste0("DGP", i, ".xlsx"))
  
  for(n in N) for(t0 in T0){
    getPlot(result, N = n, T0 = t0, T1 = T1)
    ggsave(paste0("photo/DGP", i, ",N=", n, ",T0=", t0, ".png"))}
}
